<#
# verify_raw_files.ps1
#
# Verify the NCBI download set in data/raw:
#   - gzip stream can be read to the end (implies CRC check passes)
#   - record count and total length for FASTA files
#   - line count and feature breakdown for the GFF3 file
#
# ASCII only (PowerShell 5.1 decodes BOM-less files as ANSI).
#
# Usage:
#   powershell -NoProfile -ExecutionPolicy Bypass -File scripts\verify_raw_files.ps1
#   powershell -NoProfile -ExecutionPolicy Bypass -File scripts\verify_raw_files.ps1 -RawDir D:\data\raw
#>

[CmdletBinding()]
param(
    [string]$RawDir,
    [string]$OutFile
)

$ErrorActionPreference = 'Stop'

if (-not $RawDir) {
    $repoRoot = Split-Path -Parent $PSScriptRoot
    $RawDir = Join-Path $repoRoot 'data\raw'
}
if (-not (Test-Path -LiteralPath $RawDir)) { throw "RawDir not found: $RawDir" }

$report = New-Object System.Collections.Generic.List[string]
$report.Add('== verify_raw_files ==')
$report.Add('dir: ' + $RawDir)
$report.Add('')

function Measure-GzFasta {
    param([string]$Path, [string]$Tag)

    $fs = [System.IO.File]::OpenRead($Path)
    try {
        $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
        $sr = New-Object System.IO.StreamReader($gz)
        $n = 0
        $nt = 0
        $first = ''
        while (($l = $sr.ReadLine()) -ne $null) {
            if ($l.StartsWith('>')) {
                $n++
                if ($n -eq 1) { $first = $l.Substring(0, [Math]::Min(70, $l.Length)) }
            }
            else { $nt += $l.Trim().Length }
        }
        $sr.Close()
        $gz.Close()
    }
    finally { $fs.Close() }

    $report.Add(('{0,-28} records={1,-8} total_nt={2,-14} {3}' -f $Tag, $n, $nt, $first))
}

function Measure-GzGff {
    param([string]$Path, [string]$Tag)

    $fs = [System.IO.File]::OpenRead($Path)
    try {
        $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
        $sr = New-Object System.IO.StreamReader($gz)
        $lines = 0
        $feat = @{}
        $seqids = @{}
        while (($l = $sr.ReadLine()) -ne $null) {
            $lines++
            if ($l.StartsWith('#')) { continue }
            $c = $l.Split("`t")
            if ($c.Length -lt 8) { continue }
            if ($feat.ContainsKey($c[2])) { $feat[$c[2]]++ } else { $feat[$c[2]] = 1 }
            if ($seqids.ContainsKey($c[0])) { $seqids[$c[0]]++ } else { $seqids[$c[0]] = 1 }
        }
        $sr.Close()
        $gz.Close()
    }
    finally { $fs.Close() }

    $report.Add(('{0,-28} lines={1,-8} seqids={2}' -f $Tag, $lines, $seqids.Count))
    foreach ($k in ($feat.Keys | Sort-Object { -$feat[$_] })) {
        $report.Add(('    {0,-14} {1}' -f $k, $feat[$k]))
    }
}

$fastaFiles = Get-ChildItem -LiteralPath $RawDir -Filter '*_from_genomic.fna.gz' -File -ErrorAction SilentlyContinue
foreach ($f in (Get-ChildItem -LiteralPath $RawDir -Filter '*.fna.gz' -File)) { Measure-GzFasta $f.FullName $f.Name }
foreach ($f in (Get-ChildItem -LiteralPath $RawDir -Filter '*.faa.gz' -File)) { Measure-GzFasta $f.FullName $f.Name }
foreach ($f in (Get-ChildItem -LiteralPath $RawDir -Filter '*.gff.gz' -File)) { Measure-GzGff   $f.FullName $f.Name }

if ($OutFile) {
    [System.IO.File]::WriteAllLines($OutFile, $report, (New-Object System.Text.UTF8Encoding($false)))
    'report written to ' + $OutFile
}
$report | ForEach-Object { $_ }
