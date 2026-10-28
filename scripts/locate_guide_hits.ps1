<#
  locate_guide_hits.ps1

  Find every exact occurrence of a 20 nt guide in the genome (both strands)
  and report chromosome + coordinate + strand.

  Usage:
    powershell -NoProfile -ExecutionPolicy Bypass -File scripts\locate_guide_hits.ps1 -Guide GTTTCCGCCACATTGATGGA
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Guide,
    [string]$GenomeGz
)

$ErrorActionPreference = 'Stop'
if (-not $GenomeGz) {
    $repoRoot = Split-Path -Parent $PSScriptRoot
    $GenomeGz = Join-Path $repoRoot 'data\raw\GCF_014633365.1_ASM1463336v1_genomic.fna.gz'
}

$g = $Guide.ToUpper()
$map = @{ 'A' = [char]'T'; 'C' = [char]'G'; 'G' = [char]'C'; 'T' = [char]'A'; 'N' = [char]'N' }
$rcChars = New-Object 'char[]' $g.Length
for ($i = 0; $i -lt $g.Length; $i++) {
    $ch = [string]$g[$g.Length - 1 - $i]
    if ($map.ContainsKey($ch)) { $rcChars[$i] = $map[$ch] } else { $rcChars[$i] = [char]'N' }
}
$rc = -join $rcChars

$hits = New-Object System.Collections.Generic.List[string]

function Invoke-Search {
    $seq = $sb.ToString().ToUpper()
    $len = $seq.Length
    $p = 0
    while (($p = $seq.IndexOf($g, $p, [System.StringComparison]::Ordinal)) -ge 0) {
        $hits.Add(('{0}  plus   {1}  {2}' -f $name, ($p + 1), $seq.Substring($p, [Math]::Min(30, $len - $p))))
        $p++
    }
    $p = 0
    while (($p = $seq.IndexOf($rc, $p, [System.StringComparison]::Ordinal)) -ge 0) {
        $hits.Add(('{0}  minus  {1}  {2}' -f $name, ($p + $g.Length), $seq.Substring($p, [Math]::Min(30, $len - $p))))
        $p++
    }
}

$fs = [System.IO.File]::OpenRead($GenomeGz)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$name = ''
$sb = New-Object System.Text.StringBuilder
try {
    while (($line = $sr.ReadLine()) -ne $null) {
        if ($line.StartsWith('>')) {
            if ($sb.Length -gt 0) { Invoke-Search }
            $name = ($line.Substring(1) -split ' ')[0]
            $sb = New-Object System.Text.StringBuilder
            continue
        }
        [void]$sb.Append($line.Trim())
    }
    Invoke-Search
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

'guide    : ' + $g
'revcomp  : ' + $rc
'hits     : ' + $hits.Count
$hits | ForEach-Object { $_ }
