<#
  design_grna_single_gene.ps1

  Single-gene guide design with a strict uniqueness filter.

  Difference from design_grna_candidates.ps1:
    - scans ALL CDS segments, not just the first two
    - keeps only guides whose 20 nt occurs EXACTLY ONCE in the genome
    - also reports the number of sites with <= MaxMismatch mismatches

  Made for genes sitting inside duplicated clusters (e.g. BAM1 / LOC127799073).

  Usage:
    powershell -NoProfile -ExecutionPolicy Bypass -File scripts\design_grna_single_gene.ps1 -Gene LOC127799073
#>

[CmdletBinding()]
param(
    [Parameter(Mandatory = $true)][string]$Gene,
    [string]$GenomeGz,
    [string]$GffGz,
    [string]$OutTsv,
    [int]$MaxMismatch = 3,
    [int]$MaxPerExon = 15
)

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot

if (-not $GenomeGz) { $GenomeGz = Join-Path $repoRoot 'data\raw\GCF_014633365.1_ASM1463336v1_genomic.fna.gz' }
if (-not $GffGz)    { $GffGz    = Join-Path $repoRoot 'data\raw\GCF_014633365.1_ASM1463336v1_genomic.gff.gz' }
if (-not $OutTsv)   { $OutTsv   = Join-Path $repoRoot ("data\targets\grna_{0}.tsv" -f $Gene) }

trap {
    'ERROR line ' + $_.InvocationInfo.ScriptLineNumber + ' : ' + $_.Exception.Message
    exit 1
}

function Get-RevComp {
    param([string]$s)
    $map = @{ 'A' = [char]'T'; 'C' = [char]'G'; 'G' = [char]'C'; 'T' = [char]'A'; 'N' = [char]'N' }
    $src = $s.ToUpper()
    $n = $src.Length
    $out = New-Object 'char[]' $n
    for ($i = 0; $i -lt $n; $i++) {
        $ch = [string]$src[$n - 1 - $i]
        if ($map.ContainsKey($ch)) { $out[$i] = $map[$ch] } else { $out[$i] = [char]'N' }
    }
    return -join $out
}

function Get-GC {
    param([string]$s)
    $gc = 0
    foreach ($ch in $s.ToCharArray()) { if ($ch -eq 'G' -or $ch -eq 'C') { $gc++ } }
    return [math]::Round(100.0 * $gc / $s.Length, 1)
}

function Test-GuideOk {
    param([string]$g)
    if (-not ($g -cmatch '^[ACGT]{20}$')) { return $false }
    $gc = Get-GC $g
    if ($gc -lt 40 -or $gc -gt 70) { return $false }
    if ($g -cmatch 'TTTT') { return $false }
    if ($g -cmatch 'GGGG') { return $false }
    if ($g -cmatch 'GGTCTC' -or $g -cmatch 'CGTCTC') { return $false }
    return $true
}

# ---------------------------------------------------------------- GFF: all CDS of the gene
$rna2gene = @{}
$cdsMap = @{}
$geneInfo = $null
$geneField = 'gene-' + $Gene

$fs = [System.IO.File]::OpenRead($GffGz)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
try {
    while (($line = $sr.ReadLine()) -ne $null) {
        if ($line.Length -eq 0 -or $line[0] -eq '#') { continue }
        $c = $line.Split("`t")
        if ($c.Length -lt 9) { continue }
        $type = $c[2]
        if ($type -ne 'gene' -and $type -ne 'mRNA' -and $type -ne 'CDS') { continue }

        $attr = @{}
        foreach ($kv in $c[8].Split(';')) { $p = $kv.Split('=', 2); if ($p.Length -eq 2) { $attr[$p[0]] = $p[1] } }

        if ($type -eq 'gene') {
            if ($attr['ID'] -eq $geneField) {
                $geneInfo = [pscustomobject]@{ SeqId = $c[0]; Strand = $c[6]; Start = [int64]$c[3]; End = [int64]$c[4] }
            }
        }
        elseif ($type -eq 'mRNA') {
            if ($attr['Parent'] -eq $geneField) {
                $rna2gene[$attr['ID']] = $Gene
                $cdsMap[$attr['ID']] = New-Object System.Collections.Generic.List[object]
            }
        }
        elseif ($type -eq 'CDS') {
            $par = $attr['Parent']
            if ($cdsMap.ContainsKey($par)) {
                $cdsMap[$par].Add([pscustomobject]@{ Start = [int64]$c[3]; End = [int64]$c[4]; Strand = $c[6]; SeqId = $c[0] })
            }
        }
    }
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

if ($null -eq $geneInfo) { throw "gene not found in GFF: $Gene" }

$best = $null; $bestLen = 0
foreach ($rid in $cdsMap.Keys) {
    $len = 0
    foreach ($s in $cdsMap[$rid]) { $len += ($s.End - $s.Start + 1) }
    if ($len -gt $bestLen) { $bestLen = $len; $best = $rid }
}
if ($null -eq $best) { throw "no CDS for $Gene" }

$segArr = @()
foreach ($s in $cdsMap[$best]) { $segArr += $s }
if ($geneInfo.Strand -eq '+') { $segArr = @($segArr | Sort-Object -Property Start) }
else                          { $segArr = @($segArr | Sort-Object -Property Start -Descending) }

# ---------------------------------------------------------------- genomic sequence of the gene region
$needSeq = $geneInfo.SeqId
$chrom = ''

$fs = [System.IO.File]::OpenRead($GenomeGz)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$curId = ''; $sb = $null
try {
    while (($line = $sr.ReadLine()) -ne $null) {
        if ($line.StartsWith('>')) {
            if ($null -ne $sb -and $curId -eq $needSeq) { $chrom = $sb.ToString().ToUpper() }
            $curId = ($line.Substring(1) -split ' ')[0]
            $sb = New-Object System.Text.StringBuilder
            continue
        }
        if ($curId -eq $needSeq) { [void]$sb.Append($line.Trim()) }
    }
    if ($null -ne $sb -and $curId -eq $needSeq) { $chrom = $sb.ToString().ToUpper() }
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

# ---------------------------------------------------------------- candidates from every CDS segment
$guides = New-Object System.Collections.Generic.List[object]
$exonNo = 0
foreach ($seg in $segArr) {
    $exonNo++
    $sub = $chrom.Substring([int]($seg.Start - 1), [int]($seg.End - $seg.Start + 1))
    $rc = Get-RevComp $sub
    $added = 0

    for ($i = 0; $i -le $sub.Length - 23; $i++) {
        if ($added -ge $MaxPerExon) { break }

        if ($sub[$i + 21] -eq 'G' -and $sub[$i + 22] -eq 'G') {
            $g = $sub.Substring($i, 20)
            if (Test-GuideOk $g) {
                $guides.Add([pscustomobject]@{
                    Gene = $Gene; Exon = $exonNo; Chrom = $seg.SeqId; Strand = '+'
                    Start = $seg.Start + $i; End = $seg.Start + $i + 19
                    Guide = $g; PAM = $sub.Substring($i + 20, 3); GC = (Get-GC $g)
                    Exact = 0; Le = 0
                })
                $added++
            }
        }

        if ($rc[$i + 21] -eq 'G' -and $rc[$i + 22] -eq 'G') {
            $g = $rc.Substring($i, 20)
            if (Test-GuideOk $g) {
                $startOnChrom = $seg.End - ($i + 22) + 1
                $guides.Add([pscustomobject]@{
                    Gene = $Gene; Exon = $exonNo; Chrom = $seg.SeqId; Strand = '-'
                    Start = $startOnChrom; End = $startOnChrom + 19
                    Guide = $g; PAM = $rc.Substring($i + 20, 3); GC = (Get-GC $g)
                    Exact = 0; Le = 0
                })
                $added++
            }
        }
    }
}

# ---------------------------------------------------------------- genome-wide scan
function Get-Mismatch {
    param([string]$Genome, [int]$Pos, [string]$Guide)
    $m = 0
    for ($k = 0; $k -lt 20; $k++) {
        if ($Genome[$Pos + $k] -ne $Guide[$k]) { $m++; if ($m -gt 3) { return 99 } }
    }
    return $m
}

$probe = New-Object System.Collections.Generic.List[object]
foreach ($g in $guides) {
    $rc = Get-RevComp $g.Guide
    $probe.Add([pscustomobject]@{
        Ref = $g
        Plus = $g.Guide; PlusSeed = $g.Guide.Substring(0, 8)
        Minus = $rc;     MinusSeed = $rc.Substring(0, 8)
    })
}

function Invoke-Scan2 {
    param([string]$Seq)
    $len = $Seq.Length
    foreach ($p in $probe) {
        $i = 0
        while (($i = $Seq.IndexOf($p.PlusSeed, $i, [System.StringComparison]::Ordinal)) -ge 0) {
            if (($i + 22) -lt $len) {
                if ($Seq[$i + 21] -eq 'G' -and $Seq[$i + 22] -eq 'G') {
                    $mm = Get-Mismatch $Seq $i $p.Plus
                    if ($mm -ne 99) {
                        if ($mm -eq 0) { $p.Ref.Exact = $p.Ref.Exact + 1 }
                        $p.Ref.Le = $p.Ref.Le + 1
                    }
                }
            }
            $i++
        }
        $i = 0
        while (($i = $Seq.IndexOf($p.MinusSeed, $i, [System.StringComparison]::Ordinal)) -ge 0) {
            if ($i -ge 3 -and ($i + 19) -lt $len) {
                if ($Seq[$i - 3] -eq 'C' -and $Seq[$i - 2] -eq 'C') {
                    $mm = Get-Mismatch $Seq $i $p.Minus
                    if ($mm -ne 99) {
                        if ($mm -eq 0) { $p.Ref.Exact = $p.Ref.Exact + 1 }
                        $p.Ref.Le = $p.Ref.Le + 1
                    }
                }
            }
            $i++
        }
    }
}

$fs = [System.IO.File]::OpenRead($GenomeGz)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$chrName = ''; $sb2 = New-Object System.Text.StringBuilder
try {
    while ($true) {
        $line = $sr.ReadLine()
        if ($null -eq $line) { break }
        if ($line.StartsWith('>')) {
            if ($sb2.Length -gt 0) { Invoke-Scan2 $sb2.ToString().ToUpper() }
            $chrName = ($line.Substring(1) -split ' ')[0]
            $sb2 = New-Object System.Text.StringBuilder
            continue
        }
        [void]$sb2.Append($line.Trim())
    }
    if ($sb2.Length -gt 0) { Invoke-Scan2 $sb2.ToString().ToUpper() }
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

$guides | Select-Object Gene, Exon, Chrom, Start, End, Strand, Guide, PAM, GC, Exact, Le |
    Export-Csv -LiteralPath $OutTsv -Delimiter "`t" -NoTypeInformation -Encoding UTF8

''
'== design_grna_single_gene =='
'gene        : ' + $Gene + '  (' + $geneInfo.SeqId + ':' + $geneInfo.Start + '-' + $geneInfo.End + ' ' + $geneInfo.Strand + ')'
'transcript  : ' + $best + '  CDS=' + $bestLen + '  exons(scanned)=' + $segArr.Count
'candidates  : ' + $guides.Count
'unique only : ' + (@($guides | Where-Object { $_.Exact -eq 1 }).Count)
'output      : ' + $OutTsv
''
'-- candidates with exactly one genomic match --'
$guides | Where-Object { $_.Exact -eq 1 } | Sort-Object Exon, Le, Start | ForEach-Object {
    'exon{0,-2} {1} {2,10} {3} PAM={4} GC={5} le={6}' -f $_.Exon, $_.Strand, $_.Start, $_.Guide, $_.PAM, $_.GC, $_.Le
}
''
'-- per-exon uniqueness summary --'
$guides | Group-Object Exon | Sort-Object { [int]$_.Name } | ForEach-Object {
    $u = (@($_.Group | Where-Object { $_.Exact -eq 1 })).Count
    'exon{0,-2} total={1,-3} unique={2}' -f $_.Name, $_.Count, $u
}
