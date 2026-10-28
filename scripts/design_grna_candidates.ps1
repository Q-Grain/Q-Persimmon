<#
  design_grna_candidates.ps1

  Stage 1 (default):
    - take the fixed target gene list
    - locate gene / mRNA / CDS from the GFF3
    - pick the longest transcript per gene
    - take the first two CDS segments (in transcript direction)
    - extract their genomic sequence
    - generate guide candidates (20 nt + NGG) on both strands
    - apply basic filters (GC, poly-T, poly-G, BsaI/BsmBI sites)
    - write TSV

  Stage 2 (-Offtarget):
    - stream the whole genome
    - for each candidate count sites with <= MaxMismatch mismatches
      (8-nt seed on both strands, then full 20-nt comparison + PAM check)
    - write TSV with offtarget columns

  ASCII only (PowerShell 5.1 decodes BOM-less files as ANSI).

  Usage:
    powershell -NoProfile -ExecutionPolicy Bypass -File scripts\design_grna_candidates.ps1
    powershell -NoProfile -ExecutionPolicy Bypass -File scripts\design_grna_candidates.ps1 -Offtarget
#>

[CmdletBinding()]
param(
    [string]$GenomeGz,
    [string]$GffGz,
    [string]$OutTsv,
    [string]$OutOffTsv,
    [switch]$Offtarget,
    [int]$MaxMismatch = 3,
    [int]$MaxGuidesPerGene = 12
)

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot

trap {
    'ERROR line ' + $_.InvocationInfo.ScriptLineNumber + ' : ' + $_.Exception.GetType().Name + ' : ' + $_.Exception.Message
    '  at: ' + ($_.InvocationInfo.Line).Trim()
    exit 1
}

if (-not $GenomeGz) { $GenomeGz = Join-Path $repoRoot 'data\raw\GCF_014633365.1_ASM1463336v1_genomic.fna.gz' }
if (-not $GffGz)    { $GffGz    = Join-Path $repoRoot 'data\raw\GCF_014633365.1_ASM1463336v1_genomic.gff.gz' }
if (-not $OutTsv)   { $OutTsv   = Join-Path $repoRoot 'data\targets\grna_candidates.tsv' }
if (-not $OutOffTsv){ $OutOffTsv= Join-Path $repoRoot 'data\targets\grna_offtarget.tsv' }

# ---------------------------------------------------------------- target genes
# Chosen from the 36 core targets: starch-degradation side only (knock-out makes sense),
# annotated as active enzymes located in the plastid/amyloplast.
$targets = @(
    [pscustomobject]@{ Gene = 'LOC127799073'; Role = 'BAM1';   Note = 'beta-amylase 1, chloroplastic' }
    [pscustomobject]@{ Gene = 'LOC127786821'; Role = 'BAM3';   Note = 'beta-amylase 3, chloroplastic' }
    [pscustomobject]@{ Gene = 'LOC127807677'; Role = 'ISA3';   Note = 'isoamylase 3, chloroplastic' }
    [pscustomobject]@{ Gene = 'LOC127803838'; Role = 'PHS-L1'; Note = 'glucan phosphorylase L-1, plastidic' }
    [pscustomobject]@{ Gene = 'LOC127806752'; Role = 'PUL1';   Note = 'pullulanase 1, chloroplastic' }
)
$targetIds = @{}
foreach ($t in $targets) { $targetIds[$t.Gene] = $t }

# ---------------------------------------------------------------- helpers
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

function Get-Attr {
    param([string]$field)
    $h = @{}
    foreach ($kv in $field.Split(';')) {
        $p = $kv.Split('=', 2)
        if ($p.Length -eq 2) { $h[$p[0]] = $p[1] }
    }
    return $h
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
    if ($g -cmatch 'TTTT') { return $false }          # Pol III terminator
    if ($g -cmatch 'GGGG') { return $false }          # G-quadruplex risk
    if ($g -cmatch 'GGTCTC' -or $g -cmatch 'CGTCTC') { return $false }  # BsaI / BsmBI
    return $true
}

# ---------------------------------------------------------------- parse GFF
$lines = @{}
$mrnaCount = @{}    # rnaId -> cds count placeholder
$cdsMap = @{}       # rnaId -> list of segment objects
$gene2rna = @{}     # geneId -> list of rnaId
$rna2gene = @{}     # rnaId -> geneId

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

        $a = Get-Attr $c[8]
        if ($type -eq 'gene') {
            if ($a.ContainsKey('gene') -and $targetIds.ContainsKey($a['gene'])) {
                $id = $a['ID']
                $lines[$id] = [pscustomobject]@{ Gene = $a['gene']; SeqId = $c[0]; Strand = $c[6] }
            }
        }
        elseif ($type -eq 'mRNA') {
            if (-not $a.ContainsKey('Parent')) { continue }
            if (-not $lines.ContainsKey($a['Parent'])) { continue }
            $rid = $a['ID']
            $rna2gene[$rid] = $a['Parent']
            if (-not $gene2rna.ContainsKey($a['Parent'])) { $gene2rna[$a['Parent']] = New-Object System.Collections.Generic.List[string] }
            $gene2rna[$a['Parent']].Add($rid)
            $cdsMap[$rid] = New-Object System.Collections.Generic.List[object]
        }
        elseif ($type -eq 'CDS') {
            if (-not $a.ContainsKey('Parent')) { continue }
            if (-not $cdsMap.ContainsKey($a['Parent'])) { continue }
            $cdsMap[$a['Parent']].Add([pscustomobject]@{
                Start  = [int64]$c[3]
                End    = [int64]$c[4]
                Strand = $c[6]
                SeqId  = $c[0]
            })
        }
    }
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

# pick longest transcript per gene
$plan = New-Object System.Collections.Generic.List[object]
foreach ($geneId in $lines.Keys) {
    $info = $lines[$geneId]
    if (-not $gene2rna.ContainsKey($geneId)) { continue }
    $best = $null; $bestLen = 0
    foreach ($rid in $gene2rna[$geneId]) {
        $segs = $cdsMap[$rid]
        if ($segs.Count -eq 0) { continue }
        $len = 0
        foreach ($s in $segs) { $len += ($s.End - $s.Start + 1) }
        if ($len -gt $bestLen) { $bestLen = $len; $best = $rid }
    }
    if ($null -eq $best) { continue }

    $segSrc = $cdsMap[$best]
    if ($null -eq $segSrc) { continue }

    $segArr = @()
    foreach ($s in $segSrc) { $segArr += $s }
    if ($segArr.Count -eq 0) { continue }

    if ($info.Strand -eq '+') { $segArr = @($segArr | Sort-Object -Property Start) }
    else                      { $segArr = @($segArr | Sort-Object -Property Start -Descending) }

    # first two CDS segments in transcript direction
    $use = @($segArr | Select-Object -First 2)
    $plan.Add([pscustomobject]@{
        Gene   = $info.Gene
        Role   = $targetIds[$info.Gene].Role
        Note   = $targetIds[$info.Gene].Note
        Rna    = $best
        SeqId  = $info.SeqId
        Strand = $info.Strand
        Segs   = $use
        CdsLen = $bestLen
    })
}

# ---------------------------------------------------------------- extract genomic sequence for needed seqIds
$need = @{}
foreach ($p in $plan) { foreach ($s in $p.Segs) { $need[$s.SeqId] = $true } }
$seqs = @{}

$fs = [System.IO.File]::OpenRead($GenomeGz)
$gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
$sr = New-Object System.IO.StreamReader($gz)
$curId = ''; $sb = $null
try {
    while (($line = $sr.ReadLine()) -ne $null) {
        if ($line.StartsWith('>')) {
            if ($null -ne $sb -and $need.ContainsKey($curId)) { $seqs[$curId] = $sb.ToString() }
            $curId = ($line.Substring(1) -split ' ')[0]
            $sb = New-Object System.Text.StringBuilder
        }
        else { [void]$sb.Append($line.Trim()) }
    }
    if ($null -ne $sb -and $need.ContainsKey($curId)) { $seqs[$curId] = $sb.ToString() }
}
finally { $sr.Close(); $gz.Close(); $fs.Close() }

# ---------------------------------------------------------------- generate candidates
$guides = New-Object System.Collections.Generic.List[object]

foreach ($p in $plan) {
    if (-not $seqs.ContainsKey($p.SeqId)) { continue }
    $chrom = $seqs[$p.SeqId]
    $idx = 0
    $added = 0

    foreach ($seg in $p.Segs) {
        if ($added -ge $MaxGuidesPerGene) { break }
        $idx++
        $sub = $chrom.Substring([int]($seg.Start - 1), [int]($seg.End - $seg.Start + 1)).ToUpper()
        $rc  = Get-RevComp $sub

        for ($i = 0; $i -le $sub.Length - 23; $i++) {
            if ($added -ge $MaxGuidesPerGene) { break }

            # sense strand: 20 nt then NGG
            if ($sub[$i + 21] -eq 'G' -and $sub[$i + 22] -eq 'G') {
                $g = $sub.Substring($i, 20)
                if (Test-GuideOk $g) {
                    $guides.Add([pscustomobject]@{
                        Gene = $p.Gene; Role = $p.Role; Exon = $idx
                        Chrom = $p.SeqId; Strand = '+'
                        Start = $seg.Start + $i; End = $seg.Start + $i + 19
                        Guide = $g; PAM = $sub.Substring($i + 20, 3)
                        GC = (Get-GC $g)
                    })
                    $added++
                }
            }

            # antisense strand: GG on the reverse complement
            if ($rc[$i + 21] -eq 'G' -and $rc[$i + 22] -eq 'G') {
                $g = $rc.Substring($i, 20)
                if (Test-GuideOk $g) {
                    $startOnChrom = $seg.End - ($i + 22) + 1
                    $guides.Add([pscustomobject]@{
                        Gene = $p.Gene; Role = $p.Role; Exon = $idx
                        Chrom = $p.SeqId; Strand = '-'
                        Start = $startOnChrom; End = $startOnChrom + 19
                        Guide = $g; PAM = $rc.Substring($i + 20, 3)
                        GC = (Get-GC $g)
                    })
                    $added++
                }
            }
        }
    }
}

$guides | Select-Object Gene, Role, Exon, Chrom, Start, End, Strand, Guide, PAM, GC |
    Export-Csv -LiteralPath $OutTsv -Delimiter "`t" -NoTypeInformation -Encoding UTF8

''
'== design_grna_candidates (stage 1) =='
'genes planned : ' + $plan.Count
'guides        : ' + $guides.Count
'output        : ' + $OutTsv
foreach ($p in $plan) {
    '{0,-12} {1,-6} {2,-14} cds={3,-6} segs={4}' -f $p.Gene, $p.Role, $p.SeqId, $p.CdsLen, $p.Segs.Count
}

# ---------------------------------------------------------------- stage 2: offtarget
if ($Offtarget) {

    $rows = @(Import-Csv -LiteralPath $OutTsv -Delimiter "`t")

    $forms = New-Object System.Collections.Generic.List[object]
    $i = 0
    foreach ($r in $rows) {
        $plus = $r.Guide
        $minus = Get-RevComp $r.Guide
        $forms.Add([pscustomobject]@{
            Idx = $i; Gene = $r.Gene; Role = $r.Role; Exon = $r.Exon
            Chrom = $r.Chrom; Start = $r.Start; Strand = $r.Strand
            Guide = $r.Guide; PAM = $r.PAM; GC = $r.GC
            Plus = $plus; PlusSeed = $plus.Substring(0, 8)
            Minus = $minus; MinusSeed = $minus.Substring(0, 8)
            End = $r.End
            Hit0 = 0; HitLe = 0
        })
        $i++
    }

    function Get-Mismatch {
        param([string]$Genome, [int]$Pos, [string]$Guide)
        $m = 0
        for ($k = 0; $k -lt 20; $k++) {
            if ($Genome[$Pos + $k] -ne $Guide[$k]) {
                $m++
                if ($m -gt 3) { return 99 }
            }
        }
        return $m
    }

    function Invoke-Scan {
        param([string]$Name, [string]$Seq)
        $len = $Seq.Length
        foreach ($f in $forms) {
            # site on the plus strand: 20 nt then NGG
            $p = 0
            while (($p = $Seq.IndexOf($f.PlusSeed, $p, [System.StringComparison]::Ordinal)) -ge 0) {
                if (($p + 22) -lt $len) {
                    if ($Seq[$p + 21] -eq 'G' -and $Seq[$p + 22] -eq 'G') {
                        $mm = Get-Mismatch $Seq $p $f.Plus
                        if ($mm -ne 99) {
                            if ($mm -eq 0) { $f.Hit0 = $f.Hit0 + 1 }
                            $f.HitLe = $f.HitLe + 1
                        }
                    }
                }
                $p++
            }

            # site on the minus strand: CCN then 20 nt (revcomp form)
            $p = 0
            while (($p = $Seq.IndexOf($f.MinusSeed, $p, [System.StringComparison]::Ordinal)) -ge 0) {
                if ($p -ge 3 -and ($p + 19) -lt $len) {
                    if ($Seq[$p - 3] -eq 'C' -and $Seq[$p - 2] -eq 'C') {
                        $mm = Get-Mismatch $Seq $p $f.Minus
                        if ($mm -ne 99) {
                            if ($mm -eq 0) { $f.Hit0 = $f.Hit0 + 1 }
                            $f.HitLe = $f.HitLe + 1
                        }
                    }
                }
                $p++
            }
        }
    }

    $fs = [System.IO.File]::OpenRead($GenomeGz)
    $gz = New-Object System.IO.Compression.GZipStream($fs, [System.IO.Compression.CompressionMode]::Decompress)
    $sr = New-Object System.IO.StreamReader($gz)
    $chromName = ''; $sb2 = $null; $chromCount = 0
    try {
        while ($true) {
            $line = $sr.ReadLine()
            if ($null -ne $line -and $line.StartsWith('>')) {
                if ($null -ne $sb2) { Invoke-Scan $chromName $sb2.ToString().ToUpper() }
                $chromName = ($line.Substring(1) -split ' ')[0]
                $sb2 = New-Object System.Text.StringBuilder
                $chromCount++
                continue
            }
            if ($null -eq $line) { break }
            [void]$sb2.Append($line.Trim())
        }
        if ($null -ne $sb2) { Invoke-Scan $chromName $sb2.ToString().ToUpper() }
    }
    finally { $sr.Close(); $gz.Close(); $fs.Close() }

    $forms | Select-Object Gene, Role, Exon, Chrom, Start, End, Strand, Guide, PAM, GC, Hit0, HitLe |
        Export-Csv -LiteralPath $OutOffTsv -Delimiter "`t" -NoTypeInformation -Encoding UTF8

    ''
    '== stage 2 (offtarget, mismatch <= ' + $MaxMismatch + ') =='
    'chromosomes scanned : ' + $chromCount
    'output              : ' + $OutOffTsv

    # ------------------------------------------------ stage 3: short list
    # de-duplicate within a gene/strand (start positions closer than 10 nt),
    # prefer fewer offtargets and GC near 50%%, keep at most 3 per gene
    $all = @(Import-Csv -LiteralPath $OutOffTsv -Delimiter "`t")
    $shortList = New-Object System.Collections.Generic.List[object]
    foreach ($grp in ($all | Group-Object -Property Role)) {
        $kept = @()
        $sorted = $grp.Group | Sort-Object @{ Expression = { [int]$_.HitLe } }, @{ Expression = { [math]::Abs([int]$_.GC - 50) } }, @{ Expression = { [int]$_.Start } }
        foreach ($r in $sorted) {
            $ok = $true
            foreach ($k in $kept) {
                if ($k.Strand -eq $r.Strand -and $k.Chrom -eq $r.Chrom -and
                    [math]::Abs([int]$k.Start - [int]$r.Start) -lt 10) { $ok = $false; break }
            }
            if ($ok) { $kept += $r; $shortList.Add($r) }
            if ($kept.Count -ge 3) { break }
        }
    }
    $shortPath = Join-Path (Split-Path -Parent $OutOffTsv) 'grna_shortlist.tsv'
    $shortList | Select-Object Gene, Role, Exon, Chrom, Start, End, Strand, Guide, PAM, GC, Hit0, HitLe |
        Export-Csv -LiteralPath $shortPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8

    ''
    '== stage 3 (short list) =='
    'entries : ' + $shortList.Count
    'output  : ' + $shortPath
    foreach ($r in $shortList) {
        '{0,-6} exon{1} {2} {3,10} {4} PAM={5} GC={6} hit0={7} hit<=3={8}' -f `
            $r.Role, $r.Exon, $r.Strand, $r.Start, $r.Guide, $r.PAM, $r.GC, $r.Hit0, $r.HitLe
    }
}
