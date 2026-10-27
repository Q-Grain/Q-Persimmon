<#
  extract_starch_genes.ps1

  Pull starch / sugar metabolism gene models out of an NCBI
  *_cds_from_genomic.fna file and emit:
     1) data/targets/starch_genes.tsv       one row per gene (longest CDS as representative)
     2) data/targets/starch_genes.cds.fna   representative CDS nucleotide sequences

  ASCII only on purpose: PowerShell 5.1 decodes BOM-less files as ANSI.

  Usage:
     powershell -NoProfile -ExecutionPolicy Bypass -File scripts\extract_starch_genes.ps1
     powershell -NoProfile -ExecutionPolicy Bypass -File scripts\extract_starch_genes.ps1 -Fasta <path> -OutDir <dir>
#>

[CmdletBinding()]
param(
    [string]$Fasta,
    [string]$OutDir
)

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot

if (-not $Fasta) {
    $Fasta = Join-Path $repoRoot 'GCF_014633365.1_ASM1463336v1_cds_from_genomic.fna\GCF_014633365.1_ASM1463336v1_cds_from_genomic.fna'
}
if (-not $OutDir) {
    $OutDir = Join-Path $repoRoot 'data\targets'
}

if (-not (Test-Path -LiteralPath $Fasta)) { throw "FASTA not found: $Fasta" }
if (-not (Test-Path -LiteralPath $OutDir)) { [void](New-Item -ItemType Directory -Path $OutDir -Force) }

# ---------------------------------------------------------------- classification
# Order matters: first match wins, so put the specific patterns first.
#   Group = Inhibit  -> starch breakdown, candidate knock-out targets
#   Group = Enhance  -> starch synthesis, candidate over-expression / de-repression targets
#   Group = Source   -> sucrose supply / transport, secondary support targets
$rules = @(
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'BAM';    Pattern = 'beta-amylase';                                                              Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'AMY';    Pattern = 'alpha-amylase';                                                             Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'ISA';    Pattern = 'isoamylase';                                                                Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'PUL';    Pattern = 'pullulanase';                                                               Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'PHS';    Pattern = 'alpha-1,4 glucan phosphorylase|alpha-glucan phosphorylase|starch phosphorylase'; Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'AGD';    Pattern = 'alpha-glucosidase';                                                         Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'DPE';    Pattern = 'disproportionating enzyme';                                                 Exclude = '' }
    [pscustomobject]@{ Group = 'Inhibit'; Role = 'MEX';    Pattern = 'maltose exporter|maltose-transport';                                        Exclude = '' }
    [pscustomobject]@{ Group = 'Enhance'; Role = 'GBSS';   Pattern = 'granule-bound starch synthase';                                             Exclude = '' }
    [pscustomobject]@{ Group = 'Enhance'; Role = 'SSS';    Pattern = 'starch synthase';                                                           Exclude = '' }
    [pscustomobject]@{ Group = 'Enhance'; Role = 'SBE';    Pattern = '1,4-alpha-glucan-branching|alpha-glucan-branching|(?<!de)branching enzyme'; Exclude = '' }
    [pscustomobject]@{ Group = 'Enhance'; Role = 'AGPase'; Pattern = 'glucose-1-phosphate adenylyltransferase';                                   Exclude = '' }
    [pscustomobject]@{ Group = 'Source';  Role = 'SPS';    Pattern = 'sucrose-phosphate synthase';                                                Exclude = '' }
    [pscustomobject]@{ Group = 'Source';  Role = 'SUS';    Pattern = 'sucrose synthase';                                                          Exclude = '' }
    [pscustomobject]@{ Group = 'Source';  Role = 'UGP';    Pattern = 'UDP-glucose pyrophosphorylase';                                             Exclude = '' }
    [pscustomobject]@{ Group = 'Source';  Role = 'INV';    Pattern = 'invertase|beta-fructofuranosidase';                                         Exclude = 'inhibitor' }
    [pscustomobject]@{ Group = 'Source';  Role = 'STP';    Pattern = 'sugar transporter|sucrose transporter|SWEET|monosaccharide transporter|tonoplast monosaccharide'; Exclude = 'SWEETIE' }
)

# RefSeq accession -> readable chromosome name (from the NCBI assembly report)
$chromNames = @{
    'NC_068338.1' = 'chr1';  'NC_068339.1' = 'chr2';  'NC_068340.1' = 'chr3'
    'NC_068341.1' = 'chr4';  'NC_068342.1' = 'chr5';  'NC_068343.1' = 'chr6'
    'NC_068344.1' = 'chr7';  'NC_068345.1' = 'chr8';  'NC_068346.1' = 'chr9'
    'NC_068347.1' = 'chr10'; 'NC_068348.1' = 'chr11'; 'NC_068349.1' = 'chr12'
    'NC_068350.1' = 'chr13'; 'NC_068351.1' = 'chr14'; 'NC_068352.1' = 'chr15'
    'NC_030786.1' = 'Pltd'
}

function Get-Tag {
    param([string]$Header, [string]$Name)
    $m = [regex]::Match($Header, '\[' + [regex]::Escape($Name) + '=([^\]]*)\]')
    if ($m.Success) { return $m.Groups[1].Value }
    return ''
}

# ---------------------------------------------------------------- single pass over the FASTA
$hits     = New-Object 'System.Collections.Generic.List[object]'
$totalRec = 0
$reader   = New-Object System.IO.StreamReader($Fasta)
$cur      = $null
$sb       = $null

try {
    while ($null -ne ($line = $reader.ReadLine())) {
        if ($line.StartsWith('>')) {
            if ($null -ne $cur) { $cur.Seq = $sb.ToString(); $hits.Add($cur) }
            $totalRec++
            $cur = $null
            $sb  = $null

            $desc = Get-Tag $line 'protein'
            $rule = $null
            foreach ($r in $rules) {
                if ($desc -notmatch $r.Pattern) { continue }
                if ($r.Exclude -and $desc -match $r.Exclude) { continue }
                $rule = $r
                break
            }
            if ($null -eq $rule) { continue }

            $loc    = Get-Tag $line 'location'
            $strand = '+'
            if ($loc -match 'complement') { $strand = '-' }

            $min = [int64]::MaxValue
            $max = [int64]0
            foreach ($n in [regex]::Matches($loc, '\d+')) {
                $v = [int64]$n.Value
                if ($v -lt $min) { $min = $v }
                if ($v -gt $max) { $max = $v }
            }

            $firstTok = ($line.Substring(1) -split ' ')[0]
            $chrom    = (($firstTok -split '\|')[1]) -replace '_cds_.*$', ''

            $cur = [pscustomobject]@{
                Group      = $rule.Group
                Role       = $rule.Role
                Gene       = (Get-Tag $line 'gene')
                Protein    = $desc
                Chrom      = $chrom
                Start      = $min
                End        = $max
                Strand     = $strand
                Transcript = (Get-Tag $line 'protein_id')
                Inactive   = $(if ($desc -match 'inactive') { 'yes' } else { '' })
                LowQuality = $(if ($desc -match 'LOW QUALITY') { 'yes' } else { '' })
                Partial    = $(if ($loc -match '[<>]') { 'yes' } else { '' })
                Exception  = (Get-Tag $line 'exception')
                Seq        = ''
            }
            $sb = New-Object System.Text.StringBuilder
        }
        elseif ($null -ne $cur) {
            [void]$sb.Append($line.Trim())
        }
    }
    if ($null -ne $cur) { $cur.Seq = $sb.ToString(); $hits.Add($cur) }
}
finally { $reader.Close() }

# ---------------------------------------------------------------- collapse to gene level
$genes = foreach ($g in ($hits | Group-Object Gene)) {
    $rep = $g.Group | Sort-Object @{ Expression = { $_.Seq.Length }; Descending = $true } | Select-Object -First 1
    $chromName = $rep.Chrom
    if ($chromNames.ContainsKey($rep.Chrom)) { $chromName = $chromNames[$rep.Chrom] }

    [pscustomobject]@{
        Group       = $rep.Group
        Role        = $rep.Role
        Gene        = $rep.Gene
        Protein     = $rep.Protein
        ChromName   = $chromName
        Chrom       = $rep.Chrom
        Start       = $rep.Start
        End         = $rep.End
        Strand      = $rep.Strand
        CDS_len     = $rep.Seq.Length
        Transcripts = $g.Count
        Plastid     = $(if ($rep.Chrom -eq 'NC_030786.1') { 'yes' } else { '' })
        Inactive    = $rep.Inactive
        LowQuality  = $rep.LowQuality
        Partial     = $rep.Partial
        Exception   = $rep.Exception
        Transcript  = $rep.Transcript
        Seq         = $rep.Seq
    }
}

$groupOrder = @{ 'Inhibit' = 1; 'Enhance' = 2; 'Source' = 3 }
$sorted = $genes | Sort-Object `
    @{ Expression = { $groupOrder[$_.Group] } }, `
    @{ Expression = { $_.Role } }, `
    @{ Expression = { $_.Gene } }

# ---------------------------------------------------------------- output
$tsvPath = Join-Path $OutDir 'starch_genes.tsv'
$sorted |
    Select-Object Group, Role, Gene, Protein, ChromName, Chrom, Start, End, Strand,
                  CDS_len, Transcripts, Plastid, Inactive, LowQuality, Partial, Exception, Transcript |
    Export-Csv -LiteralPath $tsvPath -Delimiter "`t" -NoTypeInformation -Encoding UTF8

$fnaPath = Join-Path $OutDir 'starch_genes.cds.fna'
$out = New-Object 'System.Collections.Generic.List[string]'
foreach ($g in $sorted) {
    $out.Add('>' + $g.Gene + '|' + $g.Role + '|' + $g.ChromName + ':' + $g.Start + '-' + $g.End + '(' + $g.Strand + ')|' + $g.Transcript)
    for ($i = 0; $i -lt $g.Seq.Length; $i += 60) {
        $out.Add($g.Seq.Substring($i, [Math]::Min(60, $g.Seq.Length - $i)))
    }
}
[System.IO.File]::WriteAllLines($fnaPath, $out, (New-Object System.Text.UTF8Encoding($false)))

# ---------------------------------------------------------------- report
''
'== extract_starch_genes =='
'input            : ' + $Fasta
'CDS records read : ' + $totalRec
'hits (records)   : ' + $hits.Count
'genes (unique)   : ' + $sorted.Count
'tsv              : ' + $tsvPath
'fasta            : ' + $fnaPath
''
'-- per group / role --'
$sorted | Group-Object Group, Role | Sort-Object Name | ForEach-Object {
    '{0,4}  {1}' -f $_.Count, $_.Name
}
''
'-- suspicious models --'
$sorted | Where-Object { $_.Inactive -or $_.LowQuality -or $_.Partial -or $_.CDS_len -lt 600 } | ForEach-Object {
    '{0,-9} {1,-6} {2,-14} len={3,-5} {4}' -f $_.Group, $_.Role, $_.Gene, $_.CDS_len, $_.Protein
}
