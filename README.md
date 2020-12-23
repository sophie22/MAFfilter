# MAFfilter
script to filter vcf for variants with MAF&lt;1% and save transcripts used in a tsv

## Step 1:
filter down to only homozygous variants, using bcftools built-in filtering functionality

## Step 2: pycellbase
use CellBase to annotate variants and only keep those with refAlleleFreq < 1(in GNOMAD_EXOMES, ALL)
1. convert variant description into Chrom:Pos:Ref:Alt format
2. query CellBase API, parse response
For variants with MAF < 1%:

## Step 3: pyvcf
3. grab transcript options from X111923_annot_MAF.vcf, check if any are in the GEMINI genes2transcripts file:
* if yes, then save that
* if no, then save all options

## Step 4: pandas
4. save overall output as a csv with columns:
chromosome, position, ref, alt, MAF, HGVS, HGNC, RefSeq, consequence, other
