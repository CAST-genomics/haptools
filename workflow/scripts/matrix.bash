#!/usr/bin/env bash

# arg1: the superpopulation to sample from
# arg2: the locus from which to obtain the genotypes
# arg3: VCF file from which to obtain genotypes
# arg4: samples file from 1000G website (contains columns: Sample name|Sex|Biosample ID|Population code|Population name|Superpopulation code|Superpopulation name|Population elastic ID|Data collections)
# ex: workflow/scripts/gt_matrix.bash EUR '1:98001984-99001984' /projects/ps-gymreklab/resources/datasets/snpstr/1kg.snp.str.chr1.vcf.gz data/igsr_samples.tsv



samp="$1"
loc="$2"
vcf_1000g="$3"
samps_file="$4"


comm -12 <(
    bcftools query -l "$vcf_1000g" | \
    sort
) <(
    awk -F $'\t' -v 'OFS=\t' '$6 == "'"$samp"'" {print $1;}' "$samps_file" | \
    sort
) | {
    echo -ne "POS\tID\tMAF\talleles\t"
    paste -s -d $'\t'
} | \
tee >(
    bcftools view --samples-file <(tr $'\t' '\n' | tail -n+5) "$vcf_1000g" "$loc" | \
    bcftools +fill-tags - -- -t MAF | \
    bcftools query -f '%POS\t%ID\t%INFO/MAF\t%REF,%ALT\t[%GT\t]\n' | \
    sed 's/\t$//'
)
