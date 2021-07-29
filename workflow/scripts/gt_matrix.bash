#!/usr/bin/env bash

# arg1: the superpopulation to sample from
# arg2: the minimum minor allele frequency
# arg3: the locus from which to obtain the genotypes
# arg4: VCF file from which to obtain genotypes
# arg5: samples file from 1000G website (contains columns: Sample name|Sex|Biosample ID|Population code|Population name|Superpopulation code|Superpopulation name|Population elastic ID|Data collections)
# ex: workflow/scripts/make_genotype_matrix.bash EUR 0.1 /projects/ps-gymreklab/mousavi/results/1000genomes/outs/merged/all_merged.sorted.vcf.gz data/igsr_samples.tsv


samp="$1"
min_maf="$2"
loc="$3"
vcf_1000g="$4"
samps_file="$5"


comm -12 <(
    bcftools query -l "$vcf_1000g" | \
    sort
) <(
    awk -F $'\t' -v 'OFS=\t' '$6 == "'"$samp"'" {print $1;}' "$samps_file" | \
    sort
) | {
    echo -ne "POS,ID\talleles\t"
    paste -s -d $'\t'
} | \
tee >(
    bcftools view --min-af "$min_maf":minor --samples-file <(tr $'\t' '\n' | tail -n+2) "$vcf_1000g" "$loc" | \
    bcftools query -f '%POS,%ID\t%REF,%ALT\t[%GT\t]\n' | \
    sed 's/\t$//'
)
