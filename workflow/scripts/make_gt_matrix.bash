#!/usr/bin/env bash

# arg1: the superpopulation to sample from
# arg2: the minimum minor allele frequency
# arg3: VCF file from which to obtain genotypes
# arg4: samples file from 1000G website (contains columns: Sample name|Sex|Biosample ID|Population code|Population name|Superpopulation code|Superpopulation name|Population elastic ID|Data collections)
# ex: workflow/scripts/make_genotype_matrix.bash EUR 0.1 /projects/ps-gymreklab/resources/datasets/1000Genomes/phase3VCFs/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz data/igsr_samples.tsv




samp="$1"
min_maf="$2"
vcf_1000g="$3"
samps_file="$4"

# zcat "$vcf_1000g" 2>/dev/null | \
#     grep -E -m1 '^#CHROM' | \
#     cut -f10- | \
#     tr '\t' '\n' | \
#     sort | \
#     join -j1 -t $'\t' - <(
#         sort -k1,1 "$samps_file"
#     ) > "$samps_file"

comm -12 <(
    bcftools query -l "$vcf_1000g" | \
    sort
) <(
    awk -F $'\t' -v 'OFS=\t' '$6 == "'"$samp"'" {print $1;}' "$samps_file" | \
    sort
) | \
paste -s -d $'\t' | \
tee >(
    bcftools view --min-af "$min_maf":minor --samples-file <(tr $'\t' '\n') "$vcf_1000g" | \
    bcftools query -f '[%GT\t]\n' | \
    sed 's/\t$//; s/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g'
)
