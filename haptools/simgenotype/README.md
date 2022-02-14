# Haptools simgenotype

`haptools simgenotype` takes as input a reference set of haplotypes in VCF format and a user-specified admixture model. It outputs a VCF file with simulated genotype information for admixed genotypes, as well as a breakpoints file that can be used for visualization.

## Basic usage

```
haptools simgenotype \
  --invcf REFVCF \
  --sample_info SAMPLEINFOFILE \
  --model MODELFILE \
  --map GENETICMAP \
  --out OUTPREFIX
```

Detailed information about each option, and example commands using publicly available files, are shown below.

## Detailed usage

## File formats

## Examples