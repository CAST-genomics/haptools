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

Model Format

Structure of model.dat file
```
{num_samples} Admixed Pop_label1 Pop_label2 ... Pop_labeln
{num_generations} {admixed_freq} {pop_label1_freq} {pop_label2_freq} ... {pop_labeln_freq}
```

Example model.dat file

```
40   Admixed   CEU    YRI
6    0         0.2    0.8
```

Map Format

Outfile Format


## Examples
