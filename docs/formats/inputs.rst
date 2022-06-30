.. _formats-inputs:


Inputs
=========

Genotype file format
--------------------
Genotype files must be specified as VCF or BCF files.

Phenotype file format
---------------------
Phenotypes are expected to be in a tab-separated format with two columns:

1. sample ID, which must match those from the genotype file
2. the phenotype value

There should be no header line in the file.

See ``tests/data/simple.tsv`` for an example of a phenotype file.

Covariate file format
---------------------
Covariates are expected to be in a tab-separated format with the first column
corresponding to the sample ID. Subsequent columns should contain each of your
covariates.

The first line of the file will be treated as a header. The column with sample IDs
should be named "sample" and subsequent columns can be labeled however you'd like.

See ``tests/data/covars.tsv`` for an example of a covariate file.

1000 Genomes sample_info file format
------------------------------------
Within the subcommand ``haptools simgenotype`` we use a file to map samples in the 
reference to their population listed in the model file. This code is expected to work
out of the box with 1000 genomes data and we have pre-computed this mapping file as 
well as given the command to how to compute it if you desire another as well.

``cut -f 1,4 igsr-1000\ genomes\ on\ grch38.tsv | sed '1d' | sed -e 's/ /\t/g' > 1000genomes_sampleinfo.tsv``

See ``example-files/1000genomes_sampleinfo.tsv`` for an example of the 1000genomes 
GRCh38 samples mapped to their subpopulations.
