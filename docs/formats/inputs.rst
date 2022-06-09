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
