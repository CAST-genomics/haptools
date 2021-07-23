# `config.yml` - config input and options for the pipeline
In order to run the pipeline on your own data, you must specify inputs and options for the pipeline within the `config.yml` file.
## Required input
1. `snp_vcf` (str) - The path to a VCF for SNPs in 1000G
    If there are multiple, per-chromosome VCFs, just specify one path but replace the contig name in the file name with `{chr}`.
2. `samples` (str) - The path to a 1000G samples TSV file
    You can get this file by clicking "Download the list" [here](https://www.internationalgenome.org/data-portal/sample).
3. `loci` (list of str) - A list of loci to simulate
	Each locus can be specified as a str composed of a contig name, a colon, the start position, a dash, and the end position.
## Options
1. `superpop` (str) - The 1000G superpopulation code that should be used
    Please provide a three letter code, corresponding to the symbol in the samples file. This will default to EUR (aka european ancestry), if not specified.
2. `min_maf` (float) - The pipeline will discard rare SNPs with an MAF below this number
	This will default to 0 if not specified, but we recommend a value of at least 0.1
4. `out` (str) - The path to a directory in which to place all of the output of the pipeline
	This will default to "out" in the current working directory if not specified.
