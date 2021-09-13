# `config.yml` - config input and options for the pipeline
In order to run the pipeline on your own data, you must specify inputs and options for the pipeline within the `config.yml` file.
## Required input
1. `ref_panel` (str) - The path to a SNP-STR haplotype reference panel for 1000G (downloaded from [here](http://gymreklab.com/2018/03/05/snpstr_imputation.html))

    If there are multiple, per-chromosome VCFs, just specify one path but replace the contig name in the file name with `{chr}`.

    The VCF(s) must be sorted and indexed (with a `.tbi` file in the same directory). In additition, the varaint IDs of each STR in the VCF must have 'STR_' prepended.
2. `samples` (str) - The path to a 1000G samples TSV file

    You can get this file by clicking "Download the list" [here](https://www.internationalgenome.org/data-portal/sample).
3. `loci` (list of dict) - A list of loci to simulate

	Each locus dictionary has the following key-value pairs for SNP loci, a chosen STR, and a chosen SNP (only if this is preferred over a random SNP):
	  - `locus` (str) - composed of a contig name, a colon, the start position, a dash, and the end position of the locus
	  - `str` (str) - composed of a contig name, a colon, and the start position of the STR
	  - `snp` (str - _optional_) - composed of a contig name, a colon, and the start position of the SNP
## Options
1. `mode` (str) - Whether to run the simulation framework for a SNP ("snp") or an STR ("str")

    You can also provide a list of strings intead of a single string, if you want the pipeline to do multiple modes at once.
2. `superpop` (str) - The 1000G superpopulation code that should be used

    Please provide a three letter code, corresponding to the symbol in the samples file. This will default to EUR (aka european ancestry), if not specified.
3. `min_maf` (float) - The pipeline will discard rare SNPs with an MAF below this number

	This will default to 0 if not specified, but we recommend a value of at least 0.1
4. `beta` (float) - The strength of association between the chosen STR or SNP and the simulated phenotype

	This will default to 0.1 if not specified. You can also provide a list of floats intead of a single float, if you want the pipeline to try a number of different values.
