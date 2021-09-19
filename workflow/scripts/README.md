# scripts
This directory contains various scripts used by the pipeline. However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the --help argument. For bash and R scripts, you can run `head <script>` to read about their usage.

## [finemapping_methods.R](finemapping_methods.R)
An R script to execute FINEMAP and SuSie and generate Manhattan plots against the true causal variants from the simulation framework.

## [generate_phenotypes.py](generate_phenotypes.py)
Uses pre-existing gentoype matrices to create a file of phenotypes.

## [matrix.bash](matrix.bash)
Extract fields (POS, ID, MAF, GT, and alleles) from a VCF of STRs and SNPs. Each column is a sample and the rows are variants.

## [gt_matrix.py](gt_matrix.py)
Transform the TSV from `matrix.bash` into a genotype matrix, with entries 0, 1, and 2. Columns are variants and rows are samples. Column names have the variant type appended.

## [ld_heatmap.py](ld_heatmap.py)
Creates a heatmap of the LD pattern among the variants in a genotype matrix.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [str_gt.py](str_gt.py)
Compute the sum of the differences of the allele lengths from the reference allele length.

