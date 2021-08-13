# scripts
This directory contains various scripts used by the pipeline. However, you can use most of these scripts on their own, too. Some may even be helpful in day-to-day use.

All python scripts implement the --help argument. For bash and R scripts, you can run `head <script>` to read about their usage.

## [finemapping_methods.R](finemapping_methods.R)
An R script to execute FINEMAP and SuSie and generate Manhattan plots against the true causal variants from the simulation framework.

## [generate_phenotypes.py](generate_phenotypes.py)
Uses pre-existing gentoype matrices to create a file of phenotypes.

## [gt_matrix.bash](gt_matrix.bash)
Creates a genotype matrix (with entries 0, 1, or 2) from a VCF of SNPs. Each column is a sample and the rows are SNPs.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [str_gt.py](str_gt.py)
Compute the sum of the differences of the allele lengths from the reference allele length.

