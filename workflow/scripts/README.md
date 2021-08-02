# scripts

## [generate_phenotypes.py](generate_phenotypes.py)
Uses pre-existing gentoype matrices to create a file of phenotypes.

## [plot_phenotypes.py](plot_phenotypes.py)
Performs a linear regression and outputs a plot for phenotypes from `generate_phenotypes.py`.

## [gt_matrix.bash](gt_matrix.bash)
Creates a genotype matrix (with entries 0, 1, or 2) from a VCF of SNPs. Each column is a sample and the rows are SNPs.

## [str_gt.py](str_gt.py)
Compute the sum of the differences of the allele lengths from the reference allele length.

