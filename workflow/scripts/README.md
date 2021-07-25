# scripts

## [generate_phenotypes.py](generate_phenotypes.py)
Uses pre-existing gentoype matrices to create a file of phenotypes.

## [snp_gt_matrix.bash](snp_gt_matrix.bash)
Creates a genotype matrix (with entries 0, 1, or 2) from a VCF of SNPs. Each column is a sample and the rows are SNPs.

## [str_gt_matrix.bash](str_gt_matrix.bash)
Creates a genotype matrix from a VCF of STRs. Each column is a sample.

## [str_gt.py](str_gt.py)
Compute the sum of the differences of the allele lengths from the reference allele length.

