#!/usr/bin/env Rscript

# This R script runs the fine-mapping methods SusieR and FINEMAP

# param1: The path to a TSV containing the genotype data.
# param2: The path to a TSV containing the phenotype data.
# param3: The path to a directory in which to write output
#         This will be created if it doesn't exist.
# param4: 1 if the causal variant should be removed from the genotype matrix and
#         0 otherwise


library(abind)


# first, we import the phenotype file into X and Y vectors
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
out = args[3]
exclude_causal = as.logical(as.integer(args[4]))

# gt = "out/1_98001984-99001984/gt_matrix.tsv.gz"
# phen = "out/1_98001984-99001984/phens.tsv.gz"
# out = "out/1_98001984-99001984/susieR"
dir.create(out, showWarnings = FALSE)


# import genotype matrices as proper matrices
gt = read.csv(gt, sep="\t", header=T)
phen = read.csv(phen, sep="\t", header=T)
# create matrices without unecessary columns
X = as.matrix(gt[,-1])
# the number of samples and the number of variants:
n = nrow(X)
p = ncol(X)
storage.mode(X) = 'double'
y = as.matrix(phen[,ncol(phen)])
# what is the column name of the causal variant?
causal_variant = colnames(phen)[2]


# remove the causal variant if requested
if (exclude_causal) {
  X = X[,!(colnames(X) %in% c(causal_variant))]
}


# compute summary statistics for FINEMAP
mm_regression = function(X, Y, Z=NULL) {
  if (!is.null(Z)) {
      Z = as.matrix(Z)
  }
  reg = lapply(seq_len(ncol(Y)), function (i) simplify2array(susieR:::univariate_regression(X, Y[,i], Z)))
  reg = do.call(abind, c(reg, list(along=0)))
  # return array:
  #   out[1,,] is beta hat (the least-squares estimates of the coefficients)
  #   out[2,,] is se betahat (the standard errors of the beta hats)
  return(aperm(reg, c(3,2,1)))
}
sumstats = mm_regression(as.matrix(X), as.matrix(y))
dat = list(X=X,Y=as.matrix(y))
input = paste0(out,'/sumstats.rds')
saveRDS(list(data=dat, sumstats=sumstats), input)


# run FINEMAP
# and set an output path; the results will be written to an RDS file with this basename
output = paste0(out, "/finemap")
args = "--n-causal-snps 1"
commandArgs = function(...) 1
source(paste0(.libPaths(), '/susieR/code/finemap.R'))


# run SuSiE
# write the output to an RDS file
fitted = susieR::susie(X, y, L=1)

# when writing the output, also include information about which variant is causal
# and whether it was included in the simulation
saveRDS(
  list(causal_var=causal_variant, causal_excluded=exclude_causal, fitted=fitted),
  paste0(out, '/susie.rds')
)
