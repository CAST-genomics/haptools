#!/usr/bin/env Rscript

# This R script runs the fine-mapping methods SusieR and FINEMAP

# param1: The path to a TSV containing the genotype data.
# param2: The path to a TSV containing the phenotype data.
# param3: The path to a directory in which to write output
#         This will be created if it doesn't exist.


library(abind)


# first, we import the phenotype file into X and Y vectors
args = commandArgs(trailingOnly = TRUE)
gt = args[1]
phen = args[2]
out = args[3]

# gt = "out/1_98001984-99001984/gt_matrix.tsv.gz"
# phen = "out/1_98001984-99001984/phens.tsv.gz"
# out = "out/1_98001984-99001984/susieR"
dir.create(out, showWarnings = FALSE)


# import genotype matrices as proper matrices
gt = read.csv(gt, sep="\t", header=T)
phen = read.csv(phen, sep="\t", header=T)
# the number of samples and the number of variants:
n = nrow(gt)
p = ncol(gt)
# create matrices without unecessary columns
X = as.matrix(gt[,-1])
storage.mode(X) = 'double'
y = as.matrix(phen[,ncol(phen)])


# create vector with causal status
# this vector indicates which variant is truly causal
b = rep(0,p)
names(b) = colnames(X)
b[colnames(phen)[2]] = 1


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


output = paste0(out, "/finemap")
args = "--n-causal-snps 2"
commandArgs = function(...) 1
source(paste0(.libPaths(), '/susieR/code/finemap.R'))


finemap = readRDS(paste0(out,"/finemap.rds"))[[1]]
snp = finemap$snp
pip = snp[order(as.numeric(snp$snp)),]$snp_prob


pdf(paste0(out,'/finemap.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(pip, y='PIP', b=b, main = 'Bayesian sparse regression')
dev.off()


fitted = susieR::susie(X, y, L=5,
               estimate_residual_variance=TRUE, 
               scaled_prior_variance=0.2,
               tol=1e-3, track_fit=TRUE, min_abs_corr=0.02)

pdf(paste0(out,'/susie.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(fitted, y='PIP', b=b, max_cs=0.4, main = paste('SuSiE, ', length(fitted$sets$cs), 'CS identified'))
dev.off()


bhat = coef(fitted)[-1]
pdf(paste0(out,'/susie_eff.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(bhat, y='bhat', b=b, main = 'SuSiE, effect size estimate') 
dev.off()

saveRDS(fitted, paste0(out, '/susie.rds'))
