#!/usr/bin/env Rscript

# This R script summarizes the output of several executions of FINEMAP and SuSiE
# It cannot be executed independently of Snakemake




# handle error messages properly:
# they should all get written to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type='message')

# import the finemap and susie results
# and the path to an output directory
finemap_results = readRDS(snakemake@input[["finemap"]])[[1]]
susie_results = readRDS(snakemake@input[["susie"]])
out = snakemake@params[["outdir"]]

dir.create(out, showWarnings = FALSE)

# parse the susie data:
# 1) the truly causal variant, as defined in the simulation
# 2) whether the causal variant was provide in the genotypes
# 3) the list provided by susie as output
causal_variant = susie_results$causal_var
exclude_causal = susie_results$causal_excluded
fitted = susie_results$fitted
p = length(fitted$pip) # the number of variants tested

# first, we must create a vector with the causal status
# this vector indicates which variant is truly causal
b = rep(0,p)
names(b) = names(fitted$pip)
if (!exclude_causal) {
	b[causal_variant] = 1
}

# parse the finemap data
snp = finemap_results$snp
pip = snp[order(as.numeric(snp$snp)),]$snp_prob

# plot the results of susie
pdf(paste0(out,'/susie.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(fitted, y='PIP', b=b, max_cs=0.4, main = paste('SuSiE, ', length(fitted$sets$cs), 'CS identified'))
dev.off()

bhat = coef(fitted)[-1]
pdf(paste0(out,'/susie_eff.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(bhat, y='bhat', b=b, main = 'SuSiE, effect size estimate') 
dev.off()

# plot the results of running FINEMAP
pdf(paste0(out,'/finemap.pdf'), width =5, height = 5, pointsize=16)
susieR::susie_plot(pip, y='PIP', b=b, main = 'Bayesian sparse regression')
dev.off()

close(log)
