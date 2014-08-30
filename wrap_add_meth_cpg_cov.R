#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to be run by HTCondor in parallel 
# that wraps add_meth_cpg_cov.R

# Usage:
# Rscript ./wrap_add_meth_cpg_cov.R
#   <feat_gr_and_all_bs.dat>  
#   <out_dir>
#   <number of cores>

.libPaths('/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.1/')

library(foreach)
library(doMC)
library(bsseq)
source('/mnt/lustre1/users/lazar/GIBBONS/gibbon_meth/R_meth_functions.R')

args <- commandArgs(TRUE)
#Load  feat.gr, all.bs
load(args[1])
cores <- as.numeric(args[2])
registerDoMC(cores)
outdir <- args[3]

feat.gr <- add_meth_cpg_cov(feat.gr, all.bs)

# Save data to file
save(feat.gr, file=paste0(outdir, 'feat.dat'))

print(1)