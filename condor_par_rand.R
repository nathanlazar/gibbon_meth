#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# R script to be run by HTCondor in parallel 
# that wraps par_rand.R

library(foreach)
library(doMC)

# Usage:
# Rscript ./condor_par_rand.R
#   <file_of_R_data.dat>
#   <type = all, gene, etc.>
#   <size of ends of chroms to be excluded>
#   <number of reps to be run on each node>
#   <number of cores>

# Example:
# Rscript ./par_rand.R
#   par_permute.dat
#   all
#   1000
#   20
#   24

args <- commandArgs(TRUE)
#Load  feat.gr, bp.lr.gr, all.bs, breaks, sizes, lengths
load(args[1])

type <- args[2]
n <- as.numeric(args[3])
end.exclude <- as.numeric(args[4])
registerDoMC(as.numeric(args[5]))

sizes <- bp.lr.gr$size

rand <- par_rand(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, reps)

# Save data to file (Condor will append numbers)
save(rand, file='rand.dat')

