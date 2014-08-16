# Functions used in analysis of bisulfite reads using BSmooth
# nathan dot lazar at gmail dot com

library(bsseq)

dir <- '~/gibbon_meth/'

source(paste0(dir, 'make_seqinfo.R'))
source(paste0(dir, 'make_bs_all.R'))
source(paste0(dir, 'read_bp.R'))
source(paste0(dir, 'read_cpg_island.R'))
source(paste0(dir, 'make_cpg_shore.R'))
source(paste0(dir, 'make_lr_bp.R'))
source(paste0(dir, 'make_rep_gr.R'))
source(paste0(dir, 'gtf2GRanges.R'))
source(paste0(dir, 'add_meth_cpg_cov.R'))
source(paste0(dir, 'get_region.R'))
source(paste0(dir, 'permute.R'))