#!/usr/bin/env Rscript

#nathan dot lazar at gmail dot com

#R script to colate data on CpG methylation from cpg10 evidence files

# Usage:
# Rscript ./bin/full_R_analysis.R 
#   <file of chromosome lengths>
#   <sorted and tabulated cpg evidence file>
#   <breakpoint region file>
#   <gene_file>
#   <repmask file>
#   <cpg_island_file>

# Example: 
# Rscript ./bin/full_R_analysis.R 
#   ~/VOK_BS_GENOME/NomLeu1.0_lengths.txt
#   ~/VOK_BS_GENOME/ALIGN/ALL/cpg10gzip
#   ~/VOK_BS_GENOME/NomLeu1.0_bp_clean.txt
#   ~/VOK_BS_GENOME/Nomascus_leucogenys.Nleu1.0.70.fixed.gtf
#   ~/VOK_BS_GENOME/NomLeu1.0_repmask.txt
#   ~/GIBBON_METH/NOMLEU1/features/gibbon_cpgislands.gff

args <- commandArgs(TRUE)
len_file <- args[1]
drive <- args[2]
bp_file <- args[3]
gene_file <- args[4]
rep_file <- args[5]
cpg_isl_file <- args[6]

#bac_file <- 'BAC_amp_on_NomLeu1.0.txt'

source('~/gibbon_meth/R_meth_functions.R')

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, 'NomLeu1.0')

# Read in CpG data and make BSeq object with coverage over 4
############################################################
all.bs <- make_all_bs(drive, 'NLE_Vok', seqinfo, 4)

#Measure methylation of CpGs with coverage of 4 or more
########################################################
all.cpg.meth <- mean(mcgetMeth(all.bs, type='raw', what='perBase'),
                     na.rm=T)

############################
# Breakpoint region analysis
############################

# Read in Breakpoint region data and make GRanges object
bp.gr <- read_bp(bp_file, seqinfo)

# Make GRanges object of 10kb regions on each side of breaks
bp.lr.gr <- make_lr_bp(bp.gr, 10000)

#Write breakpoint regions to bed file
bp.bed <- data.frame(seqnames(bp.lr.gr), start(bp.lr.gr), end(bp.lr.gr))
write.table(bp.bed, file='bp_regions.bed', quote=F, sep='\t',
            row.names=F, col.names=F)

# Add mean methylation, number of CpGs and coverage
# of sides of bp regions to GRanges object
bp.lr.gr <- add_meth_cpg_cov(bp.lr.gr, all.bs)

# Run permutation analysis to determine whether the breakpoint
# regions have lower methylation or coverage than would be seen
# in random regions
bp.permute <- permute(bp.lr.gr, bp.lr.gr, all.bs, n=1000, type='all', 
                      min.chr.size=12000, end.exclude=1000)

###############
# Gene analysis
###############

gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)

# Add methylation, number of CpGs and mean CpG coverage to each
# of these GRanges objects
gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs)

# Run permutation analysis
gene.permute.list <- lapply(gene.gr.list, permute, bp.lr.gr, all.bs, n=1000, 
                            type='gene', min.chr.size=12000,
                            end.exclude=1000)

#################
# Repeat analysis
#################

# Read in repmask data
rep.gr <- make_rep_gr(rep_file, seqinfo(all.bs))

# Add mean methylation, number of CpGs and mean CpG
# coverage to repeat GRanges object
rep.gr <- add_meth_cpg_cov(rep.gr, all.bs)

# All repeat permutation analysis
rep.permute <- permute(rep.gr, bp.lr.gr, all.bs, n=1000,
                       type='repeat', min.chr.size=12000,
                       end.exclude=1000)

# Run permutation analysis on major repeat classes 
# LINE, SINE, DNA, LTR, Satellite

rep.gr.list <- list()
rep.gr.list$LINE <- rep.gr[grepl('LINE', rep.gr$rep_family)]
rep.gr.list$SINE <- rep.gr[grepl('SINE', rep.gr$rep_family)]
rep.gr.list$DNA <- rep.gr[grepl('DNA', rep.gr$rep_family)]
rep.gr.list$LTR <- rep.gr[grepl('LTR', rep.gr$rep_family)]
rep.gr.list$Satellite <- rep.gr[grepl('Satellite', rep.gr$rep_family)]

rep.permute.list <- lapply(rep.gr.list, permute, bp.lr.gr, all.bs, n=1000,
                           type='repeat', min.chr.size=12000,
                           end.exclude=1000)


#####################
# CpG island analysis
#####################

# Read in CpG island data
cpg_island.gr <- gff2GRanges(cpg_island_file, seqinfo)

cpg_shore.gr <- make_cpg_shore(cpg_island.gr, 1000)

# Add mean methylation, number of CpGs and mean CpG
# coverage to CpGisland GRanges object
cpg_island.gr <- add_meth_cpg_cov(cpg_island.gr, all.bs)
cpg_shore.gr <- add_meth_cpg_cov(cpg_shore.gr, all.bs)

# Look at differences between classes of breakpoints
####################################################

#bp_Class1.gr <- bp_region.gr[bp_region.gr$class=='Class_I']
#bp_Class2a.gr <- bp_region.gr[bp_region.gr$class=='Class_II-a']
#bp_Class2b.gr <- bp_region.gr[bp_region.gr$class=='Class_II-b']














