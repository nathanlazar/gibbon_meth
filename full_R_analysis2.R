#!/usr/bin/env Rscript

#nathan dot lazar at gmail dot com

#R script to colate data on CpG methylation from cpg10 evidence files

# Usage:
# Rscript ./bin/full_R_analysis.R
#   <bin directory>
#   <output directory> 
#   <file of chromosome lengths>
#   <sorted and tabulated cpg evidence file>
#   <breakpoint region file>
#   <gene_file>
#   <repmask file>
#   <cpg_island_file>

# Example: 
# Rscript ./bin/full_R_analysis.R 
#   gibbon_meth/
#   VOK_GENOME/all_BP_full_analysis/
#   VOK_GENOME/NomLeu1.0_lengths.txt
#   VOK_GENOME/ALIGN/ALL/cpg10gzip
#   VOK_GENOME/NomLeu1.0_bp_clean.txt
#   VOK_GENOME/NomLeu1.0.70.genes.gtf
#   VOK_GENOME/NomLeu1.0_repmask.txt
#   VOK_GENOME/NomLeu1.0_cpg_islands.gff

library(bsseq)

args <- commandArgs(TRUE)
bindir <- args[1]
outdir <- args[2]
len_file <- args[3]
cpg_drive <- args[4]
bp_file <- args[5]
gene_file <- args[6]
rep_file <- args[7]
cpg_isl_file <- args[8]

#bac_file <- 'BAC_amp_on_NomLeu1.0.txt'

source(paste0(bindir, 'R_meth_functions.R'))

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file, 'NomLeu1.0')

# Read in CpG data and make BSeq object with coverage over 4
############################################################
all.bs <- make_all_bs(cpg_drive, 'NLE_Vok', seqinfo, 4)

# Write bed files of CpG coverage and methylation
#################################################
make_tracks(all.bs, outdir, 'meth.bedgraph', 'cov.bedgraph')

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
# unless bp_file has 'BAC' in it (hack)
if(grepl('BAC', bp_file)) {
  bp.lr.gr <- bp.gr
} else {
  bp.lr.gr <- make_lr_bp(bp.gr, 10000)
}

# Write breakpoint regions to bed file
bp.bed <- data.frame(seqnames(bp.lr.gr), start(bp.lr.gr), end(bp.lr.gr))
write.table(bp.bed, file=paste0(outdir,'bp_regions.bed'), quote=F, 
            sep='\t', row.names=F, col.names=F)

# Add mean methylation, number of CpGs and coverage
# of sides of bp regions to GRanges object
bp.lr.gr <- add_meth_cpg_cov(bp.lr.gr, all.bs)

# Write out coverage of targeted regions
write.table(data.frame(chr=as.vector(seqnames(bp.lr.gr)), 
                       start=start(bp.lr.gr), end=end(bp.lr.gr), 
                       cov=bp.lr.gr$cov), 
            'amp.cov.txt', quote=F, sep='\t', row.names=F,
            col.names=T)

# Run permutation analysis to determine whether the breakpoint
# regions have lower methylation or coverage than would be seen
# in random regions
min.chr.size <- max(bp.lr.gr$size)+2000 # Choose only from chroms that
                                        # are big enough to 
                                        # contain regions
bp.permute <- par_permute(outdir, bindir, bp.lr.gr, bp.lr.gr, all.bs, 
                          n=1000, type='all', 
                          min.chr.size=min.chr.size, end.exclude=1000)

###############
# Gene analysis
###############

gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)

# Add methylation, number of CpGs and mean CpG coverage to each
# of these GRanges objects
gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs)

# Run permutation analyses
gene.permute <- par_permute(outdir, bindir, gene.gr.list$gene, 
                            bp.lr.gr, all.bs, n=1000, type='gene', 
                            min.chr.size=min.chr.size, end.exclude=1000)
exon.permute <- par_permute(outdir, bindir, gene.gr.list$exon, 
                            bp.lr.gr, all.bs, n=1000, type='exon', 
                            min.chr.size=min.chr.size, end.exclude=1000)
intron.permute <- par_permute(outdir, bindir, gene.gr.list$intron, 
                              bp.lr.gr, all.bs, n=1000, type='intron',
                              min.chr.size=min.chr.size, end.exclude=1000)
promoter.permute <- par_permute(outdir, bindir, gene.gr.list$promoter, 
                                bp.lr.gr, all.bs, n=1000, type='promoter',
                                min.chr.size=min.chr.size, end.exclude=1000)
gene.permute.list <- list(gene=gene.permute, exon=exon.permute, 
                          intron=intron.permute, promoter=promoter.permute)

#################
# Repeat analysis
#################

# Read in repmask data
rep.gr <- make_rep_gr(rep_file, seqinfo(all.bs))

# Add mean methylation, number of CpGs and mean CpG
# coverage to repeat GRanges object
# rep.gr <- add_meth_cpg_cov(rep.gr, all.bs)

rep.gr <- condor_add_meth_cpg_cov(rep.gr, all.bs, 
            paste0(outdir, 'rep_add'))

rep.permute <- par_permute(outdir, bindir, rep.gr, bp.lr.gr, 
                           all.bs, n=1000,
                           type='repeat', min.chr.size=min.chr.size,
                           end.exclude=1000)

# Run permutation analysis on major repeat classes 
# LINE, SINE, DNA, LTR, Satellite

rep.gr.list <- list()
rep.gr.list$LINE <- rep.gr[grepl('LINE', rep.gr$rep_family)]
rep.gr.list$SINE <- rep.gr[grepl('SINE', rep.gr$rep_family)]
rep.gr.list$DNA <- rep.gr[grepl('DNA', rep.gr$rep_family)]
rep.gr.list$LTR <- rep.gr[grepl('LTR', rep.gr$rep_family)]
rep.gr.list$Satellite <- rep.gr[grepl('Satellite', rep.gr$rep_family)]

LINE.permute <- par_permute(outdir, bindir, rep.gr.list$LINE, bp.lr.gr,
                            all.bs, n=1000, type='LINE',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
SINE.permute <- par_permute(outdir, bindir, rep.gr.list$SINE, bp.lr.gr,
                            all.bs, n=1000, type='SINE',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
DNA.permute <- par_permute(outdir, bindir, rep.gr.list$DNA, bp.lr.gr,
                            all.bs, n=1000, type='DNA',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
LTR.permute <- par_permute(outdir, bindir, rep.gr.list$LTR, bp.lr.gr,
                            all.bs, n=1000, type='LTR',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
Satellite.permute <- par_permute(outdir, bindir, rep.gr.list$Satellite, 
                            bp.lr.gr, all.bs, n=1000, type='Satellite',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
rep.permute.list <- list(LINE=LINE.permute, SINE=SINE.permute,
                         DNA=DNA.permute, LTR=LTR.permute,
                         Satellite=Satellite.permute)

# Run permutation analysis on major classes of SINEs
# Alu, MIR, AluS, AluJ, AluY

SINE.gr.list <- list()
SINE.gr.list$Alu <- rep.gr[grepl('Alu', rep.gr$rep_class)]
SINE.gr.list$MIR <- rep.gr[grepl('MIR', rep.gr$rep_class)]
SINE.gr.list$AluS <- rep.gr[grepl('AluS', rep.gr$rep_class)]
SINE.gr.list$AluJ <- rep.gr[grepl('AluJ', rep.gr$rep_class)]
SINE.gr.list$AluY <- rep.gr[grepl('AluY', rep.gr$rep_class)] 

Alu.permute <-  par_permute(outdir, bindir, SINE.gr.list$Alu, bp.lr.gr,
                            all.bs, n=1000, type='Alu',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
MIR.permute <-  par_permute(outdir, bindir, SINE.gr.list$MIR, bp.lr.gr,
                            all.bs, n=1000, type='MIR',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluJ.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluJ, bp.lr.gr,
                            all.bs, n=1000, type='AluJ',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluS.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluS, bp.lr.gr,
                            all.bs, n=1000, type='AluS',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
AluY.permute <-  par_permute(outdir, bindir, SINE.gr.list$AluY, bp.lr.gr,
                            all.bs, n=1000, type='AluY',
                            min.chr.size=min.chr.size,
                            end.exclude=1000)
SINE.permute.list <- list(Alu=Alu.permute, MIR=MIR.permute,
                          AluJ=AluJ.permute, AluS=AluS.permute,
                          AluY=AluY.permute)

#####################
# CpG island analysis
#####################

# Read in CpG island data
cpg_island.gr <- gff2GRanges(cpg_isl_file, seqinfo)

cpg_shore.gr <- make_cpg_shore(cpg_island.gr, 1000)

# Add mean methylation, number of CpGs and mean CpG
# coverage to CpGisland GRanges object

cpg_island.gr <- condor_add_meth_cpg_cov(cpg_island.gr, all.bs,
            paste0(outdir, 'cpg_island_add'))
cpg_shore.gr <- condor_add_meth_cpg_cov(cpg_shore.gr, all.bs,
            paste0(outdir, 'cpg_shore_add'))

#cpg_island.gr <- add_meth_cpg_cov(cpg_island.gr, all.bs)
#cpg_shore.gr <- add_meth_cpg_cov(cpg_shore.gr, all.bs)

# Run permutation analysis
cpg_island.permute <- par_permute(outdir, bindir, cpg_island.gr, bp.lr.gr,
                                  all.bs, n=1000, type='CpGisl',
                                  min.chr.size=min.chr.size, end.exclude=1000)
cpg_shore.permute <- par_permute(outdir, bindir, cpg_shore.gr, bp.lr.gr,
                                 all.bs, n=1000, type='CpGshore',
                                 min.chr.size=min.chr.size, end.exclude=1000)

cpg.permute.list <- list(cpg_island=cpg_island.permute, 
                         cpg_shore=cpg_shore.permute)

#################################################
# Combine permutation results into one data frame
#################################################

per.results <- combine_per(bp.permute, gene.permute.list, rep.permute, 
  rep.permute.list, SINE.permute.list, cpg.permute.list)

# Write out data
#********************************************************************
save.image(file=paste0(outdir, 'R.dat'))
#load(paste0(outdir, 'R.dat'))
#********************************************************************

# Look at differences between classes of breakpoints
####################################################

#bp_Class1.gr <- bp_region.gr[bp_region.gr$class=='Class_I']
#bp_Class2a.gr <- bp_region.gr[bp_region.gr$class=='Class_II-a']
#bp_Class2b.gr <- bp_region.gr[bp_region.gr$class=='Class_II-b']
