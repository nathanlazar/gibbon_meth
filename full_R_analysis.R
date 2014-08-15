#!/usr/bin/env Rscript

#nathan dot lazar at gmail dot com

#R script to colate data on CpG methylation from cpg10 evidence files

#Usage Rscript ./bin/full_R_analysis.R 
#  <file of chromosome lengths>
#  <sorted and tabulated cpg evidence file>
#  <breakpoint region file>
#  <gene_file>
#  <repmask file>
#  <cpg_island_file>

# Example: 
# Rscript ./bin/full_R_analysis.R 
#   ~/VOK_BS_GENOME/NomLeu1.0_lengths.txt
#   <> 
#   <> 
#   ~/VOK_BS_GENOME/Nomascus_leucogenys.Nleu1.0.70.fixed.gtf
#   ~/GIBBON_METH/NOMLEU1/features/gibbon_cpgislands.gff

args <- commandArgs(TRUE)
len_file <- args[1]
drive <- args[2]
bp_file <- args[3]
gene_file <- args[4]
rep_file <- args[5]

bac_file <- 'BAC_amp_on_NomLeu1.0.txt'

library(bsseq)
source('~/gibbon_meth/R_meth_functions.R')

# Read in lengths of chromsomes and make seqinfo object
#######################################################
seqinfo <- make_seqinfo(len_file)

# Read in CpG data and make BSeq object with coverage over 4
############################################################
make_bs_all(drive, 'NLE_Vok', seqinfo, 4)

#Measure methylation of CpGs with coverage of 4 or more
########################################################
all.cpg.meth <- mean(getMeth(bs.all, type='raw', what='perBase'),
                     na.rm=T)

############################
# Breakpoint region analysis
############################

# Read in Breakpoint region data and make GRanges object
bp.gr <- read_bp(bp_file)

# Make GRanges object of 10kb regions on each side of breaks
bp.lr.gr <- make_lr_bp(bp.gr, 10000)

#Write breakpoint regions to bed file
bp.bed <- data.frame(seqnames(bp.lr.gr), start(bp.lr.gr), end(bp.lr.gr))
write.table(bp.bed, file='bp_regions.bed', quote=F, sep='\t',
            row.names=F, col.names=F)

# Add mean methylation, number of CpGs and coverage
# of sides of bp regions to GRanges object
bp.lr.gr <- add_meth_cpg_cov(bp.lr.gr, bs.all)

# Run permutation analysis to determine whether the breakpoint
# regions have lower methylation or coverage than would be seen
# in random regions
bp.permute <- permute(bp.lr.gr, bp.lr.gr, bs.all, n=1000, type='all', 
                      min.chr.size=12000, end.exclude=1000)

###############
# Gene analysis
###############

gene.gr.list <- gtf2GRanges(gene_file, seqinfo, prom_size=1000)

# Add methylation, number of CpGs and mean CpG coverage to each
# of these GRanges objects
gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, bs.all)

# Run permutation analysis


#################
# Repeat analysis
#################

# Read in repmask data
rep.gr <- make_rep_gr(rep_file, seqinfo(bs.all))

# Add mean methylation, number of CpGs and mean CpG
# coverage to repeat GRanges object
reps.gr <- add_meth_cpg_cov(reps.gr, bs.all)

#####################
# CpG island analysis
#####################

# Read in CpG island data
cpg_island.gr <- read_cpg_island(cpg_island_file, seqinfo)

cpg_shore.gr <- make_cpg_shore(cpg_island.gr, shore_size)


# Add mean methylation, number of CpGs and mean CpG
# coverage to CpGisland GRanges object



#Try some graphical model inference
###################################
library(bnlearn)

bp_meta <- mcols(bp.lr.gr)[,c('cpgs', 'mean_meth',
                              'sd_meth', 'median_meth',
                              'mean_cov', 'sd_cov',
                              'median_cov')]
#mark small regions
bp_meta$small <- 0
bp_meta$small[mcols(bp.lr.gr)$size < 8000] <- 1
bp_meta$classI <- 0
bp_meta$classI[mcols(bp.lr.gr)$class=='Class_I'] <- 1
bp_meta$classII <- 0
bp_meta$classII[grepl('Class_II', mcols(bp.lr.gr)$class)] <- 1
bp_meta$classIIa_b <- 0
bp_meta$classIIa_b[mcols(bp.lr.gr)$class=='Class_II-b'] <- 1

#Remove two sides of BP with 0 or 1 CpG
bp_meta <- bp_meta[!is.na(bp_meta$sd_cov),]
bp_meta <- as.data.frame(as.matrix(bp_meta))

bnet1 <- gs(data.frame(bp_meta))
bnet2 <- iamb(bp_meta)
bnet3 <- hc(bp_meta)

png('Bayesian_graphs.png', width=1440)
par(mfrow=c(1,3))
plot(bnet1)
plot(bnet2)
plot(bnet3)
dev.off()

m <- mean(mcols(bp.lr.gr)$mean_meth, na.rm=T)
wm <- weighted.mean(mcols(bp.lr.gr)$mean_meth,
                    mcols(bp.lr.gr)$cpgs, na.rm=T)

for(type in c('Class_I', 'Class_II-a', 'Class_II-b')) {
  m <- mean(mcols(bp.lr.gr)$mean_meth[mcols(bp.lr.gr)$class==type], na.rm=T)
  wm <- weighted.mean(mcols(bp.lr.gr)$mean_meth[mcols(bp.lr.gr)$class==type],
                      mcols(bp.lr.gr)$cpgs[mcols(bp.lr.gr)$class==type], na.rm=T)
  print(type)
  print(m)
  print(wm)
}

#######################################################
# Region permutation type 1
# Compare methylation of CpGs near breaks to the same 
# number of random CpGs in the genome (permutation analysis)

rand.mean.meths <- rep(0,10000)
rand.sd.meths <- rep(0,10000)
rand.mean.covs <- rep(0,10000)
rand.sd.covs <- rep(0,10000)
for(i in 1:10000) {
  rands <- sample(1:20832153, 15032)
  rand.meth <- getMeth(bs.all, type='raw')[rands]
  rand.cov <- getCoverage(bs.all)[rands]
  rand.mean.meths[i] <- mean(rand.meth)
  rand.sd.meths[i] <- sd(rand.meth)
  rand.mean.covs[i] <- mean(rand.cov)
}

#See how many of these have a value as extreme as CpGs near BP
# regions (get a p-value)
bp.w.av.meth <- weighted.mean(mcols(bp.lr.gr)$mean_meth,
                    mcols(bp.lr.gr)$cpgs)
bp.w.av.cov <- weighted.mean(mcols(bp.lr.gr)$mean_cov,
                    mcols(bp.lr.gr)$cpgs)
sum(rand.mean.meths < bp.w.av.meth)/10000
sum(rand.mean.covs > bp.w.av.cov)/10000

##########################################################
# Region permutation type 2
# find random regions matching the sizes of the BP regions
# and measure the methylation of all CpGs in these regions

#Make a list for sampling according to chrom lengths
# only sample from chroms at least 12kb
lengths <- seqlengths(bs.all)[seqlengths(bs.all) >= 12000]
breaks <- cumsum(as.numeric(lengths)) / 
            sum(as.numeric(lengths))
names(breaks) <- names(lengths)

sizes <- bp.lr.gr$size

rand_regions <- data.frame(chr=character(),
                           start=numeric(0),
                           end=numeric(0), 
                           stringsAsFactors=F)

for(i in 1:1000) {
  regions <- lapply(sizes, get_region, breaks, lengths)
  regions <- data.frame(do.call(rbind, regions))
  names(regions) <- c('chr', 'start', 'end')
  rand_regions <- rbind(rand_regions, regions)
}

rand.gr <- GRanges(seqnames=rand_regions[,1],
                     ranges=IRanges(start=as.numeric(matrix(rand_regions[,2])),
                                    end=as.numeric(matrix(rand_regions[,3]))))
meth <- getMeth(bs.all, regions=rand.gr, type='raw', what='perBase')
cov <- getCoverage(bs.all, regions=rand.gr, what='perBase')

#Get means in grouped by the number of BP regions
rand_mean_meth <- rep(0,1000)
rand_mean_cov <- rep(0,1000)
n <- length(sizes)
for(i in 1:1000) {
  rand_mean_meth[i] <- mean(unlist(meth[((i-1)*n+1):(i*n)]))
  rand_mean_cov[i] <- mean(unlist(cov[((i-1)*n+1):(i*n)]))
}

#See how many of these permutations have methylation as low as
# the breakpoint regions
sum(rand_mean_meth < bp.w.av.meth)/1000
sum(rand_mean_cov > bp.w.av.cov)/1000

######################################

#CLASS I
bp_class1.gr <- bp.lr.gr[bp.lr.gr$class=="Class_I"]
bp_class1.gr$meth <- getMeth(bs.all, regions=bp_class1.gr,
		             type='raw', what='perRegion')
#CLASS I region meth
bp_class1.gr <- bp.lr.gr[bp.lr.gr$class=="Class_I"]
bp_class1.gr$meth <- getMeth(bs.all, regions=bp_class1.gr,
		             type='raw', what='perRegion')

######################################
# Feature permutation
# Select the same number of features randomly throughout the 
# genome and calculate the average methylation, CpG count and 
# coverage for these. Compare values to features within 10kb
# of breakpoints

#-----------------
permute(bp.lr.gr, bp.lr.gr, all.bs, n=1000, type='all')
#-----------------

#SINE
#####
sines <- reps.gr[grepl('SINE', reps.gr$family)]
sines$meth <- getMeth(bs.all, regions=sines,
	              type='raw', what='perRegion')
sine_cpgs <- getCoverage(bs.all, regions=sines,
	                 what='perBase')
sines$cpgs <- unlist(lapply(sine_cpgs, length))
sines$cov <- unlist(lapply(sine_cpgs, mean, na.rm=T))

sines.bp.gr <- subsetByOverlaps(sines, bp.lr.gr)

sines_in_classI <- subsetByOverlaps(sines, bp_class1.gr)
mean(sines_in_classI$meth, na.rm=T)

#Alu analysis
#############

#Mesure methylation in all Alu elements
alu.gr <- reps.gr[reps.gr$family=='SINE/Alu']
alu.gr$meth <- getMeth(bs.all, regions=alu.gr,
	               type='raw', what='perRegion')
alu.gr$cpgs <- sapply(getMeth(bs.all, regions=alu.gr,
	                      type='raw', what='perBase'),
                      length)
mean(alu.gr$meth, na.rm=T)

#Measure methylation of Alu in all BP regions
alu_in_bp <- subsetByOverlaps(alu.gr, bp.lr.gr)
mean(alu_in_bp$meth, na.rm=T)
aluY_in_bp <- alu_in_bp[grepl('AluY', alu_in_bp$class)]
mean(aluY_in_bp$meth, na.rm=T)

#Methylation of Alu in class I
reps_in_classI <- subsetByOverlaps(reps.gr, bp_class1.gr)
alu_in_classI <- reps_in_classI[reps_in_classI$family=='SINE/Alu']
alu_in_classI$meth <- getMeth(bs.all, regions=alu_in_classI,
		              type='raw', what='perRegion')
mean(alu_in_classI$meth, na.rm=T)

mean(alu_in_classI$meth[grepl('AluY', alu_in_classI$class)], na.rm=T)


#Measure methylation in AluY elements
aluY.gr <- alu.gr[grepl('AluY', alu.gr$class)]
mean(aluY.gr$meth, na.rm=T)

aluY_in_bp <- alu_in_bp[grepl('AluY', alu_in_bp$class)]
mean(aluY_in_bp$meth, na.rm=T)

wilcox.test(alu.gr$meth, alu_in_bp$meth)
wilcox.test(aluY.gr$meth, aluY_in_bp$meth)

#Permutation (call python function)
###################################

#Write data to a file
#####################
write_amp_meth <- function(gr, meth, file) {
  out <- data.frame(names(gr),
                    seqnames(gr), 
                    start(gr),
                    end(gr),
		    gr$size,
  		    gr$bp)
  names(out) <- c('name', 'chrom', 'start', 'end',
  	          'size', 'bp')
  out <- cbind(out, meth)
  write.table(out, file=file, sep='\t', quote=F)
}

write_amp_meth(amp.gr, amp.meth, 'amp_meth.txt')

write_rep_meth <- function(gr, meth, file) {
  out <- data.frame(names(gr),
                    seqnames(gr),
                    start(gr),
                    end(gr),
                    gr$bp,
                    gr$notes,
                    gr$class,
                    gr$meth,
                    gr$left_meth,
                    gr$right_meth)
  names(out) <- c('name', 'chrom', 'start', 'end',
                  'size', 'notes', 'class', 'meth',
                  'left_meth', 'right_meth')
  write.table(out, file=file, sep='\t', quote=F)
}

write_meth(bp_region.gr, 'bp_meth.txt')

#Make bedgraph of methylation for viewing in IGV
################################################
cpgs.out <- data.frame(seqnames(bs.all),
                       start(bs.all),
		       end(bs.all), 1)
writeLines("track type=bedGraph name=Vok_CpGs description=center_label visibility=display_mode color=0,0,255 graphType=bar viewLimits=0:1 yLineOnOff=off", 'Vok_cpgs.bedgraph')
write.table(cpgs.out, file='Vok_cpgs.bedgraph', append=T,
            row.names=F, col.names=F, quote=F, sep='\t')

meth.out <- data.frame(seqnames(bs.all),
                       start(bs.all),
		       end(bs.all),
		       getMeth(bs.all, type='raw'))
writeLines("track type=bedGraph name=Vok_Meth description=Vok_Meth visibility=display_mode color=0,0,255 graphType=bar viewLimits=0:1 yLineOnOff=off", 'Vok_meth.bedgraph')

write.table(meth.out, file='Vok_meth.bedgraph', append=T,
            row.names=F, col.names=F, quote=F, sep='\t')

cov.out <- data.frame(seqnames(bs.all),
                      start(bs.all),
                      end(bs.all),
	              getCoverage(bs.all))
writeLines("track type=bedGraph name=Vok_CpG_Cov description=CpG_Cov visibility=display_mode color=255,0,0 graphType=bar yLineOnOff=off", 'Vok_cov.bedgraph')

write.table(cov.out, file='Vok_cov.bedgraph', append=T,
            row.names=F, col.names=F, quote=F, sep='\t')


#Read in amplified regions used on the BACs to compare here
###########################################################
BAC_amp <- read.table(bac_file, sep="\t", header=T,
  blank.lines.skip=T, strip.white=T)
names(BAC_amp)=c('BAC', 'BAC_start', 'BAC_end', 'BP', 'chr', 'strand',
                 'start', 'end', 'size')

BAC_amp.gr <- GRanges(
  seqnames = Rle(BAC_amp$chr),
  ranges = IRanges(start=BAC_amp$start,
    end=BAC_amp$end, 
    names=paste0(BAC_amp$BAC, ":", BAC_amp$BAC_start, "-", BAC_amp$BAC_end)),
  strand = Rle(BAC_amp$strand),
  BP = BAC_amp$BP,
  size = BAC_amp$size,
  seqinfo = seqinfo(bs.all))
BAC_amp.gr <- sort(BAC_amp.gr)

BAC_amp.gr$meth <- getMeth(bs.all, regions=BAC_amp.gr,
                           type='raw', what='perRegion')

mean(BAC_amp.gr$meth[BAC_amp.gr$BP==0], na.rm=T)
mean(BAC_amp.gr$meth[BAC_amp.gr$BP==1], na.rm=T)
wilcox.test(BAC_amp.gr$meth[BAC_amp.gr$BP==0], BAC_amp.gr$meth[BAC_amp.gr$BP==1])
  #This doesn't show significance

BAC_meth.per_cpg <- getMeth(bs.all, regions=BAC_amp.gr,
		            type='raw', what='perBase')
wilcox.test(unlist(BAC_meth.per_cpg[BAC_amp.gr$BP==0]),
	unlist(BAC_meth.per_cpg[BAC_amp.gr$BP==1]))
  #This is highly significant

#Methylation of Alu for BP and non from BAC analysis
alu_in_BAC_bp <- subsetByOverlaps(alu.gr, 
  BAC_amp.gr[BAC_amp.gr$BP==1])
alu_in_BAC_non <- subsetByOverlaps(alu.gr, 
  BAC_amp.gr[BAC_amp.gr$BP==0])
mean(alu_in_BAC_bp$meth, na.rm=T)
mean(alu_in_BAC_non$meth, na.rm=T)

aluY_in_BAC_bp <- alu_in_BAC_bp[grepl('AluY', alu_in_BAC_bp$class)]
aluY_in_BAC_non <- alu_in_BAC_non[grepl('AluY', alu_in_BAC_non$class)]
mean(aluY_in_BAC_bp$meth, na.rm=T)
mean(aluY_in_BAC_non$meth, na.rm=T)
wilcox.test(aluY_in_BAC_bp$meth, aluY_in_BAC_non$meth)

#Same but not taking the average of meth in each Alu
alu_in_BAC_bp.per_cpg <- getMeth(bs.all,
  alu_in_BAC_bp, type='raw', what='perBase')
alu_in_BAC_non.per_cpg <- getMeth(bs.all,
  alu_in_BAC_non, type='raw', what='perBase')
mean(unlist(alu_in_BAC_bp.per_cpg))
mean(unlist(alu_in_BAC_non.per_cpg))
wilcox.test(unlist(alu_in_BAC_bp.per_cpg),
            unlist(alu_in_BAC_non.per_cpg))

aluY_in_BAC_bp.per_cpg <- getMeth(bs.all,
  aluY_in_BAC_bp, type='raw', what='perBase')
aluY_in_BAC_non.per_cpg <- getMeth(bs.all,
  aluY_in_BAC_non, type='raw', what='perBase')
mean(unlist(aluY_in_BAC_bp.per_cpg))
mean(unlist(aluY_in_BAC_non.per_cpg))
wilcox.test(unlist(aluY_in_BAC_bp.per_cpg),
            unlist(aluY_in_BAC_non.per_cpg))



##
#Add column to alu.gr counting the number of CpGs
##

##################################
#Plot mean methylation by coverage
##################################
cov_meth <- data.frame(cov=4:1000, meth=NA, cpgs=0)
cnt <- function(i) length(getMeth(bs.all, type='raw')[getCoverage(bs.all)==i])
meth <- function(i) mean(getMeth(bs.all, type='raw')[getCoverage(bs.all)==i])
cov_meth$cpgs <- unlist(lapply(cov_meth$cov, cnt))
cov_meth$meth <- unlist(lapply(cov_meth$cov, meth))
write.table(cov_meth, file='cov_meth.txt', quote=F,
  sep='\t', col.names=T, row.names=F)
tiff('Vok_cov_meth.tiff', width=1440)
par(mfrow=c(1,3))
plot(c(0,30), c(0,1), type='n', xlab='Coverage',
  ylab='Methylation')
points(cov_meth$cov, cov_meth$meth)
plot(c(0,30), c(0,2500000), type='n', xlab='Coverage',
  ylab='Number of CpGs')
points(cov_meth$cov, cov_meth$cpgs)
plot(c(0,2200000), c(0,1), type='n', xlab='Number of CpGs',
  ylab='Methylation')
points(cov_meth$cpgs, cov_meth$meth)
dev.off()


cov_meth2 <- data.frame(cov=4:1000, meth=NA, cpgs=0)
cnt <- function(i) length(getMeth(bs.all, type='raw')[getCoverage(bs.all)>=i])
meth <- function(i) mean(getMeth(bs.all, type='raw')[getCoverage(bs.all)>=i])
cov_meth2$cpgs <- unlist(lapply(cov_meth2$cov, cnt))
cov_meth2$meth <- unlist(lapply(cov_meth2$cov, meth))
write.table(cov_meth2, file='cov_meth2.txt', quote=F,
  sep='\t', col.names=T, row.names=F)
tiff('Vok_cov_meth2.tiff', width=1440)
par(mfrow=c(1,3))
plot(c(0,30), c(0,1), type='n', 
  xlab='Coverage', ylab='Methylation')
points(cov_meth2$cov, cov_meth2$meth)
plot(c(0,20), c(0,21000000), type='n', 
  xlab='Coverage', ylab='Number of CpGs')
points(cov_meth2$cov, cov_meth2$cpgs)
plot(c(0,21000000), c(0,1), type='n', xlab='Number of CpGs',
  ylab='Methylation')
points(cov_meth2$cpgs, cov_meth2$meth)
dev.off()


##########
#Find CpGs in breakpoint regions
################################
findOverlaps(bp_region.gr, bs.all,ignore.strand=T)
cpg_in_bp <- bs.all[bs.all %within% bp_region.gr]






bp_Class1.gr <- bp_region.gr[bp_region.gr$class=='Class_I']
bp_Class2a.gr <- bp_region.gr[bp_region.gr$class=='Class_II-a']
bp_Class2b.gr <- bp_region.gr[bp_region.gr$class=='Class_II-b']














