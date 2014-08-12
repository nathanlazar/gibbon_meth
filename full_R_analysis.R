#!/usr/bin/env Rscript

#nathan dot lazar at gmail dot com

#R script to colate data on CpG methylation from cpg10 evidence files

#Usage Rscript ./bin/full_R_analysis.R 
#  <file of chromosome lengths>
#  <sorted and tabulated cpg evidence file>
#  <breakpoint region file>
#  <repmask file>

args <- commandArgs(TRUE)
len_file <- args[1]
drive <- args[2]
bp_file <- args[3]
rep_file <- args[4]
bac_file <- 'BAC_amp_on_NomLeu1.0.txt'

library(bsseq)
library(rpython)

#Read in lengths
################
lengths <- read.table(len_file, sep="\t")
lengths <- lengths[with(lengths, order(V1)), ]
  #ordering
seqinfo <- Seqinfo(as.character(lengths[,1]), 
  lengths[,2], NA, "NomLeu1.0")

#Read in CpG data and make BSeq object
######################################
bs.all <- read.bsmooth(drive)
bs.all <- orderBSseq(bs.all)
sampleNames(bs.all) <- "NLE_Vok"
seqlevels(bs.all) <- seqnames(seqinfo)
seqlengths(bs.all) <- seqlengths(seqinfo)

#Subset CpG data with coverage of at least 4 reads
##################################################
bs_cov4 <- bs.all[getCoverage(bs.all) >= 4]

#Read in Breakpoint region data
###############################
bps <- read.table(bp_file, sep="\t", header=T, 
  blank.lines.skip=T, strip.white=T)
bps <- bps[,1:7]
names(bps)=c('BP_name', 'chr', 'start', 'end', 'size', 'notes', 'class')

bp.gr <- GRanges(
  seqnames = Rle(bps$chr),
  ranges = IRanges(start=bps$start,
    end=bps$end, names=bps$BP_name),
  strand = Rle(strand(c("+")), c(nrow(bps))),
  size = bps$size,
  notes = bps$notes,
  class = bps$class,
  seqinfo = seqinfo(bs.all))
bp.gr <- sort(bp.gr)

#Make GRanges objects for regions around breakpoints
#bp_region.gr <- bp.gr
#start(bp_region.gr) <- start(bp.gr) -10000
#end(bp_region.gr) <- end(bp.gr)+10000
#bp_region.gr$size <- end(bp_region.gr) - start(bp_region.gr)

bp_left.gr <- bp.gr[!grepl('Start', bp.gr$notes)]
end(bp_left.gr) <- start(bp_left.gr)
start(bp_left.gr) <- start(bp_left.gr) - 10000
bp_left.gr$size <- end(bp_left.gr) - start(bp_left.gr)
bp_left.gr$side <- 'left'

bp_right.gr <- bp.gr[!grepl('End', bp.gr$notes)]
start(bp_right.gr) <- end(bp_right.gr)
end(bp_right.gr) <- end(bp_right.gr) + 10000
bp_right.gr$size <- end(bp_right.gr) - start(bp_right.gr)
bp_right.gr$side <- 'right'

bp_lr.gr <- c(bp_left.gr, bp_right.gr)
bp_lr.gr <- bp_lr.gr[bp_lr.gr$size> 200]
bp_lr.gr <- sort(bp_lr.gr)
tmp <- reduce(bp_lr.gr, with.mapping=T, min.gapwidth=0)
rows <- unlist(lapply(tmp$mapping, '[', 1))
names(tmp) <- names(bp_lr.gr)[rows]
tmp$size <- end(tmp) - start(tmp)
tmp$notes <- bp_lr.gr$notes[rows]
tmp$class <- bp_lr.gr$class[rows]
tmp$side <- bp_lr.gr$side[rows]
bp_lr.gr <- sort(tmp)

#Write breakpoint regions to bed file
#####################################
bp.bed <- data.frame(seqnames(bp_lr.gr), start(bp_lr.gr), end(bp_lr.gr))
write.table(bp.bed, file='bp_regions.bed', quote=F, sep='\t',
            row.names=F, col.names=F)

#Read in repmask data
#####################
reps <- read.table(rep_file, skip=3)
names(reps) <- c('SW_score', 'perc_div', 'per_del',
	         'per_ins', 'query', 'query_start',
		 'query_end', 'rep_left', 'matching_rep',
		 'rep_class', 'rep_family', 'rep_start', 
		 'rep_end', 'rep_left', 'ID')
#Remove 5 repeats on "random" chromosomes
reps <- reps[grepl('chr', reps$query)==F, ]

reps.gr <- GRanges(
  seqnames = Rle(reps$query),
  ranges = IRanges(start=reps$query_start,
                   end=reps$query_end),
  strand = Rle(strand(c("+")), c(nrow(reps))),
  size = reps$query_end - reps$query_start,
  class = reps$rep_class,
  family = reps$rep_family,
  seqinfo = seqinfo(bs.all))
reps.gr <- sort(reps.gr)

####too slow 
#reps.gr$meth <- getMeth(bs_cov4, 
#	                regions=reps.gr,
#			type='raw',
#			what='perRegion')


#Measure methylation of CpGs with coverage of 4 or more
########################################################

#Full genome
mean(getMeth(bs_cov4, type='raw', what='perBase'))


#bp_region.gr$meth <- getMeth(bs_cov4, regions=bp_region.gr, 
#                             type="raw", what="perRegion")
#bp_region.gr$left_meth <- getMeth(bs_cov4, regions=bp_left.gr,
#                                  type="raw", what="perRegion")
#bp_region.gr$right_meth <- getMeth(bs_cov4, regions=bp_right.gr,
#                                   type="raw", what="perRegion")

#Add methylation of sides of amp regions to gr
##############################################
bp_meth <- getMeth(bs_cov4, regions=bp_lr.gr,
                   type='raw', what='perBase')
bp_cov <- getCoverage(bs_cov4, regions=bp_lr.gr,
                      what='perBase')
bp_lr.gr$cpgs <- unlist(lapply(bp_meth, length))
bp_lr.gr$mean_meth <- unlist(lapply(bp_meth, mean))
bp_lr.gr$sd_meth <- unlist(lapply(bp_meth, sd))
bp_lr.gr$median_meth <- unlist(as.numeric(as.character(
                               lapply(bp_meth, median))))
bp_lr.gr$mean_cov <- unlist(lapply(bp_cov, mean))
bp_lr.gr$sd_cov <- unlist(lapply(bp_cov, sd))
bp_lr.gr$median_cov <- unlist(as.numeric(as.character(
                              lapply(bp_cov, median))))

#Try some graphical model inference
###################################
library(bnlearn)

bp_meta <- mcols(bp_lr.gr)[,c('cpgs', 'mean_meth',
                              'sd_meth', 'median_meth',
                              'mean_cov', 'sd_cov',
                              'median_cov')]
#mark small regions
bp_meta$small <- 0
bp_meta$small[mcols(bp_lr.gr)$size < 8000] <- 1
bp_meta$classI <- 0
bp_meta$classI[mcols(bp_lr.gr)$class=='Class_I'] <- 1
bp_meta$classII <- 0
bp_meta$classII[grepl('Class_II', mcols(bp_lr.gr)$class)] <- 1
bp_meta$classIIa_b <- 0
bp_meta$classIIa_b[mcols(bp_lr.gr)$class=='Class_II-b'] <- 1

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

m <- mean(mcols(bp_lr.gr)$mean_meth, na.rm=T)
wm <- weighted.mean(mcols(bp_lr.gr)$mean_meth,
                    mcols(bp_lr.gr)$cpgs, na.rm=T)

for(type in c('Class_I', 'Class_II-a', 'Class_II-b')) {
  m <- mean(mcols(bp_lr.gr)$mean_meth[mcols(bp_lr.gr)$class==type], na.rm=T)
  wm <- weighted.mean(mcols(bp_lr.gr)$mean_meth[mcols(bp_lr.gr)$class==type],
                      mcols(bp_lr.gr)$cpgs[mcols(bp_lr.gr)$class==type], na.rm=T)
  print(type)
  print(m)
  print(wm)
}

#######################################################
# Permutation type 1
# Compare methylation of CpGs near breaks to the same 
# number of random CpGs in the genome (permutation analysis)

rand.mean.meths <- rep(0,10000)
rand.sd.meths <- rep(0,10000)
rand.mean.covs <- rep(0,10000)
rand.sd.covs <- rep(0,10000)
for(i in 1:10000) {
  rands <- sample(1:20832153, 15032)
  rand.meth <- getMeth(bs_cov4, type='raw')[rands]
  rand.cov <- getCoverage(bs_cov4)[rands]
  rand.mean.meths[i] <- mean(rand.meth)
  rand.sd.meths[i] <- sd(rand.meth)
  rand.mean.covs[i] <- mean(rand.cov)
}

#See how many of these have a value as extreme as CpGs near BP
# regions (get a p-value)
bp.w.av.meth <- weighted.mean(mcols(bp_lr.gr)$mean_meth,
                    mcols(bp_lr.gr)$cpgs)
bp.w.av.cov <- weighted.mean(mcols(bp_lr.gr)$mean_cov,
                    mcols(bp_lr.gr)$cpgs)
sum(rand.mean.meths < bp.w.av.meth)/10000
sum(rand.mean.covs > bp.w.av.cov)/10000

##########################################################
#Permutation type 2
# find random regions matching the sizes of the BP regions
# and measure the methylation of all CpGs in these regions

get_region <- function(size, breaks, seqlen) {
  len <- Inf
  start <- 0
  end <- 0
  while(start < 1 | end > len) {
    #choose chromosome according to length
    idx <- min(which(runif(1) < breaks))
    len <- seqlen[idx]
    #Find random start
    start <- round(runif(1) * (len-size))
    end <- start+size
  }
  c(names(seqlen[idx]), start, end)
}

#Make a list for sampling according to chrom lengths
breaks <- cumsum(as.numeric(seqlengths(bs.all))) / 
            sum(as.numeric(seqlengths(bs.all)))

sizes <- bp_lr.gr$size
seqlen <- seqlengths(bs.all)

rand_mean_meth <- rep(0,1000)
rand_mean_cov <- rep(0,1000)

for(i in 1:1000) {
  regions <- lapply(sizes, get_region, breaks, seqlen)
  regions <- data.frame(do.call(rbind, regions))
  rand.gr <- GRanges(seqnames=regions[,1],
                     ranges=IRanges(start=as.numeric(matrix(regions[,2])),
                                    end=as.numeric(matrix(regions[,3]))))
  meth <- getMeth(bs_cov4, regions=rand.gr, type='raw', what='perBase')
  cov <- getCoverage(bs_cov4, regions=rand.gr, what='perBase')
  rand_mean_meth[i] <- mean(unlist(meth))
  rand_mean_cov[i] <- mean(unlist(cov))
}

#See how many of these permutations have methylation as low as
# the breakpoint regions
sum(rand_mean_meth < bp.w.av.meth)/1000
sum(rand_mean_cov > bp.w.av.cov)/1000

######################################

#CLASS I
bp_class1.gr <- bp_lr.gr[bp_lr.gr$class=="Class_I"]
bp_class1.gr$meth <- getMeth(bs_cov4, regions=bp_class1.gr,
		             type='raw', what='perRegion')
#CLASS I region meth
bp_class1.gr <- bp_lr.gr[bp_lr.gr$class=="Class_I"]
bp_class1.gr$meth <- getMeth(bs_cov4, regions=bp_class1.gr,
		             type='raw', what='perRegion')



#SINE
#####
sines <- reps.gr[grepl('SINE', reps.gr$family)]
sines$meth <- getMeth(bs_cov4, regions=sines,
	              type='raw', what='perRegion')
sines_in_classI <- subsetByOverlaps(sines, bp_class1.gr)
mean(sines_in_classI$meth, na.rm=T)

#Alu analysis
#############

#Mesure methylation in all Alu elements
alu.gr <- reps.gr[reps.gr$family=='SINE/Alu']
alu.gr$meth <- getMeth(bs_cov4, regions=alu.gr,
	               type='raw', what='perRegion')
alu.gr$cpgs <- sapply(getMeth(bs_cov4, regions=alu.gr,
	                      type='raw', what='perBase'),
                      length)
mean(alu.gr$meth, na.rm=T)

#Measure methylation of Alu in all BP regions
alu_in_bp <- subsetByOverlaps(alu.gr, bp_lr.gr)
mean(alu_in_bp$meth, na.rm=T)
aluY_in_bp <- alu_in_bp[grepl('AluY', alu_in_bp$class)]
mean(aluY_in_bp$meth, na.rm=T)

#Methylation of Alu in class I
reps_in_classI <- subsetByOverlaps(reps.gr, bp_class1.gr)
alu_in_classI <- reps_in_classI[reps_in_classI$family=='SINE/Alu']
alu_in_classI$meth <- getMeth(bs_cov4, regions=alu_in_classI,
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

meth.out <- data.frame(seqnames(bs_cov4),
                       start(bs_cov4),
		       end(bs_cov4),
		       getMeth(bs_cov4, type='raw'))
writeLines("track type=bedGraph name=Vok_Meth description=Vok_Meth visibility=display_mode color=0,0,255 graphType=bar viewLimits=0:1 yLineOnOff=off", 'Vok_meth.bedgraph')

write.table(meth.out, file='Vok_meth.bedgraph', append=T,
            row.names=F, col.names=F, quote=F, sep='\t')

cov.out <- data.frame(seqnames(bs_cov4),
                      start(bs_cov4),
                      end(bs_cov4),
	              getCoverage(bs_cov4))
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

BAC_amp.gr$meth <- getMeth(bs_cov4, regions=BAC_amp.gr,
                           type='raw', what='perRegion')

mean(BAC_amp.gr$meth[BAC_amp.gr$BP==0], na.rm=T)
mean(BAC_amp.gr$meth[BAC_amp.gr$BP==1], na.rm=T)
wilcox.test(BAC_amp.gr$meth[BAC_amp.gr$BP==0], BAC_amp.gr$meth[BAC_amp.gr$BP==1])
  #This doesn't show significance

BAC_meth.per_cpg <- getMeth(bs_cov4, regions=BAC_amp.gr,
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
alu_in_BAC_bp.per_cpg <- getMeth(bs_cov4,
  alu_in_BAC_bp, type='raw', what='perBase')
alu_in_BAC_non.per_cpg <- getMeth(bs_cov4,
  alu_in_BAC_non, type='raw', what='perBase')
mean(unlist(alu_in_BAC_bp.per_cpg))
mean(unlist(alu_in_BAC_non.per_cpg))
wilcox.test(unlist(alu_in_BAC_bp.per_cpg),
            unlist(alu_in_BAC_non.per_cpg))

aluY_in_BAC_bp.per_cpg <- getMeth(bs_cov4,
  aluY_in_BAC_bp, type='raw', what='perBase')
aluY_in_BAC_non.per_cpg <- getMeth(bs_cov4,
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
cnt <- function(i) length(getMeth(bs_cov4, type='raw')[getCoverage(bs_cov4)==i])
meth <- function(i) mean(getMeth(bs_cov4, type='raw')[getCoverage(bs_cov4)==i])
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
cnt <- function(i) length(getMeth(bs_cov4, type='raw')[getCoverage(bs_cov4)>=i])
meth <- function(i) mean(getMeth(bs_cov4, type='raw')[getCoverage(bs_cov4)>=i])
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














