# Functions used in analysis of bisulfite reads using BSmooth
# nathan dot lazar at gmail dot com

make_seqinfo <- function(len_file) {
#########################################
# Read in lengths and make seqinfo object
#########################################
  lengths <- read.table(len_file, sep="\t")
  #order
  lengths <- lengths[with(lengths, order(V1)), ]
  Seqinfo(as.character(lengths[,1]),
    lengths[,2], NA, "NomLeu1.0")
}


make_bs_all <- function(drive, names='', seqinfo, cov=4) {
#########################################
# Meke BSseq object with methylation data
# for CpGs with at least <cov> coverage
#########################################
  bs.all <- read.bsmooth(drive)
  bs.all <- orderBSseq(bs.all)
  sampleNames(bs.all) <- names
  seqinfo(bs.all) <- seqinfo
  bs.all[getCoverage(bs.all) >= cov]
}


read_bp <- function(bp_file) {
########################################################
# Read in Breakpoint region data and make GRanges object
########################################################
  bps <- read.table(bp_file, sep="\t", header=T,
    blank.lines.skip=T, strip.white=T)
  bps <- bps[,1:7]
  names(bps)=c('BP_name', 'chr', 'start', 'end', 'size', 'notes', 'class')

  bp.gr <- makeGRangesFromDataFrame(bps, keep.extra.columns=T,
             seqinfo=seqlengths(bs.all))
  sort(bp.gr)
}


make_lr_bp <- function(bp.gr, span) {
##############################################################
# Make GRanges object of <span> regions on each side of breaks
##############################################################
  bp.left.gr <- bp.gr[!grepl('Start', bp.gr$notes)]
  end(bp.left.gr) <- start(bp.left.gr)
  start(bp.left.gr) <- start(bp.left.gr) - span
  bp.left.gr$size <- end(bp.left.gr) - start(bp.left.gr)
  bp.left.gr$side <- 'left'

  bp.right.gr <- bp.gr[!grepl('End', bp.gr$notes)]
  start(bp.right.gr) <- end(bp.right.gr)
  end(bp.right.gr) <- end(bp.right.gr) + span
  bp.right.gr$size <- end(bp.right.gr) - start(bp.right.gr)
  bp.right.gr$side <- 'right'

  bp.lr.gr <- c(bp.left.gr, bp.right.gr)
  bp.lr.gr <- sort(bp.lr.gr)

  # Combine ranges that overlap
  tmp <- reduce(bp.lr.gr, with.mapping=T, min.gapwidth=0)

  # Keep metadata for first of overlapping regions
  rows <- unlist(lapply(tmp$mapping, '[', 1))
  mcols(tmp) <- cbind(mcols(tmp), mcols(bp.lr.gr)[rows,])

  # If regions were combined, concatenate metadata BP_name,
  # notes, class, side with a '-'
  cats <- which(lapply(tmp$mapping, length) > 1)
  tmp.names <- as.character(tmp$BP_name)
  tmp.names[cats] <-  paste(as.character(
    bp.lr.gr$BP_name[unlist(tmp[cats]$mapping)]), collapse='-')
  tmp$BP_name <- as.factor(tmp.names)
  tmp.notes <- as.character(tmp$notes)
  tmp.notes[cats] <-  paste(as.character(
    bp.lr.gr$notes[unlist(tmp[cats]$mapping)]), collapse='-')
  tmp$notes <- as.factor(tmp.notes)
  tmp.notes <- as.character(tmp$notes)
  tmp.notes[cats] <-  paste(as.character(
    bp.lr.gr$notes[unlist(tmp[cats]$mapping)]), collapse='-')
  tmp$notes <- as.factor(tmp.notes)
  tmp.side <- as.character(tmp$side)
  tmp.side[cats] <-  paste(as.character(
    bp.lr.gr$side[unlist(tmp[cats]$mapping)]), collapse='-')
  tmp$side <- as.factor(tmp.side)

  tmp
}

make_rep_gr <- function(rep_file, seqinfo) {
##################################################
# Reads repeat data from repeat masker output file
# <rep_file> and makes a GRanges object
##################################################
  reps <- read.table(rep_file, skip=3)
  names(reps) <- c('SW_score', 'perc_div', 'per_del',
                   'per_ins', 'chr', 'start',
                   'end', 'rep_left', 'matching_rep',
                   'rep_class', 'rep_family', 'rep_start',
                   'rep_end', 'rep_left', 'ID')
  #Remove 5 repeats on "random" chromosomes
  reps <- reps[grepl('chr', reps$chr)==F, ]

  reps.gr <- makeGRangesFromDataFrame(reps,
               keep.extra.columns=T,
               seqinfo=seqinfo)
  sort(reps.gr)
}

add_meth_cpg_cov <- function(gr, bs.all) {
########################################################
# Add metadata column of mean methylation per region,
# number of CpGs per region and mean coverage per region
########################################################
  gr$meth <- getMeth(bs.all, regions=gr, type='raw',
                     what='perRegion')
  cov <- getCoverage(bs.all, regions=gr,
                     what='perBase')
  gr$cpgs <- unlist(lapply(cov, length))
  gr$cov <- unlist(lapply(cov, mean, na.rm=T))
  gr
}


get_region <- function(size, breaks, lengths) {
#############################################
# Choose a random region from the genome
# Chromsomes are chosen by length and regions
# exclude 1kb at the start and end of each 
# chromosome
#############################################
  idx <- min(which(runif(1) < breaks))
  #Find random start (excluding first & last 1kb)
  len <- lengths[idx]
  start <- round(runif(1, 1000, len-size-1000))
  end <- start + size
  c(names(lengths)[idx], start, end)
}

permute <- function(feats.gr, bp.gr, all.bs, n=1000,
                    type=c('gene', 'exon', 'intron', 'promoter', '3UTR', 
                           '5UTR', 'CpGisl', 'CpGshore', 'repeat', 'LINE', 
                           'SINE', 'DNA', 'LTR', 'SINE', 'Alu', 'AluS', 
                           'AluJ', 'AluY', 'MIR', 'all'),
                    min.chr.size=12000, end.exclude=1000)
# feats.gr is a GRanges object of features of the given type with
# metadata columns meth, cpgs and cov telling the average methylation,
# the number of cpgs in each range and the average coverage of those cpgs.
# bp.gr is a GRanges object of BP regions
#
# Random regions of the same size and number as in bp.gr are chosen
# from the genome (excluding 1kb on the end of scaffolds).
# Methylation of features in these regions is calculated as the average
# of the feature methylation weighted by how many CpGs there are in each
# feature. This is equivalent to an average of all CpG methylation for
# CpGs in the given set of features.
# 
# The percentage of the regions covered by the type of feature is also
# recorded and compared to the BP regions.
#
# If type=='all' then all CpGs in the random regions are measured
# Prints out p-values and returns dataframe of mean methylation,
# mean coverage and area covered (if type is a feature)

  # Make a list for sampling according to chrom lengths
  # only sample from chroms at least <min.chr.size>
  lengths <- seqlengths(bp.gr)[seqlengths(bp.gr) >=
                               min.chr.size]
  breaks <- cumsum(as.numeric(lengths-(end.exclude*2))) /
              sum(as.numeric(lengths-(end.exclude*2)))
  names(breaks) <- names(lengths)

  sizes <- bp.gr$size
  num <- length(sizes)

  rand_regions <- data.frame(chr=character(),
                             start=numeric(0),
                             end=numeric(0),
                             stringsAsFactors=F)

  for(i in 1:n) {
    # Choose a set of regions randomly with probability
    # proportional to chromosome lengths
    regions <- lapply(sizes, get_region, breaks, lengths)
    regions <- data.frame(do.call(rbind, regions),
                          stringsAsFactors=F)
    names(regions) <- c('chr', 'start', 'end')
    rand_regions <- rbind(rand_regions, regions)
  }

  rand_regions$start <- as.numeric(rand_regions$start)
  rand_regions$end <- as.numeric(rand_regions$end)
  rand.gr <- makeGRangesFromDataFrame(rand_regions)

  rand <- data.frame(mean.meth=rep(0,n),
                     mean.cov=rep(0,n),
	             tot.cpgs=rep(0,n),
                     per.cov=rep(0,n))

  if(type=='all') {
    meth <- getMeth(all.bs, regions=rand.gr, type='raw', what='perBase')
    cov <- getCoverage(all.bs, regions=rand.gr, what='perBase')

    # Get means in groups by the number of regions in each group
    for(i in 1:n) {
      rand$mean.meth[i] <- mean(unlist(meth[((i-1)*num+1):(i*num)]))
      rand$mean.cov[i] <- mean(unlist(cov[((i-1)*num+1):(i*num)]))
      rand$tot.cpgs[i] <- length(unlist(cov[((i-1)*num+1):(i*num)]))
    }

    # See how many of these permutations have methylation as low as
    # the breakpoint regions
    bp.w.av.meth <- weighted.mean(bp.gr$meth, bp.gr$cpgs)
    bp.w.av.cov <- weighted.mean(bp.gr$cov, bp.gr$cpgs)

    ######Report p-values############
    cat('Permutation p-values:\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.gr$cpgs)), '\n')

  } else {

    tot.size <- sum(sizes)
    feat.in.bp <- subsetByOverlaps(feat.gr, bp.gr)

    bp.w.av.meth <- weighted.mean(feat.in.bp$meth, feat.in.bp$cpgs)
    bp.w.av.cov <- weighted.mean(feat.in.bp$cov, feat.in.bp$cpgs)
    bp.cpgs <- sum(feat.in.bp$cpgs)
    bp.cvg <- coverage(c(bp.gr, feat.in.bp))
    bp.per.cov <- sum(unlist(lapply(bp.cvg[bp.cvg>1], runLength)))/tot.size

    for(i in 1:n) {
      feat.in.rand <- subsetByOverlaps(feat.gr, rand.gr[((i-1)*num+1):(i*num)])

      rand$mean.meth[i] <- weighted.mean(feat.in.rand$meth, feat.in.rand$cpgs)
      rand$mean.cov[i] <- weighted.mean(feat.in.rand$cov, feat.in.rand$cpgs)
      rand$tot.cpgs[i] <- sum(feat.in.rand$cpgs)

      #get percentage of regions covered by features
      cvg <- coverage(c(rand.gr[((i-1)*num+1):(i*num)], feat.in.rand))
      rand$per.cov[i] <- sum(unlist(lapply(cvg[cvg>1], runLength)))/tot.size
    }

    ######Report p-values############
    cat('Permutation p-values:\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.gr$cpgs))/n, '\n')
    cat('Percent region covered: \t', sum(rand$per.cov < bp.per.cov)/n, '\n')
  }
  ####return dataframe of permutation values#########
  rand
}

