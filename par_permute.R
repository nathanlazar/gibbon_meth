# nathan dot lazar at gmail dot com

par_permute <- function(feat.gr, bp.lr.gr, all.bs, n=1000,
                        type=c('gene', 'exon', 'intron', 'promoter', '3UTR',
                               '5UTR', 'CpGisl', 'CpGshore', 'repeat', 'LINE',
                               'SINE', 'DNA', 'LTR', 'SINE', 'Alu', 'AluS',
                               'AluJ', 'AluY', 'MIR', 'all'),
                        min.chr.size=12000, end.exclude=1000) {
# This function performs the same function as permute, but utilizes parallel 
# processing through HTCondor
#
# feat.gr is a GRanges object of features of the given type with
# metadata columns meth, cpgs and cov telling the average methylation,
# the number of cpgs in each range and the average coverage of those cpgs.
# bp.lr.gr is a GRanges object of BP regions
#
# Random regions of the same size and number as in bp.lr.gr are chosen
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
  lengths <- seqlengths(bp.lr.gr)[seqlengths(bp.lr.gr) >=
                               min.chr.size]
  breaks <- cumsum(as.numeric(lengths-(end.exclude*2))) /
              sum(as.numeric(lengths-(end.exclude*2)))
  names(breaks) <- names(lengths)
  sizes <- bp.lr.gr$size

  # Save the objects feat.gr, bp.lr.gr and all.bs, breaks and lengths
  #  to a file that can be read by all workers
  save(feat.gr, bp.lr.gr, all.bs, breaks, sizes, lengths, file='par_permute.dat')

  if(type=='all') {

    # Create HTCondor script and call it
    # Call Condor script

    # Read in files written by HTCondor
    rand <- data.frame()
    for( i in 1:n/reps) {
      fil <- paste0('rand', i, '.dat')
      load(paste0('rand', i, '.dat')
      rand <- rbind(rand, 
    #par_rand(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, reps)

    # See how many of these permutations have methylation as low as
    # the breakpoint regions
    bp.w.av.meth <- weighted.mean(bp.lr.gr$meth, bp.lr.gr$cpgs)
    bp.w.av.cov <- weighted.mean(bp.lr.gr$cov, bp.lr.gr$cpgs)

    ######Report p-values############
    n <- length(!is.na(rand$mean.cov))
    cat('Permutation p-values (random < observed):\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.lr.gr$cpgs))/n, '\n')

  } else {

    tot.size <- sum(sizes)
    feat.in.bp <- subsetByOverlaps(feat.gr, bp.lr.gr)

    bp.w.av.meth <- weighted.mean(feat.in.bp$meth, feat.in.bp$cpgs,
                                  na.rm=T)
    bp.w.av.cov <- weighted.mean(feat.in.bp$cov, feat.in.bp$cpgs)
    bp.cpgs <- sum(feat.in.bp$cpgs)
    overlap <- intersect(bp.lr.gr, feat.in.bp, ignore.strand=T)
    bp.per.cov <- sum(width(overlap))/tot.size

    par_rand(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, n)

    ######Report p-values############
    n <- length(!is.na(rand$mean.cov))
    cat('Permutation p-values (random < observed):\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/n, '\n')
    cat('CpG count: \t', sum(rand$tot.cpgs < sum(bp.lr.gr$cpgs, na.rm=T))/n, '\n')
    cat('Percent region covered: \t', sum(rand$per.cov < bp.per.cov, na.rm=T)/n, '\n')
  }
  ####return dataframe of permutation values#########
  rand
}
