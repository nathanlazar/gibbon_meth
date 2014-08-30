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

  # Make directory to store output
  wdir <- paste0('/mnt/lustre1/users/lazar/GIBBONS/VOK_GENOME/', type)
  if(!file.exists(wdir)) {
    dir.create(wdir)
    dir.create(paste0(wdir, '/logs')) 
  }

  # Make condor submit script
  cores <- 16
  jobs <- 63
  make_per_submit('/mnt/lustre1/users/lazar/GIBBONS',
    '/gibbon_meth/wrap_par_rand.R',
    c('$(dir)/VOK_GENOME/par_permute.dat', type, '1000', '$$(Cpus)'),
    wdir, cores, '2 GB', '2 GB', jobs, 'condor.submit')

  # Run condor script
  system("condor_submit condor.submit")

  # If writing over output files, wait a minute so condor can create empty files
  if(file.exists(paste0(wdir, '/permute.', as.character(jobs-1), '.txt')))
    Sys.sleep(60) 

  # Wait for these to be done (there's probably a better way to do this)
  written.files <- rep(0,jobs)
  while(sum(written.files) < jobs) {
    for(i in 1:jobs) {
      if(written.files[i] < 1) {
        f <- paste0(wdir, '/permute.', as.character(i-1), '.txt')
        if (length(readLines(f)) > 0)
          written.files[i] <- 1
      }
    }
    Sys.sleep(5) #check every 5 seconds
  }

  # Read in files written by HTCondor
  rand <- data.frame()
  for(i in 0:(jobs-1)) {
    file <- paste0(wdir, '/permute.', i, '.txt')
    results <- read.table(file, header=T)
    rand <- rbind(rand, data.frame(results))
  }

  tot.size <- sum(sizes)

  if(type=='all') {

    # See how many of these permutations have methylation as low as
    # the breakpoint regions
    bp.w.av.meth <- weighted.mean(bp.lr.gr$meth, bp.lr.gr$cpgs)
    bp.w.av.cov <- weighted.mean(bp.lr.gr$cov, bp.lr.gr$cpgs)
    bp.cpgs.per.kb <- sum(bp.lr.gr$cpgs)/tot.size*1000

    ######Report p-values############
    n <- length(!is.na(rand$mean.cov))
    cat('Permutation p-values (random < observed):\n')
    cat('Methylation:\t',  sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/n, '\n')
    cat('Coverage: \t', sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/n, '\n')
    cat('CpGs per Kb: \t', sum(rand$cpgs.per.kb < bp.cpgs.per.kb)/n, '\n')

  } else {

    feat.in.bp <- subsetByOverlaps(feat.gr, bp.lr.gr)

    bp.w.av.meth <- weighted.mean(feat.in.bp$meth, feat.in.bp$cpgs,
                                  na.rm=T)
    bp.w.av.cov <- weighted.mean(feat.in.bp$cov, feat.in.bp$cpgs)
    bp.cpgs <- sum(feat.in.bp$cpgs)
    bp.cpgs.per.kb <- bp.cpgs/sum(width(feat.in.bp))*1000
    overlap <- intersect(bp.lr.gr, feat.in.bp, ignore.strand=T)
    bp.per.cov <- sum(width(overlap))/tot.size

    ######Report p-values############
    cat('Permutation p-values (random < observed):\n')
    meth.n <- sum(!is.nan(rand$mean.meth))
    meth.p <- sum(rand$mean.meth < bp.w.av.meth, na.rm=T)/meth.n
    cat('Methylation ( n=', meth.n, ') :\t', meth.p, '\n')

    cov.n <- sum(!is.nan(rand$mean.cov))
    cov.p <- sum(rand$mean.cov < bp.w.av.cov, na.rm=T)/cov.n
    cat('Coverage ( n=', cov.n, ') :\t', cov.p, '\n')

    cpg.n <- sum(!is.nan(rand$cpgs.per.kb))
    cpg.p <- sum(rand$cpgs.per.kb < bp.cpgs.per.kb, na.rm=T)/cpg.n
    cat('CpGs per Kb ( n=', cpg.n, ') :\t', cpg.p, '\n')

    per.n <- sum(!is.nan(rand$per.cov))
    per.p <- sum(rand$per.cov < bp.per.cov, na.rm=T)/per.n
    cat('% region covered: ( n=', per.n, ') :\t', per.p, '\n')
  }
  ####return dataframe of permutation values#########
  rand
}
