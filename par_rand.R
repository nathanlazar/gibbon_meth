# nathan dot lazar at gmail dot com

par_rand <- function(all.bs, feat.gr, breaks, sizes, lengths, end.exclude, type, reps) {
# Finds random regions in a given genome,
# overlaps with features (if given) and gets the
# mean methylation in features, mean coverage of 
# CpGs in features, total count of CpGs in features
# and percentage of the random regions covered by
# the features

  num <- length(sizes)
  nn <- reps*num

  rand_regions <- data.frame(chr=rep('', nn),
                             start=rep(0, nn),
                             end=rep(0, nn),
                             stringsAsFactors=F)

  for(i in 1:reps) {
    # Choose a set of regions randomly with probability
    # proportional to chromosome lengths
    regions <- mclapply(sizes, get_region, breaks, lengths,
                      end.exclude)
    for(j in 1:num) {
      rand_regions[(i-1)*num + j,] <- regions[[j]]
    }
  }

  rand_regions$start <- as.numeric(rand_regions$start)
  rand_regions$end <- as.numeric(rand_regions$end)
  rand.gr <- makeGRangesFromDataFrame(rand_regions)

  rand <- data.frame(mean.meth=rep(0,n),
                     mean.cov=rep(0,n),
                     tot.cpgs=rep(0,n),
                     per.cov=rep(0,n))

  if(type=='all') {
    # Get methylation and coverage for all random regions at once
    # to minimize overhead
    meth <- mcgetMeth(all.bs, regions=rand.gr, type='raw', what='perBase')
    cov <- mcgetCoverage(all.bs, regions=rand.gr, what='perBase')

    # Get means in groups by the number of regions in each group
    rand$mean.meth <- foreach(i=1:reps) %dopar%
      mean(unlist(meth[((i-1)*num+1):(i*num)]), na.rm=T)
    rand$mean.cov <- foreach(i=1:reps) %dopar% 
      mean(unlist(cov[((i-1)*num+1):(i*num)]), na.rm=T)
    rand$tot.cpgs <- foreach(i=1:reps) %dopar%
      length(unlist(cov[((i-1)*num+1):(i*num)]))
    
  } else {

    foreach(i=1:reps) %dopar% {
      feat.in.rand <- subsetByOverlaps(feat.gr, rand.gr[((i-1)*num+1):(i*num)])

      rand$mean.meth[i] <- weighted.mean(feat.in.rand$meth, feat.in.rand$cpgs,
                                         na.rm=T)
      rand$mean.cov[i] <- weighted.mean(feat.in.rand$cov, feat.in.rand$cpgs,
                                        na.rm=T)
      rand$tot.cpgs[i] <- sum(feat.in.rand$cpgs)

      #get percentage of regions covered by features
      rand.lap <- intersect(rand.gr[((i-1)*num+1):(i*num)], feat.in.rand,
                            ignore.strand=T)
      rand$per.cov[i] <- sum(width(rand.lap))/tot.size
    }
  }
  rand
}
