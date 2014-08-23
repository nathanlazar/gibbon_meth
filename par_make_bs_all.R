# nathan dot lazar at gmail dot com

make_bs_all <- function(drive, name='', seqinfo, cov=4) {
###########################################
# Use HTCondor to schedule parallel jobs to
# meke BSseq object with methylation data
# for CpGs with at least <cov> coverage
###########################################

  # Each worker will read one chromosome file from the 
  # given drive and write out a BSseq object
  # These will all be read in and combined
  # into one big bs.all object



  bs.all <- read.bsmooth(drive)
  bs.all <- orderBSseq(bs.all)
  sampleNames(bs.all) <- name
  seqlevels(bs.all) <- seqlevels(seqinfo)
  seqlengths(bs.all) <- seqlengths(seqinfo)
  genome(bs.all) <- genome(seqinfo)
#  seqinfo(bs.all) <- seqinfo
  bs.all[getCoverage(bs.all) >= cov]
}
