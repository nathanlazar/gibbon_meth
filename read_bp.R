# nathan dot lazar at gmail dot com

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

