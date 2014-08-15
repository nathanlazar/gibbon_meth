# nathan dot lazar at gmail dot com

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
