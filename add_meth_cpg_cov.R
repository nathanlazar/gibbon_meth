# nathan dot lazar at gmail dot com

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
