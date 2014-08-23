source('~/gibbon_meth/R_meth_functions.R')

test1 <- function()
{
  seqinfo <- make_seqinfo('test.lengths', 'test')
  all.bs <- make_bs_all('cpg', 'test', seqinfo, 4)
  bp.gr <- read_bp('test.bp.txt', seqinfo)
  bp.lr.gr <- make_lr_bp(bp.gr, 500)
  bp.lr.gr <- add_meth_cpg_cov(bp.lr.gr, all.bs)
  bp.permute <- permute(bp.lr.gr, bp.lr.gr, all.bs, n=100,
    type='all', min.chr.size=100, end.exclude=10)

  # Gene analysis
  gene.gr.list <- gtf2GRanges('test.gene.gtf', seqinfo, prom_size=20)
  gene.gr.list <- lapply(gene.gr.list, add_meth_cpg_cov, all.bs)

  gene.permute.list <- lapply(gene.gr.list, permute, bp.lr.gr, all.bs, n=20,
                              type='gene', min.chr.size=100,
                              end.exclude=10)




}