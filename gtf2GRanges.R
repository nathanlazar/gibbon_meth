# Modified from Jiang (River) Li

library(IRanges)
library(GenomicRanges)

#Return GRanges object of genes or exon information in gtf
gtf2GRanges <- function(myfile="my.gff", seqinfo) {
  gtf <- read.delim(myfile, header=FALSE)
  colnames(gtf) <- c("seqname", "source", "type", "start", "end", "score", "strand", "frame",      
                     "attributes")

  len <- nrow(gtf)

  gene_id <- rep('', len)                                      #get gene_id from attributes column
  idx <- grepl('gene_id', gtf$attributes)
  gene_id[idx] <- 
    gsub(".*gene_id (.*?);.*", "\\1", gtf$attributes[idx])

  transcript_id <- rep('', len)                                #get transcript_id
  idx <- grepl('transcript_id', gtf$attributes)
  transcript_id[idx] <- 
    gsub(".*transcript_id (.*?);.*", "\\1", gtf$attributes[idx])

  exon_number <- rep('', len)                                  #get exon_number
  idx <- grepl('exon_number', gtf$attributes)
  exon_number[idx] <- 
    gsub(".*exon_number (.*?);.*", "\\1", gtf$attributes[idx])

  gene_name <- rep('', len)                                    #get gene_name
  idx <- grepl('gene_name', gtf$attributes)
  gene_name[idx] <- 
    gsub(".*gene_name (.*?);.*", "\\1", gtf$attributes[idx])

  gene_biotype <- rep('', len)                                 #get gene_biotype
  idx <- grepl('gene_biotype', gtf$attributes)
  gene_biotype[idx] <- 
    gsub(".*gene_biotype (.*?);.*", "\\1", gtf$attributes[idx])

  transcript_name <- rep('', len)                              #get transcript_name
  idx <- grepl('transcript_name', gtf$attributes)
  transcript_name[idx] <- 
    gsub(".*transcript_name (.*?);.*", "\\1", gtf$attributes[idx])

  exon_id <- rep('', len)                                      #get exon_id
  idx <- grepl('exon_id', gtf$attributes)
  exon_id[idx] <- 
    gsub(".*exon_id (.*?);.*", "\\1", gtf$attributes[idx])

  gene.gr<-GRanges(seqnames=gtf$seqname,
                   ranges=IRanges(gtf$start,gtf$end),
                   strand=gtf$strand,
		   source=gtf$source,
		   type=gtf$type,
                   gene_id=gene_id,
                   transcript_id=transcript_id,
                   exon_number=as.numeric(exon_number),
		   gene_name=gene_name,
		   gene_biotype=gene_biotype,
		   transcript_name=transcript_name,
		   exon_id=exon_id)

  seqlevels(gene.gr) <- seqlevels(seqinfo)
  seqlengths(gene.gr) <- seqlengths(seqinfo)
  gene.gr
}