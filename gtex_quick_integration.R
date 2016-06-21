# quick integration of GTEx data
# using public available data
library(dplyr)
library(data.table)

gtex_median <- fread('~/Downloads/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct')
######## convert rpkm (fpkm) to tpm
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# first create new data table
gtex_tpm <- gtex_median[,Name:Description,with=FALSE]
fpkm_to_tpm <- function(fpkm) {exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
gtex_tpm <- 
  data.table(apply(gtex_median[,`Adipose - Subcutaneous`:`Whole Blood`,with=FALSE], 2, function(x) fpkm_to_tpm(x)))

gtex_tpm <- gtex_tpm * 2
gtex_tpm$Name <- gtex_median$Name


# convert IDs
hgnc <- fread('~/git/unified_gene_expression/gencode.v24.metadata.HGNC',header=F)
tx2gene <- data.frame(hgnc)
colnames(tx2gene) <- c("TXNAME","GENEID")
# add column without # of transcripts (e.g. .2)
tx2gene$TXNAME2 <- sapply(tx2gene$TXNAME,function(x) strsplit(x,'\\.')[[1]][1])


# bring in conversion info from gene_id to hgnc name
# source("https://bioconductor.org/biocLite.R")
# biocLite("EnsDb.Hsapiens.v79")
library(EnsDb.Hsapiens.v79)

ensdb <-transcripts(EnsDb.Hsapiens.v79, return.type='data.frame')

# bring together (ENSG_ ENST_ GENENAME)

trio <- data.table(left_join(tx2gene,ensdb,by=c("TXNAME2"="tx_id")))


# ok, now remove .2 from gtex_median Name
gtex_tpm$Name2 <- sapply(gtex_tpm$Name,function(x) strsplit(x,'\\.')[[1]][1])

# join gtex data with GENEID
together <- left_join(data.frame(trio),gtex_tpm,by=c("gene_id"="Name2"))

