library(data.table)
library(tximport)
library(dplyr)


example_kallisto <- fread("/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/E-MTAB-513/ERR030876_kallisto/abundance.tsv")
# add column without the .number enst
example_kallisto$target_id.less <- sapply(example_kallisto$target_id,function(x) strsplit(x,'\\.')[[1]][1])

                                          
gtex_rpkm <- fread('gzcat /Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/GTEx/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz')
gtex_counts <- fread('~/Downloads/Flux_SetSummary_Merged.txt')
gtex_rpkm[,1:10,with=FALSE]
gtex_counts[,1:10,with=FALSE]
fpkm_to_tpm <- function(fpkm) {exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
gtex_tpm <- 
  data.table(apply(gtex_rpkm[,3:8557,with=FALSE], 2, function(x) fpkm_to_tpm(x)))
samples <- colnames(gtex_tpm)

gtex_enst <- sapply(gtex_rpkm$TargetID,function(x) strsplit(x,'\\.')[[1]][1])
for(i in samples) {
  temp <- cbind(gtex_enst,gtex_tpm[,i,with=FALSE],gtex_counts[,i,with=FALSE])
  temp <- left_join(temp,example_kallisto,by=c("gtex_enst"="target_id.less"))
  temp <- 
  colnames(temp) <- c("TargetID","TPM","Counts")
  
  #output_name <- paste('/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/GTEx/',i,"_tpm_fluxCounts.Rdata",collapse='')
  #save(temp,file=output_name)
}