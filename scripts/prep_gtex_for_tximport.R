library(data.table)
library(tximport)
library(dplyr)


example_kallisto <- fread("/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/E-MTAB-513/ERR030876_kallisto/abundance.tsv")
# add column without the .number enst
example_kallisto$target_id.less <- sapply(example_kallisto$target_id,function(x) strsplit(x,'\\.')[[1]][1])

# gtex RPKM (will be converted to TPM)
gtex_rpkm <- fread('gzcat /Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/GTEx/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_reads.txt.gz')
# gtex counts (calculated from "Flux", which no one else uses)
gtex_counts <- fread('~/Downloads/Flux_SetSummary_Merged.txt')
# print out some info
gtex_rpkm[,1:10,with=FALSE]
gtex_counts[,1:10,with=FALSE]
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
fpkm_to_tpm <- function(fpkm) {exp(log(fpkm) - log(sum(fpkm)) + log(1e6))}
gtex_samples <- colnames(gtex_rpkm)[-c(1,2,3,4)]
# http://stackoverflow.com/questions/9236438/how-do-i-run-apply-on-a-data-table
# replace rpkm with tpm
for (sample in gtex_samples){
  gtex_rpkm[, sample := fpkm_to_tpm(gtex_rpkm[[sample]]), with=FALSE ]
}
# pull out enst and remove the .number notation for left_join
gtex_enst <- sapply(gtex_rpkm$TargetID,function(x) strsplit(x,'\\.')[[1]][1])
# loop through each gtex sample and create individual file for calculate_lenghtScaledTPM.R
for(i in gtex_samples) {
  temp <- cbind(gtex_enst,gtex_rpkm[,i,with=FALSE],gtex_counts[,i,with=FALSE])
  colnames(temp) <- c("target_id.less","tpm","est_counts")
  temp <- plyr::join(example_kallisto[,c("target_id","target_id.less","length","eff_length"),with=FALSE],temp,type="left")
  #temp <- temp[,-1,with=FALSE]
  temp <- temp %>% select(target_id,length,eff_length,est_counts,tpm)
  temp[is.na(temp)] <- 0
  output_name <- paste(c('/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/GTEx/',i,"_tpm_fluxCounts.tsv"),collapse='')
  write.table(temp,file=output_name,sep="\t",quote=F,row.names=F)
}