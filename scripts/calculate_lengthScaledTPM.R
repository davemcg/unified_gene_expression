source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')
# low median count files removed, from 'analysis/QC_and_qsmooth.Rmd'
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed_01_27_2017.Rdata')

library(tidyverse)
library(stringr)
library(edgeR)

# clean up metadata, removing run accession, and then removing dups. Some samples have multiple runs
core_tight <- core_tight %>% dplyr::select(-run_accession)
core_tight <- core_tight[!duplicated(core_tight),]
core_tight$Tissue = trimws(core_tight$Tissue)

# biowulf Salmon counts 
working_dir <- '/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/salmon_counts_bootstrap50_txUsed/'
# eyeMac
working_dir <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/salmon_counts_bootstrap50_txUsed/'
setwd(working_dir)
files <- list.files(path=working_dir,recursive=TRUE,pattern='quant.sf')

# remove outlier files
# and ENCODE
eye_and_gtex_samples <- core_tight %>% 
  filter(sample_accession %in% colnames(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)) %>% # low median count outliers
  filter(!Tissue %in% 'ENCODE Cell Line') %>% 
  filter(!sample_accession %in% c('SRS523795','SRS360124','SRS360123')) %>% # retina in RPE
  filter(!sample_accession %in% c('SRS332999', 'SRS389531', 'SRS626623', 'SRS627133', 'SRS623923', 'SRS629483')) %>% # few gtex samples that outliers (in weird clusters)
  .[['sample_accession']]

# update files to import
eye_and_gtex_files <- files[files %in% paste(eye_and_gtex_samples,'/quant.sf',sep='')]

# Gene TX to name conversion
load('~/git/unified_gene_expression/data/gencode_v25_annotation.Rdata')
anno <- gencode_v25_annotation %>% dplyr::select(Transcript.ID, Gene.Name)

# tximport to aggregate counts to gene level, counts
txi.eye_and_gtex.counts <- tximport(eye_and_gtex_files, type = "salmon", tx2gene = anno, reader = read_tsv)

# tximport to aggregate counts to gene level, lengthScaledTPM
txi.eye_and_gtex.lsTPM <- tximport(eye_and_gtex_files, type = "salmon", tx2gene = anno, reader = read_tsv, countsFromAbundance = c("lengthScaledTPM"))

# output counts
counts <- data.frame(txi.eye_and_gtex.counts$counts)
names <- sapply(eye_and_gtex_files, function(x) strsplit(x,"\\/")[[1]][1])
colnames(counts) <- names
#save(counts, file='~/git/unified_gene_expression/data/eye_and_gtex_counts_2017_02.Rdata')

# output lsTPM
lengthScaledTPM <- data.frame(txi.eye_and_gtex.lsTPM$counts)
names <- sapply(eye_and_gtex_files, function(x) strsplit(x,"\\/")[[1]][1])
colnames(lengthScaledTPM) <- names
#save(lengthScaledTPM, file='~/git/unified_gene_expression/data/lengthScaledTPM_eye_gtex_2017_02.Rdata')

# low gene count and qsmooth

QC <- function(data, metadata, qsmoothFactor){
  # density plot, pre low gene count removal and qsmooth(?)
  gather_lST<-gather(data, sample_accession) %>% left_join(.,metadata)
  print(ggplot(gather_lST,aes(x=log2(value+1),group=sample_accession)) +
          geom_density() +
          facet_wrap(~Tissue) + 
          ggtitle('No normalization') +
          theme_Publication())
  
  # remove genes with an average lsTPM < 1
  table(rowSums(data)<(ncol(data)))
  data_geneRemoval <- data[(rowSums(data)>ncol(data)),]
  
  # re-do density plot
  gather_lST<-gather(data_geneRemoval,sample_accession) %>% left_join(.,metadata)
  print(ggplot(gather_lST,aes(x=log2(value+1),group=sample_accession)) +
          geom_density() +
          facet_wrap(~Tissue) +
          ggtitle('Low gene count removal') +
          theme_Publication())
  
  # library size normalize
  norm <- DGEList(data_geneRemoval)
  norm <- calcNormFactors(norm)
  tpm <- norm$counts
  correction <- norm$samples %>% data.frame() %>% .[['norm.factors']]
  lsTPM_librarySize <- tpm %*% diag(correction)
  colnames(lsTPM_librarySize) <- colnames(data_geneRemoval)
  
  # re-redo density plot
  gather_lST<-gather(data.table(lsTPM_librarySize),sample_accession) %>% left_join(.,metadata)
  print(ggplot(gather_lST,aes(x=log2(value+1),group=sample_accession)) +
          geom_density() +
          facet_wrap(~Tissue) +
          ggtitle('Low gene count removal AND library size normalize') +
          theme_Publication())
  
  # qsmooth
  qs <- qsmooth(object = lsTPM_librarySize,groupFactor = as.factor(metadata[,qsmoothFactor]))
  lsTPM_qsmooth <- qsmoothData(qs)
  
  # re-re-redo density plot
  gather_lST<-gather(data.table(lsTPM_qsmooth),sample_accession) %>% left_join(.,metadata)
  print(ggplot(gather_lST,aes(x=log2(value+1),group=sample_accession)) +
          geom_density() +
          facet_wrap(~Tissue) +
          ggtitle('Low gene count removal AND library size normalize AND qsmooth') +
          theme_Publication())

  
  return(lsTPM_qsmooth)
}
metadata <- core_tight %>% filter(sample_accession %in% colnames(lengthScaledTPM))
metadata <- metadata[match(colnames(lengthScaledTPM), metadata$sample_accession),]
lengthScaledTPM_processed <- QC(lengthScaledTPM, metadata, 'Tissue')

save(lengthScaledTPM_processed, file='~/git/unified_gene_expression/data/lengthScaledTPM_processed_2017_02.Rdata')



