#### t-sne (PCA-like)
library(Rtsne)
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed_2017_02.Rdata')
source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')

core_tight <- core_tight %>% dplyr::select(-run_accession)
core_tight <- core_tight[!duplicated(core_tight),]
core_tight$Tissue = trimws(core_tight$Tissue)

######## first run without fetal and cell line eye tissues###########
tsne_list = list()
eye_and_gtex_samples <- core_tight %>% 
  filter(sample_accession %in% colnames(lengthScaledTPM_processed)) %>% 
  filter(!Tissue %in% 'ENCODE Cell Line') %>% 
  filter(!sample_accession %in% c('SRS523795','SRS360124','SRS360123')) %>% 
  filter(Origin %in% c('Tissue','Adult Tissue')) %>% 
  .[['sample_accession']]
eye_and_gtex_TPM <- lengthScaledTPM_processed[,eye_and_gtex_samples]
for (n in seq(5,50)) {
  set.seed(935489)
  tsne_out <- Rtsne(as.matrix(log2(t(eye_and_gtex_TPM )+1)),perplexity = n, check_duplicates = FALSE, theta=0.0 )
  # Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
  tsne_plot <- data.frame(tsne_out$Y)
  tsne_plot$sample_accession <- colnames(eye_and_gtex_TPM )
  tsne_plot$perplexity <- n
  tsne_list[[n]]<-tsne_plot
}

long_tsne_plot <- do.call(rbind, tsne_list)
save(long_tsne_plot, file='~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity_2017-02.Rdata')

########### second run with fetal, cell line and adult eye tissues with gtex #############

tsne_list = list()
eye_and_gtex_samples <- core_tight %>% 
  filter(sample_accession %in% colnames(lengthScaledTPM_processed)) %>% 
  filter(!Tissue %in% 'ENCODE Cell Line') %>% 
  filter(!sample_accession %in% c('SRS523795','SRS360124','SRS360123')) %>% 
  #filter(Origin %in% c('Tissue','Adult Tissue')) %>% 
  .[['sample_accession']]
eye_and_gtex_TPM <- lengthScaledTPM_processed[,eye_and_gtex_samples]
for (n in seq(5,50)) {
  set.seed(935489)
  tsne_out <- Rtsne(as.matrix(log2(t(eye_and_gtex_TPM )+1)),perplexity = n, check_duplicates = FALSE, theta=0.0 )
  # Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
  tsne_plot <- data.frame(tsne_out$Y)
  tsne_plot$sample_accession <- colnames(eye_and_gtex_TPM )
  tsne_plot$perplexity <- n
  tsne_list[[n]]<-tsne_plot
}

long_tsne_plot <- do.call(rbind, tsne_list)
save(long_tsne_plot, file='~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity_2017-02__all_eye.Rdata')