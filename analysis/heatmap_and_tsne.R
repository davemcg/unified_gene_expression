#### t-sne (PCA-like)
library(ggplot2)
library(tidyverse)
library(Rtsne)
library(gganimate)
library(formattable)
source('~/git/scripts/theme_Publication.R')
load('~/git/unified_gene_expression/data/lengthScaledTPM_eye_gtex.Rdata')
load('~/git/unified_gene_expression/data/gencode_v25_gtf_annotation.Rdata')


tsne_list = list()
# remove NA-filled samples
lengthScaledTPM <- lengthScaledTPM[,!(is.na(lengthScaledTPM[1,]))]
# skip X, Y, mito chromosomes
gene_names_to_keep <- gtf_info %>% filter(!chr %in% c('chrM','chrX','chrY')) %>% .[['gene.Name']]
lengthScaledTPM<-lengthScaledTPM[row.names(lengthScaledTPM) %in% gene_names_to_keep,]
for (n in seq(5,50)) {
  set.seed(935489)
  tsne_out <- Rtsne(as.matrix(log2(t(lengthScaledTPM)+1)),perplexity = n, check_duplicates = FALSE, theta=0.0 )
  # Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
  tsne_plot <- data.frame(tsne_out$Y)
  tsne_plot$sample_accession <- colnames(lengthScaledTPM)
  tsne_plot$perplexity <- n
  tsne_list[[n]]<-tsne_plot
}
long_tsne_plot <- do.call(rbind, tsne_list)
save(long_tsne_plot, file='~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity.Rdata')
