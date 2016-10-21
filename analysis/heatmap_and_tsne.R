#### t-sne (PCA-like)
library(Rtsne)
library(ggrepel)

source('~/git/scripts/theme_Publication.R')
source('~/git/unified_gene_expression/scripts/calculate_lengthScaledTPM.R')
source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')

for (n in seq(5,50)) {
  set.seed(935489)
  tsne_out <- Rtsne(as.matrix(log2(t(lengthScaledTPM)+1)),perplexity = n, check_duplicates = FALSE, theta=0.0 )
  # Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
  tsne_plot <- data.frame(tsne_out$Y)
  tsne_plot$sample_accession <- colnames(lengthScaledTPM)
  tsne_plot2<-left_join(tsne_plot,core_eye_info)
  plot <- ggplot(tsne_plot2,aes(x=X1,y=X2,colour=Eye_Structure,shape=study_accession)) + 
    geom_point(size=4) + 
    scale_shape_manual(values=c(0:20)) + 
    ggtitle(paste0("t-sne. Perplexity = ",n)) +
    theme_Publication()
  
  plot_name <- paste0('~/git/unified_gene_expression/analysis/tsne/tsne_perplexity_',str_pad(n,width=2,pad=0),'.png')
  ggsave(filename = plot_name, plot = plot, width=20, height=20, units='cm' )
}
# on-off plot
set.seed(935489)
n=30
tsne_out <- Rtsne(as.matrix(log2(t(lengthScaledTPM)+1)),perplexity = n, check_duplicates = FALSE, theta=0.0 )
# Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$sample_accession <- colnames(lengthScaledTPM)
tsne_plot2<-left_join(tsne_plot,core_eye_info)
ggplot(tsne_plot2,aes(x=X1,y=X2,colour=Eye_Structure,shape=study_accession)) + 
  geom_point(size=4) + 
  scale_shape_manual(values=c(0:20)) + 
  ggtitle(paste0("t-sne. Perplexity = ",n)) +
  theme_Publication()



tsne_plot <- left_join(tsne_plot,core_data,by=c("run_accession"="run_accession"))
# calculate number of studies (for setting num of shapes)
shape_max <- length(table((tsne_plot$ArrayExpressAccession)))
tsne_plot$Tissue <- tolower(tsne_plot$Tissue)
ggplot(tsne_plot,aes(x=X1,y=X2,label=Tissue,shape=ArrayExpressAccession,colour=Tissue)) + 
  geom_text_repel() + geom_point() + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering") +
  scale_shape_manual(values=1:shape_max)
# svg(file='tsne_2016-06-29.svg',width=20,height=15)
# ggplot(tsne_plot,aes(x=X1,y=X2,label=Tissue,shape=ArrayExpressAccession,colour=Tissue)) + 
#   geom_text_repel() + geom_point() + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering") +
#   scale_shape_manual(values=1:shape_max)
# dev.off()
