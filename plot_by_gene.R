# plotting
library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RSQLite)

# grab lengthScaledTPM data from calculate_lengthScaledTPM.R
load("lengthScaledTPM.Rdata")

### remove outlier
#lengthScaledTPMc <- lengthScaledTPM %>% select(-contains("E-MTAB-4377.RNA11"))
experiments <- ncol(lengthScaledTPM)

# grab metadata (massaged/created code in metadata_to_sqlite.R)
con <- dbConnect(RSQLite::SQLite(), "metaData.sqlite")
core_data <- dbReadTable(con,"core_data")

# begin plotting reshaping for individual gene boxplots
geneTPM <- lengthScaledTPM
geneTPM$Gene <- row.names(geneTPM)
melt_lSTPM <- melt(geneTPM)
colnames(melt_lSTPM) <- c("Gene","run_accession","lsTPM")
metaData_lSTPM <- full_join(melt_lSTPM,core_data[,c("ArrayExpressAccession","run_accession","scientific_name","Name","Tissue","Tissue_Source")])
metaData_lSTPM$Tissue <- tolower(metaData_lSTPM$Tissue)
metaData_lSTPM$Tissue_Source <- tolower(metaData_lSTPM$Tissue_Source)
# manually mark eye tissues in red
metaData_lSTPM <- metaData_lSTPM %>% mutate(EyeMarker = ifelse(Tissue %in% c('retina','eye','rpe'),"Eye","Not Eye"))


shape_max <- length(table((metaData_lSTPM$ArrayExpressAccession)))
gene <- 'ABCA4'
ggplot(data=metaData_lSTPM %>% filter(Gene %in% gene),aes(x=Tissue,y=log2(lsTPM+1),colour=EyeMarker)) + 
  geom_jitter(aes(shape=ArrayExpressAccession),size=2) + geom_boxplot(alpha=0.5) + facet_wrap(~Gene,ncol=1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Gene Expression | log2(lengthScaledTPM+1) ") +
  scale_shape_manual(values=1:shape_max)



### make a heatmap
library(pheatmap)
library(NMF) #aheatmap
library(RColorBrewer)
data_heatmap <- log2(t(lengthScaledTPM[,1:experiments]+1))
ann_heatmap <- data.frame(row.names(data_heatmap))
colnames(ann_heatmap) <- c("run_accession")
ann_heatmap <- left_join(ann_heatmap,core_data,by=c("run_accession"="run_accession"))

data_heatmap <- as.matrix(dist(data_heatmap))
aheatmap(data_heatmap,annCol=ann_heatmap$Tissue_Source)

# left_join with dplyr
data_f <- left_join(data.frame(test),data.frame(sra_info[,c("Run_s","Library_Name_s")]),by=c("SRR"="Run_s"))
# move Library names to row.names
row.names(data_f)<-data_f$Library_Name_s
# remove extraneous columns
data_f <- data_f[,2:27238]

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(as.matrix(dist(log10(data_f+1))),col=colors)

#### t-sne (PCA-like)
library(Rtsne)
library(ggrepel)

set.seed(935489)
tsne_out <- Rtsne(as.matrix(log2(t(lengthScaledTPM[,1:experiments])+1)),perplexity = 30, check_duplicates = FALSE)
# Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$run_accession <- colnames(lengthScaledTPM[,1:experiments])
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



##### pca 
# first get top 500 genes by variance
# hand-pull out outlier "E-MTAB-4377.RNA11"

vars<-apply(lengthScaledTPM,1,function(x) var(x))
topX<-names(head(sort(-vars),n=50000))
lengthScaledTPM_topVar<-lengthScaledTPM[topX,]
pca <- prcomp(log2(t(lengthScaledTPM_topVar[,1:experiments])+1))
pca_data<-data.frame(pca$x)
pca_data$run_accession <- row.names(pca_data)
pca_data <- left_join(pca_data,core_data,by=c("run_accession"="run_accession"))
# show stdev for each PC
plot(pca)
# plot pca mapping
ggplot(data=pca_data,aes(x=PC1,y=PC2, colour=ArrayExpressAccession, label=Tissue, shape=ArrayExpressAccession)) + 
  geom_text_repel(size=2) + geom_point() + theme_bw() + scale_shape_manual(values=1:shape_max)
