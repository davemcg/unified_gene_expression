# plotting
library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)
library(RSQLite)

# grab lengthScaledTPM data from calculate_lengthScaledTPM.R
load("lengthScaledTPM.Rdata")

# grab metadata (massaged/created code in metadata_to_sqlite.R)
con <- dbConnect(RSQLite::SQLite(), "metaData.sqlite")
core_data <- dbReadTable(con,"core_data")

# begin plotting reshaping
lengthScaledTPM$Gene <- row.names(lengthScaledTPM)
melt_lSTPM <- melt(lengthScaledTPM)
colnames(melt_lSTPM) <- c("Gene","run_accession","lsTPM")
metaData_lSTPM <- full_join(melt_lSTPM,core_data[,c("ArrayExpressAccession","run_accession","scientific_name","Name","Tissue","Tissue_Source")])

ggplot(data=subset(metaData_lSTPM,Gene=='A2MP1'),aes(x=Name,y=log2(lsTPM+1),colour=Tissue_Source)) + 
  geom_point() + facet_grid(~ArrayExpressAccession,space='free',scales='free') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("Gene Expression | log2(lengthScaledTPM+1) ")

### make a heatmap
library(pheatmap)
library(RColorBrewer)
## rename SRR to library names
# first move row.names to column
test<- data.frame(t(lengthScaledTPM)) %>% dplyr::add_rownames(var="SRR") 
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
experiments <- ncol(lengthScaledTPM)-1
tsne_out <- Rtsne(as.matrix(log(t(lengthScaledTPM[,1:experiments])+1)),perplexity = 5)
# Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$run_accession <- colnames(lengthScaledTPM[,1:experiments])
tsne_plot <- left_join(tsne_plot,core_data,by=c("run_accession"="run_accession"))
ggplot(tsne_plot,aes(x=X1,y=X2,label=Name,colour=Tissue_Source,shape=ArrayExpressAccession)) + 
  geom_text_repel() + geom_point(size=3) + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")



##### pca 
# first get top 500 genes by variance
# hand-pull out outlier "E-MTAB-4377.RNA11"
lengthScaledTPMc <- lengthScaledTPM %>% select(-contains("E-MTAB-4377.RNA11"))
vars<-apply(lengthScaledTPMc,1,function(x) var(x))
topX<-names(head(sort(-vars),n=500))
lengthScaledTPM_topVar<-lengthScaledTPMc[topX,]
pca <- prcomp(log(t(lengthScaledTPM_topVar[,1:experiments])+1))
pca_data<-data.frame(pca$x)
pca_data$run_accession <- row.names(pca_data)
pca_data <- left_join(pca_data,core_data,by=c("run_accession"="run_accession"))
# show stdev for each PC
plot(pca)
# plot pca mapping
ggplot(data=pca_data,aes(x=PC1,y=PC2, colour=Tissue, label=run_accession, shape=ArrayExpressAccession)) + 
  geom_text_repel() + geom_point() + theme_bw()
