# plotting
library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)

load("lengthScaledTPM.Rdata")
load('fastq_info.Rdata')
sra_info <- fread('SraRunTable.txt')

# hand-fix missing info
sra_info <- data.frame(sra_info)
sra_info[13,'Library_Name_s'] <- 'GSM1099813: hfRPESeq072611'
sra_info[14,'Library_Name_s'] <- "GSM1099814: H1RPE.ksr.061511"
sra_info[15,'Library_Name_s'] <- "GSM1099815: H9RPE.ksr.041511"
sra_info[16,'Library_Name_s'] <- "GSM1099816: H9RPE.ksr.060311"

# hand-edit tissue
sra_info$tissue_s <- c("Retina",rep("RPE",8),rep("Retina",3),rep("RPE",4))

# begin plotting reshaping
lengthScaledTPM$Gene <- row.names(lengthScaledTPM)
melt_lSTPM <- melt(lengthScaledTPM)
metaData_lSTPM <- merge(melt_lSTPM,sra_info[,c("Run_s","SRA_Study_s","Library_Name_s","tissue_s")],by.x="variable",by.y="Run_s")

ggplot(data=subset(metaData_lSTPM,Gene=='A2MP1'),aes(x=Library_Name_s,y=log2(value+1),colour=tissue_s)) + 
  geom_point() + facet_grid(~SRA_Study_s,space='free',scales='free') + 
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
tsne_out <- Rtsne(as.matrix(log(t(lengthScaledTPM[,1:16])+1)),perplexity = 3)
# Perplexity is a measure for information that is defined as 2 to the power of the Shannon entropy. The perplexity of a fair die with k sides is equal to k. In t-SNE, the perplexity may be viewed as a knob that sets the number of effective nearest neighbors. It is comparable with the number of nearest neighbors k that is employed in many manifold learners.
tsne_plot <- data.frame(tsne_out$Y)
tsne_plot$SRR <- colnames(lengthScaledTPM[,1:16])
tsne_plot <- left_join(tsne_plot,sra_info,by=c("SRR"="Run_s"))
ggplot(tsne_plot,aes(x=X1,y=X2,label=SRR,shape=tissue_s,colour=SRA_Study_s)) + 
  geom_text_repel() + geom_point() + theme_bw() + xlab("") + ylab("") + ggtitle("t-sne Clustering")



##### pca 
# first get top 500 genes by variance
vars<-apply(lengthScaledTPM,1,function(x) var(x))
topX<-names(head(sort(-vars),n=5000))
lengthScaledTPM_topVar<-lengthScaledTPM[topX,]
pca <- prcomp(log(t(lengthScaledTPM_topVar[,1:16])+1))
pca_data<-data.frame(pca$x)
pca_data$SRR <- row.names(pca_data)
pca_data <- left_join(pca_data,sra_info,by=c("SRR"="Run_s"))
# show stdev for each PC
plot(pca)
# plot pca mapping
ggplot(data=pca_data,aes(x=PC1,y=PC2, colour=SRA_Study_s, label=Library_Name_s, shape=tissue_s)) + 
  geom_text_repel() + geom_point() + theme_bw()
