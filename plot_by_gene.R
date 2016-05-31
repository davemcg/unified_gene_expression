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
library("RColorBrewer")
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
