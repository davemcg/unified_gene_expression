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

# begin plotting reshaping
lengthScaledTPM$Gene <- row.names(lengthScaledTPM)
melt_lSTPM <- melt(lengthScaledTPM)
metaData_lSTPM <- merge(melt_lSTPM,sra_info[,c("Run_s","SRA_Study_s","Library_Name_s")],by.x="variable",by.y="Run_s")

ggplot(data=subset(metaData_lSTPM,Gene=='A2MP1'),aes(x=Library_Name_s,y=log2(value+1))) + 
  geom_point() + facet_grid(~SRA_Study_s,space='free',scales='free') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("log2(lengthScaledTPM+1) Gene Expression")
