library(tximport)
library(data.table)
library(dplyr)
library(readr)

# sonnesson 2016 inspired ID of low-usage tx
# http://biorxiv.org/content/early/2015/08/24/025387

# working dir, biowulf2
working_dir <- '/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/salmon_counts'
# eyeMac
working_dir <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/salmon_counts'

setwd(working_dir)

# pull in salmon files
files <- list.files(path=working_dir,recursive=TRUE,pattern='quant.sf')

# Gene TX to name conversion
load('~/git/unified_gene_expression/data/gencode_v25_annotation.Rdata')
anno <- gencode_v25_annotation %>% select(Transcript.ID, Gene.Name)

# pull counts
txi <- tximport(files, type = "salmon", tx2gene = anno, reader = read_tsv, txOut = T)

# get counts
tx_c <- data.frame(txi$counts)

# remove samples with low median counts
sample_medians <- apply(tx_c[,2:ncol(tx_c)],2,function(x) median(x))
# remove all with 0 median
samples_to_keep <- names(sample_medians[sample_medians>0])
tx_temp <- tx_c %>% select_(.dots = samples_to_keep)
tx_c <- cbind(tx_c[,1:2],tx_temp)

# get gene name added
tx_c <- merge(data.frame(anno),tx_c,by.x='Transcript.ID',by.y='row.names')
# group by gene and sum counts by gene
gene_sums <- tx_c %>% group_by(Gene.Name) %>% summarise_each(funs(sum), 3:ncol(tx_c))
# add tx name back
gene_sum_tx <- tx_c %>% select(Transcript.ID, Gene.Name) %>% left_join(., gene_sums)
# matrix divide to get ratio
all_ratios <- tx_c[,3:ncol(tx_c)]/gene_sum_tx[,3:ncol(gene_sum_tx)]
# find number of samples for each transcripts which are < 10% of the total
low_usage <- apply(all_ratios,1, function(x) sum(x<0.1,na.rm=T))
# summary
summary(low_usage)
# num with all samples have < 10% tx usage
table(low_usage[low_usage>(ncol(tx_c)-3)])

