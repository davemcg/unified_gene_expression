library(tximport)
library(data.table)
library(dplyr)
library(readr)

# thank our noodly lord for Love
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

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

# merge tx specific counts to gene level and scale TPM
txi <- tximport(files, type = "salmon", tx2gene = anno, reader = read_tsv, countsFromAbundance = c("lengthScaledTPM"))

lengthScaledTPM <- data.frame(txi$counts)
names <- sapply(files, function(x) strsplit(x,"\\/")[[1]][1])
colnames(lengthScaledTPM) <- names


