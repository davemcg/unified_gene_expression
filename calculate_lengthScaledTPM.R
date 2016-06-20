library(tximport)
library(data.table)
library(readr)

# thank our noodly lord for Love
# http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html

# get file paths of kallisto counts
#files <- file.path('~/Desktop/kallisto',paste(fastq_info$run_accession,'_kallisto',sep=''),'abundance.tsv')
files <- paste(list.files(path='/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression', pattern=".*.*_kallisto", full.names=TRUE, include.dirs=TRUE, recursive=TRUE), '/abundance.tsv', sep='')
# convert ensembl tx to gene names
hgnc <- fread('~/git/unified_gene_expression/gencode.v24.metadata.HGNC')
tx2gene <- data.frame(hgnc)
colnames(tx2gene) <- c("TXNAME","GENEID")

#"scaledTPM", or "lengthScaledTPM", for whether to generate estimated counts using 
# abundance estimates scaled up to library size (scaledTPM) or additionally scaled 
# using the average transcript length over samples and the library size (lengthScaledTPM)
txi.lengthScaledTPM <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader = read_tsv, countsFromAbundance = c("lengthScaledTPM"))
lengthScaledTPM <- data.frame(txi.lengthScaledTPM$counts)
# now we have a matrix with each experiment/run in a column and the gene counts (lengthScaledTPM)
# for each row

# but we need to add the names back to the columns, taken from the _kallisto folder names
names <- sapply(files, function(x) strsplit(x,"\\/|_kallisto")[[1]][8])
colnames(lengthScaledTPM) <- names

# save
# head over to plot_by_gene.R
save(lengthScaledTPM,file="lengthScaledTPM.Rdata")
