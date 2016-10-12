# bring in roadmap epigenomics data

library(RSQLite)
library(SRAdb)
library(tidyverse)
library(stringr)
sqlfile <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/SRAmetadb.sqlite'
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)


# load up roadmap_epigenomics metadata
# http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/
load('roadmap_epigenomics_metadata.Rdata')

# grab only mRNA experiments with an run accession num
grab <- 
  roadmap_epigenomics %>% filter(grepl("mrna",tolower(Experiment))) %>% 
  filter(!grepl("smRNA",Experiment)) %>% 
  filter(grepl('ftp',`SRA FTP`))
# get run accession from ftp link
grab$experiment_accession <- 
  sapply(grab$`SRA FTP`,function(x) strsplit(x,'\\/')[[1]][11])


# use experiment accession to match up with full sra db

sra <- dbGetQuery(sra_con, 'SELECT * FROM sra WHERE experiment_accession IN' )
roadmap_epigenomics_experiments <- sra %>% filter(experiment_accession %in% grab$experiment_accession)
rm(sra)

save(roadmap_epigenomics_experiments,file='data/roadmap_epigenomics_sraMetadata.Rdata')




