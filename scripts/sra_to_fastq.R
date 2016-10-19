library(RSQLite)
library(SRAdb)
library(tidyverse)
library(data.table)


sqlfile <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/SRAmetadb.sqlite'
sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)

################
# eye RNA-seq search
# generated from sraDB_search_select.R
#################
load('data/eye_rnaSeq_experiments_sraMetadata.Rdata')
runs <- eye_rnaseq_experiments %>% filter(study_accession!='SRP080886') %>% .[['run_accession']]
###
# testing streaming is going poorly. EBI ENA appears to be only sporadically working.
# let's proceed with using NCBI's sra
# create a call for my biowulf-hosted script which downloads the sra, extracts fastq, then run salmon
swarm_call <- 
  eye_rnaseq_experiments %>% 
  select(sample_accession, library_layout) %>% distinct() %>% 
  left_join(.,listSRAfile(in_acc = runs, sra_con),by=c('sample_accession'='sample')) %>% 
  group_by(sample_accession, library_layout) %>% summarise(ftps=paste(ftp,collapse=',')) %>% 
  mutate(library=ifelse(grepl('PAIR', library_layout, ignore.case = T),'paired','single')) %>% 
  select(sample_accession, ftps, library) %>% ungroup %>% 
  mutate(swarm = paste('~/git/unified_gene_expression/scripts/./sra_to_salmon.py', sample_accession, ftps, library, sep=' ')) %>% 
  select(swarm)

####################
# GTEx
####################
# run table from my dbGaP availables
run_table <- fread('data/gtex_sra_run_table_from_dbGaP.txt')
gtex_table <- run_table %>% filter(biospecimen_repository_s=='GTEx')
table(gtex_table$body_site_s)
table(gtex_table$histological_type_s)
# Let's take a GTEx subset for now
# 5 each from body_site_s, skipping female or male only parts (e.g testis, uterus)
skip_sites <- c('Breast - Mammary Tissue', 'Cervix - Ectocervix', 'Cervix - Endocervix', 'Fallopian Tube', 'Ovary', 'Prostate', 'Testis', 'Uterus', 'Vagina')
set.seed(90243)
gtex_select_runs <- gtex_table %>% filter(!body_site_s %in% skip_sites, analyte_type_s=='RNA:Total RNA') %>% group_by(body_site_s,sex_s) %>% sample_n(5) %>% .[['Run_s']]
# download links
listSRAfile(in_acc = gtex_select_runs, sra_con) %>% select(run)
