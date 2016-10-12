library(RSQLite)
library(SRAdb)
library(tidyverse)
library(stringr)
sqlfile <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/SRAmetadb.sqlite'
sqlfile <- getSRAdbFile()
sra_con <- dbConnect(SQLite(),sqlfile)
# list all tables
dbListTables(sra_con)
# list fields in a table
dbListFields(sra_con, 'sra')
# list fields in all tables
sapply(dbListTables(sra_con), function(x) dbListFields(sra_con,x))
# grab human transcriptome (no miRNA) against keyword function
human_transcriptome_sra_info <- function(keyword) {
  sql_query <- 
    'select *  from sra WHERE library_source="TRANSCRIPTOMIC" AND 
  (study_abstract LIKE keyword OR
  experiment_name LIKE keyword OR
  study_name LIKE keyword OR 
  sample_ID LIKE keyword OR
  sample_name LIKE keyword OR
  study_title LIKE keyword OR
  study_description LIKE keyword OR
  description LIKE keyword) AND 
  library_strategy!="miRNA-Seq" AND
  taxon_id="9606"'
  new_query <- gsub('keyword',keyword, sql_query)
  dbGetQuery(sra_con, new_query)
}

human_rpe <- human_transcriptome_sra_info("\"%RPE%\"") 

human_tx_studies <- rbind(
  human_transcriptome_sra_info("\"%RPE%\""),
  human_transcriptome_sra_info("\"%macula%\""),
  human_transcriptome_sra_info("\"%fovea%\""),
  human_transcriptome_sra_info("\"%retina%\""),
  human_transcriptome_sra_info("\"%choroid%\""),
  human_transcriptome_sra_info("\"%sclera%\""),
  human_transcriptome_sra_info("\"%iris%\""),
  human_transcriptome_sra_info("\"%lens%\""),
  human_transcriptome_sra_info("\"%cornea%\""),
  human_transcriptome_sra_info("\"%eye%\"")) %>% distinct()

# check study titles and abstracts by hand to figure out what to keep
human_tx_studies %>% select(study_accession, study_title, study_abstract) %>% distinct()
# these rows. get the accessions for them
# 4, 5, 10, 13, 17, 18, 20, 21, 22, 24, 25, 26, 28, 32, 33, 43, 44, 45, 46 
select_studies <- human_tx_studies %>% select(study_accession) %>% distinct() %>% 
  slice(c(4, 5, 10, 13, 17, 18, 20, 21, 22, 24, 25, 26, 28, 32, 33, 43, 44, 45, 46)) %>% .[['study_accession']]
# OK, now need to select runs. Don't want gene KO / treated samples
# 202 to scan
human_tx_studies %>% filter(study_accession %in% select_studies) %>% nrow()
human_tx_studies %>% filter(study_accession %in% select_studies) %>% select(study_accession, study_title, experiment_title, sample_name, sample_attribute)
# rows to toss because they involved chemical treatment
# 56,57,58,61,63,64,65,66,67
eye_rnaseq_experiments <- human_tx_studies %>% filter(study_accession %in% select_studies) %>% slice(-c(56,57,58,61,63,64,65,66,67))

# check which have functioning ENA fastq links
runs<- eye_rnaseq_experiments$run_accession
fastq_status <- getFASTQinfo(in_acc = runs, sra_con,srcType='ftp')
missing_studies <- fastq_status %>% select(study, ftp) %>% filter(is.na(ftp)) %>% select(study) %>% distinct() %>% .[['study']]
missing_studies
# of these, only 'SRP080886' needs dbGaP permission. I downloaded it through dbGap
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=dbgap_use
# https://hpc.nih.gov/apps/sratoolkit.html
# biowulf2:/data/mcgaugheyd/projects/nei/mcgaughey/unified_gene_expression/dbGap/11588/sra
# the other missing will have to have their reads downloaded in sra format

# the one's that are NOT missing can be streamed via curl in salmon
# find samples with multiple runs (will need to aggregate the fastqs when making the salmon call)
eye_rnaseq_experiments %>% select(sample_accession) %>% group_by(sample_accession) %>% arrange(sample_accession) %>% filter(n()>1)

save(eye_rnaseq_experiments, file='data/eye_rnaSeq_experiments_sraMetadata.Rdata')



fastq_status %>% filter(!study %in% missing_studies, !grepl('_2.fastq.gz',ftp)) %>% group_by(sample) %>% summarise(ftp=paste(ftp,collapse=',')) %>% data.frame()
###
# testing streaming is going poorly. EBI ENA appears to be only sporadically working.
# let's proceed with using NCBI's sra
# listSRAfile(in_acc = runs, sra_con) %>% group_by(sample) %>% summarise(ftp=paste(ftp,collapse=',')) %>% data.frame()
listSRAfile(in_acc = runs, sra_con) %>% filter(study!='SRP080886') %>% mutate(wget_call=paste0('wget ', ftp)) %>% select(wget_call)
# next is to move them into sample_specific folders
listSRAfile(in_acc = runs, sra_con) %>% filter(study!='SRP080886') %>% select(sample) %>% distinct() %>% mutate(mkdir=paste0('mkdir ',sample)) %>% select(mkdir)
listSRAfile(in_acc = runs, sra_con) %>% filter(study!='SRP080886') %>% mutate(mkdir=paste0('mv ', run, '.sra ',sample)) %>% select(mkdir)
# above three commands copied to ~/git/unified_gene_expression/scripts/download_eye_sra_files.sh

