library(RSQLite)
library(SRAdb)
library(tidyverse)


sqlfile <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/SRAmetadb.sqlite'
sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)

load('data/eye_rnaSeq_experiments_sraMetadata.Rdata')
runs <- eye_rnaseq_experiments %>% filter(study_accession!='SRP080886') %>% .[['run_accession']]
###
# testing streaming is going poorly. EBI ENA appears to be only sporadically working.
# let's proceed with using NCBI's sra
listSRAfile(in_acc = runs, sra_con) %>% mutate(wget_call=paste0('wget ', ftp)) %>% select(wget_call)
# above copied to ~/git/unified_gene_expression/scripts/download_eye-related_sra.sh
# next is to move them into sample_specific folders
listSRAfile(in_acc = runs, sra_con) %>% select(sample) %>% distinct() %>% mutate(mkdir=paste0('mkdir ',sample)) %>% select(mkdir)
listSRAfile(in_acc = runs, sra_con) %>% mutate(mkdir=paste0('mv ', run, '.sra ',sample)) %>% select(mkdir)
# above two commands copied to ~/git/unified_gene_expression/scripts/download_eye_sra_files.sh
