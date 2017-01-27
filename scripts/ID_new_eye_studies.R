library(RSQLite)
# Checks for eye-related studies
getSRAdbFile(destdir='/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/',destfile='SRAmetadb_2017-01-19.sqlite.gz')
sqlfile <- '/Volumes/ThunderBay/PROJECTS/mcgaughey/unified_gene_expression/SRAmetadb_2017-01-19.sqlite'
sra_con <- dbConnect(RSQLite::SQLite(),sqlfile)


# grab latest possible studies
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

#human_rpe <- human_transcriptome_sra_info("\"%RPE%\"") 

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


# load previous search results
load('data/eye_studies_considered_2016-10-12.Rdata')


# ID studies not in the last search
old_studies <- eye_studies_considered %>% .[['study_accession']]
current_studies <- human_tx_studies %>% select(study_accession, study_title, study_abstract) %>% distinct() %>% .[['study_accession']]
new_studies <- setdiff(current_studies, old_studies)

# hand check to see which of the new studies will be added
human_tx_studies %>% filter(study_accession %in% new_studies) %>% select(study_accession,study_title,study_abstract) %>% distinct()
# on 2017-01-19, keeping these two studies: SRP091605 and SRP080002
# now check if any individual samples should be tossed:
human_tx_studies %>% filter(study_accession %in% c('SRP091605','SRP080002')) %>% select(sample_attribute)
# add to existing eye studies
load('data/eye_rnaSeq_experiments_sraMetadata.Rdata')
eye_rnaseq_experiments <- bind_rows(eye_rnaseq_experiments,human_tx_studies %>% filter(study_accession %in% c('SRP091605','SRP080002')) )

# save new eye metadata
save(eye_rnaseq_experiments, file='data/eye_rnaSeq_experiments_sraMetadata.Rdata')
# save new studies considered

save()