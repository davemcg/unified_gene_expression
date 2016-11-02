library(jsonlite)
library(tidyverse)
library(data.table)

# https://www.encodeproject.org/report/?type=Experiment&assay_title=ChIP-seq&assay_title=RNA-seq&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&replicates.library.biosample.biosample_type=primary+cell&replicates.library.biosample.biosample_type=immortalized+cell+line&replicates.library.biosample.biosample_type=in+vitro+differentiated+cells&replicates.library.biosample.biosample_type=stem+cell&replicates.library.biosample.biosample_type=induced+pluripotent+stem+cell+line&files.file_type=fastq&files.run_type=paired-ended&field=%40id&field=accession&field=assay_term_name&field=assay_title&field=target.label&field=target.gene_name&field=biosample_summary&field=biosample_term_name&field=description&field=lab.title&field=award.project&field=status&field=replicates.biological_replicate_number&field=replicates.technical_replicate_number&field=replicates.antibody.accession&field=replicates.library.biosample.organism.scientific_name&field=replicates.library.biosample.life_stage&field=replicates.library.biosample.age&field=replicates.library.biosample.age_units&field=replicates.library.biosample.treatments.treatment_term_name&field=replicates.library.biosample.treatments.treatment_term_id&field=replicates.library.biosample.treatments.concentration&field=replicates.library.biosample.treatments.concentration_units&field=replicates.library.biosample.treatments.duration&field=replicates.library.biosample.treatments.duration_units&field=replicates.library.biosample.synchronization&field=dbxrefs
encode_metaData <- fread('~/git/unified_gene_expression/data/encode_pairedEnd_RNA-seq_cellLines.tsv')

json_fastq_base_url <- 'https://www.encodeproject.org/search/?type=file&dataset=/experiments/__REPLACE_ME__/&file_format=fastq&format=json&frame=object&limit=all'

encode_base_url <- 'https://www.encodeproject.org'

create_fastq <- function(encode_id){
  json_url <- gsub('__REPLACE_ME__',encode_id, json_fastq_base_url)
  json_metaData <- fromJSON(json_url)
  json_metaData$`@graph` %>% 
    arrange(as.character(biological_replicates), paired_end) %>% 
    filter(!is.na(paired_with)) %>% 
    select(biological_replicates, paired_end, href) %>%
    mutate(href_complete=paste0(encode_base_url,href)) %>% 
    group_by(as.character(biological_replicates)) %>% 
    summarise(reads=paste(as.character(href_complete),collapse=',')) %>% 
    mutate(all_fastq=paste(reads,collapse="__")) %>% 
    select(all_fastq) %>% 
    distinct() %>% 
    .[['all_fastq']]
}

bash_commands_for__encode_to_salmon.py <-apply(encode_metaData,1,function(x) paste(x['Accession'], create_fastq(x['Accession']), 'paired', paste = ' '))
write.table(bash_commands_for__encode_to_salmon.py,file='~/git/unified_gene_expression/scripts/encode_to_salmon.sh',quote=F,row.name=F,col.names = F)
