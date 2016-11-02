library(jsonlite)
library(tidyverse)
library(data.table)

hand_selected_encode <- fread('~/git/unified_gene_expression/data/encode_web_search_paired_end_RNA-seq.tsv')

search_url_base <- 'https://www.encodeproject.org/search/?type=file&dataset=/experiments/REPLACE_ME/&file_format=fastq&format=json&frame=object&limit=all'
encode_sample <- 'ENCSR971GPJ'
json_url <- gsub('REPLACE_ME',encode_sample,search_url_base)

metadata<-fromJSON(json_url, flatten=T)
download_base <- 'https://www.encodeproject.org'

read_pairs <- metadata$`@graph` %>% arrange(as.character(technical_replicates)) %>% 
  mutate(nhref=paste0(download_base, href)) %>% select(technical_replicates, paired_end,nhref)


