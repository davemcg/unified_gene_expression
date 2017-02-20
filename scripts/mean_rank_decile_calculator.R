#calculate mean_rank_deciles for genes by sub_tissue

source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed_2017_02.Rdata')
lengthScaledTPM_processed <- data.frame(lengthScaledTPM_processed)
lengthScaledTPM_processed$Gene.Name <- row.names(lengthScaledTPM_processed)
core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
core_tight$Sub_Tissue <- gsub('_',' - ',core_tight$Sub_Tissue)

mean_rank_decile <- lengthScaledTPM_processed %>% 
  gather(sample_accession, value, -Gene.Name) %>% 
  left_join(.,core_tight) %>% 
  arrange(-value, Sub_Tissue) %>% 
  dplyr::select(Sub_Tissue, sample_accession, Gene.Name, value ) %>% 
  group_by(Sub_Tissue, Gene.Name) %>% summarise(meanlsTPM = mean(value)) %>% 
  mutate(Rank = rank(-meanlsTPM, ties.method='first'), Decile = ntile(meanlsTPM, n = 10)) %>% 
  arrange(Sub_Tissue, Rank)

save(mean_rank_decile, file='~/git/unified_gene_expression/data/mean_rank_decile.Rdata')
