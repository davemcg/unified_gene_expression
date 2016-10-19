library(stringr)
library(tidyverse)

load('data/eye_rnaSeq_experiments_sraMetadata.Rdata')

# get a feel for how insanely messy this is
eye_rnaseq_experiments %>% select(sample_attribute) %>% sample_n(10)

grab_attribute <- function(full_attribute, keyword, delimiter){
  attribute_vector <- sapply(full_attribute, function(x) list(str_split(x, delimiter)[[1]]))
  attribute_vector
  attribute <- sapply(attribute_vector, function(x) grep(pattern = keyword,x = x,value = T))
  sapply(attribute, function(x) ifelse(length(x)>0, strsplit(x,':')[[1]][2], NA))
  #attribute <- attribute_vector[grepl(x = attribute_vector, pattern = keyword)]
  #attribute
}

test <- c('isolate: not applicable || age: not applicable || biomaterial_provider: not applicable || sex: not applicable || tissue: retinal pigment epithelial || cell_line: H1 || cell_subtype: human embryonic stem cells || treatment: induced by three-dimensional cultures methods || treatement_time: 100 days || BioSampleModel: Human',
          'source_name: hES derived neural retina || cell type: GFP positive neural retina cells || day post induction: Day 90',
          'source_name: retina ,  disease status: normal ,  tissue: retina',
          'source_name: Retinal Pigmented Epithelial Cells (hTERT-RPE1) || cell line: RPE1 WT || cell type: Retinal Pigmented Epithelial Cells (hTERT-RPE1) || protocol: RNA-seq')

grab_attribute(test,'tissue','\\|\\|')

eye_rnaseq_experiments %>% mutate(Tissue=grab_attribute(sample_attribute,'tissue','\\|\\|'),
                                  Cell=grab_attribute(sample_attribute,'cell type','\\|\\|'),
                                  Source=grab_attribute(sample_attribute,'source_name','\\|\\|'),
                                  Histological=grab_attribute(sample_attribute,'histological type','\\|\\|'),
                                  Origin=grab_attribute(sample_attribute,'origen','\\|\\|')) %>% 
  select(Tissue, Cell, Source, Histological, Origin) %>% sample_n(10)
# types
eye_structure <- c('RPE','Retina', 'Cornea','EyeLid')
orig_source <- c('Cell','Tissue')

# ugly logic to fill out eye structure and sample origin (cell line or human tissue)
eye_rnaseq_experiments %>% mutate(Tissue=grab_attribute(sample_attribute,'tissue','\\|\\|'),
                                        Cell=grab_attribute(sample_attribute,'cell type','\\|\\|'),
                                        Source=grab_attribute(sample_attribute,'source_name|Origen','\\|\\|'),
                                        Histological=grab_attribute(sample_attribute,'histological type','\\|\\|')) %>% 
                            mutate(Eye_Structure=ifelse(grepl(x= sample_attribute, pattern = 'neural|Retina_|tissue: retina', ignore.case=T),'Retina',NA)) %>% 
                            mutate(Eye_Structure=ifelse(grepl(x = sample_attribute,pattern = 'RPE|pigment',ignore.case=T),'RPE',Eye_Structure)) %>% 
                            mutate(Eye_Structure=ifelse(grepl(x= sample_attribute, pattern = 'cornea|limbus', ignore.case=T),'Cornea',Eye_Structure)) %>%
                            mutate(Eye_Structure=ifelse(grepl(x= sample_attribute, pattern = 'lid', ignore.case=T),'EyeLid',Eye_Structure)) %>%
                            mutate(Eye_Structure=ifelse(is.na(Eye_Structure),'ESC',Eye_Structure)) %>% 
                            mutate(Origin=ifelse(grepl('TERT|ATCC|hES|ESC|H9', sample_attribute),'Cell_Line','Tissue')) %>% 
                            select(sample_attribute, Eye_Structure, Origin)
