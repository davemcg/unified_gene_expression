library(data.table)
library(tidyverse)
library(stringr)

load('~/git/unified_gene_expression/data/eye_rnaSeq_experiments_sraMetadata.Rdata')
eye_rnaseq_experiments <- data.frame(eye_rnaseq_experiments
                                     )
# get a feel for how insanely messy this is
eye_rnaseq_experiments %>% dplyr::select(sample_attribute) %>% sample_n(10)

grab_attribute <- function(full_attribute, keyword, delimiter){
  attribute_vector <- sapply(full_attribute, function(x) list(stringr::str_split(x, delimiter)[[1]]))
  attribute_vector
  attribute <- sapply(attribute_vector, function(x) grep(pattern = keyword,x = x,value = T))
  sapply(attribute, function(x) ifelse(length(x)>0, base::strsplit(x,':')[[1]][2], NA))
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
  dplyr::select(Tissue, Cell, Source, Histological, Origin) %>% sample_n(10)
# types
eye_structure <- c('RPE','Retina', 'Cornea','EyeLid')
orig_source <- c('Cell','Tissue')



# ugly logic to fill out eye structure and sample origin (cell line or human tissue)
# found at this point I accidentally brought along a few samples that are melanoma
# filtered out at this point
eye_rnaseq_experiments_extra <-  
  eye_rnaseq_experiments %>% 
  filter(!grepl('melanoma',sample_attribute)) %>% 
  mutate(Tissue=grab_attribute(sample_attribute,'tissue','\\|\\|'),
         Cell=grab_attribute(sample_attribute,'cell type','\\|\\|'),
         Source=grab_attribute(sample_attribute,'source_name|Origen','\\|\\|'),
         Histological=grab_attribute(sample_attribute,'histological type','\\|\\|')) %>% 
  mutate(Tissue=ifelse(grepl(x = sample_attribute, pattern = 'neural|Retina_|tissue: retina', ignore.case=T),'Retina',NA)) %>% 
  mutate(Tissue=ifelse(grepl(x = sample_attribute,pattern = 'RPE|pigment',ignore.case=T),'RPE',Tissue)) %>% 
  mutate(Tissue=ifelse(grepl(x = sample_attribute,pattern = 'histological type: neural retina',ignore.case=T),'Retina',Tissue)) %>%
  mutate(Tissue=ifelse(grepl(x = sample_attribute, pattern = 'cornea|limbus', ignore.case=T),'Cornea',Tissue)) %>%
  mutate(Tissue=ifelse(grepl(x = sample_attribute, pattern = 'lid', ignore.case=T),'EyeLid',Tissue)) %>%
  mutate(Tissue=ifelse(grepl(x = sample_attribute, pattern = 'noggin', ignore.case=T),'Lens', Tissue)) %>% 
  mutate(Tissue=ifelse(is.na(Tissue),'ESC',Tissue)) %>% 
  mutate(Origin=ifelse(grepl('TERT|ATCC|hES|ESC|H9|hfRPE|H1|hiPS2|BG01|HSF1', sample_attribute),'Cell Line','Adult Tissue')) %>%
  mutate(Origin=ifelse(grepl('fetal|fetus|wk',sample_attribute, ignore.case = T),'Fetal Tissue',Origin)) %>% 
  dplyr::select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute, Tissue, Origin)

                          
# hack in E-MTAB-4377
e_mtab_4377 <- fread('~/git/unified_gene_expression/data/E-MTAB-4377.sdrf.txt')
nums <- sapply(e_mtab_4377$`Source Name`, function(x) strsplit(x, '\\s')[[1]][2])
e_mtab_4377$sample_accession <- paste0('E-MTAB-4377.RNA',nums)
e_mtab_4377$run_accession <- paste0('E-MTAB-4377.RNA',nums)
e_mtab_4377$study_accession <- 'E_MTAB_4377'
e_mtab_4377$study_title <- 'An atlas of gene expression and gene co-regulation in the human retina'
e_mtab_4377$study_abstract <- 'The human retina is a specialized tissue involved in light stimulus transduction. Despite its unique biology, an accurate reference transcriptome is still missing. Here, we performed gene expression analysis (RNA-seq) of 50 retinal samples from non-visually impaired post-mortem donors. We identified novel transcripts with high confidence (Observed Transcriptome (ObsT)) and quantified the expression level of known transcripts (Reference Transcriptome (RefT)). The ObsT included 77 623 transcripts (23 960 genes) covering 137 Mb (35 Mb new transcribed genome). Most of the transcripts (92%) were multi-exonic: 81% with known isoforms, 16% with new isoforms and 3% belonging to new genes. The RefT included 13 792 genes across 94 521 known transcripts. Mitochondrial genes were among the most highly expressed, accounting for about 10% of the reads. Of all the protein-coding genes in Gencode, 65% are expressed in the retina. We exploited inter-individual variability in gene expression to infer a gene co-expression network and to identify genes specifically expressed in photoreceptor cells. We experimentally validated the photoreceptors localization of three genes in human retina that had not been previously reported. RNA-seq data and the gene co-expression network are available online (http://retina.tigem.it).'
e_mtab_4377 <- data.frame(e_mtab_4377)
e_mtab_4377 <- e_mtab_4377 %>% mutate(sample_attribute=paste(Characteristics.organism.part.,' || gender: ', Characteristics.sex., ' || age: ', Characteristics.age., ' || post-mortem time: ', Characteristics.total.post.mortem.time., sep='')) %>% 
  mutate(Origin='Adult Tissue', Tissue='Retina') %>% dplyr::select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute, Tissue, Origin)

core_eye_info <- bind_rows(eye_rnaseq_experiments_extra, e_mtab_4377)
core_eye_info <- core_eye_info %>% mutate(Sub_Tissue=paste(Tissue,Origin,sep='_'))

##################################
# gtex
load('~/git/unified_gene_expression/data/gtex_sraMetadata.Rdata')
core_info <- 
  gtex %>% 
  mutate(Tissue=grab_attribute(sample_attribute,'histological type:','\\|\\|')) %>% 
  mutate(Sub_Tissue=grab_attribute(sample_attribute,'body site:','\\|\\|')) %>% 
  mutate(Origin='Tissue') %>% 
  dplyr::select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute, Tissue, Sub_Tissue, Origin) %>% 
  bind_rows(.,core_eye_info) 

# brought along quite a few gender-specific tissues (prostrate, vagina, etc.)
# also taking the opportunity to drop the eye lid tissue
#keepers <- c('RPE','Retina','Cornea',' Adipose Tissue ',' Adrenal Gland ',' Blood ',' Blood Vessel ',' Brain ',' Breast ',' Colon ',' Esophagus ',' Heart ',' Liver ',' Lung ',' Muscle ',' Nerve ',' Pancreas ',' Pituitary ',' Salivary Gland ',' Skin ',' Small Intestine ',' Spleen ',' Stomach ',' Thyroid ')

#core_tight <- core_info %>% filter(Tissue %in% keepers)
#save(core_tight, file='interactive_page/metaData.Rdata')


#########################################
#encode in now
# https://www.encodeproject.org/report/?type=Experiment&assay_title=ChIP-seq&assay_title=RNA-seq&assay_slims=Transcription&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&replicates.library.biosample.biosample_type=primary+cell&replicates.library.biosample.biosample_type=immortalized+cell+line&replicates.library.biosample.biosample_type=in+vitro+differentiated+cells&replicates.library.biosample.biosample_type=stem+cell&replicates.library.biosample.biosample_type=induced+pluripotent+stem+cell+line&files.file_type=fastq&files.run_type=paired-ended&field=%40id&field=accession&field=assay_term_name&field=assay_title&field=target.label&field=target.gene_name&field=biosample_summary&field=biosample_term_name&field=description&field=lab.title&field=award.project&field=status&field=replicates.biological_replicate_number&field=replicates.technical_replicate_number&field=replicates.antibody.accession&field=replicates.library.biosample.organism.scientific_name&field=replicates.library.biosample.life_stage&field=replicates.library.biosample.age&field=replicates.library.biosample.age_units&field=replicates.library.biosample.treatments.treatment_term_name&field=replicates.library.biosample.treatments.treatment_term_id&field=replicates.library.biosample.treatments.concentration&field=replicates.library.biosample.treatments.concentration_units&field=replicates.library.biosample.treatments.duration&field=replicates.library.biosample.treatments.duration_units&field=replicates.library.biosample.synchronization&field=dbxrefs
encode_metaData <- fread('~/git/unified_gene_expression/data/encode_pairedEnd_RNA-seq_cellLines.tsv')
to_join <- encode_metaData %>% 
  mutate(study_accession='ENCODE',study_title='ENCODE',study_abstract='',sample_accession=Accession, run_accession='',sample_attribute=Description, Tissue='ENCODE Cell Line', Sub_Tissue=Biosample, Origin='Cell Line') %>% 
  dplyr::select(study_accession, study_title, study_abstract, sample_accession, run_accession, sample_attribute, Tissue, Sub_Tissue, Origin)

core_tight <- rbind(core_info,to_join)