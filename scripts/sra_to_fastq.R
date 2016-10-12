# get you some fastq from NCBI's sra

load('data/eye_rnaSeq_experiments_sraMetadata.Rdata')

# create fastq-dump call, differentiating for single and paired end
# passing on SRP080886 because it's in dbGaP (protected)
# paired
# first check all are labeled paired or single
table(eye_rnaseq_experiments$library_layout)
eye_rnaseq_experiments %>% filter(study_accession!='SRP080886',grepl('PAIR',library_layout)) %>% mutate(fastqdump_call=paste0('fastq-dump --split-files --gzip ', run_accession)) %>% select(fastqdump_call)
eye_rnaseq_experiments %>% filter(study_accession!='SRP080886',grepl('SINGLE',library_layout)) %>% mutate(fastqdump_call=paste0('fastq-dump --gzip ', run_accession)) %>% select(fastqdump_call)
# copy the above two lines to ~/git/unified_gene_expression/data/eye_search_fastq-dump_2016-10-12.swarm
# invoke with swarm -f ~/git/unified_gene_expression/data/eye_search_fastq-dump_2016-10-12.swarm --module sratoolkit --time=16:00:00