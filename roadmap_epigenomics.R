# bring in roadmap epigenomics data

# use download_SRAdb.R to get sqlite database for sra

library(RSQLite)
library(stringr)
library(httr)
# connect to full SRA db
# this HUGE file is created with download_SRAdb.R
# change your sqlite filename to match your creation (date will be the day you run download_SRAdb.R)
sra_con<- dbConnect(SQLite(), '2016-06-27.SRAmetadb.sqlite')

# connect to my db
# uge_con <- dbConnect(SQLite(), 'metaData.sqlite')

# load up roadmap_epigenomics metadata
# http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/
load('roadmap_epigenomics_metadata.Rdata')

# grab only mRNA experiments with an run accession num
grab <- 
  roadmap_epigenomics %>% filter(grepl("mrna",tolower(Experiment))) %>% 
  filter(!grepl("smRNA",Experiment)) %>% 
  filter(grepl('ftp',`SRA FTP`))
# get run accession from ftp link
grab$experiment_accession <- 
  sapply(grab$`SRA FTP`,function(x) strsplit(x,'\\/')[[1]][11])


# grab SRAdb metadata to merge into main table in sqlite3 database
ALL_sra_fields <- colnames(dbGetQuery(sra_con,"SELECT * FROM sra WHERE experiment_accession='SRX023777'"))
main_sra_fields <- c('experiment_accession','sample_alias', 'submission_accession', 'run_accession','library_layout','instrument_model')
# initialize empty data frame
main_run_info <- read.table(text='',col.names = main_sra_fields)
run_info <- read.table(text='',col.names = ALL_sra_fields)
for (srx in grab$experiment_accession) {
  fields <- 
  # sql SELECT WHERE query that grab select sra columns for the roadmap epigenetics experiment accessions
  query <- paste(c("SELECT experiment_accession, sample_alias, submission_accession, run_accession, library_layout, instrument_model FROM sra WHERE experiment_accession='",srx,"'"),collapse='')
  data <- dbGetQuery(sra_con, query)
  main_run_info <- rbind(main_run_info,data)
  # now redo and grab all info
  query <- paste(c("SELECT * FROM sra WHERE experiment_accession='",srx,"'"),collapse='')
  data <- dbGetQuery(sra_con, query)
  run_info <- rbind(run_info, data)
}

roadmap_main <- left_join(main_run_info[,!names(main_run_info) %in% c("instrument_model")],grab[,c("experiment_accession","Sample Name")])


# build fastq link from accession
# http://www.ebi.ac.uk/ena/browse/read-download
# explains how to build EBI fastq links from SRR/ERR
# Have SRX with the roadmap_epignomics
ebi_fastq_PE_link_builder <- function(srr) {
  front <- rep('ftp.sra.ebi.ac.uk/vol1/fastq/')
  end_1 <- paste(c(srr,'/',srr,'_1.fastq.gz'),collapse='')
  #end_2 <- paste(c(srr,'/',srr,'_2.fastq.gz'),collapse='')
  if (nchar(srr) == 9) {
    middle <- paste(substr(srr,1,6),'/',sep='')
    link1 <- paste(c(front,middle,end_1),collapse='')
    #link2 <- paste(c(front,middle,end_2),collapse='')
  }
  else {
    middle <- paste(substr(srr,1,6),'/',sep='')
    srr_end_to_grab <- nchar(srr) - 10 
    build_next_middle <- substr(srr,nchar(srr)-srr_end_to_grab,nchar(srr))
    build_next_middle <- str_sub(paste('00000',build_next_middle,sep=''),-3)
    link1 <- paste(c(front,middle,build_next_middle,'/',end_1),collapse='')
    #link2 <- paste(c(front,middle,build_next_middle,'/',end_2),collapse='')
  }
  #out <- paste(link1,link2,sep=';')
  return(link1)
}

fastq_urls <- ''
for (srr in roadmap_main$run_accession){
  if (try(!http_error(ebi_fastq_PE_link_builder(srr)),TRUE) == TRUE) {
    fastq_urls <- c(fastq_urls, ebi_fastq_PE_link_builder(srr))
    print(c(srr,ebi_fastq_PE_link_builder(srr)))
  }
  else {
    link <- ebi_fastq_PE_link_builder(srr)
    link <- gsub("_1.fastq.gz",".fastq.gz",link)
    print(c(srr, link))
    fastq_urls <- c(fastq_urls,link)
  }
}

fastq_urls <- fastq_urls[-1]
roadmap_main$comment_fastq_uri <- fastq_urls



# rename columns to match "main" from grabbing arrayExpress info (main_data_db_pull.R)
colnames(roadmap_main)[which(names(roadmap_main) == "experiment_accession")] <- "comment_ena_experiment"
colnames(roadmap_main)[which(names(roadmap_main) == "sample_alias")] <- "source_name"
colnames(roadmap_main)[which(names(roadmap_main) == "submission_accession")] <- "comment_ena_run"



# only keep the paired-end runs. 
roadmap_core <- roadmap_core %>% filter(grepl("PAIRED",library_layout))



roadmap_core$fastq_ftp <- sapply(roadmap_core$run_accession, function(x) ebi_fastq_PE_link_builder(x))
# high chance of this going pear-shaped next. make a copy
hold <- roadmap_core



















# now need to reshape roadmap_core to match core_data from my database
roadmap_core$ArrayExpressAccession <- 'Roadmap_Epigenomics'
roadmap_core$study_accession <- roadmap_core$submission_accession
roadmap_core$secondary_study_accession <- ''
roadmap_core$sample_accession <- ''
roadmap_core$secondary_sample_accession <- ''
roadmap_core$tax_id <- 9606
roadmap_core$scientific_name <- "Homo sapiens"
roadmap_core$fastq_galaxy <- ''
roadmap_core$submitted_ftp <- ''
roadmap_core$submitted_galaxy <- ''
roadmap_core$cram_index_ftp <- ''
roadmap_core$cram_index_galaxy <- ''
roadmap_core$notes <- 'roadmap_epigenomics.R used to bring data in'
roadmap_core$Tissue <- sapply(roadmap_core$`Sample Name`,function(x) strsplit(x,',')[[1]][1])
roadmap_core$Tissue_Source <- roadmap_core$`Sample Name`
colnames(roadmap_core)[which(names(roadmap_core) == "# GEO Accession")] <- "geo"
roadmap_core$Name <- roadmap_core$sample_alias
colnames_of_core_data <- colnames(dbGetQuery(uge_con,"SELECT * FROM core_data"))

roadmap_core <- roadmap_core[,colnames_of_core_data]

# hand check the above, then can append to core_data
core_data <- dbGetQuery(uge_con,"SELECT * FROM core_data")
core_data <- rbind(core_data,roadmap_core)

# ok, now write to sqlite
dbWriteTable(con,"core_data",data.frame(core_data),overwrite=TRUE)
