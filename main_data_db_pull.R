# Bring in metadata from ArrayExpress
# and get into sqlite3

# https://cran.r-project.org/web/packages/RSQLite/README.html

library(data.table)
library(RSQLite)
library(dplyr)
library(stringr)

# add to this vector as you bring in more experiments
arrayExpress_accessions <- c(
  "E-GEOD-22765", 
  "E-GEOD-36695",
  "E-GEOD-40524",
  "E-MTAB-513",
  "E-MTAB-3716",
  "E-GEOD-26284",
  #"E-MTAB-3871", # Small portion of Roadmap Epigenomics Data, will hand add larger set
  #"E-MTAB-4344", # Small portion of ENCODE, will hand add larger set
  "E-MTAB-2836",
  "E-MTAB-4377")


uge_metadata <- data.frame('')

for(i in arrayExpress_accessions){
  print(i)
  base_url <- "https://www.ebi.ac.uk/arrayexpress/files/"
  end_url <- c(i,"/",i,".sdrf.txt")
  end_url <- paste(end_url,collapse="")
  url <- paste(base_url,end_url,sep="")
  # import and force columns to be unique with base make.unique()
  ae_metadata <- data.frame(fread(url,check.names=TRUE))
  # replace default accession dash with underscore
  db_name <- toupper(gsub("-","_",i))
  
  # clean for sqlite
  # replace special characters with _
  new_colnames <- str_replace_all(colnames(ae_metadata), "[^[:alnum:]]", "_")
  # remove trailing _
  new_colnames <- gsub("_+$","",new_colnames)
  # remove any runs of __ with _
  new_colnames <- gsub("_+","_",new_colnames)
  colnames(ae_metadata) <- new_colnames
  
  # lowercase all column names to normalize any cAmElCaSe situations
  colnames(ae_metadata) <- tolower(colnames(ae_metadata))
  
  # add accession (in this case all arrayExpress, but custom may be added)
  ae_metadata$project_accession <- i

  # use dplyr's bind_rows to bind all together, automatically adding null info
  # if the data frame doesn't have a column in the full-set (guaranteed yes)
  uge_metadata <- bind_rows(list(uge_metadata,ae_metadata))
}

# remove dummy first row and column
uge_metadata <- uge_metadata[-1,-1]

# make sure no dup rows are present
uge_metadata <- uge_metadata[!(duplicated(uge_metadata)),]

# convert to data.table
uge_metadata <- data.table(uge_metadata)

# create main table which holds the crucial info
# this table is a good start for analysis as it holds
# run accession, the ftp link for the ebi/ena fastq, 
# and the tissue
main <- uge_metadata %>% select(project_accession, source_name, comment_ena_run, comment_ena_experiment, characteristics_organism, characteristics_organism_part, factor_value_cell_line, factor_value_sex, comment_library_layout, comment_fastq_uri)

#############
# Up to this point, this is a fully automated process
# However, larger consortium projects tend to maintain their own
# metadata information and are imcompletely shared with 
# large umbrella orgs like NCBI/EBI
# ENCODE, GTEx*, Roadmap Epigenomics are offenders; 
# but they contain huge amounts of useful info
# so will be hand-added
# *GTEx access requires permission granted through dbGaP 
##############

# Roadmap Epigenomics
# http://www.ncbi.nlm.nih.gov/geo/roadmap/epigenomics/
# Accessed June 2016, "Export" button used
source('roadmap_epigenomics.R')








# open sql database
con <- dbConnect(RSQLite::SQLite(), "UGE_metadata.sqlite")