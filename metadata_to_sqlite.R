# Bring in metadata from ArrayExpress
# and get into sqlite3

# https://cran.r-project.org/web/packages/RSQLite/README.html

library(data.table)
library(RSQLite)
library(dplyr)
library(stringr)
# open sql database
con <- dbConnect(RSQLite::SQLite(), "metaData.sqlite")

arrayExpress_accessions <- c("E-GEOD-22765","E-GEOD-36695","E-GEOD-40524")
  # E-MTAB-4377 not here because: 
    # this one had to be modified slightly to make fake "Source Name" match up with the main table
    # I hope this to be a rare case, and thus can fully automate the import 
    # manually add E-MTAB-4377
    open("E_MTAB_4377_metaData.Rdata")
    dbWriteTable(con,"E_MTAB_4377",e.mtab.4377)
for(i in arrayExpress_accessions){
  base_url <- "https://www.ebi.ac.uk/arrayexpress/files/"
  end_url <- c(i,"/",i,".sdrf.txt")
  end_url <- paste(end_url,collapse="")
  url <- paste(base_url,end_url,sep="")
  # import and force columns to be unique with base make.unique()
  ae_metadata <- fread(url,check.names=TRUE)
  # replace default accession dash with underscore
  db_name <- toupper(gsub("-","_",i))
  # replace special characters with _
  new_colnames <- str_replace_all(colnames(ae_metadata), "[^[:alnum:]]", "_")
  # remove trailing _
  new_colnames <- gsub("_+$","",new_colnames)
  # remove any runs of __ with _
  new_colnames <- gsub("_+","_",new_colnames)
  colnames(ae_metadata) <- new_colnames
  
  current_tables <- dbListTables(con)
  if (!(db_name %in% current_tables)) {
    dbWriteTable(con,db_name,ae_metadata)
    print(paste(db_name,"added"))
  } else{print(paste(db_name,"already present, skipping!"))}
}


open("ena_metadata.Rdata")
# the above is a hand massaged info pulled from EBI ENA Study "Read File" info.
  # see notes for my comments
  # I also hand add the GEO ID and ArrayExpress Accession
# Example for e.geod.36695: 
# http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=SRP011895&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,cram_index_ftp,cram_index_galaxy&download=txt
# can also use fread to pull in
# not sure how to quickly do this, since the SRPXXXXXXX is not in the array expressed metadata, unfortunately
# so for now, just have to use the website

# manually add ena_metadata table
# can either drop and re-add in the future as I add more experiments
# which isn't ideal, but easier
# or do sql queries
dbWriteTable(con,"ena_metadata",ena_metadata)


# example queries
# print out ena_metadata, which holds the core info (ena run info) for the R plotting
# more granular info comes from the individual tables for each experiment
rs <- dbSendQuery(con,"SELECT * FROM ena_metadata")
dbFetch(rs)
# just prints out the ArrayExpress experiments
rs <- dbSendQuery(con,"SELECT ArrayExpressAccession FROM ena_metadata")
unique(dbFetch(rs))
# show metadata available for E_GEOD_40524
rs <- dbSendQuery(con, "SELECT * FROM E_GEOD_40524")
dbFetch(rs)
# show ena run info for E_GEOD_40524 joined with ena_metadata
rs <- dbSendQuery(con, "SELECT * FROM E_GEOD_40524
                  JOIN ena_metadata
                  ON ena_metadata.experiment_accession =
                  E_GEOD_40524.`Comment [ENA_EXPERIMENT]` ")
# print first row for each arrayexpress metadata table to see column names
# easier and more readable to just run sqlite metaData.sqlite '.schema' in bash
rs <- dbSendQuery(con,"SELECT ArrayExpressAccession FROM ena_metadata")
for(table in unique(dbFetch(rs)$ArrayExpressAccession)){
  print(table)
  sql_query <- c("SELECT * FROM ", table, " JOIN ena_metadata ON ena_metadata.experiment_accession = ", 
                 table,".Comment_ENA_EXPERIMENT LIMIT 1")
  sql_query <- paste(sql_query,collapse="")
  rs <- dbSendQuery(con, sql_query)
  print(dbFetch(rs))
  print("")
  }


###########

# close connection
dbDisconnect(con)
