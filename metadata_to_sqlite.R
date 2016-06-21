# Bring in metadata from ArrayExpress
# and get into sqlite3

# https://cran.r-project.org/web/packages/RSQLite/README.html

library(data.table)
library(RSQLite)
library(dplyr)
library(stringr)
# open sql database
con <- dbConnect(RSQLite::SQLite(), "metaData.sqlite")

# add to this vector as you bring in more experiments
arrayExpress_accessions <- c("E-GEOD-22765","E-GEOD-36695","E-GEOD-40524","E-MTAB-513","E-MTAB-3716")
  # E-MTAB-4377 not here because: 
    # this one had to be modified slightly to make fake "Source Name" match up with the main table
    # I hope this to be a rare case, and thus can fully automate the import 
    # manually add E-MTAB-4377
    # open("E_MTAB_4377_metaData.Rdata")
    # dbWriteTable(con,"E_MTAB_4377",e.mtab.4377)
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

# example code to semi-automate editing of new core_data info 
# first get current core_data
core_data <- dbReadTable(con,"core_data")
# then handfill in accesions for ArrayExpress, GEO, and ENA
ae_accession <- "E_MTAB_513"  # hand-fill
geo_accession <- "GSE30611"   # hand-fill
erp_accession <- "ERP000546"  # hand-fill
# Pull out arrayExpress info in to data frame
arrayExpressTable <- dbReadTable(con,ae_accession)
# build URL to get ENA info that "core_data" is built off of
url <- paste(c("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",erp_accession,"&result=read_run&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,fastq_ftp,fastq_galaxy,submitted_ftp,submitted_galaxy,cram_index_ftp,cram_index_galaxy&download=txt"),collapse='')
new_core_data <- fread(url)
new_core_data$ArrayExpressAccession <- ae_accession
new_core_data$geo <- geo_accession
# move some arrayExpressTable info over to core_data
########## will have to tweak depending on how scientists filled in metadata
arrayExpressTable <- arrayExpressTable %>% select(Comment_ENA_EXPERIMENT, Characteristics_organism_part)
colnames(arrayExpressTable) <- c("experiment_accession","Tissue")
arrayExpressTable$Tissue_Source <- arrayExpressTable$Tissue
# do a (left) join to ensure rows match up
new_core_data <-left_join(data.table(arrayExpressTable),new_core_data,by = setNames("experiment_accession","experiment_accession"))
new_core_data$notes <- NA # don't have notes for this one
new_core_data$Name <- new_core_data$experiment_accession
# CUSTOM FOR THIS EXPERIMENT. DON'T WANT SINGLE-END DATA
new_core_data <- new_core_data %>% filter(library_layout == "PAIRED")
# OK, now get columns in order
new_core_data <- new_core_data %>% select(run_accession,ArrayExpressAccession,geo,study_accession:secondary_sample_accession,experiment_accession,tax_id:cram_index_galaxy,notes,Name,Tissue,Tissue_Source)
new_core_data <- rbind(core_data,new_core_data)
# OK, ready to update the table
# FIRST HAND CHECK
dbWriteTable(con,"core_data",data.frame(new_core_data),overwrite=TRUE)

#################
# example queries
# print out core_data, which holds the core info (ena run info) for the R plotting
# more granular info comes from the individual tables for each experiment
rs <- dbSendQuery(con,"SELECT * FROM core_data")
dbFetch(rs)
# just prints out the ArrayExpress experiments
rs <- dbSendQuery(con,"SELECT ArrayExpressAccession FROM core_data")
unique(dbFetch(rs))
# show metadata available for E_GEOD_40524
rs <- dbSendQuery(con, "SELECT * FROM E_GEOD_40524")
dbFetch(rs)
# show ena run info for E_GEOD_40524 joined with core_data
rs <- dbSendQuery(con, "SELECT * FROM E_GEOD_40524
                  JOIN core_data
                  ON core_data.experiment_accession =
                  E_GEOD_40524.`Comment [ENA_EXPERIMENT]` ")
# print first row for each arrayexpress metadata table to see column names
# BUT easier and more readable to just run sqlite metaData.sqlite '.schema' in bash
rs <- dbSendQuery(con,"SELECT ArrayExpressAccession FROM core_data")
for(table in unique(dbFetch(rs)$ArrayExpressAccession)){
  print(table)
  sql_query <- c("SELECT * FROM ", table, " JOIN core_data ON core_data.experiment_accession = ", 
                 table,".Comment_ENA_EXPERIMENT LIMIT 1")
  sql_query <- paste(sql_query,collapse="")
  rs <- dbSendQuery(con, sql_query)
  print(dbFetch(rs))
  print("")
  }

#prints all across each
rs <- dbSendQuery(con,"SELECT ArrayExpressAccession FROM core_data")
for(table in unique(dbFetch(rs)$ArrayExpressAccession)){
  print(table)
  sql_query <- c("SELECT * FROM ", table)
  sql_query <- paste(sql_query,collapse="")
  rs <- dbSendQuery(con, sql_query)
  print(dbFetch(rs))
  print("")
}



###########

# close connection
dbDisconnect(con)
