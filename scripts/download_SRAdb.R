library(SRAdb)
# run only when you want the latest SRA database
#sqlfile <- getSRAdbFile(destfile = paste(Sys.Date(), "SRAmetadb.sqlite.gz", sep= '.'))
# > sqlfile
# [1] "/Users/mcgaugheyd/git/unified_gene_expression/2016-06-27.SRAmetadb.sqlite"
sra_con<- dbConnect(SQLite(), sqlfile)
# example query
dbGetQuery(con,"select run_accession, experiment_accession from sra where experiment_accession = 'SRX263865'")
