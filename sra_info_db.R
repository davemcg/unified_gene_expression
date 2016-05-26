# import SRA info into sqlite

library(SRAdb)
library(DBI)
srafile = getSRAdbFile()
con = dbConnect(RSQLite::SQLite(), srafile)
# Give up and just use EBI's web page for now. Can use GSE to search and export easily.
# getFASTQinfo not pulling enough useful info