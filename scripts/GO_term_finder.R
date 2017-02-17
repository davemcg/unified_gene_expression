library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(GOstats)

# returns gene symbols that have a user-given GO ID

# initializing biomaRt reference
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)

GO_term_finder <- function(genes,goID) {
  
  # getting entrez IDs for genes of interest (enriched set)
  names(genes) = "Gene"
  geneSymbol_entrez = getBM(attributes=c("external_gene_name", "entrezgene"),filters="external_gene_name",values=genes, mart=mart)
  entrez_ids <- geneSymbol_entrez[,'entrezgene']
  entrez_ids_matched_to_goID <- entrez_ids[entrez_ids %in% get(goID, org.Hs.egGO2ALLEGS)]
  return(geneSymbol_entrez %>% data.frame() %>% filter(entrezgene %in% entrez_ids_matched_to_goID) %>% .[['external_gene_name']])
}
