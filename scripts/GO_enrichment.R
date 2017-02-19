library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)
library(GOstats)

# initializing biomaRt reference
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)

GO_enrichment <- function(genes, background_genes, ontology) {
  # getting entrez IDs for background genes
  names(background_genes) = "Gene"
  background.entrez = getBM(attributes=c("external_gene_name", "entrezgene"),filters="external_gene_name",values=background_genes, mart=mart)$entrezgene

  # getting entrez IDs for genes of interest (enriched set)
  names(genes) = "Gene"
  entrez_ids = getBM(attributes=c("external_gene_name", "entrezgene"),filters="external_gene_name",values=genes, mart=mart)$entrezgene
  
  # Testing over-enrichment
  # genes given must be unique or GOstats won't work
  params.over <- new('GOHyperGParams',
                   geneIds=unique(entrez_ids),
                   universeGeneIds=unique(background.entrez),
                   ontology=ontology, #BP, CC, MF
                   pvalueCutoff=1,
                   conditional=F,
                   testDirection='over',
                   annotation="org.Hs.eg.db"
  )
  
  hgOver <- hyperGTest(params.over)
  result.over <- summary(hgOver)
  result.over <- result.over %>% 
    mutate(`Expected Count` = ExpCount, 
           `P value (FDR)`= p.adjust(Pvalue,method='fdr'), 
           `P value` = Pvalue, 
           `Odds Ratio` = OddsRatio, 
           `Expected Count` = round(`Expected Count`, 1)) %>% 
    dplyr::select(1, `P value`, `P value (FDR)`, `Odds Ratio`, `Expected Count`, Count, Size, Term)
  result.over$`P value` <- format(result.over$`P value`, digits=3)
  result.over$`P value (FDR)` <- format(result.over$`P value (FDR)`, digits=3)
  result.over$`Odds Ratio` <- round(result.over$`Odds Ratio`, digits=2)
  colnames(result.over)[1] = paste('GO', ontology, 'ID',sep=' ')
  return(result.over)
}
