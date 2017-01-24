library(broom)
library(tidyverse)
library(plotly)
library(shiny)
library(ggplot2)
library(tidyverse)

source('~/git/scripts/theme_Publication.R')
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed.Rdata')
load('~/git/unified_gene_expression/interactive_page/metaData.Rdata')
load('~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity.Rdata')
lengthScaledTPM_qsmooth_highExp_remove_lowGenes <- data.frame(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
lengthScaledTPM_qsmooth_highExp_remove_lowGenes$Gene.Name <- row.names(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
shiny_data <- lengthScaledTPM_qsmooth_highExp_remove_lowGenes
core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
long_tsne_plot$sample_accession<-gsub('E-MTAB-','E.MTAB.',long_tsne_plot$sample_accession)


gene <- c('RABGAP1','TYR')
tissue<- c(" Whole Blood ",
" Pancreas ",
" Cells - EBV-transformed lymphocytes ",
" Cells - Transformed fibroblasts ",
" Liver ",
" Lung ",
"Cornea",
"fetalRetina",
"fetalRPE",
"RPE",
"Retina")
col_num <- 1
plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>%
gather(sample_accession, value, -Gene.Name) %>%
left_join(.,core_tight)
plot_data <- plot_data %>% filter(Sub_Tissue %in% tissue)


bench <- c('RPE')
# calculates log2 Fold Change against a user-defined 'bench' (reference)
plot_data %>% 
  group_by(Gene.Name) %>% 
  mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>% 
  group_by(Gene.Name, Sub_Tissue) %>% 
  summarise(log2FC=mean(log2(value+1)) - mean(BenchValue))

# does t.test against a user-defined reference
tissue_subset <- plot_data %>% filter(Sub_Tissue %in% bench)
plot_data %>% 
  group_by(Gene.Name, Sub_Tissue) %>%
  do(tidy(t.test(.$value, tissue_subset$value)))

# calculates mean, standard deviation, variance                                                                                                                                                                 
plot_data %>% 
  group_by(Gene.Name, Sub_Tissue) %>% 
  summarise(mean=mean(value),median=median(value),variance=var(value),sd=sd(value))