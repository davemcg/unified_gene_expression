library(ggvis)
library(tidyverse)
library(Cairo)

load('~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity.Rdata')
source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')

dynamic_plot<- long_tsne_plot%>% left_join(.,core_info) %>% filter(perplexity==40, Tissue %in% keepers)
dynamic_plot$id <- 1:nrow(dynamic_plot)

all_values <- function(x) {
  if(is.null(x)) return(NULL)
  row <- dynamic_plot[dynamic_plot$id == x$id, c('Tissue','Sub_Tissue','Origin','study_title','sample_attribute','sample_accession', 'study_accession')]
  paste0(names(row), ": ", format(row), collapse = "<br />")
}

dynamic_plot%>% ggvis(~X1,~X2,fill = ~ Tissue,shape = ~Tissue, key := ~id) %>% 
  layer_points() %>% add_tooltip(all_values, c('hover','click'))



p <- 
  ggplot(dynamic_plot,aes(x=X1,y=X2, shape=Sub_Tissue, colour=Tissue, label=sample_attribute)) + geom_point() + scale_shape_manual(values=c(0:24,35:75))
plot_ly(p)
