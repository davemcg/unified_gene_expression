# creates plot_ly heatmap for web app for GO enrichment heatmap

load('~/git/unified_gene_expression/data/go_enrichment_all_vs_all.Rdata')
library(superheat)
library(tidyverse)

# select top go ids for up
up_go_ids <- all_vs_all_go %>% 
  filter(Ontology=='BP', grepl('Body', Set),Test=='Up') %>% 
  arrange(as.numeric(`P value (FDR)`)) %>% 
  filter(as.numeric(`P value (FDR)`)<0.01) %>% 
  dplyr::select(`GO ID`) %>% unique() %>%
  .[['GO ID']]

.[['GO ID']]
# and down
down_go_ids <- all_vs_all_go %>% 
  filter(Ontology=='BP', grepl('Body', Set),Test=='Down') %>% 
  arrange(as.numeric(`P value (FDR)`)) %>% 
  filter(as.numeric(`P value (FDR)`)<0.01) %>% 
  dplyr::select(`GO ID`) %>% unique() %>% 
  .[['GO ID']]

go_ids <- c(up_go_ids, down_go_ids)
wide_data <- all_vs_all_go %>% 
  filter(Ontology=='BP', grepl('Body', Set), `GO ID` %in% go_ids) %>% 
  separate(Set, c('Comparison','Base'), sep='_vs_') %>% 
  mutate(Comparison=gsub('_',' (',Comparison), 
         Comparison=gsub('Adult.Tissue','adult',Comparison),
         Comparison=gsub('Fetal.Tissue','fetal',Comparison),
         Comparison=gsub('Stem.Cell.Line','stem cell',Comparison),
         Comparison=gsub('Cell.Line','immortalized cell',Comparison),
         Comparison=gsub('$',')', Comparison), 
         Comparison=gsub('cell',' cell', Comparison)) %>% 
  mutate(Base=gsub('_',' (',Base), 
         Base=gsub('Adult.Tissue','adult',Base),
         Base=gsub('Fetal.Tissue','fetal',Base),
         Base=gsub('Stem.Cell.Line','stem cell',Base),
         Base=gsub('Cell.Line','immortalized cell',Base),
         Base=gsub('$',')', Base), 
         Base=gsub('cell',' cell', Base), 
         Base=gsub('Body \\(Tissue\\)','Body (adult)',Base)) %>% 
  mutate(Comparisons = ifelse(Test=='Up', paste(Comparison, Base, sep=' > '), paste(Comparison, Base, sep=' < '))) %>% 
  mutate(`-log10(FDR)` = -log(as.numeric(`P value (FDR)`), base=10)) %>% 
  mutate(GO = paste(`GO ID`, Term, sep=' ')) %>% 
  dplyr::select(Comparisons, GO, `-log10(FDR)`) %>% 
  spread(Comparisons, `-log10(FDR)`, fill=1)
head(wide_data)
row.names(wide_data) <- wide_data$GO

# manually cluster the rows
m <- wide_data[,2:ncol(wide_data)]
go_order <-hclust(dist(m))$order
# manually cluster the cols
sample_order <-hclust(dist(t(m)))$order

# now actually reorder
m <- m[go_order, sample_order]
# get colnames and rownames as ordered factors
order_columns <- factor(colnames(m), levels = colnames(m))
order_rows <-  factor(row.names(m), levels =row.names(m))


save(m, file='~/git/human_eyeIntegration_web_app/www/go_heatmap.Rdata')


        
# example plot_ly below
options <- list(
  l = 500,
  b = 300,
  tickangle = 90,
  type = 'category')
options2 <- list(
  l = 500,
  b = 300,
  type = 'category')
# example plot
plot_ly(z=as.matrix(m),
                       y = order_rows,
                       x = order_columns, width=1500, height=20000,
                       type='heatmap', zmax=20, zauto=FALSE) %>%
                 layout(margin=options, xaxis=options, yaxis=options2)
