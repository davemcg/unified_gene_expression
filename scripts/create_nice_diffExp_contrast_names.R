# little script to create list of all vs all DiffExp contrast names for web app
load('~/git/unified_gene_expression/data/limma_voom_DE_all_by_all.Rdata')
library(tidyverse)



all_by_all_names <- colnames(efit_all) %>% data.frame()
colnames(all_by_all_names) <- 'Comparisons'

nice_names <- all_by_all_names %>% data.frame() %>% 
  separate(Comparisons, c('Comparison','Base'), sep='_vs_') %>% 
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
  mutate(Comparisons = paste(Comparison, Base, sep=' vs ')) %>% 
  mutate(Comparisons = gsub('  ',' ',Comparisons))

de_comparison_contrast_names <- sort(as.character(all_by_all_names$Comparisons))
de_comparison_contrast_names <- set_names(de_comparison_contrast_names, as.character(nice_names$Comparisons))


save(de_comparison_contrast_names, file='data/de_comparison_name_list.Rdata')
