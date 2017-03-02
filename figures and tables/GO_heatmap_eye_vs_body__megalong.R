#mega long heatmap code 
# too slow to have in rmarkdown

load('../data/go_enrichment_all_vs_all.Rdata')
library(superheat)
library(tidyverse)

# select top go ids for up
up_go_ids <- all_vs_all_go %>% 
  filter(Ontology=='BP', grepl('Body', Set),Test=='Up') %>% 
  arrange(as.numeric(`P value (FDR)`)) %>% 
  dplyr::select(`GO ID`) %>% unique() %>%
  head(n=250) %>% 
  .[['GO ID']]
# and down
down_go_ids <- all_vs_all_go %>% 
  filter(Ontology=='BP', grepl('Body', Set),Test=='Down') %>% 
  arrange(as.numeric(`P value (FDR)`)) %>% 
  dplyr::select(`GO ID`) %>% unique() %>% 
  head(n=250) %>%
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
  mutate(`-log10(FDR)` = -log10(as.numeric(`P value (FDR)`))) %>% 
  mutate(GO = paste(`GO ID`, Term, sep=' ')) %>% 
  dplyr::select(Comparisons, GO, `-log10(FDR)`) %>% 
  spread(Comparisons, `-log10(FDR)`, fill=1)
head(wide_data)
row.names(wide_data) <- wide_data$GO
superheat((wide_data[,2:ncol(wide_data)]),
          #heat.pal = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
          #heat.pal=c("#984ea3","#47039FFF","#7301A8FF","#9C179EFF","#BD3786FF","#D8576BFF","#ED7953FF","#FA9E3BFF","#FDC926FF","#F0F921FF",'white'),
          #heat.pal = viridis(n=10, option = 'plasma'),
          #heat.col.scheme = 
          pretty.order.cols = T, force.left.label = T,
          grid.hline.col = 'white', grid.vline.col = 'white',
          pretty.order.rows = T, force.grid.hline = T,
          scale = F, 
          bottom.label.text.angle = 90,
          left.label.size = 0.85,
          bottom.label.size = 0.025, 
          col.dendrogram = F,
          left.label.text.alignment = "left",
          bottom.label.text.alignment = "right")

# save 1500 pixels wide and 20000 pixels high