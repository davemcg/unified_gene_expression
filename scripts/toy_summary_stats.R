library(broom)
library(tidyverse)

gene <- c('ABCA4')
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


bench <- c('RPE','Retina')
# calculates log2 Fold Change against a user-defined 'bench' (reference)
plot_data %>% 
  mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>% 
  group_by(Sub_Tissue) %>% 
  summarise(log2FC=mean(log2(value+1)) - mean(BenchValue))

# does t.test against a user-defined reference
tissue_subset <- plot_data %>% filter(Sub_Tissue %in% bench)
plot_data %>% 
  group_by(Sub_Tissue) %>%
  do(tidy(t.test(.$value, tissue_subset$value)))

# calculates mean, standard deviation, variance                                                                                                                                                                 
