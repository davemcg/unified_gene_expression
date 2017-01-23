library(broom)
library(tidyverse)

bench <- c('RPE','Retina')
# calculates log2 Fold Change against a user-defined 'bench' (reference)
plot_data %>% 
  mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>% 
  group_by(Sub_Tissue) %>% 
  summarise(log2FC=mean(log2(value+1)) - mean(BenchValue))

# does t.test against a user-defined reference
plot_data %>% 
  group_by(Sub_Tissue) %>%
  do(tidy(t.test(Sub_Tissue ~ subset(plot_data, Sub_Tissue %in% bench)$value, data=.)))

pairwise.t.test(plot_data$value,plot_data$Sub_Tissue,p.adjust.method = 'none')