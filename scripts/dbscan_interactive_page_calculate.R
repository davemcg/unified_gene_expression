####### script to calculate dbscan clusters and build labels for interactive pagea

library(dbscan)
load('~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity_2017-02__all_eye.Rdata')



all_tsne_plot_prepped <- data.frame(matrix(ncol = 16, nrow = 1))
colnames(all_tsne_plot_prepped) <- c("X1","X2","sample_accession","perplexity","Cluster","study_accession","study_title","study_abstract","run_accession","sample_attribute","Tissue","Sub_Tissue","Origin","Study","Cluster_Tissues","Label")
for (i in 5:50){
  tsne_plot <- long_tsne_plot %>% filter(perplexity==i)
  # dbscan
  dbscan_cluster <- dbscan(tsne_plot %>% dplyr::select(X1,X2), eps=2)
  tsne_plot$Cluster <- dbscan_cluster$cluster
  
  # create label points for each cluster
  cluster_centers <- tsne_plot %>% 
    left_join(.,core_tight) %>% 
    mutate(Tissue = ifelse(grepl('fibroblasts',Sub_Tissue),'Fibroblasts',Tissue)) %>% 
    group_by(Cluster) %>%
    summarise(C1=mean(X1),C2=mean(X2),Tissue=paste(unique(Tissue),collapse=','))
  
  # samples closest to center in each cluster
  center_samples <- tsne_plot %>% left_join(.,core_tight)  %>%
    left_join(.,cluster_centers, by=c('Cluster')) %>% 
    mutate(Distance = abs((X1-C1)-(X2-C2))) %>% 
    group_by(Cluster) %>% 
    dplyr::slice(which.min(Distance)) %>% 
    .[['sample_accession']]
  
  # cluster stats
  cluster_stats <- tsne_plot %>% left_join(.,core_tight)  %>% 
    mutate(Tissue = ifelse(grepl('fibroblasts',Sub_Tissue),' Fibroblasts',Tissue)) %>% 
    mutate(Cluster = as.factor(Cluster)) %>% 
    group_by(Cluster, Tissue) %>% 
    summarise(Cluster_Counts = n(), rTissue = ifelse(Cluster_Counts <= 3, '__', Tissue)) %>%  #Remove outlier tissues from label
    group_by(Cluster) %>% 
    summarise(Cluster_Tissues = paste(unique(rTissue), collapse=', ')) %>% 
    mutate(Cluster_Tissues = gsub('__, ','',Cluster_Tissues))
  
  
  # set up for ggplot
  tsne_plot_prep <- tsne_plot %>% left_join(.,core_tight)  %>% 
    mutate(Tissue = ifelse(grepl('fibroblasts',Sub_Tissue),' Fibroblasts',Tissue)) %>% 
    mutate(Study = study_accession) %>% 
    mutate(Cluster = as.factor(Cluster)) %>%
    left_join(., cluster_stats, by=c('Cluster')) %>% 
    mutate(Cluster = paste(formatC(Cluster,width=2,format='d',flag='0') ,': ', Tissue,sep='')) %>% 
    mutate(Label = paste('Cluster: ', as.factor(Cluster), '<br>', 'Sub-Tissue: ', Sub_Tissue, '<br>', 
                         'SRA: ', sample_accession, '<br>', 
                         gsub('\\|\\|', '<br>', sample_attribute, '<br>', study_title), sep =''))
  tsne_plot_prep$perplexity <- i
  all_tsne_plot_prepped <- bind_rows(all_tsne_plot_prepped, tsne_plot_prep)
}
all_tsne_plot_prepped <- all_tsne_plot_prepped %>% dplyr::slice(-1) %>% select(X1,X2, Cluster, Tissue, Sub_Tissue, Label, perplexity)
save(all_tsne_plot_prepped, file='~/git/unified_gene_expression/interactive_page/all_tsne_plot_prepped__2017_02.Rdata')

the_plot<-ggplot(tsne_plot_prep,aes(x=X1,y=X2, label=Label)) +
  ggtitle(paste0('Pan tissue t-sne')) +
  geom_point(size=4, alpha=0.2, aes(colour=Cluster)) +
  geom_point(size=1, alpha=0.2)   +
  theme_Publication() + guides(colour=guide_legend(nrow=4), shape=guide_legend(nrow=4))

