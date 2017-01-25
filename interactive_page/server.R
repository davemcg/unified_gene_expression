# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

# load stuff for server
library(plotly)
library(shiny)
library(ggplot2)
library(tidyverse)
library(broom)
library(DT)
source('~/git/scripts/theme_Publication.R')
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed.Rdata')
load('~/git/unified_gene_expression/interactive_page/metaData.Rdata')
load('~/git/unified_gene_expression/data/tsne_plotting_5_50_perplexity.Rdata')
load('~/git/unified_gene_expression/data/mean_rank_decile.Rdata')
lengthScaledTPM_qsmooth_highExp_remove_lowGenes <- data.frame(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
lengthScaledTPM_qsmooth_highExp_remove_lowGenes$Gene.Name <- row.names(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
shiny_data <- lengthScaledTPM_qsmooth_highExp_remove_lowGenes
core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
long_tsne_plot$sample_accession<-gsub('E-MTAB-','E.MTAB.',long_tsne_plot$sample_accession)
# responsive stuff!
shinyServer(function(input, output, session) {
  # tissues
  observe({
    selected_tissue <- input$Tissue
    updateSelectInput(session, "Bench",
      label='Select Reference Tissue(s): ',
      choices = selected_tissue,
      selected = selected_tissue)
  })
  
  #########
  # pan - tissue boxplot
  #########
  output$boxPlot <- renderPlot({
    gene <- input$Gene
    tissue <- input$Tissue
    col_num <- input$num
    plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>% 
      gather(sample_accession, value, -Gene.Name) %>% 
      left_join(.,core_tight)
    plot_data <- plot_data %>% filter(Sub_Tissue %in% tissue)
    p<-ggplot(data=data.frame(plot_data),aes(x=Sub_Tissue,y=log2(value+1),colour=Tissue)) + 
      geom_jitter(size=2) + geom_boxplot(alpha=0.5) + xlab('') + facet_wrap(~Gene.Name, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      ggtitle('Box Plot of Pan-Human Gene Expression') +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") 
    p
  }, height=function(){(500*length(input$Gene))/min(input$num,length(input$Gene))})
  
  ############
  # fold change
  ############
  output$FC <- renderPlot({
    gene <- input$Gene
    tissue <- input$Tissue
    col_num <- input$num
    bench <- input$Bench
    plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>% 
      gather(sample_accession, value, -Gene.Name) %>% 
      left_join(.,core_tight)
    plot_data <- plot_data %>% 
      filter(Sub_Tissue %in% tissue) %>% 
      group_by(Gene.Name) %>% 
      mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>% 
      group_by(Gene.Name, Sub_Tissue) %>% 
      summarise(log2FC=mean(log2(value+1)) - mean(BenchValue))
    p<-ggplot(data=data.frame(plot_data),aes(x=Sub_Tissue,y=log2FC,fill=Sub_Tissue)) + 
      geom_bar(stat = 'identity') + xlab('') + facet_wrap(~Gene.Name, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      geom_hline(aes(yintercept=0,colour='Red')) + 
      ggtitle('Fold Change (log2) of pan-human gene expression') +
      ylab("log2 Fold Change of Gene Expression") 
    p
  }, height=function(){(500*length(input$Gene))/min(input$num,length(input$Gene))})
  
  ##########
  # boxplot stats
  ##########
  output$rankStats <- DT::renderDataTable(({
    gene <- input$Gene
    tissue <- input$Tissue
    mean_rank_decile %>% 
      filter(Gene.Name %in% gene, Sub_Tissue %in% tissue) %>% 
      mutate(`Gene Name` = Gene.Name, Tissue = `Sub_Tissue`) %>% 
      ungroup() %>% 
      select(`Gene Name`, Tissue, Rank, Decile) %>% 
      arrange(`Gene Name`, Tissue) %>% 
      DT::datatable(options = list(pageLength = 20)) 
  }))
  
  ########
  # FC table stats
  ########
  output$basicStats <- DT::renderDataTable({
    gene <- input$Gene
    tissue <- input$Tissue
    bench <- input$Bench
    plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>% 
      gather(sample_accession, value, -Gene.Name) %>% 
      left_join(.,core_tight)
    base_stats <- plot_data %>% 
      filter(Sub_Tissue %in% tissue) %>% 
      group_by(Gene.Name) %>% 
      mutate(Bench=ifelse(Sub_Tissue %in% bench, 1, 0), BenchValue=mean(log2(value[Bench==1]+1))) %>% 
      group_by(Gene.Name, Sub_Tissue) %>% 
      summarise(log2DeltaFC=mean(log2(value+1)) - mean(BenchValue), mean=mean(log2(value+1)))
  
    # does t.test against a user-defined reference
    # corrects for number of tests
    tissue_subset <- plot_data %>% filter(Sub_Tissue %in% bench)
    pvals <- plot_data %>% 
      group_by(Sub_Tissue, Gene.Name) %>%
      do(tidy(t.test(.$value, tissue_subset$value))) %>%
      # multiple test correction
      mutate(`t test p` = signif(min(1,p.value * length(unique(plot_data$Sub_Tissue)))),3) %>%
      select(Gene.Name, Sub_Tissue, `t test p`)
    stat_join <- left_join(base_stats, pvals) %>% 
      mutate(`Gene Name` = Gene.Name, Tissue = `Sub_Tissue`, `log2 Fold Change` = log2DeltaFC, `Fold Change` = 2^log2DeltaFC, `Mean Expression` = mean) %>% 
      ungroup() %>% 
      select(`Gene Name`, Tissue, `log2 Fold Change`, `Fold Change`, `Mean Expression`, `t test p`) %>% 
      arrange(`Gene Name`, Tissue) %>% 
      DT::datatable(options = list(pageLength = 20))  %>% 
      DT::formatRound(c('log2 Fold Change','Mean Expression'), digits=2) %>% 
      DT::formatRound('Fold Change', digits=6)
    stat_join})
  
  #########
  # eye-only boxplot
  #########
  output$eyeBoxPlot <- renderPlotly({
    gene <- input$eyeGene
    col_num <- input$eyeNum
    eye_plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>% 
      gather(sample_accession, value, -Gene.Name) %>% 
      left_join(.,core_tight) %>% 
      mutate(Info = paste('<br>','Sub-Tissue: ', Sub_Tissue, '<br>', 'SRA: ', sample_accession, '<br>', gsub('\\|\\|', '<br>', sample_attribute), sep =''))
    eye_plot_data <- eye_plot_data %>% filter(Tissue %in% c('Retina','RPE','Cornea'))
    eye_p<-ggplot(data=data.frame(eye_plot_data),aes(x=Sub_Tissue,y=log2(value+1),colour=Tissue, label=Info)) + 
      geom_jitter(size=1) + xlab('')  +  facet_wrap(~Gene.Name, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") +
      theme(text = element_text(size=12), axis.title.y=element_text(margin=margin(0,40,0,0)))
    ggplotly(eye_p, width = 800, height=600) %>% layout(margin=list(b=100))
  })

  
  ##############
  # tsne
  ##############
  
  output$tsne <- renderPlotly({
    perplexity_level <- input$perplexity
    tsne_plot<- long_tsne_plot%>% left_join(.,core_tight) %>% filter(perplexity==perplexity_level)
    
    p <- tsne_plot %>% 
      mutate(Info = paste('<br>','Sub-Tissue: ', Sub_Tissue, '<br>', 'SRA: ', sample_accession, '<br>', gsub('\\|\\|', '<br>', sample_attribute), sep ='')) %>% 
      ggplot(aes(X1,X2,colour = Tissue, shape = Tissue, label=Info ))  +
      geom_point(size=4) + scale_shape_manual(values=c(0:24,35:50)) +
      theme_bw()
    ggplotly(p)
  })  
  ##########
  # data table
  ##########
  output$table = DT::renderDataTable({
      shiny_data %>% filter(Gene.Name == input$table_gene) %>% 
        gather(sample_accession, value, -Gene.Name) %>% 
        left_join(.,core_tight) %>% filter(Tissue == input$table_tissue) %>% 
        select(one_of(input$table_columns)) %>% 
        DT::datatable() %>% DT::formatRound(c('value'), digits=2)
  })
})