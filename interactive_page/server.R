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
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed_2017_02.Rdata')
#load('~/git/unified_gene_expression/interactive_page/metaData.Rdata')
source('~/git/unified_gene_expression/scripts/parse_sample_attribute.R')

load('~/git/unified_gene_expression/interactive_page/all_tsne_plot_prepped.Rdata')
load('~/git/unified_gene_expression/data/mean_rank_decile.Rdata')
lengthScaledTPM_processed <- data.frame(lengthScaledTPM_processed)
lengthScaledTPM_processed$Gene.Name <- row.names(lengthScaledTPM_processed)
shiny_data <- lengthScaledTPM_processed
core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
core_tight$Sub_Tissue <- gsub('_',' - ',core_tight$Sub_Tissue)
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
      theme_Publication() + theme(axis.text.x = element_text(angle = 90)) +
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
      theme_Publication() + theme(axis.text.x = element_text(angle = 90)) +
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
    eye_p<-ggplot(data=data.frame(eye_plot_data),aes(x=Tissue,y=log2(value+1),colour=Tissue, label=Info, shape=Origin)) + 
      geom_jitter(size=2, alpha=0.7) + xlab('')  +  facet_wrap(~Gene.Name, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") +
      ggtitle('Interactive scatter plot of eye-tissue gene expression') +
      theme(text = element_text(size=11), axis.title.y=element_text(margin=margin(0,0,0,100))) +
      scale_colour_manual(values=c('#ff6a6e', '#00c1c1','#6aad27'))
    ggplotly(eye_p, width = 800, height=500) %>% layout(margin=list(b=150))
  })

  
  ##############
  # tsne
  ##############
  
  output$tsne <- renderPlotly({
    perplexity_level <- input$perplexity
    tsne_plot<- all_tsne_plot_prepped %>% filter(perplexity==perplexity_level)
    
    p <- ggplot(tsne_plot) +
        ggtitle('Pan tissue t-sne') +
        geom_point(size=4, alpha=0.2, aes(x=X1,y=X2,colour=Cluster, label=Label)) +
        geom_point(data=tsne_plot %>% select(X1,X2), size=0.5, alpha=0.3, colour='black', aes(x=X1,y=X2)) +
        theme_Publication()
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