# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

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
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$boxPlot <- renderPlot({
    gene <- input$Gene
    tissue <- input$Tissue
    col_num <- input$num
    plot_data <- shiny_data %>% filter(Gene.Name %in% gene) %>% 
      gather(sample_accession, value, -Gene.Name) %>% 
      left_join(.,core_tight)
    plot_data <- plot_data %>% filter(Sub_Tissue %in% tissue)
    # draw the histogram with the specified number of bins
    p<-ggplot(data=data.frame(plot_data),aes(x=Sub_Tissue,y=log2(value+1),colour=Tissue)) + 
      geom_jitter(size=2) + geom_boxplot(alpha=0.5) + xlab('') + facet_wrap(~Gene.Name, ncol=col_num) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 75, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") 
    p
  }, height=function(){(500*length(input$Gene))/(input$num)})
  
  output$tsne <- renderPlotly({
    tsne_plot<- long_tsne_plot%>% left_join(.,core_tight) %>% filter(perplexity==40)
    p <- tsne_plot %>% ggplot(aes(X1,X2,colour = Tissue, shape = Tissue, label=paste(sample_accession, gsub('\\|\\|', '<br>', sample_attribute), sep='<br>')))  +
      geom_point(size=4) + scale_shape_manual(values=c(0:24,35:45)) +
      theme_bw()
    ggplotly(p)
  })
})