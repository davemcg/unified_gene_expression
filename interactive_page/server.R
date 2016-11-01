# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)
library(tidyverse)
source('~/git/scripts/theme_Publication.R')
load('~/git/unified_gene_expression/data/lengthScaledTPM_processed.Rdata')
load('~/git/unified_gene_expression/interactive_page/metaData.Rdata')
lengthScaledTPM_qsmooth_highExp_remove_lowGenes <- data.frame(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
lengthScaledTPM_qsmooth_highExp_remove_lowGenes$Gene.Name <- row.names(lengthScaledTPM_qsmooth_highExp_remove_lowGenes)
shiny_data <- lengthScaledTPM_qsmooth_highExp_remove_lowGenes
core_tight$sample_accession<-gsub('E-MTAB-','E.MTAB.',core_tight$sample_accession)
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$boxPlot <- renderPlot({
    gene <- input$Gene
    data <- shiny_data %>% filter(Gene.Name == gene) %>% data.frame()
    
    plot_data <- t(data) %>% data.frame(stringsAsFactors = F) %>% 
      rownames_to_column(var='sample_accession') %>% left_join(.,core_tight)
    # draw the histogram with the specified number of bins
    colnames(plot_data)[2]<-'lsTPM'
    p<-ggplot(data=data.frame(plot_data),aes(x=Sub_Tissue,y=log2(as.numeric(lsTPM)+1),colour=Tissue)) + 
      geom_jitter(size=2) + geom_boxplot(alpha=0.5) + ggtitle(gene) + facet_wrap(~gene) +
      theme_Publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") 
    
    print(p)
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
  }, height=700)
})