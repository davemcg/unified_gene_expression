# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)
library(dplyr)

load('gene_plotting_data_shiny.Rdata')
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
    data <- metaData_lSTPM %>% filter(Gene %in% gene)  
    # draw the histogram with the specified number of bins
    p<-ggplot(data=data,aes(x=Tissue,y=log2(lsTPM+1),colour=EyeMarker)) + 
      geom_jitter(aes(shape=ProjectID),size=2) + geom_boxplot(alpha=0.5) + facet_wrap(~Gene,ncol=1) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ") +
      scale_shape_manual(values=1:8)
    
    print(p)
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
  }, height=700)
})