# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)

load('metaData_with_lengthScaledTPM.Rdata')
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Expression that generates a histogram. The expression is
  # wrapped in a call to renderPlot to indicate that:
  #
  #  1) It is "reactive" and therefore should re-execute automatically
  #     when inputs change
  #  2) Its output type is a plot
  
  output$distPlot <- renderPlot({
    
    gene <- input$Gene
    data <- subset(metaData_lSTPM,Gene==gene)  
    # draw the histogram with the specified number of bins
    p<-ggplot(data=data,aes(x=Library_Name_s,y=log2(value+1),colour=tissue_s)) + 
      geom_point() + facet_grid(~SRA_Study_s,space='free',scales='free') + 
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ylab("Gene Expression | log2(lengthScaledTPM+1) ")
    
    print(p)
    #hist(x, breaks = bins, col = 'darkgray', border = 'white')
  })
})