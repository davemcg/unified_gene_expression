# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)
load('~/git/unified_gene_expression/data/tx_genes.Rdata')

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Gene Expression for Human Tissues and Cells"),

  # Sidebar with a slider input for the number of bins
  sidebarLayout(
    sidebarPanel(width=0),
    # Show a plot of the generated distribution
    mainPanel(
      h3('Interactive boxplot of pan-human gene expression', align="center"),
      selectInput("Gene","Genes:", choices=tx_genes$gene.Name, 
                  selected='RP1',multiple=TRUE),
      plotOutput("boxPlot",height=700) 
      #h1(''),
      #h3('Distance between each RNA-seq experiment', align="center"),
      #h5('Closer points are more related',align="center"),
      #img(src='tsne_2016-06-29.svg'), 
      #width='90%'
    )
  )
))