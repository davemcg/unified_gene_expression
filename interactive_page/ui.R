# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)
load('~/git/unified_gene_expression/data/tx_genes.Rdata')
load('~/git/unified_gene_expression/interactive_page/metaData.Rdata')

# Define UI for application that draws a histogram
shinyUI(
  navbarPage('eyeIntegration',
    tabPanel('BoxPlot',
      fluidPage(
       # Application title
       titlePanel("Gene Expression for Human Tissues and Cells"),

       # Show a plot of the generated distribution
       sidebarLayout(
        sidebarPanel(width=3,
          selectInput("Gene","Select Genes:", choices=unique(sort(tx_genes$gene.Name)), 
            selected=c('ABCA4','TYRP1'),multiple=TRUE),
          selectInput("Tissue","Select Tissues:", choices=unique(sort(core_tight$Sub_Tissue)), 
            selected=c(" Whole Blood ",
                       " Pancreas ",
                       " Cells - EBV-transformed lymphocytes ",
                       " Cells - Transformed fibroblasts ",
                       " Liver ",
                       " Lung ",
                       "Cornea",
                       "fetalRetina",
                       "fetalRPE",
                       "RPE",
                       "Retina"),multiple=TRUE)
        ),
      mainPanel(
        h3('Interactive boxplot of pan-human gene expression', align="center"),
        plotOutput("boxPlot")
      )
    )
  )),
    tabPanel('2D Tissue Clustering',plotOutput("tsne"))
      )
    )

