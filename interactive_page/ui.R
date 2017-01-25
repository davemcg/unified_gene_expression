# https://bioconductor.org/help/course-materials/2015/CSAMA2015/lab/shiny.html
# http://shiny.rstudio.com/gallery/telephones-by-region.html
# http://davetang.org/muse/2014/01/03/using-shiny/

library(shiny)
library(ggplot2)
library(plotly)
load('~/git/unified_gene_expression/interactive_page/gene_names.Rdata')
load('~/git/unified_gene_expression/interactive_page/tissue_info.RData')
print('UI Start')
print(Sys.time())

# Define UI for application that draws a histogram
shinyUI(
  navbarPage('eyeIntegration',
    tabPanel('Pan-Tissue Expression',
      fluidPage(
        fluidRow(
          column(2,
            radioButtons('plot_type','Visualization:',
                         choices = c('Box Plot','Fold Change'),
                         selected = 'Box Plot'),
            selectInput('Gene','Genes:', 
              choices=unique(sort(gene_names$Gene.Name)), 
              selected=c('ABCA4','TYRP1'),multiple=TRUE),
            selectInput('Tissue','Tissues:', 
              choices=unique(sort(tissue_info$Sub_Tissue)), 
              selected=c(' Whole Blood ', ' Pancreas ',
                       ' Cells - EBV-transformed lymphocytes ',
                       ' Cells - Transformed fibroblasts ',
                       ' Liver ', ' Lung ', 'Cornea',
                       'fetalRetina', 'fetalRPE',
                       'RPE', 'Retina'),multiple=TRUE),
            numericInput('num', 'Number of columns:', 
              value = 2, min = 1)
          ),
        
        column(6,
            conditionalPanel(condition = "input.plot_type == 'Box Plot'",
              plotOutput('boxPlot')
            ),
            conditionalPanel(condition = "input.plot_type == 'Fold Change'",
              selectInput('Bench','Select Reference Tissue(s):', 
                          unique(sort(tissue_info$Sub_Tissue)),multiple = TRUE),
              plotOutput('FC')
            )
        ),
        column(4,
          conditionalPanel(condition = "input.plot_type == 'Fold Change'",
            div(DT::dataTableOutput('basicStats'),style='font-size:75%')),
          conditionalPanel(condition = "input.plot_type == 'Box Plot'",
            div(DT::dataTableOutput('rankStats'),style='font-size:75%'))
        )
      )
    )
  ),
    tabPanel('Eye Plot',
      fluidPage(
        fluidRow(
          column(2,
            selectInput('eyeGene','Genes:', 
              choices=unique(sort(gene_names$Gene.Name)), 
              selected=c('ABCA4','RPE65','TYRP1'),multiple=TRUE),
            numericInput('eyeNum', 'Number of columns:', 
              value = 3, min = 1)
        ),
          column(10,
            mainPanel(
              h3('Interactive scatter plot of eye-tissue gene expression', align='center'),
              plotlyOutput('eyeBoxPlot')
            )
          )
        
      ))),
    tabPanel('2D Tissue Clustering',
      fluidPage(
        fluidRow(plotlyOutput('tsne',height = '800px')),
        fluidRow(numericInput('perplexity','Perplexity (5 - 50):', value=40, min=5, max=50))
      )
    ),
    tabPanel('Data Table',
      fluidPage(
        fluidRow(
          column(2,
            selectInput('table_tissue',
              'Tissue:',
              choices = unique(as.character(tissue_info$Tissue)),
              selected = 'Retina')),
          column(2,
            selectInput('table_gene',
              'Gene:',
              choices = unique(as.character(gene_names$Gene.Name)),
              selected = 'CFH')),
          column(6,
              selectInput('table_columns',
                'Columns: ', 
                multiple = TRUE,
                choices = c('Gene.Name', 'sample_accession', 'value', colnames(core_tight)),
                selected = c('Gene.Name', 'sample_accession', 'value', 'study_title', 'sample_attribute', 'Tissue', 'Sub_Tissue','Origin'))
          )
        ), 
        fluidRow(DT::dataTableOutput('table')
        )
      )
    )
  )
)
print('UI End')
print(Sys.time())
