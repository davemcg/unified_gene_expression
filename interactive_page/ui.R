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
  navbarPage('eyeIntegration', theme='bootstrap.css', 
    tabPanel('Pan-Tissue Expression',
      fluidPage(
        fluidRow(
          column(2,
            img(src='NIH_NEI_Vertical_Logo_Black90.png',align='left'), br(),br(),br(),br(),
            radioButtons('plot_type',strong('Visualization:'),
                         choices = c('Box Plot','Fold Change'),
                         selected = 'Box Plot'),
            selectInput('Gene',strong('Genes:'), 
              choices=unique(sort(gene_names$Gene.Name)), 
              selected=c('ABCA4','TYRP1'),multiple=TRUE),
            selectInput('Tissue',strong('Tissues:'), 
              choices=unique(sort(tissue_info$Sub_Tissue)), 
              selected=c(' Whole Blood ', ' Pancreas ',
                       ' Cells - EBV-transformed lymphocytes ',
                       ' Cells - Transformed fibroblasts ',
                       ' Liver ', ' Lung ', 'Cornea',
                       'Retina - Adult Tissue', 'RPE - Fetal Tissue',
                       'Cornea - Cell Line'),multiple=TRUE),
            numericInput('num', strong('Number of columns:'), 
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
            img(src='NIH_NEI_Vertical_Logo_Black90.png',align='left'), br(),br(),br(),br(),
            selectInput('eyeGene',strong('Genes:'), 
              choices=unique(sort(gene_names$Gene.Name)), 
              selected=c('ABCA4','RPE65','TYRP1'),multiple=TRUE),
            numericInput('eyeNum', strong('Number of columns:'), 
              value = 3, min = 1)
        ),
          column(10,
            mainPanel(
              plotlyOutput('eyeBoxPlot')
            )
          )
        
      ))),
    tabPanel('2D Tissue Clustering',
      fluidPage(
        fluidRow(column(10,
          img(src='NIH_NEI_Vertical_Logo_Black90.png',align='left'), br(),br(),br(),br(),
          plotlyOutput('tsne',height = '600px')),
        fluidRow(column(10,
          numericInput('perplexity','Perplexity (5 - 50):', value=40, min=5, max=50)))
      ))
    ),
    tabPanel('Data Table',
      fluidPage(
        fluidRow(
          column(12,img(src='NIH_NEI_Vertical_Logo_Black90.png',align='left'))), br(),
        fluidRow(
          column(2,
              selectInput('table_tissue',
              strong('Tissue:'),
              choices = unique(as.character(tissue_info$Tissue)),
              selected = 'Retina')),
          column(2,
            selectInput('table_gene',
              strong('Gene:'),
              choices = unique(as.character(gene_names$Gene.Name)),
              selected = 'CFH')),
          column(6,
              checkboxGroupInput('table_columns',
                strong('Columns: '), 
                inline = T,
                choices = c('Gene.Name', 'sample_accession', 'value', 'study_accession', 'study_title', 'study_abstract', 'sample_attribute', 'Tissue', 'Sub_Tissue','Origin'),
                selected = c('Gene.Name', 'sample_accession', 'value', 'study_title', 'sample_attribute', 'Tissue', 'Sub_Tissue','Origin'))
          )
        ), 
        fluidRow(DT::dataTableOutput('table')
        )
      )
    ),
  tabPanel('Information',
           fluidPage(
             fluidRow(column(12,img(src='NIH_NEI_Vertical_Logo_Black90.png',align='left'))), br(),
             fluidRow(column(width = 8, offset = 1, h2('Mission and Brief Methods'))),
             fluidRow(column(width = 8, offset = 1, 
                    'The goal of this project is to provide a high quality source collating all publical available human eye RNA expression 
                    data-sets, allowing researchers to have a single source to find gene expression in eye tissues and to compare their e
                    xpression against the tissues from the ', a(href='http://www.gtexportal.org','GTEx Project'), br(), br(),
                    'The raw sequencing reads were obtained and processed with an identical bioinformatic workflow. Extensive QC was done
                     to remove samples with both lower sequencing quality and unexpected behavior during undirected clustering. Read scoring 
                     was scaled for transcript size, library depth, and quantile normalized by tissue type. Transcript level quantitation 
                     was merged to the gene level to improve power and lower variance. Finally genes with low expression across all tissues 
                     were removed from analysis.')),
             fluidRow(column(width = 8, offset = 1, h2('Attribution'))),
             fluidRow(column(width = 8, offset = 1, 'This project was conceived and implemented by', 
              tags$a(href = "mailto:mcgaugheyd@mail.nih.gov?Subject=eyeIntegration%20Comment", "David McGaughey"), 
              'in the', tags$a(href='https://nei.nih.gov/intramural/ogcsb','OGVFB'), '/',
              tags$a(href='https://nei.nih.gov','NEI'), '/',
              tags$a(href='https://www.nih.gov','NIH'), 
              'The gene network analysis, visualization, and page design was done by ', tags$a(href='mailto:john.bryan@nih.gov', 'John Bryan'))),
             fluidRow(column(width = 8, offset = 1, h2('Source Code'))),
             fluidRow(column(width = 8, offset = 1, h2('Manuscript')))
             )
             )
           )
  )

print('UI End')
print(Sys.time())
