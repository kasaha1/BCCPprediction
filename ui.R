library(tidyverse)
library(data.table)
library(shiny)
library(plotly)
library(heatmaply)
library(shinycssloaders)
library(shiny)
library(kasaBasicFunctions)
library(classpredict)
library(shinythemes)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  theme = shinytheme("journal"),
  # Application title
  # tags$img(src="PICSicon.jpg",width = 100,height = 50),
  titlePanel("Bayesian Compound Covariate Prediction (BCCP)"),
  # titlePanel( div(column(width = 3, tags$img(src = "PICSicon.jpg",width = 80,height = 40)),column(width = 8, h3("PICS100 prediction by Ju-Seog's Lab"))),windowTitle="PICS100" ),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      tags$img(
        src = "icon.png",
        width = 220,
        height = 100
      ),
      br(),
      br(),
      fileInput(
        'TrainingGeneExprfile',
        h4('The mRNA expression file for training(txt/csv)'),
        accept = c('text/csv',
                   'text/comma-separated-values,text/plain',
                   '.csv')
      ),
      fileInput(
        'ClassOfTraining',
        h4('The file for class of the training(txt/csv)'),
        accept = c('text/csv',
                   'text/comma-separated-values,text/plain',
                   '.csv')
      ),
      fileInput(
        'PredictedGeneExprfile',
        h4('The mRNA expression file for prediction(txt/csv)'),
        accept = c('text/csv',
                   'text/comma-separated-values,text/plain',
                   '.csv')
      ),
      radioButtons(
        'standardizationType',
        'Standardization',
        c(
          'Non-standardization' = 'NoneS',
          # 'Mendian-centering only' = 'medianCenter',
          'Median-centering and dividing by SD' = 'devidedBySD',
          'Robust Z score' = 'RobustZ',
          'Z score'='Zscore'
        )
      ),
      actionButton(inputId = "doPrediction", label = "Prediction", class = "btn-primary"),
      br(),
      br(),
      downloadButton('downloadResults', 'Download result table')
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        type = "tabs",
        
        tabPanel(
          "How to use",
          HTML(
            "&nbsp; <p>Hi. This is the prediction tool for BCCP using mRNA expression data.</p><p>Just upload your dataset. And press the prediction button. That's all. You can download example dataset from"
          ),
          tags$a(href = "exampleSet.zip", " here"),
          HTML("or the below button."),
          br(),
          HTML("</p><p>&nbsp;</p><p><b>Step 1. Prepare of dataset.</b></p>"),
          HTML("<b>Step 1.1 Training dataset (txt/csv file).</b>"),
          img(
            src = "1signature.png",
            width = 400,
            height = 150
          ),
          br(),
          HTML("</p><p>&nbsp;</p><p><b>Step 1.2 Class(group) of training dataset (txt/csv).</b></p>"),
          img(
            src = "2class.png",
            width = 150,
            height = 150
          ),
          br(),
          HTML("</p><p>&nbsp;</p><p><b>Step 1.3 Test(for prediction) dataset (txt/csv).</b></p>"),
          img(
            src = "3test.png",
            width = 600,
            height = 200
          ),
          
          HTML(
            "<p> &nbsp;</p><p> The first line contains the labels Name(<em>HUGO Gene Nomenclature</em>) followed by the identifiers for each sample in the dataset.The dataset is the gene-level transcription estimates, as in log2(x+1) transformed normalized count.&nbsp; </br>* The symbols of test dataset are automatically matched to training dataset symbols .&nbsp;</p><p>&nbsp;</p><p><b>Step 2. Standardization. </b> &nbsp;</p><p> Select the data standardization method. &nbsp;</p><p><b>Step 3. Prediction.</b> &nbsp;</p><p> Press the prediction button. &nbsp;</p><p><b>Step 4. Check out the results.</b> &nbsp;</p><p>After analysis, You can find the results at the result tab. The results of dataset could be downloaded using the download button.</p>"
          )
          
        )
      
      )
    )
  )
))