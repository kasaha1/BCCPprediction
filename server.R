

library(tidyverse)
library(data.table)
library(shiny)
library(plotly)
library(heatmaply)
library(shinycssloaders)
library(shiny)
library(kasaBasicFunctions)
library(classpredict)
library(zip)

library(BiocManager)
options(repos = BiocManager::repositories())
library(impute)
library(ROC)

options(shiny.maxRequestSize=300*1024^2)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
 
  trainingDataset_r <- reactive({
    
    inFile <- input$TrainingGeneExprfile
    req(inFile)
    f <- fread(inFile$datapath) %>% as.data.frame()
    colnames(f)[1] <- "gene"
    return(f)
  })
  trainingClass_r <- reactive({
    
    inFile <- input$ClassOfTraining
    
    req(inFile)
    f <- fread(inFile$datapath) %>% as.data.frame()
    colnames(f) <- c("sample","class")
    return(f)
  })
  predictDataset_r <- reactive({
    inFile <- input$PredictedGeneExprfile
    req(inFile)
    f <- fread(inFile$datapath) %>% as.data.frame()
    colnames(f)[1] <- "gene"
   
    return(f)
  })

    DoPrediction <- observeEvent(input$doPrediction, {
   
    # duplication removal
    trainingDataset <-trainingDataset_r() %>% kasa.duplicationRemovalBySD()
    colnames(trainingDataset)[1] <- "gene"
    predictDataset <- predictDataset_r() %>% kasa.duplicationRemovalBySD()
    colnames(predictDataset)[1] <- "gene"
    
 
    # sample matching
    
    tmp <- trainingDataset %>% kasa.transposeMatrix()
    tmp.1 <- kasa.matchingRow(dataframeX = tmp,dataframeY = trainingClass_r(), keycolX = "sample",keycolY = "sample")
    trainingDataset <- tmp.1$dataframeX %>% kasa.transposeMatrix()
    colnames(trainingDataset)[1] <- "gene"
    trainingClass <- tmp.1$dataframeY
    colnames(trainingClass) <- c("sample","class")
    
    ## unmatched samples
    matchedSamples <- trainingClass$sample %>% as.vector()
    unmatchedSamples.t <- tmp$sample[!(tmp$sample %in% matchedSamples)]
    unmatchedSamples.c <- trainingClass$sample[!(trainingClass$sample%in% matchedSamples)]
    
    
    # gene matching
    tmp <- kasa.matchingRow(dataframeX = trainingDataset,dataframeY = predictDataset,keycolX = "gene",keycolY = "gene")
    trainingDataset <- tmp$dataframeX
    predictDataset <- tmp$dataframeY
    colnames(trainingDataset)[1] <- "gene"
    colnames(predictDataset)[1] <- "gene"
    
    
    ## unmatched gene
    matchedgene <- trainingDataset$gene %>% as.vector()
    unmatchedGene.t <- trainingDataset_r()$gene[!(trainingDataset_r()$gene %in% matchedgene)]
    unmatchedGene.p <- predictDataset_r()$gene[!(predictDataset_r()$gene %in% matchedgene)]
    
   
    # Standardization by STD_method ---- "none", "STD", "Z_Score", "Robust_Z_score"
    switch (
      input$standardizationType,
      devidedBySD = {
        trainingDataset <-
          trainingDataset %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
        predictDataset <-
          predictDataset %>% kasa.geneMedianCentering() %>% kasa.geneStandardization()
      },
      Zscore = {
        trainingDataset <- trainingDataset %>% kasa.geneZscoring()
        predictDataset <- predictDataset %>% kasa.geneZscoring()
      },
      RobustZ = {
        trainingDataset <- trainingDataset %>% kasa.geneRobustZscoring()
        predictDataset <- predictDataset %>% kasa.geneRobustZscoring()
      }
    )
    
    
    # data preparing for class prediction ----
    geneId_x_ <- trainingDataset$gene %>% as.data.frame()
    colnames(geneId_x_) <- c("UniqueID")
    x_x_ <- trainingDataset[-1]
    filter_x_ <- rep(1,nrow(trainingDataset))
    expdesign_x_ <- trainingClass
    
    exprTrain_x_ <- x_x_
    exprTest_x_ <- predictDataset[-1]
    
    
    k <- count(expdesign_x_,class)
    class_t.l1 <- k$n[1] %>% as.numeric()
    class_t.l2 <- k$n[2] %>% as.numeric()
    
    # projectPath <-"OutputBRBarray"
    projectPath <- paste0(tempdir(),"/OutputBRBarray")
    outputName <- "classPrediction"
    generateHTML <- TRUE
    prevalence <- c(class_t.l1/(class_t.l1+class_t.l2),class_t.l2/(class_t.l1+class_t.l2))
    names(prevalence) <- c(k$class[1], k$class[2])
    geneId_m <- geneId_x_[c("UniqueID")]
    cls_t <- expdesign_x_$class %>% as.vector()
    resList <- classPredict(exprTrain = exprTrain_x_, exprTest = exprTest_x_, isPaired = FALSE, 
                            pairVar.train = NULL, pairVar.test = NULL, geneId_m,
                            cls = cls_t,
                            pmethod = c("ccp", "bcc", "dlda", "knn", "nc", "svm"), 
                            geneSelect = "igenes.univAlpha",
                            univAlpha = 0.001, univMcr = 0, foldDiff = 0, rvm = TRUE, filter = filter_x_, 
                            ngenePairs = 25, nfrvm = 10, cvMethod = 1, kfoldValue = 10, bccPrior = 1, 
                            bccThresh = 0.8, nperm = 0, svmCost = 1, svmWeight =1, fixseed = 1, 
                            prevalence = prevalence, projectPath = projectPath, 
                            outputName = outputName, generateHTML = generateHTML)
    if (generateHTML)
      browseURL(file.path(projectPath, "Output", outputName,
                          paste0(outputName, ".html")))
    print(projectPath)
    
    output$downloadResults <- downloadHandler(
      filename = 'AllReports.zip',
      content = function(fname) {
        
        setwd(projectPath)
        print(projectPath)
        
        fs <- c("unmatchedSamples_training.txt", "unmatchedSamples_class.txt", "unmatchedGene_training.txt","unmatchedGene_prediction.txt","1_trainingDataset.txt","2_trainingClass.txt","3_predictDataset.txt","4_PredictionResults.txt","5_probability_BCCP.txt")
        write.table(unmatchedSamples.t,file = "unmatchedSamples_training.txt",quote=F,row.names = F)
        write.table(unmatchedSamples.c,file = "unmatchedSamples_class.txt",quote=F,row.names = F)
        write.table(unmatchedGene.t,file = "unmatchedGene_training.txt",quote=F,row.names = F)
        write.table(unmatchedGene.p,file = "unmatchedGene_prediction.txt",quote=F,row.names = F)
        
        write_delim(x = trainingDataset,file = "1_trainingDataset.txt",delim = "\t")
        write_delim(x = predictDataset,file = "3_predictDataset.txt",delim = "\t")
        write_delim(x = trainingClass,file = "2_trainingClass.txt",delim = "\t")
        write_delim(x=resList$predNewSamples,file = "4_PredictionResults.txt",delim = "\t")
        write_delim(x=resList$probNew,file = "5_probability_BCCP.txt",delim = "\t")
        
        zip(zipfile=fname, files=fs)
        if(file.exists(paste0(fname, ".zip"))) {file.rename(paste0(fname, ".zip"), fname)}
      },
      contentType = "application/zip"
    )
    output$downloadExample <- downloadHandler(
      filename = function() {
        "exampleSet.zip"
      },
      content = function(file) {
        file.copy("exampleSet.zip", file)
      },
      contentType = "application/zip"
    )
   
  })
 
  
})
