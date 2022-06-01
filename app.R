#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)
library(shinydashboard)    # for Dashboard
library(shinyWidgets)      # for radio button widgets
library(shinyalert)        # for alert message very nice format
library(dplyr)             # select functions are covered in the library
library(plyr)              # empty() function is from this package
library(DT)                # for using %>% which works as a pipe in R code for filter etc
library(shinyjs)           # to perform common useful JavaScript operations in Shiny apps
library(plotly)            # to prepare interactive graphs/ charts
library(tidyverse)         # used to make pivot from dataset
library(ggplot2)           # to draw or take advantage of ggplot functions
library(InformationValue)  # is to generate KS plot and ks stat
library(shinyBS)           # use shiny alert function
library(sjPlot)            # Data Visualization for Statistics
library(GEOquery)          # get GEO datasets 
library(preprocessCore)    # Normalize the datasets
source("GSE57387model.R")  # Gather GSE57387's data to build prediction models for integrated risk predictor
source("GSE76882model.R")  # Gather GSE76882's data to build prediction models for integrated risk predictor
source("GSE22459model.R")  # Gather GSE22459's data to build prediction models for integrated risk predictor
##########################
#########Read data########
##########################

###Round function####
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
} # A function that sets all datasets round to two decimal places

### Import top gene data 
merged <- read.csv("GSE57387.csv")
gse76882 <- read.csv("GSE76882.csv")
scorecomp <- as.data.frame(read.csv("score.csv"))
gse22459 <- read.csv("GSE22459.csv")
gse25902 <- read.csv("GSE25902.csv")
## Round data & switch outcome to front position
gse57387 <- merged[,-c(14,15,16,17,18,19,20)]
gse57387_r <- round_df(gse57387, 2)
gse76882_r <- round_df(gse76882, 2)
gse22459_r <- round_df(gse22459, 2)
gse25902_r <- round_df(gse25902, 2)
gse57387_r <- gse57387_r[,c(1,15,2,3,4,5,6,7,8,9,10,11,12,13,14)]
gse76882_r <- gse76882_r[,c(1,15,2,3,4,5,6,7,8,9,10,11,12,13,14)]
gse22459_r <- gse22459_r[,c(1,15,2,3,4,5,6,7,8,9,10,11,12,13,14)]
gse25902_r <- gse25902_r[,c(1,15,16,17,2,3,4,5,6,7,8,9,10,11,12,13,14)]


###Top gene data's DE analysis
DE57387 <- read.csv("DE57387.csv")
rownames(DE57387) <- DE57387[,1]
DE76882 <- read.csv("DE76882.csv")
rownames(DE76882) <- DE76882[,1]
DE22459 <- read.csv("DE22459.csv")
rownames(DE22459) <- DE22459[,1]
DE25902 <- read.csv("DE25902.csv")
rownames(DE25902) <- DE25902[,1]
##########################
#########clean data#######
## For cadi prediction
# Eight variables in the linear regression model were selected and new tables were generated with them and outcome.
top <- c("MACC1","LAMC2","ITGB6","GABRP","DUSP6","C11ORF53","AGR2","RORA")
cadi <- merged["m12.cadi.ch1"]
rownames(cadi) <- merged[,1]
base_57387_top <- merged[,top]
rownames(base_57387_top) <- gse57387[,1] 
base_57387_top_m <- merge(base_57387_top, cadi, by = 'row.names', all = TRUE)
base_57387_top_m <- base_57387_top_m[,-1]



## For cadi vs IFTA
# Convert ifta_score and cadi_score in the dataset to character class. This is to facilitate the generation of Dencity plot
scorecomp$ifta_score <- as.character(scorecomp$ifta_score)
scorecomp$cadi_score <- as.character(scorecomp$cadi_score)

#########################
###########regression####
## use cadi_score to predict ifta_score
#This regression model are built to predict IFTA score in page 6.
lm1 = lm(ifta_score ~ as.numeric(cadi_score), scorecomp)


ui <- dashboardPage(skin = "green",
    dashboardHeader(title = "Kidney"),
### Siderbar -------------------------------------------------------------
    dashboardSidebar(
        sidebarMenu(
            menuItem(" Instruction", tabName = "instruct", icon = icon("sistrix")), 
            menuItem("Dataset Information", tabName = "Data", icon = icon("database")),
            menuItem("Dataset Links", tabName = "Links", icon = icon("compass"),
                   menuSubItem("GSE57387", href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57387', newtab = FALSE),
                   menuSubItem("GSE76882", href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76882', newtab = FALSE),
                   menuSubItem("GSE22459", href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22459', newtab = FALSE),
                   menuSubItem("GSE25902", href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25902', newtab = FALSE)),
            
            
            menuItem("DE Outcomes", tabName = "DE", icon = icon("list")),
            menuItem("CADI Prediction", tabName = "cadi", icon = icon("calendar-day")),
            menuItem("CADI vs IFTA", tabName = "vs", icon = icon("server")),
            menuItem("Prediction", icon = icon("clipboard"), tabName = "Prediction"),
            submitButton("RUN THE APP", width = '100%')
        )
    ),

### Body -----------------------------------------------------------------
    dashboardBody(
        tabItems(
          tabItem(tabName = "instruct",
                  h2("Instruction"),
                  h3("- This tab is the instruction of this shiny app, please read carefully"),
                  h3("- This shiny app focus on risk prediction of kidney trasplantation and include 5 related datasets"),
                  h3("- There are guidence buttons distributed among the tabs. Please read the guide before interacting with the app"),
                  h3("- There is an Action button on sideMenu that needs to be pressed for every execution, as well as for the use of the guide button")
                  
          ),
          
          #Tab for display outcome and 13 top genes in 4 datasets
          tabItem(tabName = "Data",
                  h2("Dataset Overview"),
            fluidRow(
              useShinyalert(),
              actionButton("guide1", "Guide", icon = icon("sistrix"), style='padding:12px; font-size:110%'),
              
              tipify(actionButton("infor6", "GSE57387",icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE57387 dataset was Public on Jul 20, 2016. The title is “Transcriptome signature in early biopsies of stably functioning kidney allografts identify patients at risk for chronic injury”. Its experiment type is expression profiling by array. During a kidney transplant, chronic damage to the kidney can occur, leading to the loss of the graft. The aim of this study was to identify a predictive genome capable of classifying renal grafts as to whether they would be at risk for progressive injury due to fibrosis.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor7", "GSE76882", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "Gene Expression in Biopsies of Acute Rejection and Interstitial Fibrosis/Tubular Atrophy Reveals Highly Shared Mechanisms that Correlate with Worse Long-term Outcomes. We conclude that undetected and/or undertreated immune rejection is leading to IFTA and graft failure.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor8", "GSE22459", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE22459 dataset was Public on Nov 01, 2010. The title is “Fibrosis with Inflammation at One Year Predicts Transplant Functional Decline”. Gene expression in kidney transplant recipients was analyzed at 1-year protocol biopsies and observed decreased graft survival in kidney transplants with interstitial fibrosis with subclinical inflammation, but not fibrosis alone. Early interventions aimed at altering rejection-like inflammation may favor improved long-term KTx survival.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor9", "GSE25902", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE25902 dataset was public on Dec 10, 2011. The title is ‘Innate and adaptive immune responses associate with progressive histological damage of renal allografts’. By comparing gene expression profiling in biopsy with biopsy gene expression in acute T cell-mediated rejection, it was found that progressive chronic histological injury after kidney transplantation was associated with significant modulation of congenital and adaptive immune responses, several months before histological impairment appeared.", placement="bottom", trigger = "hover over"),

              box(background = "blue", solidHeader = TRUE,height = 700,width = 28,
                tabBox(
                title = "Top Gene Dataset", id = "tabset 1",width = 28, height = 500, 
                tabPanel("GSE57387", mainPanel(width = 28,
                          DT::dataTableOutput("gse1"))),
              tabPanel("GSE76882",
                  mainPanel(width = 28,
                            DT::dataTableOutput("gse2"))),
              tabPanel("GSE22459",
                       mainPanel(width = 28,
                                 DT::dataTableOutput("gse3"))),
              tabPanel("GSE25902",
                       mainPanel(width = 28,
                                 DT::dataTableOutput("gse4"))),
)),
              useShinyalert()
              )),
          
          # Tab for DE analysis was performed on four original data sets to select 13 top genes
          tabItem(tabName = "DE",
                  h2("DE_Dataset Overview"),
            fluidRow(
              useShinyalert(),
              
              tipify(actionButton("infor1", "GSE57387",icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE57387 dataset was Public on Jul 20, 2016. The title is “Transcriptome signature in early biopsies of stably functioning kidney allografts identify patients at risk for chronic injury”. Its experiment type is expression profiling by array. During a kidney transplant, chronic damage to the kidney can occur, leading to the loss of the graft. The aim of this study was to identify a predictive genome capable of classifying renal grafts as to whether they would be at risk for progressive injury due to fibrosis.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor2", "GSE76882", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "Gene Expression in Biopsies of Acute Rejection and Interstitial Fibrosis/Tubular Atrophy Reveals Highly Shared Mechanisms that Correlate with Worse Long-term Outcomes. We conclude that undetected and/or undertreated immune rejection is leading to IFTA and graft failure.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor3", "GSE22459", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE22459 dataset was Public on Nov 01, 2010. The title is “Fibrosis with Inflammation at One Year Predicts Transplant Functional Decline”. Gene expression in kidney transplant recipients was analyzed at 1-year protocol biopsies and observed decreased graft survival in kidney transplants with interstitial fibrosis with subclinical inflammation, but not fibrosis alone. Early interventions aimed at altering rejection-like inflammation may favor improved long-term KTx survival.", placement="bottom", trigger = "hover over"),
              tipify(actionButton("infor4", "GSE25902", icon = icon("play-circle"), style='padding:12px; font-size:110%'), "The GSE25902 dataset was public on Dec 10, 2011. The title is ‘Innate and adaptive immune responses associate with progressive histological damage of renal allografts’. By comparing gene expression profiling in biopsy with biopsy gene expression in acute T cell-mediated rejection, it was found that progressive chronic histological injury after kidney transplantation was associated with significant modulation of congenital and adaptive immune responses, several months before histological impairment appeared.", placement="bottom", trigger = "hover over"),
       
              box(background = "green", solidHeader = TRUE,height = 700,
                tabBox(
                  title = "DE_analysis Top Gene", id = "tabset 2",width = 14, height = 500,
                  tabPanel("DE_57387", 
                    selectInput("variable1", "GENE Symbol:",
                              c(
                                "adj.P.Val" = "adj.P.Val",
                                "P.Value" = "P.Value"
                                )),
                            DT::dataTableOutput("de57387")
              ),
              tabPanel("DE_76882",
                  selectInput("variable2", "GENE Symbol:",
                              c(
                                "adj.P.Val" = "adj.P.Val",
                                "P.Value" = "P.Value"
                              )),
                  DT::dataTableOutput("de76882")
              ),
              
              tabPanel("DE_22459",
                       selectInput("variable3", "GENE Symbol:",
                                   c(
                                     "adj.P.Val" = "adj.P.Val",
                                     "P.Value" = "P.Value"
                                   )),
                       DT::dataTableOutput("de22459")),
              tabPanel("DE_25902",
                       selectInput("variable4", "GENE Symbol:",
                                   c(
                                     "adj.P.Val" = "adj.P.Val",
                                     "P.Value" = "P.Value"
                                   )),
                       DT::dataTableOutput("de25902"))
          )),
          actionButton("guide2", "Guide", icon = icon("sistrix"), style='padding:16px; font-size:160%')
          )),
          # Tab for CADI score prediction
          tabItem(tabName = "cadi",
                  h2("Using Biopsy Result to Predict CADI Value"),
                  fluidRow(
                   
                      box(
                          title = "Input Biopsy Result", background = "blue", solidHeader = TRUE,
                          column(4, numericInput("ca1", "MACC1",value = 5,min = 0,width = '100%'),numericInput("ca2", "LAMC2",value = 5,min = 0,width = '100%'), numericInput("ca3", "ITGB6",value = 5,min = 0,width = '100%'),numericInput("ca4", "GABRP",value = 5,min = 0,width = '100%')),
                          column(4, numericInput("ca5", "DUSP6",value = 5,min = 0,width = '100%'),numericInput("ca6", "C11ORF53",value = 5,min = 0,width = '100%'), numericInput("ca7", "AGR2",value = 5,min = 0,width = '100%'),numericInput("ca8", "RORA",value = 5,min = 0,width = '100%')),

                          
                          fluidRow(hr(),
                                   h3("CADI score Result"),
                            column(4, verbatimTextOutput("numeric1")))
                      ),
                  actionButton("guide3", "Guide", icon = icon("sistrix"), style='padding:16px; font-size:160%')
                  )),
          # Tab for IFTA score prediction
          tabItem(tabName = "vs",
                  h2("CADI Score VS. IFTA Score"),
                  fluidRow(
                    box(
                        title = "Input CADI score", background = "olive", solidHeader = TRUE,
                        
                        numericInput("num5", label = h3("CADI score input"), value = 5),
                        hr(),
                        h3("IFTA Score Result"),
                        fluidRow(column(6, verbatimTextOutput("CADI")))
                    ),
                    box(
                      title = "Prediction", background = "green", solidHeader = TRUE,
                      plotOutput("plot4", height = 500)
                    ),
                    useShinyalert(),
                    actionButton("guide5", "Guide", icon = icon("sistrix"), style='padding:16px; font-size:160%')
                  )
                  
            
          ),
        # Tab for Integrated risk predictor
        tabItem(tabName = "Prediction",
                    
                    tabsetPanel(
                      tabPanel("CADI Condition",
                               box(
                                 title = "Information of the patient", status = "info", solidHeader = TRUE,
                                 print("By entering the non-zero values of a patient, you will get outcomes trained by your selected genes.The outcomes will present in the form: 1. The range visualization of accuracies 2. Selected genes 3. Patient Infomation 4. The predicted outcomes."),
                                 fluidRow(
                                   column(2, numericInput("LAMC2", "LAMC2",value = 0,min = 0,width = '100%'),numericInput("CCL2", "CCL2",value = 0,min = 0,width = '100%'), numericInput("KLHL13", "KLHL13",value = 0,min = 0,width = '100%'),numericInput("CHST9", "CHST9",value = 0,min = 0,width = '100%'),numericInput("repeat57387", "CV repeated times",value = 10,min = 1,width = '100%')),
                                   column(2, numericInput("GABRP", "GABRP",value = 0,min = 0,width = '100%'),numericInput("DUSP6", "DUSP6",value = 0,min = 0,width = '100%'),numericInput("LIX1", "LIX1",value = 0,min = 0,width = '100%'),numericInput("FJX1", "FJX1",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("GCNT3", "GCNT3",value = 0,min = 0,width = '100%'),numericInput("ITGB6", "ITGB6",value = 0,min = 0,width = '100%'),numericInput("TTC39A", "TTC39A",value = 0,min = 0,width = '100%'),numericInput("NFKBIZ", "NFKBIZ",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("MACC1", "MACC1",value = 0,min = 0,width = '100%'),numericInput("PROM1", "PROM1",value = 0,min = 0,width = '100%'),numericInput("RNF182", "RNF182",value = 0,min = 0,width = '100%'),numericInput("CDC42SE2", "CDC42SE2",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("AGR2", "AGR2",value = 0,min = 0,width = '100%'),numericInput("SLC34A2", "SLC34A2",value = 0,min = 0,width = '100%'),numericInput("C11ORF53", "C11ORF53",value = 0,min = 0,width = '100%'),numericInput("CHCHD10", "CHCHD10",value = 0,min = 0,width = '100%')),
                                 )),
                               box(
                                 title = "Visualization and Prediction", status = "info", solidHeader = TRUE,
                                 plotOutput("boxplot_57387"),
                                 textOutput("predict_57387")
                               )), 
                      tabPanel("IFTA Condition", 
                               box(
                                 title = "Information of the patient", status = "info", solidHeader = TRUE,
                                 print("By entering the non-zero values of a patient, you will get outcomes trained by your selected genes.The outcomes will present in the form: 1. The range visualization of accuracies 2. Selected genes 3. Patient Infomation 4. The predicted outcomes."),
                                 fluidRow(
                                   column(2, numericInput("AGR2", "AGR2",value = 0,min = 0,width = '100%'),numericInput("AGR3", "AGR3",value = 0,min = 0,width = '100%'), numericInput("KCNT2", "KCNT2",value = 0,min = 0,width = '100%'),numericInput("SYT11", "SYT11",value = 0,min = 0,width = '100%'),numericInput("repeat76882", "CV repeated times",value = 10,min = 1,width = '100%')),
                                   column(2, numericInput("PDLIM1", "PDLIM1",value = 0,min = 0,width = '100%'),numericInput("IGK", "IGK@",value = 0,min = 0,width = '100%'),numericInput("CLU", "CLU",value = 0,min = 0,width = '100%'),numericInput("RABGGTB", "RABGGTB",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("LMBR1L", "LMBR1L",value = 0,min = 0,width = '100%'),numericInput("GABRP", "GABRP",value = 0,min = 0,width = '100%'),numericInput("MUTED", "MUTED",value = 0,min = 0,width = '100%'),numericInput("LOC652493", "LOC652493",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("IGLL1", "IGLL1",value = 0,min = 0,width = '100%'),numericInput("IGLV2", "IGLV2-23",value = 0,min = 0,width = '100%'),numericInput("TMPRSS4", "TMPRSS4",value = 0,min = 0,width = '100%'),numericInput("C20orf7", "C20orf7",value = 0,min = 0,width = '100%')),
                                   column(2, numericInput("IGH", "IGH@",value = 0,min = 0,width = '100%'),numericInput("ITM2C", "ITM2C",value = 0,min = 0,width = '100%'),numericInput("LOC729020", "LOC729020",value = 0,min = 0,width = '100%'),numericInput("LOC100130100", "LOC100130100",value = 0,min = 0,width = '100%'))
                                 )),
                               box(
                                 title = "Visualization and Prediction", status = "info", solidHeader = TRUE,
                                 plotOutput("boxplot_76882"),
                                 textOutput("predict_76882")
                               )),
                      tabPanel("Inflammation",box(
                        title = "Information of the patient", status = "info", solidHeader = TRUE,
                        print("By entering the non-zero values of a patient, you will get outcomes trained by your selected genes.The outcomes will present in the form: 1. The range visualization of accuracies 2. Selected genes 3. Patient Infomation 4. The predicted outcomes."),
                        fluidRow(
                          column(2, numericInput("RP11", "RP11",value = 0,min = 0,width = '100%'),numericInput("KRT17P5", "KRT17P5",value = 0,min = 0,width = '100%'), numericInput("AX747826", "AX747826",value = 0,min = 0,width = '100%'),numericInput("PPIF", "PPIF",value = 0,min = 0,width = '100%'),numericInput("repeat22459", "CV repeated times",value = 10,min = 1,width = '100%')),
                          column(2, numericInput("RP5", "RP5",value = 0,min = 0,width = '100%'),numericInput("DPF1", "DPF1",value = 0,min = 0,width = '100%'),numericInput("MROH6", "MROH6",value = 0,min = 0,width = '100%'),numericInput("CHRM5", "CHRM5",value = 0,min = 0,width = '100%')),
                          column(2, numericInput("SLC16A11", "SLC16A11",value = 0,min = 0,width = '100%'),numericInput("DYNLL2", "DYNLL2",value = 0,min = 0,width = '100%'),numericInput("LOC101927131", "LOC101927131",value = 0,min = 0,width = '100%'),numericInput("CEP55", "CEP55",value = 0,min = 0,width = '100%')),
                          column(2, numericInput("BSND", "BSND",value = 0,min = 0,width = '100%'),numericInput("LOC100130700", "LOC100130700",value = 0,min = 0,width = '100%'),numericInput("SDHAP1", "SDHAP1",value = 0,min = 0,width = '100%'),numericInput("SERPINE2", "SERPINE2",value = 0,min = 0,width = '100%')),
                          column(2, numericInput("C15orf27", "C15orf27",value = 0,min = 0,width = '100%'),numericInput("LA16c", "LA16c",value = 0,min = 0,width = '100%'),numericInput("DALRD3", "DALRD3",value = 0,min = 0,width = '100%'),numericInput("TMSB10", "TMSB10",value = 0,min = 0,width = '100%'))
                        )),
                        box(
                          title = "Visualization and Prediction", status = "info", solidHeader = TRUE,
                          plotOutput("boxplot_22459"),
                          textOutput("predict_22459")
                        ))
                    )

                    
                    
          )
            
)))
            



server <- function(input, output) {
    output$menu <- renderMenu({
        sidebarMenu(
                # Setting id makes input$tabs give the tabName of currently-selected tab
            id = "tabs",
            menuItem(" Instruction", tabName = "instruct", icon = icon("sistrix")),
            menuItem("Dataset Infromation", tabName = "Data", icon = icon("database")),
            menuItem("Dataset Links", tabName = "Links", icon = icon("compass")),
            menuItem("DE Outcomes", tabName = "DE", icon = icon("list")),
            menuItem("CADI Prediction", tabName = "cadi", icon = icon("calendar-day")),
            menuItem("CADI vs IFTA", tabName = "vs", icon("server")),
            menuItem("Prediction", icon = icon("clipboard"), tabName = "Prediction"),
            submitButton("RUN THE APP",width = '100%')
            )
        }) # Siderbar Menu
    
    #gse_table output
    output$gse1 <- DT::renderDataTable(
        gse57387_r, options = list(scrollX = TRUE)
        , rownames = FALSE)
    output$gse2 <- DT::renderDataTable(
        gse76882_r, options = list(scrollX = TRUE)
        , rownames = FALSE)
    output$gse3 <- DT::renderDataTable(
      gse22459_r, options = list(scrollX = TRUE)
      , rownames = FALSE)
    output$gse4 <- DT::renderDataTable(
      gse25902_r, options = list(scrollX = TRUE)
      , rownames = FALSE)

    
    ##DE table output
    output$de57387 <- DT::renderDataTable(
        DE57387[,input$variable1, drop = FALSE], options = list(scrollX = TRUE)
        , rownames = TRUE)
    output$de76882 <- DT::renderDataTable(
        DE76882[,input$variable2, drop = FALSE], options = list(scrollX = TRUE)
        , rownames = TRUE)
    output$de22459 <- DT::renderDataTable(
      DE22459[,input$variable3, drop = FALSE], options = list(scrollX = TRUE)
      , rownames = TRUE)
    output$de25902 <- DT::renderDataTable(
      DE25902[,input$variable4, drop = FALSE], options = list(scrollX = TRUE)
      , rownames = TRUE)

    
    # cadi score prediction
    output$numeric1 <- renderPrint({ cadi_predi <- round(-1.482 + 0.6143*input$ca1 + 2.0015*input$ca2 - 0.5837*input$ca3 + 0.9954*input$ca4 + 0.8147*input$ca5 - 1.0749*input$ca6 - 0.9753*input$ca7 - 1.0449*input$ca8)
    # This is the linear regression model for CADI score prediction
    if (input$ca1 < 0){
      print("MAcc1 input invalid")
    }else if (input$ca2 < 0){
      print("LACM2 input invalid")
    }else if (input$ca3 < 0){
      print("ITGB6 input invalid")
    }else if (input$ca4 < 0){
      print("GABRP input invalid")
    }else if (input$ca5 < 0){
      print("DUSP6 input invalid")
    }else if (input$ca6 < 0){
      print("C11ORF53 input invalid")
    }else if (input$ca7 < 0){
      print("AGR2 input invalid")
    }else if (input$ca8 < 0){
      print("RORA input invalid")
    }else if (cadi_predi < 0){
      print(0)
    }else if (cadi_predi > 10){
      print(10)
    }else{
      print(cadi_predi)
    }
    })# Make sure the result is between 0 to 10. And make sure the input value is greater than 0.

    
    
    # Ifta score prediction
    output$CADI <- renderPrint({ new_obs = data.frame(cadi_score = input$num5)
    x <- predict(lm1, new_obs, interval = "prediction", level = 0.90) # Use linear regression model to predict the out come.
    cadi_ifta <- round(x[1])
    if (input$num5 < 0){
      print("input invalid")
    }else if (cadi_ifta < 0){
      print(0)
    }else if (cadi_ifta > 3){
      print(3)
    }else{
      print(cadi_ifta)
    }
    })# Make sure the result is between 0 to 3. And make sure the input value is greater than 0.
    
    
    output$plot4 <- renderPlot({scorecomp %>% ggplot() + aes(x = ifta_score, fill = cadi_score) + geom_density(alpha = 0.5)})
    output$predict_22459 <- renderText({
      selected_genes = c() # Empty sets, in order to collect genes' data selected for prediction
      patientInfo = c() # Empty sets, inorder to store the input gene ID
      if(input$RP11 != 0){
        selected_genes = append(selected_genes, "RP11-672L10.6")
        patientInfo = append(patientInfo, input$RP11)
      }
      if(input$RP5 != 0){
        selected_genes = append(selected_genes, "RP5-1118D24.2")
        patientInfo = append(patientInfo, input$RP5)
      }
      if(input$SLC16A11 != 0){
        selected_genes = append(selected_genes, "SLC16A11")
        patientInfo = append(patientInfo,input$SLC16A11)
      }
      if(input$BSND != 0){
        selected_genes = append(selected_genes, "BSND")
        patientInfo = append(patientInfo,input$BSND)
      }
      if(input$C15orf27 != 0){
        selected_genes = append(selected_genes, "C15orf27")
        patientInfo = append(patientInfo,input$C15orf27)
      }
      if(input$KRT17P5 != 0){
        selected_genes = append(selected_genes, "KRT17P5")
        patientInfo = append(patientInfo,input$KRT17P5)
      }
      if(input$DPF1 != 0){
        selected_genes = append(selected_genes, "DPF1")
        patientInfo = append(patientInfo,input$DPF1)
      }
      if(input$DYNLL2 != 0){
        selected_genes = append(selected_genes, "DYNLL2")
        patientInfo = append(patientInfo,input$DYNLL2)
      }
      if(input$LOC100130700 != 0){
        selected_genes = append(selected_genes, "LOC100130700")
        patientInfo = append(patientInfo,input$LOC100130700)
      }
      if(input$LA16c != 0){
        selected_genes = append(selected_genes, "LA16c-83F12.6")
        patientInfo = append(patientInfo,input$LA16c)
      }
      if(input$AX747826 != 0){
        selected_genes = append(selected_genes, "AX747826")
        patientInfo = append(patientInfo,input$AX747826)
      }
      if(input$MROH6 != 0){
        selected_genes = append(selected_genes, "MROH6")
        patientInfo = append(patientInfo,input$MROH6)
      }
      if(input$LOC101927131 != 0){
        selected_genes = append(selected_genes, "LOC101927131")
        patientInfo = append(patientInfo,input$LOC101927131)
      }
      if(input$SDHAP1 != 0){
        selected_genes = append(selected_genes, "SDHAP1 /// SDHAP2")
        patientInfo = append(patientInfo,input$SDHAP1)
      }
      if(input$DALRD3 != 0){
        selected_genes = append(selected_genes, "DALRD3")
        patientInfo = append(patientInfo,input$DALRD3)
      }
      if(input$PPIF != 0){
        selected_genes = append(selected_genes, "PPIF")
        patientInfo = append(patientInfo,input$PPIF)
      }
      if(input$CHRM5 != 0){
        selected_genes = append(selected_genes, "CHRM5")
        patientInfo = append(patientInfo,input$CHRM5)
      }
      if(input$CEP55 != 0){
        selected_genes = append(selected_genes, "CEP55")
        patientInfo = append(patientInfo,input$CEP55)
      }
      if(input$SERPINE2 != 0){
        selected_genes = append(selected_genes, "SERPINE2")
        patientInfo = append(patientInfo,input$SERPINE2)
      }
      if(input$TMSB10 != 0){
        selected_genes = append(selected_genes, "TMSB10")
        patientInfo = append(patientInfo,input$TMSB10)
      } # Check to see if the gene is not equal to zero, which means the gene is selected for training
      paste("Selected genes: ",toString(selected_genes),".   PatientInfo: ",toString(patientInfo),".   The prediction outcome of this patient is: ",toString(GSE22459_model(selected_genes,input$repeat22459,patientInfo,TRUE)),sep="\n")
    })
    
    output$boxplot_22459 <- renderPlot({
      selected_genes = c()
      patientInfo = c()
      if(input$RP11 != 0){
        selected_genes = append(selected_genes, "RP11-672L10.6")
        patientInfo = append(patientInfo, input$RP11)
      }
      if(input$RP5 != 0){
        selected_genes = append(selected_genes, "RP5-1118D24.2")
        patientInfo = append(patientInfo, input$RP5)
      }
      if(input$SLC16A11 != 0){
        selected_genes = append(selected_genes, "SLC16A11")
        patientInfo = append(patientInfo,input$SLC16A11)
      }
      if(input$BSND != 0){
        selected_genes = append(selected_genes, "BSND")
        patientInfo = append(patientInfo,input$BSND)
      }
      if(input$C15orf27 != 0){
        selected_genes = append(selected_genes, "C15orf27")
        patientInfo = append(patientInfo,input$C15orf27)
      }
      if(input$KRT17P5 != 0){
        selected_genes = append(selected_genes, "KRT17P5")
        patientInfo = append(patientInfo,input$KRT17P5)
      }
      if(input$DPF1 != 0){
        selected_genes = append(selected_genes, "DPF1")
        patientInfo = append(patientInfo,input$DPF1)
      }
      if(input$DYNLL2 != 0){
        selected_genes = append(selected_genes, "DYNLL2")
        patientInfo = append(patientInfo,input$DYNLL2)
      }
      if(input$LOC100130700 != 0){
        selected_genes = append(selected_genes, "LOC100130700")
        patientInfo = append(patientInfo,input$LOC100130700)
      }
      if(input$LA16c != 0){
        selected_genes = append(selected_genes, "LA16c-83F12.6")
        patientInfo = append(patientInfo,input$LA16c)
      }
      if(input$AX747826 != 0){
        selected_genes = append(selected_genes, "AX747826")
        patientInfo = append(patientInfo,input$AX747826)
      }
      if(input$MROH6 != 0){
        selected_genes = append(selected_genes, "MROH6")
        patientInfo = append(patientInfo,input$MROH6)
      }
      if(input$LOC101927131 != 0){
        selected_genes = append(selected_genes, "LOC101927131")
        patientInfo = append(patientInfo,input$LOC101927131)
      }
      if(input$SDHAP1 != 0){
        selected_genes = append(selected_genes, "SDHAP1 /// SDHAP2")
        patientInfo = append(patientInfo,input$SDHAP1)
      }
      if(input$DALRD3 != 0){
        selected_genes = append(selected_genes, "DALRD3")
        patientInfo = append(patientInfo,input$DALRD3)
      }
      if(input$PPIF != 0){
        selected_genes = append(selected_genes, "PPIF")
        patientInfo = append(patientInfo,input$PPIF)
      }
      if(input$CHRM5 != 0){
        selected_genes = append(selected_genes, "CHRM5")
        patientInfo = append(patientInfo,input$CHRM5)
      }
      if(input$CEP55 != 0){
        selected_genes = append(selected_genes, "CEP55")
        patientInfo = append(patientInfo,input$CEP55)
      }
      if(input$SERPINE2 != 0){
        selected_genes = append(selected_genes, "SERPINE2")
        patientInfo = append(patientInfo,input$SERPINE2)
      }
      if(input$TMSB10 != 0){
        selected_genes = append(selected_genes, "TMSB10")
        patientInfo = append(patientInfo,input$TMSB10)
      }
      GSE22459_model(selected_genes,input$repeat22459,patientInfo,FALSE)})
    output$predict_57387 <- renderText({
      selected_genes = c()
      patientInfo = c()
      if(input$LAMC2 != 0){
        selected_genes = append(selected_genes, "LAMC2")
        patientInfo = append(patientInfo, input$LAMC2)
      }
      if(input$GABRP != 0){
        selected_genes = append(selected_genes, "GABRP")
        patientInfo = append(patientInfo, input$GABRP)
      }
      if(input$GCNT3 != 0){
        selected_genes = append(selected_genes, "GCNT3")
        patientInfo = append(patientInfo,input$GCNT3)
      }
      if(input$MACC1 != 0){
        selected_genes = append(selected_genes, "MACC1")
        patientInfo = append(patientInfo,input$MACC1)
      }
      if(input$AGR2 != 0){
        selected_genes = append(selected_genes, "AGR2")
        patientInfo = append(patientInfo,input$AGR2)
      }
      if(input$CCL2 != 0){
        selected_genes = append(selected_genes, "CCL2")
        patientInfo = append(patientInfo,input$CCL2)
      }
      if(input$DUSP6 != 0){
        selected_genes = append(selected_genes, "DUSP6")
        patientInfo = append(patientInfo,input$DUSP6)
      }
      if(input$ITGB6 != 0){
        selected_genes = append(selected_genes, "ITGB6")
        patientInfo = append(patientInfo,input$ITGB6)
      }
      if(input$PROM1 != 0){
        selected_genes = append(selected_genes, "PROM1")
        patientInfo = append(patientInfo,input$PROM1)
      }
      if(input$SLC34A2 != 0){
        selected_genes = append(selected_genes, "SLC34A2")
        patientInfo = append(patientInfo,input$SLC34A2)
      }
      if(input$KLHL13 != 0){
        selected_genes = append(selected_genes, "KLHL13")
        patientInfo = append(patientInfo,input$KLHL13)
      }
      if(input$LIX1 != 0){
        selected_genes = append(selected_genes, "LIX1")
        patientInfo = append(patientInfo,input$LIX1)
      }
      if(input$TTC39A != 0){
        selected_genes = append(selected_genes, "TTC39A")
        patientInfo = append(patientInfo,input$TTC39A)
      }
      if(input$RNF182 != 0){
        selected_genes = append(selected_genes, "RNF182")
        patientInfo = append(patientInfo,input$RNF182)
      }
      if(input$C11ORF53 != 0){
        selected_genes = append(selected_genes, "C11ORF53")
        patientInfo = append(patientInfo,input$C11ORF53)
      }
      if(input$CHST9 != 0){
        selected_genes = append(selected_genes, "CHST9")
        patientInfo = append(patientInfo,input$CHST9)
      }
      if(input$FJX1 != 0){
        selected_genes = append(selected_genes, "FJX1")
        patientInfo = append(patientInfo,input$FJX1)
      }
      if(input$NFKBIZ != 0){
        selected_genes = append(selected_genes, "NFKBIZ")
        patientInfo = append(patientInfo,input$NFKBIZ)
      }
      if(input$CDC42SE2 != 0){
        selected_genes = append(selected_genes, "CDC42SE2")
        patientInfo = append(patientInfo,input$CDC42SE2)
      }
      if(input$CHCHD10 != 0){
        selected_genes = append(selected_genes, "CHCHD10")
        patientInfo = append(patientInfo,input$CHCHD10)
      }
      paste("Selected genes: ",toString(selected_genes),".   PatientInfo: ",toString(patientInfo),".   The prediction outcome of this patient is: ",toString(GSE57387_model(selected_genes,input$repeat57387,patientInfo,TRUE)),sep="\n")
    })
    output$boxplot_57387 <- renderPlot({
      selected_genes = c()
      patientInfo = c()
      if(input$LAMC2 != 0){
        selected_genes = append(selected_genes, "LAMC2")
        patientInfo = append(patientInfo, input$LAMC2)
      }
      if(input$GABRP != 0){
        selected_genes = append(selected_genes, "GABRP")
        patientInfo = append(patientInfo, input$GABRP)
      }
      if(input$GCNT3 != 0){
        selected_genes = append(selected_genes, "GCNT3")
        patientInfo = append(patientInfo,input$GCNT3)
      }
      if(input$MACC1 != 0){
        selected_genes = append(selected_genes, "MACC1")
        patientInfo = append(patientInfo,input$MACC1)
      }
      if(input$AGR2 != 0){
        selected_genes = append(selected_genes, "AGR2")
        patientInfo = append(patientInfo,input$AGR2)
      }
      if(input$CCL2 != 0){
        selected_genes = append(selected_genes, "CCL2")
        patientInfo = append(patientInfo,input$CCL2)
      }
      if(input$DUSP6 != 0){
        selected_genes = append(selected_genes, "DUSP6")
        patientInfo = append(patientInfo,input$DUSP6)
      }
      if(input$ITGB6 != 0){
        selected_genes = append(selected_genes, "ITGB6")
        patientInfo = append(patientInfo,input$ITGB6)
      }
      if(input$PROM1 != 0){
        selected_genes = append(selected_genes, "PROM1")
        patientInfo = append(patientInfo,input$PROM1)
      }
      if(input$SLC34A2 != 0){
        selected_genes = append(selected_genes, "SLC34A2")
        patientInfo = append(patientInfo,input$SLC34A2)
      }
      if(input$KLHL13 != 0){
        selected_genes = append(selected_genes, "KLHL13")
        patientInfo = append(patientInfo,input$KLHL13)
      }
      if(input$LIX1 != 0){
        selected_genes = append(selected_genes, "LIX1")
        patientInfo = append(patientInfo,input$LIX1)
      }
      if(input$TTC39A != 0){
        selected_genes = append(selected_genes, "TTC39A")
        patientInfo = append(patientInfo,input$TTC39A)
      }
      if(input$RNF182 != 0){
        selected_genes = append(selected_genes, "RNF182")
        patientInfo = append(patientInfo,input$RNF182)
      }
      if(input$C11ORF53 != 0){
        selected_genes = append(selected_genes, "C11ORF53")
        patientInfo = append(patientInfo,input$C11ORF53)
      }
      if(input$CHST9 != 0){
        selected_genes = append(selected_genes, "CHST9")
        patientInfo = append(patientInfo,input$CHST9)
      }
      if(input$FJX1 != 0){
        selected_genes = append(selected_genes, "FJX1")
        patientInfo = append(patientInfo,input$FJX1)
      }
      if(input$NFKBIZ != 0){
        selected_genes = append(selected_genes, "NFKBIZ")
        patientInfo = append(patientInfo,input$NFKBIZ)
      }
      if(input$CDC42SE2 != 0){
        selected_genes = append(selected_genes, "CDC42SE2")
        patientInfo = append(patientInfo,input$CDC42SE2)
      }
      if(input$CHCHD10 != 0){
        selected_genes = append(selected_genes, "CHCHD10")
        patientInfo = append(patientInfo,input$CHCHD10)
      }
      GSE57387_model(selected_genes,input$repeat57387,patientInfo,FALSE)
    })
    output$predict_76882 <- renderText({
      selected_genes = c()
      patientInfo = c()
      if(input$AGR2 != 0){
        selected_genes = append(selected_genes, "AGR2")
        patientInfo = append(patientInfo, input$AGR2)
      }
      if(input$PDLIM1 != 0){
        selected_genes = append(selected_genes, "PDLIM1")
        patientInfo = append(patientInfo, input$PDLIM1)
      }
      if(input$LMBR1L != 0){
        selected_genes = append(selected_genes, "LMBR1L")
        patientInfo = append(patientInfo,input$LMBR1L)
      }
      if(input$IGLL1 != 0){
        selected_genes = append(selected_genes, "IGLL1")
        patientInfo = append(patientInfo,input$IGLL1)
      }
      if(input$IGH != 0){
        selected_genes = append(selected_genes, "IGH@")
        patientInfo = append(patientInfo,input$IGH)
      }
      if(input$AGR3 != 0){
        selected_genes = append(selected_genes, "AGR3")
        patientInfo = append(patientInfo,input$AGR3)
      }
      if(input$IGK != 0){
        selected_genes = append(selected_genes, "IGK@")
        patientInfo = append(patientInfo,input$IGK)
      }
      if(input$GABRP != 0){
        selected_genes = append(selected_genes, "GABRP")
        patientInfo = append(patientInfo,input$GABRP)
      }
      if(input$IGLV2 != 0){
        selected_genes = append(selected_genes, "IGLV2-23")
        patientInfo = append(patientInfo,input$IGLV2)
      }
      if(input$ITM2C != 0){
        selected_genes = append(selected_genes, "ITM2C")
        patientInfo = append(patientInfo,input$ITM2C)
      }
      if(input$KCNT2 != 0){
        selected_genes = append(selected_genes, "KCNT2")
        patientInfo = append(patientInfo,input$KCNT2)
      }
      if(input$CLU != 0){
        selected_genes = append(selected_genes, "CLU")
        patientInfo = append(patientInfo,input$CLU)
      }
      if(input$MUTED != 0){
        selected_genes = append(selected_genes, "MUTED")
        patientInfo = append(patientInfo,input$MUTED)
      }
      if(input$TMPRSS4 != 0){
        selected_genes = append(selected_genes, "TMPRSS4")
        patientInfo = append(patientInfo,input$TMPRSS4)
      }
      if(input$LOC729020 != 0){
        selected_genes = append(selected_genes, "LOC729020")
        patientInfo = append(patientInfo,input$LOC729020)
      }
      if(input$SYT11 != 0){
        selected_genes = append(selected_genes, "SYT11")
        patientInfo = append(patientInfo,input$SYT11)
      }
      if(input$RABGGTB != 0){
        selected_genes = append(selected_genes, "RABGGTB")
        patientInfo = append(patientInfo,input$RABGGTB)
      }
      if(input$LOC652493 != 0){
        selected_genes = append(selected_genes, "LOC652493")
        patientInfo = append(patientInfo,input$LOC652493)
      }
      if(input$C20orf7 != 0){
        selected_genes = append(selected_genes, "C20orf7")
        patientInfo = append(patientInfo,input$C20orf7)
      }
      if(input$LOC100130100 != 0){
        selected_genes = append(selected_genes, "LOC100130100")
        patientInfo = append(patientInfo,input$LOC100130100)
      }
      paste("Selected genes: ",toString(selected_genes),".   PatientInfo: ",toString(patientInfo),".   The prediction outcome of this patient is: ",toString(GSE76882_model(selected_genes,input$repeat76882,patientInfo,TRUE)),sep="\n")
      
    })
    output$boxplot_76882 <- renderPlot({
      selected_genes = c()
      patientInfo = c()
      if(input$AGR2 != 0){
        selected_genes = append(selected_genes, "AGR2")
        patientInfo = append(patientInfo, input$AGR2)
      }
      if(input$PDLIM1 != 0){
        selected_genes = append(selected_genes, "PDLIM1")
        patientInfo = append(patientInfo, input$PDLIM1)
      }
      if(input$LMBR1L != 0){
        selected_genes = append(selected_genes, "LMBR1L")
        patientInfo = append(patientInfo,input$LMBR1L)
      }
      if(input$IGLL1 != 0){
        selected_genes = append(selected_genes, "IGLL1")
        patientInfo = append(patientInfo,input$IGLL1)
      }
      if(input$IGH != 0){
        selected_genes = append(selected_genes, "IGH@")
        patientInfo = append(patientInfo,input$IGH)
      }
      if(input$AGR3 != 0){
        selected_genes = append(selected_genes, "AGR3")
        patientInfo = append(patientInfo,input$AGR3)
      }
      if(input$IGK != 0){
        selected_genes = append(selected_genes, "IGK@")
        patientInfo = append(patientInfo,input$IGK)
      }
      if(input$GABRP != 0){
        selected_genes = append(selected_genes, "GABRP")
        patientInfo = append(patientInfo,input$GABRP)
      }
      if(input$IGLV2 != 0){
        selected_genes = append(selected_genes, "IGLV2-23")
        patientInfo = append(patientInfo,input$IGLV2)
      }
      if(input$ITM2C != 0){
        selected_genes = append(selected_genes, "ITM2C")
        patientInfo = append(patientInfo,input$ITM2C)
      }
      if(input$KCNT2 != 0){
        selected_genes = append(selected_genes, "KCNT2")
        patientInfo = append(patientInfo,input$KCNT2)
      }
      if(input$CLU != 0){
        selected_genes = append(selected_genes, "CLU")
        patientInfo = append(patientInfo,input$CLU)
      }
      if(input$MUTED != 0){
        selected_genes = append(selected_genes, "MUTED")
        patientInfo = append(patientInfo,input$MUTED)
      }
      if(input$TMPRSS4 != 0){
        selected_genes = append(selected_genes, "TMPRSS4")
        patientInfo = append(patientInfo,input$TMPRSS4)
      }
      if(input$LOC729020 != 0){
        selected_genes = append(selected_genes, "LOC729020")
        patientInfo = append(patientInfo,input$LOC729020)
      }
      if(input$SYT11 != 0){
        selected_genes = append(selected_genes, "SYT11")
        patientInfo = append(patientInfo,input$SYT11)
      }
      if(input$RABGGTB != 0){
        selected_genes = append(selected_genes, "RABGGTB")
        patientInfo = append(patientInfo,input$RABGGTB)
      }
      if(input$LOC652493 != 0){
        selected_genes = append(selected_genes, "LOC652493")
        patientInfo = append(patientInfo,input$LOC652493)
      }
      if(input$C20orf7 != 0){
        selected_genes = append(selected_genes, "C20orf7")
        patientInfo = append(patientInfo,input$C20orf7)
      }
      if(input$LOC100130100 != 0){
        selected_genes = append(selected_genes, "LOC100130100")
        patientInfo = append(patientInfo,input$LOC100130100)
      }
      GSE76882_model(selected_genes,input$repeat76882,patientInfo,FALSE)
    })

    # shiny alert
    observeEvent(input$guide1, {
      shinyalert(title = "Guide", text = "This button will instruct you how to use the functions of this tab. 
There will be the same Guide button on each subsequent tab. Please enjoy using this product after reading Guide.
This TAB shows four datasets, including their Outcome and the corresponding 13 top genes. The details of the data set are available at the top buttons.
")
    })

    observeEvent(input$guide2, {
      shinyalert(title = "Guide", text = "DE analysis was performed on four original data sets to select 13 top genes.
This TAB shows the p-value & adj p-value of 13 genes in each dataset.
")
    })
    observeEvent(input$guide3, {
      shinyalert(title = "Guide", text = "Input your biopsy gene values in the corresponding Spaces to get the predicted CADI Score.
")
    })


    observeEvent(input$guide5, {
      shinyalert(title = "Guide", text = "This TAB uses the CADI score you get from the third TAB to predict the IFTA score.
IFTA score: 0, absent; 1 (mild), <25%; 2 (moderate), 25-50%; and 3 (severe), >50% of the total area.
Input your CADI score to predict the IFTA score.
The density chart on the right shows the IFTA score distribution of the CADI Score.
")
    })

    
}

shinyApp(ui, server)

