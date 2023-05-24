library(R.matlab)
library(shiny)
library(ggplot2)
library(plotly)

dist_data = readMat("data//Methods_TCRs_XY.mat")
dist_data = data.frame(dist_data)

dist_met <- data.frame(dist_func = as.factor(c("Atchley", "BLOSUM100",
                                               "Dayhoff", "Gonnet", "Hamming",
                                               "PAM50", "R Contigous")))

ui1 <- fluidPage(
  fluidRow(headerPanel('Methods: 1. Atchley, 2. BLOSUM100, 3. Dayhoff, 4. Gonnet, 5. Hamming, 6. PAM450, 7. R Contiguous
              
              |
              TCRs: 1. 868, 2. A42,  3. A6,  4. B7,  5. E7NLV,  6. G10,  7. ILA-1,  8. T5-004')
           ),
  fluidRow(
    column(width = 9,
           sidebarPanel(
              selectInput('dist_method', 'Method-name',
                          choices = dist_met$dist_func),
              #numericInput('xcol', 'Method',
              #              reactive({as.numeric(input$dist_method)})),
              numericInput('xcol', 'Method', 1, min =1, max = 7),
              numericInput('ycol', 'TCR #', 1, min = 1, max = 8),
              )
           ),
    column(width = 7,
           plotlyOutput("plot1")
           )
    ),
  
  
)

