## Amitavas code all pasted together as an app

#all needed packages
library(R.matlab)
library(shiny)
library(ggplot2)
library(plotly)

#set ui

ui <- fluidPage(
  fluidRow(headerPanel('Methods: 1. Atchley, 2. BLOSUM100, 3. Dayhoff, 4. Gonnet, 5. Hamming, 6. PAM450, 7. R Contiguous
              
              |
              TCRs: 1. 868, 2. A42,  3. A6,  4. B7,  5. E7NLV,  6. G10,  7. ILA-1,  8. T5-004')
  ),
  fluidRow(
    column(width=5,
           sidebarPanel(
             numericInput('xcol', 'Method #', 1,min = 1, max = 7),
             numericInput('ycol', 'TCR #', 1, min = 1, max = 8),
           )
    ),
    column(width = 7,
           plotlyOutput("plot1")
    )
  ),
  
)

# set server

server <- function(input, output, session) {
  
  # Combine the selected variables into a new data frame
  
  dist_data = readMat("data//Methods_TCRs_XY.mat")
  dist_data=data.frame(dist_data);
  
  nSB = c(13,21,29,35,27,30,15,13)
  nNB = c(72,58,50,51,97,2683,121,146)
  
  
  selectedData1 <- reactive({
    dist_data[[input$xcol]][[input$ycol]][1,]
  })
  
  
  selectedData2 <- reactive({
    dist_data[[input$xcol]][[input$ycol]][2,]
  })
  
  seqs = readMat("data//epitope_seqs.mat")
  s <- reactive({seqs[[1]][[input$ycol]]})
  
  
  d=reactive({data.frame(selectedData1(),selectedData2(),row.names=s())})
  
  
  c1=reactive({rgb(rep(1,nSB[input$ycol]), 0, 0)});
  c2=reactive({rgb(0,0, rep(1,nNB[input$ycol]))});
  
  
  xmin=reactive({min(dist_data[[input$xcol]][[input$ycol]][1,])})
  xmax=reactive({max(dist_data[[input$xcol]][[input$ycol]][1,])})
  ymin=reactive({min(dist_data[[input$xcol]][[input$ycol]][2,])})
  ymax=reactive({max(dist_data[[input$xcol]][[input$ycol]][2,])}) 
  
  
  output$plot1 <- renderPlotly({
    
    #ggplot(d()) + geom_point(aes(selectedData1(),selectedData2(),alpha = 0.1,color=c(c1(),c2())),show.legend = FALSE)+theme_classic()
    plot_ly(data = d(), x = ~selectedData1(), y = ~selectedData2(),
            marker = list(
              color = c(c1(),c2()),
              size = 5,
              opacity=0.5,
              line = list(
                color = 'rgb(255, 255, 255)',
                width = 0
              )
            ),
            
            type = "scatter", mode = "markers",
            hoverinfo = 'text',
            text = ~paste(row.names(d()))
    )
    
  })
  
}

# run the app

shinyApp(ui, server)
