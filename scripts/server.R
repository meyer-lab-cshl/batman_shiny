server1 <- function(input, output, session) {
  
  # Combine the selected variables into a new data frame

  dist_data = readMat("data//Methods_TCRs_XY.mat")
  dist_data=data.frame(dist_data);
  
  nSB = c(13,21,29,35,27,30,15,13)
  nNB = c(72,58,50,51,97,2683,121,146)
  
  xcol <- reactive(which(dist_met == "input$dist_method"))
  
 
  
  selectedData1 <- reactive({
   dist_data[[xcol]][[input$ycol]][1,]
  })
  
  
  selectedData2 <- reactive({
    dist_data[[xcol]][[input$ycol]][2,]
  })
  
  #seqs = readMat("data//epitope_seqs.mat")
  #s <- reactive({seqs[[1]][[input$ycol]]})

  
  #d=reactive({data.frame(selectedData1(),selectedData2(),row.names=s())})
  d=reactive({data.frame(selectedData1(),selectedData2())})
  

  c1=reactive({rgb(rep(1,nSB[input$ycol]), 0, 0)});
  c2=reactive({rgb(0,0, rep(1,nNB[input$ycol]))});
 
  
  xmin=reactive({min(dist_data[[xcol]][[input$ycol]][1,])})
  xmax=reactive({max(dist_data[[xcol]][[input$ycol]][1,])})
  ymin=reactive({min(dist_data[[xcol]][[input$ycol]][2,])})
  ymax=reactive({max(dist_data[[xcol]][[input$ycol]][2,])}) 
  
  
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

