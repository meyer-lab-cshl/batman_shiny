server1 <- function(input, output, session) {
  
  # Combine the selected variables into a new data frame. and add lists for distance functions and TCRs

  dist_data = readMat("data//Methods_TCRs_XY.mat")
  dist_data=data.frame(dist_data);
  
  dist_met <- c("Atchley", "BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM50",
                "R Contigous")
  
  TCR_name <- c("868", "A42", "A6", "B7", "E7NLV",  "G10",  "ILA-1", "T5-004")
 
  # Set which epitopes are strong and not binders
   
  nSB = c(13,21,29,35,27,30,15,13)
  nNB = c(72,58,50,51,97,2683,121,146)
  
  
  # create numeric input from chosen dist. function and TCR
  
  xcol <- reactive({
    tmp1 <- which(dist_met == input$dist_method)
  })
  
  ycol <- reactive({
    tmp2 <- which(TCR_name == input$TCR_names)
  })
    
  
  #Subset data frame to chosen element/ distances
  
  selectedData1 <- reactive({
   dist_data[[xcol()]][[ycol()]][1,]
  })
  
  
  selectedData2 <- reactive({
    dist_data[[xcol()]][[ycol()]][2,]
  })
  
  seqs = readMat("data//epitope_seqs.mat")
  s <- reactive({seqs[[1]][[ycol()]]})

  
  d=reactive({data.frame(selectedData1(),selectedData2(),row.names=s())})
  
  
  c1=reactive({rgb(rep(1,nSB[ycol()]), 0, 0)});
  c2=reactive({rgb(0,0, rep(1,nNB[ycol()]))});
 
  
  xmin=reactive({min(dist_data[[xcol()]][[ycol()]][1,])})
  xmax=reactive({max(dist_data[[xcol()]][[ycol()]][1,])})
  ymin=reactive({min(dist_data[[xcol()]][[ycol()]][2,])})
  ymax=reactive({max(dist_data[[xcol()]][[ycol()]][2,])}) 
  
  
  output$plot1 <- renderPlotly({
    
    #ggplot(d()) + geom_point(aes(selectedData1(),selectedData2(),alpha = 0.1,color=c(c1(),c2())),show.legend = FALSE)+theme_classic()
    plot_ly(data = d(), x = ~selectedData1(), y = ~selectedData2(),
            marker = list(
              color = c(c1(),c2()),
              size = 5,
              opacity=0.4,
              line = list(
                color = 'rgb(255, 255, 255)',
                width = 0
              )
            ),
            
            type = "scatter", mode = "markers",
            hoverinfo = 'text',
            text = ~paste(row.names(d()))
             )  %>%
      layout(xaxis = list(title = 'Distance in x'), 
             yaxis = list(title = 'Distance in y'), 
            legend = list(title=list(text='<b> Binding epitopes </b>'))
            ) %>%
      add_trace(x = ~selectedData1(), y = ~selectedData2(), 
                type = "scatter",
                mode = "markers",
                name = 'non binder') %>%
      add_trace(x = ~selectedData1(), y = ~selectedData2(),
                type = "scatter",
                mode = "markers",
                name = 'strong binder')
      
    })
  
}

