# Load needed packages

library(R.matlab)
library(shiny)
library(ggplot2)
library(plotly)
library(dbplyr)



# load dataframe for TCR distances and epitope binding information

load(file = "All_TCRs_epitopes.Rda")

# Create ui for shiny App

ui3 <- fluidPage(
  fluidRow(headerPanel("Strong and non binding epitopes to different TCRs,
                       using multiple distance functions")
  ),
  
  
  fluidRow(
    column(width = 9,
           sidebarPanel(
             selectizeInput('dist_method', 'Method',
                         choices = unique(all_tcrs$dist)),
             
             selectizeInput('TCR_names', 'TCR',
                         choices = unique(all_tcrs$tcr))
           )
    ),
    column(width = 9,
           plotlyOutput("plot3"), 
           # hover = hoverOpts("seq_hover"),
           # verbatimTextOutput("peptide_seq")
    )
  ),
)


# create server for shiny app

server3 <- function(input, output, session) {
  
  seqs = readMat("data//epitope_seqs.mat")
  epitope_seq <- reactive({seqs[[1]][[input$TCR_names]]})
  
  subset_tcrs <-  reactive(
    all_tcrs[all_tcrs$tcr == input$TCR_names & all_tcrs$dist == input$dist_method, ]
    )

  output$plot3 <- renderPlotly({
    
    ggplot(subset_tcrs(), aes(x , y, label = Sequence)) + 
      geom_point(alpha = 0.5, size = 2, aes(color = Binding)) +
      scale_color_manual(values = c('blue', 'red', 'darkgrey')) +
      theme_classic() + 
      theme(
        #legend.background = element_rect(colour = 'darkgrey', linetype = 1, size = 0.3),
            
            #dispose x and y axis
            axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()) +
      
      #add box arround plot
      theme(panel.border = element_rect(color = "#6C7B8B", fill = NA)) +
      
      labs(x = "relative position in X",
           y = "relative position in y")
    ggplotly(tooltip = c("color", "label"))
  })
}

