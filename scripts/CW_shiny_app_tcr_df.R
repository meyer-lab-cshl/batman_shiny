# Load needed packages

library(readr)
library(shiny)
library(ggplot2)
library(plotly)
library(dbplyr)
library(readxl)
library(igraph) #load network plotting package


# load dataframe for TCR distances and epitope binding information
all_tcrs <- read_xlsx("data/TCR_Epitope_activity_updated.xlsx")


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



## App 3.1 with new distance function calculation + df

dist_met <- c("BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM10")

ui3.1 <- fluidPage(
  fluidRow(headerPanel("Strong and non binding epitopes to different TCRs,
                       using multiple distance functions")
  ),
  
  
  fluidRow(
    column(width = 12,
           sidebarPanel(
             selectizeInput('dist_method', 'Method',
                            choices = dist_met),
             
             selectizeInput('TCR_names', 'TCR',
                            choices = unique(all_tcrs$tcr_name))
           )
    ),
    column(width = 12,
           plotlyOutput("plot3.1"), 
           
    )
  ),
)




server3.1 <- function(input, output, session) {
  
  peptide <- reactive({all_tcrs$peptide[all_tcrs$tcr_name == input$TCR_names]})

  generate_epitope_coordinates <- function(peptide_list, distance_metric){
  
  # Load Amino Acid distance matrix
  load(paste("data/", distance_metric, ".rda", sep = ""))
  AA_list <- c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')
  
  AA_distance <- matrix(aa_mat, nrow = 20, ncol = 20, dimnames = list(AA_list, AA_list))
  
  peptide_length = length(unlist(strsplit(peptide_list[1], split = ""))) 
  #extract peptide length from first peptide
  
  position_dependent_distances <- matrix(0, nrow = 1, ncol = peptide_length)
  #This will store position-dependent AA distances for each peptide pair
  
  epitope_adjacency_matrix <- matrix(0, nrow = length(peptide_list), 
                                        ncol = length(peptide_list))
  #matrix that stores epitope-epitope distances
  
  for (peptide1 in 1:length(peptide_list)) {
    for (peptide2 in 1:length(peptide_list)) {
      for (position in 1:peptide_length) {
        position_dependent_distances[1,position] <- AA_distance[
                      unlist(strsplit(peptide_list[peptide1], split = ""))[position],
                      unlist(strsplit(peptide_list[peptide2], split = ""))[position]
          ]
      }
      epitope_adjacency_matrix[peptide1,peptide2] <- sum(position_dependent_distances)
      #Add position-dependent distances to get total distance
    }
  }
  epitope_coordinates <- layout_with_fr(
    graph_from_adjacency_matrix(epitope_adjacency_matrix, weighted = TRUE)
    )
  #Output array of dims #peptides X 2
  return(epitope_coordinates)
}

coordinates <- reactive({as.data.frame(generate_epitope_coordinates(peptide(), 
                                                                    input$dist_method))})

activity <- reactive({all_tcrs$normalized_activity[all_tcrs$tcr_name == 
                                                            input$TCR_names]})

plot_df <- reactive({cbind(coordinates(), activity(), peptide())})

output$plot3.1 <- renderPlotly({
  
  ggplot(plot_df(), aes(V1 , V2)) + 
    geom_point(alpha = 0.8, size = 2, aes(color = plot_df()$activity, 
                                          label = plot_df()$peptide)) +
    scale_colour_gradient(low = "grey", high = "red") +
    labs(color = 'normalized binding activity') +
    theme_minimal() + 
    theme(
      axis.line=element_blank(),
      axis.text.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks=element_blank(),
      axis.title.x=element_blank(),
      axis.title.y=element_blank()) +
    
    #add box arround plot
    theme(panel.border = element_rect(color = "#6C7B8B", fill = NA))
    
    
  ggplotly(tooltip = c('plot_df()$activity', 'plot_df()$peptide'))
})
}
