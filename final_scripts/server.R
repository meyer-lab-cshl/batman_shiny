#.libPaths("/opt/R/4.2.1/lib/R/library")
.libPaths(c("/usr/lib64/R/library","/usr/share/R/library"))

## Libraries ####
library(plotly)
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(readxl)
library(ggplot2)
library(viridis)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(InteractiveComplexHeatmap)
library(circlize)
library(igraph)

# load dataframe for TCR distances and epitope binding information
TCR_epitope <- read_xlsx("TCR_epitope_database.xlsx")
dist_met <- c("BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM10")

##create shiny server

server <- shinyServer(function(input, output, session){
  
  
  ##Everything for the distance functions
  peptide <- reactive({TCR_epitope$peptide[TCR_epitope$tcr_name == input$TCR_names]})
  
  generate_epitope_coordinates <- function(peptide_list, distance_metric, weights){
    
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
        epitope_adjacency_matrix[peptide1,peptide2] <- 1/(1+position_dependent_distances%*%weights)
        #Add position-dependent distances to get total distance
      }
    }
    set.seed(100)
    epitope_coordinates <- layout_with_fr(
      graph_from_adjacency_matrix(epitope_adjacency_matrix, weighted = TRUE),niter=5000
    )
    #Output array of dims #peptides X 2
    return(epitope_coordinates)
  }
  
  peptide_length<-reactive({length(unlist(strsplit(peptide(), split = "")))/length(peptide())})
  
  weights <- reactive({
    array(c(input$obs1, input$obs2, input$obs3, input$obs4, input$obs5, input$obs6, 
            input$obs7, input$obs8, input$obs9, input$obs10),
          dim = c(peptide_length(),1))})
  
  coordinates <- reactive({as.data.frame(generate_epitope_coordinates(peptide(), 
                                                                      input$dist_method, weights()))})
  
  activity <- reactive({TCR_epitope$normalized_peptide_activity[TCR_epitope$tcr_name == 
                                                       input$TCR_names]})
  
  plot_df <- reactive({cbind(coordinates(), activity(), peptide())})
  
  output$plot3.1 <- renderPlotly({
    
    ggplot(plot_df(), aes(V1 , V2)) + 
      geom_point(alpha = 0.8, size = 2, aes(color = plot_df()$activity, 
                                            label = plot_df()$peptide)) +
      scale_colour_viridis_c(option = "plasma") +
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
  
  
  ##Heatmap code
  #Subset original df, only use columns of interest, make it into wide df format
  
  epitope_sub <- reactive({TCR_epitope[TCR_epitope$index_peptide == input$peptide1, ]})
  
  smaller_epitope_sub <- reactive({epitope_sub()[ , c('peptide', "tcr_name", 
                                                      "normalized_peptide_activity", "position")]})
  
  wide_epitope_sub <- reactive({pivot_wider(smaller_epitope_sub(), 
                                            id_cols = c("peptide", "position"),
                                            names_from = "tcr_name", 
                                            values_from = "normalized_peptide_activity") %>%
      column_to_rownames(., var = 'peptide')}) #set peptides as row IDs
  
  
  #Set color squeme
  col_fun = reactive({colorRamp2(c(min(smaller_epitope_sub()$normalized_peptide_activity), 
                                   0.25 * max(smaller_epitope_sub()$normalized_peptide_activity),
                                   0.5 * max(smaller_epitope_sub()$normalized_peptide_activity),
                                   0.75 * max(smaller_epitope_sub()$normalized_peptide_activity),
                                   max(smaller_epitope_sub()$normalized_peptide_activity)), 
                                 c("#0D0887FF", "#7E03A8FF", "#CC4678FF", "#F89441FF", "#F0F921FF"))})
  
  
  
  #Draw interactive Headmap
  
  row_split <- reactive({wide_epitope_sub()[ ,1]})
  sub_df <- reactive({as.matrix(wide_epitope_sub()[, 2:ncol(wide_epitope_sub())])})
  
  
  observe ({
    
    ht1 = Heatmap(
      sub_df(), 
      name = "normalized peptide activity", 
      col = col_fun(), 
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_split = row_split(), 
      row_title_rot = 0,
      column_names_rot = 45)
    
    makeInteractiveComplexHeatmap(input, output, session, ht1)
    
  })
  
  
  ##Alluvium plot code
  
  TCR_epitope_peptide <- reactive({TCR_epitope[TCR_epitope$index_peptide == input$peptide2, ] %>%
      dplyr::filter((between(normalized_peptide_activity, 
                             input$activity[1], input$activity[2])))
  })
  
  output$plot5 <- renderPlotly({
    
    if (sum(TCR_epitope_peptide()$normalized_peptide_activity) == 0) {
      ggplot() +
        theme_void() +
        theme(axis.line = element_blank()) +
        ggtitle("Oops, there is no data for epitopes within your selected activity range")
      
    } else {
      
      plot_5.1 <- ggplot(data = TCR_epitope_peptide(),
                         aes(axis1 = tcr_name, axis2 = peptide, y = normalized_peptide_activity)) +
        geom_alluvium(aes(fill = tcr_name), curve_type = "cubic") +
        geom_stratum(aes(fill = tcr_name)) +
        geom_text(stat = "stratum",
                  aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("tcr_name", "peptide"),
                         expand = c(0.25, 0.15)) +
        theme_void() +
        guides(fill = guide_legend(title = "TCRs"))
      
      ggplotly(plot_5.1)
    }
    
  })
  
})
