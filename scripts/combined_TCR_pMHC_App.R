## Libraries ####
library(shiny)
library(shinydashboard)
library(readxl)
library(ggplot2)
library(plotly)
library(ggalluvial)
library(reshape2)
library(dplyr)
library(devtools)
library(InteractiveComplexHeatmap)
library(tidyr)
library(tidyverse)
library(circlize)
library(igraph)

# load dataframe for TCR distances and epitope binding information
all_tcrs <- read_xlsx("data/TCR_Epitope_activity_updated.xlsx")
TCR_epitope <- read_csv("TCR_Epitope_activity_updated.csv")
dist_met <- c("BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM10")


## Load Browser functions and annotations ####
#source("browser.R")

## Create dashboard ####
ui <- shinyUI(
  dashboardPage(
    dashboardHeader(title = 'TCR-pMHC Binding'),
    dashboardSidebar(
      sidebarMenu(
      menuItem('About', tabName = 'about'),
      menuItem('Activity clusters - Distance functions', tabName = 'distances'),
      menuItem('Activity Heatmap', tabName = 'heatmap'),
      menuItem('TCR_pMHC mutation interactions', tabName = 'alluvium')
    )),
    dashboardBody(
      tabItems(
        tabItem(tabName = 'about'),
        
        tabItem(tabName = 'distances', 
                h1('Epitope activity by diferent distance functions'),
                fluidRow(
                  box(selectizeInput('dist_method', 'Method',
                      choices = dist_met),
              
                      selectizeInput('TCR_names', 'TCR',
                      choices = unique(all_tcrs$tcr_name)
                    )),
                  box(
                    sliderInput("obs1", "P1",
                                    min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs2", "P2",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs3", "P3",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs4", "P4",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs5", "P5",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    )),
                  box(
                    sliderInput("obs6", "P6",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs7", "P7",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs8", "P8",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs9", "P9",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    ),
                    sliderInput("obs10", "P10",
                                min = 0, max = 1, value = 0.5, step=0.1, ticks= FALSE
                    )
                   )),
                  box(plotlyOutput("plot3.1"))),
                  
      tabItem(tabName = 'heatmap', 
              h1('TCR-pMHC normalized binding activity heatmap'),
              fluidRow(
                box(selectizeInput('peptide', 'Peptide',
                                    choices = unique(TCR_epitope$index_name))
                         #choice which peptide/epitope (by trivial name) is displayed
                         #possible addidtion of choice by publication
                       )),
                box(InteractiveComplexHeatmapOutput())
      ),
      
      tabItem(tabName = 'alluvium', 
              h1('Alluvium plot showing TCR-pMHC bindings'),
              fluidRow(
                box(selectizeInput('peptide', 'Peptide',
                      choices = unique(TCR_epitope$index_name)),
                 
                 sliderInput('activity', "normalized binding activity range", value = c(0.5, 1), 
                             min = 0, max = 1, step = 0.005)
                )),
                box(plotlyOutput("plot5", inline = FALSE))
              )
    )
  )
))

##create shiny server

server <- shinyServer(function(input, output, session){
  
  
  ##Everything for the distance functions
  peptide <- reactive({all_tcrs$peptide[all_tcrs$tcr_name == input$TCR_names]})
  
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
                                                                      input$dist_method,weights()))})
  
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
  
  
  ##Heatmap code
  #Subset original df, only use columns of interest, make it into wide df format
  
  epitope_sub <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide, ]})
  
  smaller_epitope_sub <- reactive({epitope_sub()[ , c('peptide', "tcr_name", 
                                                      "normalized_activity", "position")]})
  
  wide_epitope_sub <- reactive({pivot_wider(smaller_epitope_sub(), 
                                            id_cols = c("peptide", "position"),
                                            names_from = "tcr_name", 
                                            values_from = "normalized_activity") %>%
      column_to_rownames(., var = 'peptide')}) #set peptides as row IDs
  
  
  #Set color squeme
  col_fun = reactive({colorRamp2(c(min(smaller_epitope_sub()$normalized_activity),
                                   max(smaller_epitope_sub()$normalized_activity)), 
                                 c("white", "red"))})
  
  
  
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
  
  TCR_epitope_peptide <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide, ] %>%
      dplyr::filter((between(normalized_activity, 
                             input$activity[1], input$activity[2])))
  })
  
  output$plot5 <- renderPlotly({
    
    if (nrow(TCR_epitope_peptide()) == 0) {
      ggplot() +
        theme_void() +
        theme(axis.line = element_blank()) +
        ggtitle("Oops, there is no data for epitopes within your selected activity range")
      
    } else {
      
      plot_5 <- ggplot(data = TCR_epitope_peptide(),
                       aes(axis1 = tcr_name, axis2 = peptide, y = normalized_activity)) +
        geom_alluvium(aes(fill = tcr_name), curve_type = "cubic") +
        geom_stratum(aes(fill = tcr_name)) +
        geom_text(stat = "stratum",
                  aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("tcr_name", "peptide"),
                         expand = c(0.25, 0.15)) +
        theme_void() +
        guides(fill = guide_legend(title = "TCRs"))
      
      ggplotly(plot_5, height = 750, width = 1000)
    }
    
  })
  
})
