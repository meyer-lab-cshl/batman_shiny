## Libraries ####
library(shiny)
library(shinydashboard)
library(readxl)
library(ggplot2)
library(plotly)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(InteractiveComplexHeatmap)
library(circlize)
library(igraph)

# load dataframe for TCR distances and epitope binding information
all_tcrs <- read_xlsx("data/TCR_Epitope_activity_updated.xlsx")
TCR_epitope <- read.csv("TCR_Epitope_activity_updated.csv")
dist_met <- c("BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM10")


## Load Browser functions and annotations ####
#source("browser.R")

## Create dashboard ####
ui <- shinyUI(
  dashboardPage(
    dashboardHeader(title = 'TCR-pMHC Binding'),
    dashboardSidebar( width = 250,
      sidebarMenu(
      menuItem('About', tabName = 'about'),
      menuItem('Epitope Cluster', tabName = 'distances'),
      menuItem('Activity Heatmap', tabName = 'heatmap'),
      menuItem('Mutant-TCR interaction', tabName = 'alluvium')
    )),
    dashboardBody(
      tabItems(
        tabItem(tabName = 'about', 
        "This application explores the vast diversity of T-cell receptor (TCR) interaction with antigens, based on their epitopes. TCR expression in human T-cells results from the somatic recombination of V and J genes. This process leads to a theoretically extensive repertoire of TCRs that could potentially recognize various antigens from pathogens. However, even with this theoretical diversity, the human immune system would still not be able to create an immune response against every possible pathogen. To overcome this limitation, T-cells are believed to be cross-reactive, meaning that a single T-cell can recognize many different peptides presented on MHC molecules, according to experimental extrapolations up to 10e6 different peptides.  The same antigens can also be recognized by multiple TCRs, based on distinct antigen features. This cross-reactivity is critical for providing an adequate immune response to different pathogens.
                However, while cross-reactivity is favorable in immune responses against pathogens, it can also have unintended consequences. For instance, it may play a role in the development of auto-immune diseases, or it can lead to off-target effects in cancer immunotherapies.
                The app visualizes the cross-reactivity of 48 TCRs to 5 specific index peptides and all possible single mutations of these peptides. The Epitope Cluster function clusters mutated epitopes according to their mutual sequence similarities, while showing their binding affinity. This information is based on user-defined distance functions and weights assigned to mutations at specific positions within the peptides.
                The Activity Heatmap in the app allows users to visualize the normalized binding activity of each TCR to different peptides. On the other hand, the Mutant-TCR interaction tab reveals which index peptides and their mutations are bound by multiple TCRs and which are recognized by only one TCR, emphasizing how different TCRs specific for the same index peptide recognizes different mutants of it differently.
                In summary, this application provides valuable insights into the cross-reactivity of TCRs and how it influences the immune response to various antigens."
        ),
        
        tabItem(tabName = 'distances', 
                h1('Epitope similarity by diferent distance functions'),
                fluidRow(
                  
                  box(width = 5, title = "Settings", status = 'info', solidHeader = T,
                      
                      selectizeInput('dist_method', 'Method',
                      choices = dist_met),
              
                      selectizeInput('TCR_names', 'TCR',
                      choices = unique(all_tcrs$tcr_name))
                    ),
                  
                  box(width = 3, title = "Position weights", status = 'info', solidHeader = T,
                      
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
                    )
                  ),
                  
                  box(width = 3, title = " ", status = 'info', solidHeader = T,
                      
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
                  )
                ),
                
                fluidRow(
                  box(width = 12, title = "Peptide Distances", status = 'primary', solidHeader = T,
                    plotlyOutput("plot3.1"))
                  )
                ),
                  
      tabItem(tabName = 'heatmap', 
              h1('TCR-pMHC normalized binding activity heatmap'),
              
              fluidRow(
                box(width = 5, title = "Settings", status = 'info', solidHeader = T,
                     selectizeInput('peptide1', 'Peptide',
                                    choices = unique(TCR_epitope$index_name))
                         #possible addition of choice by publication
                       ),
                
                box(width = 9, title = "TCR-pMHC binding activity", status = 'primary', solidHeader = T,
                  InteractiveComplexHeatmapOutput())
        )
      ),
      
      tabItem(tabName = 'alluvium', 
              h1('Alluvium plot showing TCR-pMHC bindings'),
              
              fluidRow(
                box(width = 5, title = "Settings", status = 'info', solidHeader = T,
                    selectizeInput('peptide2', 'Peptide',
                      choices = unique(TCR_epitope$index_name)),
                 
                 sliderInput('activity', "normalized binding activity range", value = c(0.5, 1), 
                             min = 0, max = 1, step = 0.005)
                )
              ),
                fluidRow(
                  box(width = 9, title = "TCR-pMHC binding", status = 'primary', solidHeader = T,
                  plotlyOutput("plot5"))
              )
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
  
  epitope_sub <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide1, ]})
  
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
  
  TCR_epitope_peptide <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide2, ] %>%
      dplyr::filter((between(normalized_activity, 
                             input$activity[1], input$activity[2])))
  })
  
  output$plot5 <- renderPlotly({
    
    if (sum(TCR_epitope_peptide()$normalized_activity) == 0) {
      ggplot() +
        theme_void() +
        theme(axis.line = element_blank()) +
        ggtitle("Oops, there is no data for epitopes within your selected activity range")
      
    } else {
      
      plot_5.1 <- ggplot(data = TCR_epitope_peptide(),
                       aes(axis1 = tcr_name, axis2 = peptide, y = normalized_activity)) +
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
