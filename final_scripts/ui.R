#.libPaths("/opt/R/4.2.1/lib/R/library")
.libPaths(c("/usr/lib64/R/library","/usr/share/R/library"))
## Libraries ####
library(plotly)
library(shiny)
library(shinydashboard)
library(readxl)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyverse)
library(InteractiveComplexHeatmap)
library(circlize)
library(igraph)

# load dataframe for TCR distances and epitope binding information
TCR_epitope <- read_xlsx("TCR_epitope_database.xlsx")
dist_met <- c("BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM10")


## Load Browser functions and annotations ####
#source("browser.R")

## Create dashboard ####
ui <- shinyUI(
  dashboardPage(skin = "purple",
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
                "This application explores the vast diversity of T-cell receptor (TCR) 
        interaction with antigens, based on their epitopes." ,br(),
                "TCR expression in human T-cells results from the somatic recombination of V and J genes. 
        This process leads to a theoretically extensive repertoire of TCRs that could potentially 
        recognize various antigens from pathogens. However, even with this theoretical diversity, 
        the human immune system would still not be able to create an immune response against every 
        possible pathogen. To overcome this limitation, T-cells are believed to be cross-reactive, 
        meaning that a single T-cell can recognize many different peptides presented on MHC molecules,
        according to experimental extrapolations up to 10e6 different peptides. 
        The same antigens can also be recognized by multiple TCRs, based on distinct antigen features. 
        This cross-reactivity is critical for providing an adequate immune response to different pathogens."
                , br(),
                "However, while cross-reactivity is favorable in immune responses against pathogens, it can also 
        have unintended consequences. For instance, it may play a role in the development of auto-immune 
        diseases, or it can lead to off-target effects in cancer immunotherapies.", br(), br(),
                "This App visualizes the cross-reactivity of 48 TCRs to 5 specific index peptides and all
        possible single mutations of these peptides.", br(),
                em("The Epitope Cluster"), "function clusters mutated epitopes according to their mutual sequence
        similarities, while showing their binding affinity. This information is based on user-defined
        distance functions and weights assigned to mutations at specific positions within the peptides."
                , br(),
                em("The Activity Heatmap"), "in the app allows users to visualize the normalized binding activity of 
        each TCR to different peptides.", br(),
                "On the other hand, the", em("Mutant-TCR Interaction"),  "tab reveals which index peptides and their 
        mutations are bound by multiple TCRs and which are recognized by only one TCR, emphasizing 
        how different TCRs specific for the same index peptide recognizes different mutants of it 
        differently.", br(), 
                "In summary, this application provides valuable insights into the cross-reactivity of TCRs
        and how it influences the immune response to various antigens."
        ),
        
        tabItem(tabName = 'distances', 
                h1('Epitope similarity by diferent distance functions'),
                fluidRow(
                  
                  box(width = 5, title = "Settings", status = 'info', solidHeader = T,
                      
                      selectizeInput('dist_method', 'Method',
                                     choices = dist_met),
                      
                      selectizeInput('TCR_names', 'TCR',
                                     choices = unique(TCR_epitope$tcr_name))
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
                                     choices = unique(TCR_epitope$index_peptide))
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
                                     choices = unique(TCR_epitope$index_peptide)),
                      
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
