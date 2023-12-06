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
                        menuItem('Peptide Binding Activity', tabName = 'heatmap'),
                        menuItem('Mutant-TCR Interaction', tabName = 'alluvium')
                      )),
    dashboardBody( 
      #css change of box colors
      tags$style(HTML("


                  .box.box-solid.box-info>.box-header {
                    color:#fff;
                    background:#5302a3
                                      }
                  
                  .box.box-solid.box-info{
                  border-bottom-color:#5302a3;
                  border-left-color:#5302a3;
                  border-right-color:#5302a3;
                  border-top-color:#5302a3;
                  }
                  
                  .box.box-solid.box-primary>.box-header {
                    color:#fff;
                    background:#b83289
                                      }
                  
                  .box.box-solid.box-primary{
                  border-bottom-color:#b83289;
                  border-left-color:#b83289;
                  border-right-color:#b83289;
                  border-top-color:#b83289;
                  }
                  
                                                      ")),
      tabItems(
        tabItem(tabName = 'about', 
        p("This application explores the vast diversity of T-cell receptor (TCR) 
        interaction with antigens, based on their epitopes."),
      
        p("TCR expression in human T-cells results from the somatic recombination of V and J genes. 
        This process leads to a theoretically extensive repertoire of TCRs that could potentially 
        recognize various antigens from pathogens. However, even with this theoretical diversity, 
        the human immune system would still not be able to create an immune response against every 
        existing pathogen. To overcome this limitation, T-cells are believed to be cross-reactive, 
        meaning that a single T-cell can recognize many different peptides presented on MHC molecules.
        The same antigens can also be recognized by multiple TCRs, based on distinct antigen features. 
        This cross-reactivity is critical for providing an adequate immune response to different pathogens."),
        
        p("However, while cross-reactivity is favorable in immune responses against pathogens, it can also 
        have unintended consequences. For instance, it might play a role in the development of auto-immune 
        diseases, or it can lead to off-target effects in cancer immunotherapies. Further understanding TCR 
        cross-reactivity might thus improve numerous applications, including predicting viral escape, 
        cancer neoantigen immunogenicity, autoimmunity, and off-target toxicity of T-cell-based therapies"),
  
        p("Existing computational methods are able to predict all TCRs that bind a given epitope. 
        But this application focuses on the opposite, predicting all epitopes that a given TCR binds. 
        It is build on a comprehensive experimental mutational scan database on TCR activation, visualizing
        the predictions of how peptide mutations affect TCR activation."),
      
        p("The application provides an overview over the similarity space of mutated epitopes,
        visualizes the binding activity of cognate epitopes and their single amino acid mutations and 
        illustrates the binding activity of an epitope and its mutations to one or multiple TCRs."),
    
        p("In summary, this application provides valuable insights into the cross-reactivity of TCRs
        and how it influences the immune response to various antigens."),
        ),
        
        tabItem(tabName = 'distances', 
                h1('Epitope similarity by diferent distance functions'),
                p("The", em("Epitope Cluster"), "function clusters mutated epitopes according to their mutual sequence
                  similarities, while showing their binding affinity. This information is based on user-defined
                  distance functions and weights assigned to mutations at specific positions within the peptides."),
                p(' '),
                
                fluidRow(
                  
                  box(width = 3, title = "Settings", status = 'info', solidHeader = T,
                      
                      selectizeInput('dist_method', 'Method',
                                     choices = dist_met),
                      
                      selectizeInput('TCR_names', 'TCR',
                                     choices = unique(TCR_epitope$tcr_name))
                  ),
                  
                  box(width = 9, title = "Position weights", status = 'info', solidHeader = T,
                      
                      div(style="display:inline-block", noUiSliderInput("obs1", "P1",
                                  min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                  width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs2", "P2",
                                  min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                  width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs3", "P3",
                                  min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                  width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs4", "P4",
                                  min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                  width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs5", "P5",
                                  min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                  width = '75px', height = '100px', color = '#0d0887',
                      )),
                      
                      div(style="display:inline-block", noUiSliderInput("obs6", "P6",
                                      min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                      width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs7", "P7",
                                      min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                      width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block",noUiSliderInput("obs8", "P8",
                                      min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                      width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs9", "P9",
                                      min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                      width = '75px', height = '100px', color = '#0d0887',
                      )),
                      div(style="display:inline-block", noUiSliderInput("obs10", "P10",
                                      min = 0, max = 1, value = 0.5, step=0.1, orientation = 'vertical',
                                      width = '75px', height = '100px', color = '#0d0887',
                      )),
                  )
                ),
                
                fluidRow(
                  box(width = 12, title = "Peptide Distances", status = 'primary', solidHeader = T,
                      plotlyOutput("plot3.1"))
                )
        ),
        
        tabItem(tabName = 'heatmap', 
                h1('TCR-pMHC normalized binding activity'),
                p("The", em("Peptide Binding Activity"), "tab displays a heatmapto visualize 
                  the normalized binding activity of a user-choosen petide and all its single amino acid 
                  mutations to one or multiple TCRs. All possible single amino acid mutations per peptide position
                  are displayed in alphabetical amino acid order as one block of the heatmap."),
                p(' '),
                
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
                h1('Alluvium plot displaing TCR-pMHC binding'),
                p("The", em("Mutant-TCR Interaction"),  "tab reveals which index peptides and their 
                  mutations are bound by multiple TCRs and which are recognized by only one TCR, emphasizing 
                  how different TCRs specific for the same index peptide recognize different mutants of it."),
                p(' '),
                
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

