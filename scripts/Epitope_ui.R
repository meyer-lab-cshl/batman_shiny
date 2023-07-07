#Load packages

library(shiny)
library(readxl)
library(ggplot2)
library(plotly)
library(heatmaply)
library(ggsankey)
library(ggalluvial)
library(reshape2)
library(dplyr)
library(devtools)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(tidyr)
library(tidyverse)
library(circlize)


#load in data frame
TCR_epitope <- read_csv("TCR_Epitope_activity_updated.csv")
TCR_epitope_distinct <- TCR_epitope %>% 
                       distinct(tcr_name, peptide, .keep_all = TRUE)


#create ui for Heatmap
ui4 <- fluidPage(
  fluidRow(headerPanel("TCR activity corresponding to different epitopes")
  ),
  
  
  fluidRow(
    column(width = 12,
           sidebarPanel(
             
             selectizeInput('peptide', 'Peptide',
                            choices = unique(TCR_epitope$index_name))
             #choice which peptide/epitope (by trivial name) is displayed
             #possible addidtion of choice by publication
           )
    ),
    column(width = 12,
           InteractiveComplexHeatmapOutput()
           #plotOutput('plot6')
    )
  ),
)

#create ui for connection plot

ui5 <- fluidPage(
  fluidRow(headerPanel("TCR binding to different epitopes")
  ),
  
  
  fluidRow(
    column(width = 12,
           sidebarPanel(
             selectizeInput('peptide', 'Peptide',
                            choices = unique(TCR_epitope$index_name)),
             
            sliderInput('activity', "normalized binding activity range", value = c(0.5, 1), 
                        min = 0, max = 1, step = 0.005)
           )
    ),
    column(width = 12,
           plotlyOutput("plot5")
           #uiOutput('warning')
    )
  ),
)
