#Load packages

library(shiny)
library(readxl)
library(ggplot2)
library(plotly)
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
TCR_epitope <- read_csv("TCR_Epitope_activity.csv")


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
    column(width = 9,
           InteractiveComplexHeatmapOutput()
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
             
             numericInput('activity', "threshold normalized binding activity", value = 2, min = 0, max = 15, step = 0.5)
           )
    ),
    column(width = 12,
           plotlyOutput("plot5"), 
    )
  ),
)
