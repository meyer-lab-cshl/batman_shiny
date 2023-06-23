#Load packages

library(readxl)
library(ggplot2)
library(ggsankey)
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
as.data.frame(TCR_epitope)

#create ui
ui4 <- fluidPage(
  fluidRow(headerPanel("TCR activity corresponding to different epitopes")
  ),
  
  
  fluidRow(
    column(width = 9,
           sidebarPanel(
             
             selectizeInput('peptide', 'Peptide',
                            choices = unique(TCR_epitope$index_name))
             #choice which peptide/epitope (by trivial name) is displayed
             #possible addidtion of choice by publication
           )
    ),
    column(width = 9,
           plotOutput("heatmap"), 
           # hover = hoverOpts("seq_hover"),
           # verbatimTextOutput("peptide_seq")
    )
  ),
)
