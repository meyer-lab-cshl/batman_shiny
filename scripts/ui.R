# Load needed packages

library(R.matlab)
library(shiny)
library(ggplot2)
library(plotly)

# read in data of calculated distances

dist_data = readMat("data//Methods_TCRs_XY.mat")
dist_data = data.frame(dist_data)

# generate Lists of distance functions and TCRs for Input windows

dist_met <- c("Atchley", "BLOSUM100", "Dayhoff", "Gonnet", "Hamming", "PAM50",
              "R Contigous")

TCR_name <- c("868", "A42", "A6", "B7", "E7NLV",  "G10",  "ILA-1", "T5-004")



ui1 <- fluidPage(
  fluidRow(headerPanel("Strong and non binding epitopes to different TCRs, using multiple distance functions")
           ),
  
  
  fluidRow(
    column(width = 9,
           sidebarPanel(
              selectInput('dist_method', 'Method',
                          choices = dist_met),
              
              selectInput('TCR_names', 'TCR',
                          choices = TCR_name)
              )
           ),
    column(width = 9,
           plotlyOutput("plot1")
           )
    ),
)

