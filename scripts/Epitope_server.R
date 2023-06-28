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
as.data.frame(TCR_epitope)

#create server
server4 <- function(input, output, session) {

  #Subset original df, only use columns of interest, make it into wide df format
  
  
  epitope_sub <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide, ]})
  
  smaller_epitope_sub <- reactive({epitope_sub()[ , c('peptide', "tcr_name", 
                                          "normalized_peptide_activity", "position")]})
  
  wide_epitope_sub <- reactive({pivot_wider(smaller_epitope_sub(), 
                             id_cols = c("peptide", "position"),
                             names_from = "tcr_name", 
                             values_from = "normalized_peptide_activity") %>%
    column_to_rownames(., var = 'peptide')}) #set peptides as row IDs
  
  #Set color squeme
  col_fun = reactive({colorRamp2(c(min(smaller_epitope_sub()$normalized_peptide_activity),
                         max(smaller_epitope_sub()$normalized_peptide_activity)), 
                       c("white", "red"))})
  
  #Draw interactive Headmap
 
  row_split <- reactive({wide_epitope_sub()[ ,1]})
  sub_df <- reactive({as.matrix(wide_epitope_sub()[, 2:ncol(wide_epitope_sub())])})
  
  ht1 = reactive({Heatmap(
      sub_df(), 
      name = "normalized peptide activity", 
      col = col_fun(), 
      cluster_rows = FALSE,
      show_row_names = FALSE,
      row_split = row_split(), 
      row_title_rot = 0,
      column_names_rot = 45)})
  
  makeInteractiveComplexHeatmap(input, output, session, ht1())
  
  
}


## Server for connection app
# create server for shiny app

server5 <- function(input, output, session) {
  
 TCR_epitope_peptide <- reactive({TCR_epitope[TCR_epitope$index_name == input$peptide, ] %>%
     dplyr::filter(TCR_epitope_peptide()$normalized_peptide_activity )
   })
  
  output$plot5 <- renderPlotly({
    ggplot(data = TCR_epitope_peptide(),
           aes(axis1 = tcr_name, axis2 = peptide, y = normalized_peptide_activity)) +
      geom_alluvium(aes(fill = tcr_name), curve_type = "cubic") +
      geom_stratum() +
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum))) +
      scale_x_discrete(limits = c("tcr_name", "peptide"),
                       expand = c(0.15, 0.05)) +
      theme_void()
    ggplotly()
  })
}
