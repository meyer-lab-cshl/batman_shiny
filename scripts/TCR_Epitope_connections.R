#Load packages

library(readxl)
library(ggplot2)
library(ggsankey)
library(reshape2)
library(dplyr)
library(devtools)
library(ComplexHeatmap)
library(tidyr)
library(tidyverse)
library(circlize)


#load in data frame
amitavas_df <- read_excel("data/TCR_epitope_database.xlsx")

#add column with Epitope names to df, make more efficient version!
amitavas_df$index_name <- NA
amitavas_df$index_name[amitavas_df$index_peptide == 'ALWGPDPAAA'] <- 
  'Insulin (ALWGPDPAAA)'
amitavas_df$index_name[amitavas_df$index_peptide == 'NLVPMVATV'] <- 
  'CMV pp65 (NLVPMVATV)'
amitavas_df$index_name[amitavas_df$index_peptide == 'TPQDLNTML'] <- 
  'Gag180–188 (TPQDLNTML)'
amitavas_df$index_name[amitavas_df$index_peptide == 'SLLMWITQC'] <- 
  'NY-ESO-1 (SLLMWITQC)'
amitavas_df$index_name[amitavas_df$index_peptide == 'SIINFEKL'] <-
  'OVA257-264 (SIINFEKL)'

#add column giving positions of peptide mutation (1-10)
amitavas_df$position <- mapply(function(x, y) which(x != y)[1], 
                               strsplit(amitavas_df$peptide, ""), 
                               strsplit(amitavas_df$index_peptide, "")) %>%
  replace_na(0)


#subset Amitavas df to one peptide
amitavas_df_SIN <- amitavas_df[amitavas_df$index_name == 'OVA257-264 (SIINFEKL)', ]

#make heatmap with Amitavas df
ggplot(amitavas_df_SIN, aes(tcr_name, peptide)) +
  geom_tile(aes(fill = peptide_activity), colour = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
       axis.text.y = element_blank())
 

#Try with ComplexHeatmap 

# First reduce df to columns TCR, Peptide and normalized_peptide_activity
# make a wide df with peptides as rows, tcrs as columns

smaller_SIN_df <- amitavas_df_SIN[ , c('peptide', "tcr_name", "normalized_peptide_activity", "position")]
wide_SIN_df <- pivot_wider(smaller_SIN_df, 
                           id_cols = c("peptide", "position"),
                           names_from = "tcr_name", 
                           values_from = "normalized_peptide_activity") %>%
  column_to_rownames(., var = 'peptide') %>% #set peptides as row IDs
  as.matrix() #Hatmap function needs matrix as input
  

#Set color squeme
col_fun = colorRamp2(c(min(smaller_SIN_df$normalized_peptide_activity),
                       max(smaller_SIN_df$normalized_peptide_activity)), 
                     c("white", "red"))
#Draw Headmap
ha = rowAnnotation(foo = anno_simple(1:20, 
                          pch = 1:20),
                   annotation_name_side = 'right')


Heatmap(wide_SIN_df[, 2:ncol(wide_SIN_df)], 
        name = "normalized peptide activity", 
        col = col_fun, 
        cluster_rows = FALSE,
        show_row_names = FALSE,
        row_split = wide_SIN_df[ ,1])


