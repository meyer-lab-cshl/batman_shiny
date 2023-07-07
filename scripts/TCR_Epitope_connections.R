#Load packages

library(readxl)
library(ggplot2)
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
amitavas_df <- read_excel("data/TCR_epitope_database_updated.xlsx")

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

#Safe df for further use
save(amitavas_df, file = "TCR_Epitope_activity_updated.Rda")
write.csv(amitavas_df, "TCR_Epitope_activity_updated.csv")

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
  column_to_rownames(., var = 'peptide') #set peptides as row IDs

  

#Set color squeme
col_fun = colorRamp2(c(min(smaller_SIN_df$normalized_peptide_activity),
                       max(smaller_SIN_df$normalized_peptide_activity)), 
                     c("white", "red"))
#Draw Headmap
ha = rowAnnotation(foo = anno_simple(1:20, 
                          pch = 1:20),
                   annotation_name_side = 'right')


heatmaply(as.matrix(wide_SIN_df[, 2:ncol(wide_SIN_df)]), 
        name = "normalized peptide activity", 
        col = col_fun,
        scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
          low = "white", high = "red", 
          limits = c(min(smaller_SIN_df$normalized_peptide_activity),
                    max(smaller_SIN_df$normalized_peptide_activity))),
        Rowv = FALSE,
        showticklabels = c(TRUE, FALSE),
        row_side_colors = wide_SIN_df$position,
        row_title_rot = 0,
        column_text_angle = 45
        )
        

#try to make sankey plot

tcr_sankey <- TCR_epitope[ ,c("tcr_name", "peptide", "index_name",
                              "index_peptide_activity", "normalized_peptide_activity", "peptide_activity")]

tcr_sankey_INS <- tcr_sankey[tcr_sankey$index_name == 'Insulin (ALWGPDPAAA)', ]
tcr_sankey_CMV <- tcr_sankey[tcr_sankey$index_name == 'CMV pp65 (NLVPMVATV)', ]
tcr_sankey_GAG <- tcr_sankey[tcr_sankey$index_name == 'Gag180–188 (TPQDLNTML)', ]
tcr_sankey_SIN <- tcr_sankey[tcr_sankey$index_name == 'OVA257-264 (SIINFEKL)', ]


ggplot(tcr_sankey_INS, aes(x = tcr_name, 
               next_x = peptide, 
               node = peptide_activity, 
               next_node = peptide_activity,
               fill = factor(peptide_activity),
               label = peptide_activity)) +
  geom_sankey() +
  geom_sankey_label() +
  theme_sankey(base_size = 16)

tcr_sankey_INS_small <- filter(tcr_sankey_INS, normalized_peptide_activity > 3.5)
tcr_sankey_CMV_small <- filter(tcr_sankey_CMV, normalized_peptide_activity > 0.5)
tcr_sankey_GAG_small <- filter(tcr_sankey_GAG, normalized_peptide_activity > 1.25)
tcr_sankey_SIN_small <- filter(tcr_sankey_SIN, between(normalized_peptide_activity, 5, 15))

ggplot(data = tcr_sankey_SIN_small,
       aes(axis1 = tcr_name, axis2 = peptide, y = normalized_peptide_activity)) +
  geom_alluvium(aes(fill = peptide), curve_type = "cubic") +
  geom_stratum(aes(fill = peptide)) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("tcr_name", "peptide"),
                   expand = c(0.15, 0.05)) +
  theme_void() +
  guides(fill = guide_legend(title = "Epitope"))
  
