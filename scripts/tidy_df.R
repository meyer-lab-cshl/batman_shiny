##load packages

library(R.matlab)
library(tidyr)
library(dplyr)

##load data, set names to distance functions

dist_data = readMat("data//Methods_TCRs_XY.mat")
names(dist_data) <- gsub("(.*)XY$", "\\1", names(dist_data))


# define function to extract distances for all epitopes for one TCR and distance funtion
## index corresponds to the chosen distance function (1 = Atchley, 2 = BLOSUM100, etc.)

tcr_epitope_table <- function(index, dist_data){
  picked_tcr <- dist_data[[index]]
  tcr_name <- dimnames(dist_data)[[1]][[index]]
   
  rownames(picked_tcr) <- c('x', 'y')
  colnames(picked_tcr) <- paste0(1:ncol(picked_tcr))
  
 picked_tcr <-  picked_tcr %>%
    t %>% 
    as_tibble %>%
    mutate(tcr=tcr_name)
 
 picked_tcr$epitope <- paste0(1:nrow(picked_tcr)) #add column with the number of epitopes
  
  return(picked_tcr)
}


## function to run tcr_epitope_table over all distance function lists in the dist_data

all_tcrs <- lapply(seq_along(dist_data), function(distfn_index) {
  sub_data <- dist_data[[distfn_index]] #subset for one distance function (by index given via seq_along)
  distname <- names(dist_data)[distfn_index]
  dist_tcrs <- lapply(seq_along(sub_data), tcr_epitope_table,
                         dist_data=sub_data) %>%
    bind_rows #combine all sub dfs
  dist_tcrs$dist <- distname
  return(dist_tcrs)
}) %>%
  bind_rows


#Add information about binders to new data frame

## create lists for binding information (SB = strong binder, NB = non binder)
nSB = c(13,21,29,35,27,30,15,13)
nNB = c(72,58,50,51,97,2683,121,146)

##add empty column about binding info to tcrs data frame
all_tcrs[ ,'Binding'] <- 'Weak Binder'

##loop over all_tcrs to add binding information; might not be the most efficent for bigger df's.
for (i in 1:nrow(all_tcrs)){
  if (all_tcrs$epitope[i] %in% nSB){
    all_tcrs$Binding[i] <- 'Strong Binder'
  } else if (all_tcrs$epitope[i] %in% nNB){
    all_tcrs$Binding[i] <- 'Non Binder'
  }
}

#safe new all_tcrs data frame
save(all_tcrs, file = "All_TCRs_epitopes.Rda")

##add column with epitope sequences

#load sequence list
seqs = readMat("data//epitope_seqs.mat")

# epitope_seq <- data.frame(seqs) %>%
#   t
# 
# df1 <- data.frame()
# 
# epitope_seq1 <- lapply(seq_along(epitope_seq), function(tcr_index){
#   
#   tcr <- epitope_seq[ , tcr_index] %>%
#     unlist
#   
#   df1$colnames(tcr) <- data.frame(tcr)
#   
#   
#   return(df1)
#   })
# 
# tmd_868 <- tmd[ , 1] %>%
#   unlist
# 
# tmd_test <- seqs[ , 1] %>%
#   unlist
# 
# df1 <- data.frame(tmd_868)
# 
# df1$col.name(tcr)

##new try of adding column containing sequences to all_tcrs

get_tcr_seq <- function(list_position){
  
  Epitope_Seq <- epitope_seq[ , list_position] %>%
    unlist
  Epitope_Seq <- data.frame(Epitope_Seq)
  
  Epitope_Seq <- Epitope_Seq %>%
    mutate(tcr = colnames(epitope_seq)[list_position]) %>%
    mutate(epitope = 1:n())
  
  return(Epitope_Seq)
  
}

#epitope_seqs_df <- lapply(list(1:8), get_tcr_seq)

#ugly manual version

df_tcr868 <- get_tcr_seq(1)
df_tcrA42 <- get_tcr_seq(2)
df_tcrA6 <- get_tcr_seq(3)
df_tcrB7 <- get_tcr_seq(4)
df_tcrE7NLV <- get_tcr_seq(5)
df_tcrG10 <- get_tcr_seq(6)
df_tcrILA1 <- get_tcr_seq(7)
df_tcrT5004 <- get_tcr_seq(8)

epitope_seqs_df <- rbind(df_tcr868, df_tcrA42, df_tcrA6, df_tcrB7, 
                         df_tcrE7NLV, df_tcrG10, df_tcrILA1, df_tcrT5004)

all_tcrs <- all_tcrs %>%
  left_join(epitope_seqs_df, join_by(tcr == tcr, epitope == epitope))
#Error in `left_join()`:
# ! Can't join `x$epitope` with `y$epitope` due to incompatible types.
# ℹ `x$epitope` is a <character>.
# ℹ `y$epitope` is a <integer>.
