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
all_tcrs[ ,'binder'] <- NA

##loop over all_tcrs to add binding information; might not be the most efficent for bigger df's.
for (i in 1:nrow(all_tcrs)){
  if (all_tcrs$epitope[i] %in% nSB){
    all_tcrs$binder[i] <- 'SB'
  } else if (all_tcrs$epitope[i] %in% nNB){
    all_tcrs$binder[i] <- 'NB'
  }
}
