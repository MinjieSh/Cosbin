
library("psych")
source('Cosbin_functions.R')



# Example data info:

# 48 samples separated into 4 groups, with number of samples in each group:
# 12, 12, 10, 14; 16845 genes totally.

# Cosbin conducts between-sample normalization, outputs group-specific weights,
# and normalizes the original data accordingly.



data <- data_for_normalization




#################################### Cosbin Normalization ######################################



## Setting parameters

nGroup <- 4 # number of groups
nGene <- 16845 # number of genes

nRep <- c(12,12,10,14) # number of samples in each group

data2 <- totalcount(data)$data


## Generating Super sample
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)])) 
  # data 3 will have a column number the same as your group number (4 in this example)
  start <- start + rep
}

## Cosbin functions
cosbin_out <- cosbin(data3)
# cosbin_out carries the group-specific weights

cosbin_out_full <- cosbin_convert(cosbin_out, data2, nGroup, nRep)
# cosbin_out_full carries the cosbin-normalized data


## COSBIN DONE


