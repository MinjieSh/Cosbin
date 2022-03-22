
library("psych")
source('Cosbin_functions.R')



# Example data info:

# 10 samples separated into 3 groups, with number of samples in each group:
# 4,3,3; 824 genes totally.

# Cosbin conducts between-sample normalization, outputs group-specific weights,
# and normalizes the original data accordingly.



data <- read.csv('dataPure.csv',row.names=1)




#################################### Cosbin Normalization ######################################



## Setting parameters

nGroup <- 3 # number of groups
nGene <- 824 # number of genes

nRep <- c(4,3,3) # number of samples in each group

data2 <- totalcount(data)$data # firstly use total-count method to normalize the original data


## Generating Super sample
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)])) 
  # data 3 will have a column number the same as your group number (3 in this example)
  start <- start + rep
}

## Cosbin functions
cosbin_out <- cosbin(data3)
# cosbin_out carries the group-specific weights

cosbin_out_full <- cosbin_convert(cosbin_out, data2)
# cosbin_out_full carries the cosbin-normalized data

write.csv(cosbin_out_full,'cosbin_normed.csv')

## COSBIN DONE


