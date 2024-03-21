
library("psych")
source('Cosbin_functions.R')



# Example data info:

# 10 samples separated into 3 groups, with number of samples in each group:
# 4,3,3; 824 genes totally.

# Cosbin conducts between-sample normalization, outputs group-specific weights,
# and normalizes the original data accordingly.


# Input: dataPure.csv
# Output: cosbin_normed.csv


data <- read.csv('dataPure.csv',row.names=1)




#################################### Cosbin Normalization ######################################



## Setting parameters

nGroup <- 3 # number of groups
nGene <- 824 # number of genes

nRep <- c(4,3,3) # number of samples in each group

data_tc <- totalcount(data)$data # firstly use total-count method to normalize the original data


## Generating Super Sample (basically the group mean)
data_ss <- NULL
start <- 1
for(rep in nRep){
  data_ss <- cbind(data_ss, rowMeans(data_tc[, start : (start + rep - 1)])) 
  # data_ss will have a column number the same as your group number (3 in this example)
  start <- start + rep
}

## Cosbin functions
cosbin_out <- cosbin(data_ss)
# cosbin_out carries the group-specific weights/scalers

cosbin_out_full <- cosbin_convert(cosbin_out, data_tc)
# cosbin_out_full carries the cosbin-normalized data

write.csv(cosbin_out_full,'cosbin_normed.csv')

## COSBIN DONE


