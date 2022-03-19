setwd("E:/Cosbin")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TCC")

# if (!require('devtools')) install.packages('devtools')
# install_github('ms609/Ternary')

###### Prepare - install packages - R 3.6.3

#library(ggtern) # https://rpubs.com/tskam/ternary_plot
library(Ternary)
library(hrbrthemes)
library(DirichletReg)
library(DescTools)
library(TCC)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)
library(ggplot2)

source('Cosbin_functions.R')

set.seed(7)

####################################################

###### Generate ideal simulation data
###### 1000 genes, 600 SDEGs (corner, 60 + 180 + 360), iCEGs (center) determined by cosine
nGroup <- 3
nGene <- 1000
psDEG <-  c(0.06, 0.18, 0.36) ###### absolute percentage over all genes

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}

non_DEG <- rdirichlet(400, rep(3, nGroup))

###### iDEG: corner
data_sDEG <- NULL
for (i in 1:nGroup) {
  data_sDEG <-
    rbind(data_sDEG, matrix(rep(ref[[i]], nGene * psDEG[i]), ncol = nGroup, byrow = TRUE))
}

###### sDEG (iDEG + noise)
noise <- rnorm(nGroup * nGene * sum(psDEG), mean = 0, sd = 0.04)
data_sDEG <- abs(data_sDEG + noise)


###### noise cannot be too large
# length(which(cos_iDEG(data_sDEG) < 0.99)) = 0


data <- rbind(data_sDEG, non_DEG)

####################################################

###### Add replicates, here: same number for each group, but you can do different: e.g. nRep <- c(50, 100, 1000)
nRep <- rep(10, nGroup)

data2 <- NULL
for (i in 1:nGroup) {
  data2 <-
    cbind(data2, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}

# dim(data2) ## 1000 * 30; 10 replicates for each group
noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data2 <- abs(data2 + noise)


##### modify the rowSums, this constant 1000 can be any large integer
for(i in 1:nGene){
  data2[i,] <- data2[i,] * 1000
}

#### Ground truth data (super sample, to be compared with the normalized result)
#### data2 -> super sample: ground truth
#### data2 vs data3 (input for Cosbin)
data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

###### Add label
ind_CEG <- which(cos_iCEG(data_grnd) >= 0.95)
ind_sDEG <- which(cos_iDEG(data_grnd) >= 0.99)
ind_neither <- intersect(which(cos_iCEG(data_grnd) < 0.95),
                         which(cos_iDEG(data_grnd) < 0.99))

ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(data_grnd, ref[[i]]) >= 0.99)
}

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'
lbl[ind_neither] <- '/'
lbl[ind_sDEG] <- 'sDEG'

lbl2 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
}

###### Ground truth label

CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1
sum(CEG_grnd)

sDEG_grnd <- rep(0, nGene)
sDEG_grnd[ind_sDEG] <- 1
sum(sDEG_grnd)

####################################################### 3D plot

plot_3d_data(data_grnd, lbl2, c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "/"), c("magenta", "orange", "red", "green4", "steelblue1"))

#################

#### assume same total count

factor <- 1 / totalcount(data2)$norm_factor
data2 <- totalcount(data2)$data

#### within group mean

data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

plot_3d_data(data3, lbl2, c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "/"), c("magenta", "orange", "red", "green4", "steelblue1"))

plot_ternay_data(data_grnd, lbl2,c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "/"), c("magenta", "orange", "red", "green", "blue"))
plot_ternay_data(data3, lbl2,c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "/"), c("magenta", "orange", "red", "green", "blue"))

###########################################################

length(which(cos_iCEG(data3) >= 0.99))
