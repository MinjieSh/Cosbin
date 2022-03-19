#library(ggtern)
library(Ternary)
library(hrbrthemes)
library(ggplot2)
library(DirichletReg)
library(DescTools)
library(TCC)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)

source('Cosbin_functions.R')

set.seed(17)

####################################################

nGroup <- 3
nGene <- 8 * 1000

############# 10% iCEG = center (1,1,1) + noise

piCEG <- 0.1 * nGene
data_iCEG <-
  matrix(rep(c(1, 1, 1) / 3, piCEG), ncol = nGroup, byrow = TRUE)
noise <- rnorm(nGroup * piCEG, mean = 0, sd = 0.05)
data_iCEG <- data_iCEG + noise


############# 225 SDEGs (corner, 50 + 75 + 100)

pSMG <-  c(50, 75, 100) 
ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}
data_SMG <- NULL
for (i in 1:nGroup) {
  data_SMG <-
    rbind(data_SMG, matrix(rep(ref[[i]], pSMG[i]), ncol = nGroup, byrow = TRUE))
}
noise <- rnorm(nGroup * sum(pSMG), mean = 0, sd = 0.05)
data_SMG <- abs(data_SMG + noise)

############# the rest: dirichlet

alpha <- list(c(10, 3.5, 3.5), c(3.5, 10, 3.5), c(3.5, 3.5, 10))
pDEG <- NULL
for (i in 1:nGroup) {
  pDEG[i] <- round((pSMG[[i]] / sum(pSMG)) * (nGene-piCEG-sum(pSMG)))
}
data_DEG <- NULL
for (i in 1:nGroup) {
  data_DEG <-
    rbind(data_DEG, rdirichlet(pDEG[i], alpha[[i]]))
}

############# 8000 Genes, 3 Groups

data <- rbind(data_SMG, data_DEG, data_iCEG)

####################################################

############# Add replicates
nRep <- rep(10, nGroup)

data2 <- NULL
for (i in 1:nGroup) {
  data2 <-
    cbind(data2, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}
############ dim(data2)  8000 Genes, 3 Groups, 10 replicates for each group

noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data2 <- abs(data2 + noise)


##### modify the rowSums (OPTIONAL)

### Cosbin rowsum
# for(i in 1:nGene){
#   data2[i,] <- data2[i,] * 1000
# }

### TCC rowsum
# data_TCC <- simulateReadCounts(Ngene = nGene, PDEG = 1 - (piCEG / nGene),
#                                DEG.assign = pSMG / sum(pSMG),
#                                DEG.foldchange = rep(10, nGroup),
#                                replicates = nRep)
# rs <- rowSums(data_TCC$count)
# SMG_ratio <- pSMG / sum(pSMG)
# move <- NULL
# for(i in 1:length(pSMG)){
#   if(i == 1){
#     offset <- 0
#   } else {
#     offset <- offset + SMG_ratio[i-1]
#   }
#   move <- c(move, floor(nGene * offset) + seq(1, pSMG[i]))
# }
# 
# temp <- rs[move]
# rs <- c(temp, rs[-move])
# # ind <- which(rs < 10000)
# ind <- 1:nGene
# for(i in ind){
#   data2[i,] <- data2[i,] * rs[i]
# }

############ data2 -> super sample -> ground truth

data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

################# use total count to alter the data

data2 <- totalcount(data2)$data

#### get super sample (average of each group), data3 is the final input
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}


################# get labels

lbl <- rep(0, nGene)
lbl[1:sum(pSMG)] <- 'sDEG'
lbl[(sum(pSMG) + 1) : (nGene-piCEG)] <- 'DEG (exclude sDEG)'
lbl[(nGene-piCEG+1) : nGene] <- 'iCEG (labeled)'

lbl2 <- lbl
start <- 1
for (i in 1:nGroup) {
  lbl2[start : (start + pSMG[[i]] - 1)] <- paste0('sDEG G', i)
  start <- start + pSMG[[i]]
}
for (i in 1:nGroup) {
  lbl2[start : (start + pDEG[[i]] - 1)] <- paste0('DEG (exclude sDEG) G', i)
  start <- start + pDEG[[i]]
}

# SMG_grnd <- rep(0, nGene)
# SMG_grnd[1:sum(pSMG)] <- 1
# sum(SMG_grnd) # 225
# 
# DEG_grnd <- rep(0, nGene)
# DEG_grnd[(sum(pSMG) + 1) : (nGene-piCEG)] <- 1
# sum(DEG_grnd) # 6975
# 
# CEG_grnd <- rep(0, nGene)
# CEG_grnd[(nGene-piCEG+1) : nGene] <- 1
# sum(CEG_grnd) # 800 (10%)

#############################################
group <- c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "DEG (exclude sDEG) G1", "DEG (exclude sDEG) G2", "DEG (exclude sDEG) G3")
color_options <- c("magenta", "orange", "red", "green4", "royalblue", "lightseagreen","slateblue")
plot_3d_data(data_grnd, lbl2, 
             group, color_options)

#################

plot_3d_data(data3, lbl2, 
             group, color_options)

plot_ternary_data(data_grnd, lbl2, group, color_options)
plot_ternary_data(data3, lbl2, group, color_options)

###########################################################

sum(CEG_grnd)
length(which(cos_iCEG(data3) >= 0.99))
length(which(cos_iCEG(data3) >= 0.98))

###########################################################

# ind_CEG <- which(cos_iCEG(data_grnd) >= 0.98)
# length(ind_CEG)
# length(ind_CEG) / nGene
# 
# ind_SMG <- which(cos_iDEG(data_grnd) >= 0.99)
# length(ind_SMG)
# length(ind_SMG) / nGene
# 
# ind_neither <- intersect(which(cos_iCEG(data_grnd) < 0.98),
#                          which(cos_iDEG(data_grnd) < 0.99))
# length(ind_neither)
# length(ind_neither) / nGene
# 
# ind_SMG2 <- list()
# for (i in 1:nGroup) {
#   ind_SMG2[[i]] <- which(cos_ref(data_grnd, ref[[i]]) >= 0.99)
# }
