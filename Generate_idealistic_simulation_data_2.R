library(Ternary)
library(ggtern)
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

############# 10%

piCEG <- 0.1 * nGene
data_iCEG <-
  matrix(rep(c(1, 1, 1) / 3, piCEG), ncol = nGroup, byrow = TRUE)
noise <- rnorm(nGroup * piCEG, mean = 0, sd = 0.05)
data_iCEG <- data_iCEG + noise


#############

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

#############

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

#############

data <- rbind(data_SMG, data_DEG, data_iCEG)

####################################################

# Add replicates
nRep <- rep(10, nGroup)

data2 <- NULL
for (i in 1:nGroup) {
  data2 <-
    cbind(data2, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}
dim(data2)
noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data2 <- abs(data2 + noise)


# #### modify the rowSums

# for(i in 1:nGene){
#   data2[i,] <- data2[i,] * 1000
# }

data_TCC <- simulateReadCounts(Ngene = nGene, PDEG = 1 - (piCEG / nGene),
                               DEG.assign = pSMG / sum(pSMG),
                               DEG.foldchange = rep(10, nGroup),
                               replicates = nRep)
rs <- rowSums(data_TCC$count)
move <- c(seq(1,50), 1600 + seq(1,75), 4000 + seq(1, 100))
temp <- rs[move]
rs <- c(temp, rs[-move])
ind <- which(rs < 10000)
# ind <- 1:nGene
for(i in ind){
  data2[i,] <- data2[i,] * rs[i]
}




data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

#################

data2 <- totalcount(data2)$data

#### within group mean
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}


#################


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

SMG_grnd <- rep(0, nGene)
SMG_grnd[1:sum(pSMG)] <- 1
sum(SMG_grnd)

DEG_grnd <- rep(0, nGene)
DEG_grnd[(sum(pSMG) + 1) : (nGene-piCEG)] <- 1
sum(DEG_grnd)

CEG_grnd <- rep(0, nGene)
CEG_grnd[(nGene-piCEG+1) : nGene] <- 1
sum(CEG_grnd)
#######################################################


#############################################

test <- data_grnd
coltemp <- lbl2

coltemp[coltemp == 'sDEG G1'] <- "magenta"
coltemp[coltemp == 'sDEG G2'] <- "orange"
coltemp[coltemp == 'sDEG G3'] <- "red"
coltemp[coltemp == 'iCEG (labeled)'] <- "green4"

coltemp[coltemp == 'DEG (exclude sDEG) G1'] <- "royalblue"
coltemp[coltemp == 'DEG (exclude sDEG) G2'] <- "lightseagreen"
coltemp[coltemp == 'DEG (exclude sDEG) G3'] <- "slateblue"

G1 <- test[,1]
G2 <- test[,2]
G3 <- test[,3]
plot3d( G1, G2, G3, col = coltemp, type = "p", r = 0.2)


#################

plot_data <- data_grnd
plot_data <- data3



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
