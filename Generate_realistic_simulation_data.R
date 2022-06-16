set.seed(7)

# generate realistic simulation data
nGroup <- 3
nGene <- 1000 * 1
pCEG <- 0.1
pNeither <- 0.3
psDEG <- c(60, 80, 100)
nRep <- rep(10, nGroup)

data_TCC <- simulateReadCounts(Ngene = nGene, PDEG = sum(pSMG) / nGene,
                               DEG.assign = pSMG / sum(pSMG),
                               DEG.foldchange = rep(10, nGroup),
                               replicates = nRep)

data <- data_TCC$count
sDEG_grnd <- data_TCC$simulation$trueDEG

#############################################

data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data[, start : (start + rep - 1)]))
  start <- start + rep
}


# lst <- sort(cos_iCEG(data_grnd), index.return = TRUE, decreasing = TRUE)
# ind_CEG <- lst$ix[1 : (pCEG * nGene)]

ind_CEG <- which(cos_iCEG(data_grnd) >= 0.95)

CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1

lbl <- rep(0, nGene)
lbl[sDEG_grnd == 0] <- 'non-DEG'
lbl[sDEG_grnd == 1] <- 'sDEG'
lbl[ind_CEG] <- 'iCEG (labeled)'


DEG_grnd <- data_TCC$simulation$trueDEG

#############################################

lbl2 <- lbl
start <- 1
for (i in 1:nGroup) {
  lbl2[start : (start + psDEG[i] * nGene - 1)] <- paste0('DEG G', i)
  start <- start + (psDEG[i]* nGene)
}

#############################################

# same total count assumption

norm_factor <- 1/ totalcount(data)$norm_factor
data2 <- totalcount(data)$data

# super sample
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}
