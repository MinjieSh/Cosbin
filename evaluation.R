outlis <- list(totalcount_out, 
               cosbin_out$data, 
               tmm_out, 
               DESeq2_out, 
               tcc_out)


i <- 2
plot_data <- outlis[[i]]

##################################################################

# CEG Cosine Dislocation
ind <- which(CEG_grnd == 1)
table <- matrix(0, nrow = length(outlis), ncol = 4)

for(j in 1:length(outlis)){
  plot_data <- outlis[[j]]
  
  res <- NULL
  
  nGroups <- 3
  for (i in ind) {
    ref <- data_grnd[i, ]
    res <-
      c(res, sum(ref * plot_data[i,]) / sqrt(sum(ref ^ 2) * sum(plot_data[i,] ^ 2)))
  }
  # hist(res, xlim = c(0.85, 1), breaks = seq(0.85, 1, 0.001))
  table[j, 1] <- mean(res)
  table[j, 2] <- (180 * acos(Mean(res))) / pi
  table[j, 3] <- (180 * acos(min(res))) / pi
  table[j, 4] <- (180 * acos(max(res))) / pi
}

########################################################

# cosine non-iCEG (exclude sDEG) displacement

ind <- which(DEG_grnd == 1)
table <- matrix(0, nrow = length(outlis), ncol = 4)

result <- NULL
for(j in 1:length(outlis)){
  plot_data <- outlis[[j]]
  res <- NULL
  nGroups <- 3
  for (i in ind) {
    ref <- data_grnd[i, ]
    res <-
      c(res, sum(ref * plot_data[i,]) / sqrt(sum(ref ^ 2) * sum(plot_data[i,] ^ 2)))
  }
  # hist(res, xlim = c(0.85, 1), breaks = seq(0.85, 1, 0.001))
  table[j, 1] <- mean(res)
  table[j, 2] <- (180 * acos(Mean(res))) / pi
  table[j, 3] <- (180 * acos(min(res))) / pi
  table[j, 4] <- (180 * acos(max(res))) / pi
}

#########################################################

# cosine sDEG displacement
table <- array(0, c(length(outlis), nGroup, 4))

for(k in 1:length(outlis)){
  plot_data <- outlis[[k]]
  
  for (j in 1:nGroup) {
    res <- NULL
    start <- 1
    for (i in start : (start + pSMG[j] - 1)) {
      res <-
        c(res, sum(data_grnd[i, ] * plot_data[i,]) / sqrt(sum(data_grnd[i, ]  ^ 2) * sum(plot_data[i,] ^ 2)))
    }
    start <- start + (pSMG[j]* nGene)
    
    table[k, j, 1] <- mean(res)
    table[k, j, 2] <- (180 * acos(Mean(res))) / pi
    table[k, j, 3] <- (180 * acos(min(res))) / pi
    table[k, j, 4] <- (180 * acos(max(res))) / pi
  }
}

##############################

# ROC (iCEG vs non-iCEG)

TPRlis <- list()
FPRlis <- list()
AUClis <- list()
for(i in 1:length(outlis)){
  lst <- sort(cos_iCEG(outlis[[i]]), index.return = TRUE, decreasing = TRUE)
  res <- NULL
  P <- sum(CEG_grnd)
  N <- nGene - P
  TPR <- NULL
  FPR <- NULL
  res <- NULL
  for (j in 1:nGene) {
    res <- c(res, sum(CEG_grnd[lst$ix[1:j]]))
    TP <- sum(CEG_grnd[lst$ix[1:j]] == 1)
    FP <- sum(CEG_grnd[lst$ix[1:j]] == 0)
    TPR <- c(TPR, TP / P)
    FPR <- c(FPR, FP / N)
  }
  TPRlis[[i]] <- TPR
  FPRlis[[i]] <- FPR
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR, TPR)) / 0.15
}

colorlis <- c("red", "green", "blue", "orange", "magenta")

plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlab = "FPR",
  ylab = "TPR",
  xlim = c(0, 0.15),
  ylim = c(0, 1),
  type= "l",
  lty = 2,
  col = colorlis[1]
)


for(i in 2:length(outlis)){
  lines(FPRlis[[i]], TPRlis[[i]], lty=2, col = colorlis[i])
}


legend("bottomright",
       legend = c("Total count", "Cosbin", "TMM", "DESeq2",
                  "TCC"),
       col= colorlis,
       lty=rep(2, length(outlis)), cex=0.8)

#########################################################

# ROC (DEG (exclude SMG) vs non DEG (exclude SMG))

TPRlis <- list()
FPRlis <- list()
AUClis <- list()

for(i in 1:length(outlis)){
  temp <- outlis[[i]]
  lst <- sort(cos_iDEG(temp), 
              index.return = TRUE, decreasing = TRUE)$ix
  lst <- lst[lst > sum(pSMG)]
  
  res <- NULL
  P <- sum(DEG_grnd)
  N <- nGene - sum(pSMG) - P 
  TPR <- NULL
  FPR <- NULL
  res <- NULL
  for (j in 1:(N+P)) {
    res <- c(res, sum(DEG_grnd[lst[1:j]]))
    TP <- sum(DEG_grnd[lst[1:j]] == 1)
    FP <- sum(DEG_grnd[lst[1:j]] == 0)
    TPR <- c(TPR, TP / P)
    FPR <- c(FPR, FP / N)
  }
  TPRlis[[i]] <- TPR
  FPRlis[[i]] <- FPR
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR, TPR)) / 0.8
}

colorlis <- c("red", "green", "blue", "orange", "pink")

plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlab = "FPR",
  ylab = "TPR",
  xlim = c(0, 1),
  ylim = c(0, 1),
  type= "l",
  lty = 2,
  col = colorlis[1]
)


for(i in 2:length(outlis)){
  lines(FPRlis[[i]], TPRlis[[i]], lty=2, col = colorlis[i])
}


legend("bottomright", 
       legend = c("Total count", "Cosbin", "TMM", "DESeq2",
                  "TCC"),
       col= colorlis,
       lty=rep(2, length(outlis)), cex=0.8)

########################################################



###########################################################

# average of pairwise LFC MSE

ind <- which(CEG_grnd == 1)
table <- matrix(0, nrow = length(outlis), ncol = 1)

for(k in 1:length(outlis)){
  plot_data <- outlis[[k]]
  res <- NULL
  nGroups <- 3
  for (i in 1:(nGroups - 1)) {
    for (j in (i + 1):nGroups) {
      normCounts1 <- plot_data[ind, i]
      normCounts2 <- plot_data[ind, j]
      res <- c(res, mean(log(normCounts1 / normCounts2) ^ 2))
    }
  }
  table[k,] <- mean(res)
}

###########################################################

save.image(file='642021-2.RData')
