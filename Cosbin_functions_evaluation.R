# CEG dislocation: cosine_dislocation(data_grnd, CEG_grnd, norm_data_list)
# DEG dislocation (exclude aDEG): cosine_dislocation(data_grnd, DEG_grnd, norm_data_list)
cosine_dislocation <- function(data_grnd, grnd_label, norm_data_list){
  
  number_of_methods <- length(norm_data_list)
  ind <- which(grnd_label == 1)
  table <- matrix(0, nrow = number_of_methods, ncol = 4)
  
  for(j in 1:number_of_methods){
    norm_data <- norm_data_list[[j]]
    
    res <- NULL
    
    for (i in ind) {
      ref <- data_grnd[i, ]
      res <- c(res, mycosine(ref, norm_data[i,]))
    }
    
    # hist(res, xlim = c(0.85, 1), breaks = seq(0.85, 1, 0.001))
    table[j, 1] <- mean(res)
    table[j, 2] <- (180 * acos(Mean(res))) / pi
    table[j, 3] <- (180 * acos(min(res))) / pi
    table[j, 4] <- (180 * acos(max(res))) / pi
  }
  return(table)
}

# aDEG dislocation: groupwise_cosine_dislocation(data_grnd, norm_data_list, nGroup, pSMG)
groupwise_cosine_dislocation <- function(data_grnd, norm_data_list, nGroup, pSMG){
  
  number_of_methods <- length(norm_data_list)
  table <- array(0, c(number_of_methods, nGroup, 4))
  
  for(k in 1:number_of_methods){
    norm_data <- norm_data_list[[k]]
    
    for (j in 1:nGroup) {
      res <- NULL
      
      start <- 1
      for (i in start : (start + pSMG[j] - 1)) {
        ref <- data_grnd[i, ]
        res <- c(res, mycosine(ref, norm_data[i,]))
      }
      start <- start + (pSMG[j]* nGene)
      
      table[k, j, 1] <- mean(res)
      table[k, j, 2] <- (180 * acos(Mean(res))) / pi
      table[k, j, 3] <- (180 * acos(min(res))) / pi
      table[k, j, 4] <- (180 * acos(max(res))) / pi
    }
  }
  return(table)
}

# ROC (iCEG vs non-iCEG): ROC(norm_data_list, CEG_grnd, "iCEG", c("red", "green", "blue", "orange", "magenta"))
# ROC (DEG vs non DEG) : ROC(norm_data_list, DEG_grnd, "DEG", c("red", "green", "blue", "orange", "magenta"))
ROC <- function(norm_data_list, grnd_label, option = 'iCEG', color_list) {
  number_of_methods <- length(norm_data_list)
  
  TPR_list <- list()
  FPR_list <- list()
  AUC_list <- list()
  
  for (i in 1:number_of_methods) {
    if (option == "iCEG") {
      lst <-
        sort(cos_iCEG(norm_data_list[[i]]),
             index.return = TRUE,
             decreasing = TRUE)$ix
      P <- sum(grnd_label)
      N <- nGene - P
      limit <- nGene
      threshold <- 0.15
    } else {
      
      # value <-
      #   sort(cos_iCEG(norm_data_list[[i]]),
      #        index.return = FALSE,
      #        decreasing = FALSE)

      lst <-
        sort(cos_iCEG(norm_data_list[[i]]),
             index.return = TRUE,
             decreasing = FALSE)$ix
      
      #lst <- cbind(value,lst)
      lst <- lst[which(lst > 600)]
      
      limit <- 400
      P <- sum(grnd_label)
      N <- 400 - P
      threshold <- 1
    }
    
    res <- NULL
    TPR <- NULL
    FPR <- NULL
    res <- NULL
    
    for (j in 1:limit) {
      res <- c(res, sum(grnd_label[lst[1:j]]))
      TP <- sum(grnd_label[lst[1:j]] == 1)
      FP <- sum(grnd_label[lst[1:j]] == 0)
      TPR <- c(TPR, TP / P)
      FPR <- c(FPR, FP / N)
    }
    TPR_list[[i]] <- TPR
    FPR_list[[i]] <- FPR
    AUC_list[[i]] <-
      1 - (1 - DescTools::AUC(FPR, TPR)) / threshold # depending on the results
  }
  
  plot(
    FPR_list[[1]],
    TPR_list[[1]],
    xlab = "FPR",
    ylab = "TPR",
    xlim = c(0, threshold),
    ylim = c(0, 1),
    type = "l",
    lty = 2,
    col = color_list[1]
  )
  
  
  for (i in 2:length(norm_data_list)) {
    lines(FPR_list[[i]], TPR_list[[i]], lty = 2, col = color_list[i])
  }
  
  
  legend(
    "bottomright",
    legend = c("Total count", "Cosbin", "CSS", "Qsmooth","DESeq2","TMM/edgeR","DEGES/TCC"),
    col = color_list,
    lty = rep(2, length(norm_data_list)),
    cex = 0.65
  )
  
  return(list(FPR_list, TPR_list, AUC_list))
}
# average_of_pairwise_Log_Fold_Change_Mean_Squared_Error: average_of_pairwise_LFC_MSE(norm_data_list, CEG_grnd, nGroup)
average_of_pairwise_LFC_MSE <- function(norm_data_list, grnd_label, nGroup){
  number_of_methods <- length(norm_data_list)
  
  ind <- which(grnd_label == 1)
  table <- matrix(0, nrow = number_of_methods, ncol = 1)
  
  for(k in 1:number_of_methods){
    norm_data <- norm_data_list[[k]]
    res <- NULL
    for (i in 1:(nGroup - 1)) {
      for (j in (i + 1):nGroup) {
        normCounts1 <- norm_data[ind, i]
        normCounts2 <- norm_data[ind, j]
        res <- c(res, mean((log(normCounts1 / normCounts2) - log(1/1)) ^ 2))
      }
    }
    table[k,] <- mean(res)
  }
  return(table)
}
