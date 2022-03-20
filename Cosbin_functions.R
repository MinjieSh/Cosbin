#' Normalize data (in super sample) using Cosbin algorithm
#'
#' @param data data (super sample: mean of each group).
#' @return a list contains normalized data and normalization factor.
#' @examples
#' cosbin(data)
cosbin <- function(data) {
  result <- data
  row.names(result) <- 1:nrow(result)
  
  # Step 1: identify sDEG
  # mycosine(c(1, 0,0), c(1, 1, 1)) = 0.5773503
  # solve 0.5773503 = 0.99 * x - sin(arccos(x)) * sin(arccos(0.99))
  
  #threshold <- 0.686758 #cos(arccos(0.5773503) - arccos(0.99))
  threshold <- 0.728284 #cos(arccos(0.5773503) - arccos(0.98))
  
  while (max(cos_iDEG(result)) >= threshold) {
    temp_sDEG_ind <- which.max(cos_iDEG(result))
    result <- result[-temp_sDEG_ind, ]
    result <- totalcount(result)$data
  }
  
  
  # STEP 2: identify CEG, normalize based on CEG
  # converge: # not change
  while (length(which(cos_iCEG(result) >= 0.98)) != dim(result)[1]) {
    temp_CEG_ind <- which(cos_iCEG(result) >= 0.98)
    result <- result[temp_CEG_ind,]
    result <- totalcount(result)$data
  }
  
  # Final step: normalization
  ind <- as.numeric(row.names(result))
  scalar <- colMeans(data[ind, ] / result)
  for (i in 1:ncol(data)) {
    data[, i] <- data[, i] / scalar[i]
  }
  
  return(list(data = data, norm_factor = scalar, CEG_index = ind))
}

#' convert the normalized data (in super sample) to the correct scale (with replicates)
#'
#' @param cosbin_out normalized data, output of the cosbin function
#' @param data2 original data (with replicates, not in super sample)
#' @return final normalization results
#' @examples
#' cosbin_convert(cosbin_out, original_data)
cosbin_convert <- function(cosbin_out, data2) {
  CEG <- data2[cosbin_out$CEG_index, ]
  CEG_norm <- totalcount(CEG)
  factor <- CEG_norm$norm_factor
  for (i in 1:dim(data2)[2]) {
    data2[, i] <- data2[, i] /  factor[i]
  }
  return(data2)
}


########################## helper functions
#' Calculate the cosine between vector A and vector B
mycosine <- function(A, B) {
  return(sum(A * B) / sqrt(sum(A ^ 2) * sum(B ^ 2)))
}

#' Calculate the cosine between data and reference vector
cos_ref <- function(data, ref) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}

#' Calculate the cosine between data and CEG reference, e.g. c(1, 1, 1)
cos_iCEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- rep(1, smp)
  
  res <- array(0, dim = gene)
  for (i in 1:gene) {
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}

#' Calculate the cosine between data and closest iDEG reference, e.g. c(1, 0, 0)
cos_iDEG <- function(data) {
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for (i in 1:smp) {
    ref[[i]] <- rep(0, smp)
    ref[[i]][i] <- 1
  }
  
  res <- array(0, dim = gene)
  
  for (i in 1:gene) {
    temp <- NULL
    for (j in 1:smp) {
      temp <- c(temp, mycosine(ref[[j]], data[i, ]))
    }
    res[i] <- max(temp)
  }
  return(res)
}

#' total count normalization
totalcount <- function(data) {
  scalar <- colSums(data) / mean(colSums(data))
  data <- apply(data, 2, function(x)
    x / sum(x)) * mean(colSums(data))
  
  return(list(data = data, norm_factor = scalar))
}

########################## visualization functions
# Example:
# group <- c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "DEG (exclude sDEG) G1", "DEG (exclude sDEG) G2", "DEG (exclude sDEG) G3")
# color_options <- c("magenta", "orange", "red", "green4", "royalblue", "lightseagreen","slateblue")
# plot_3d_data(data_grnd, lbl2, group, color_options)
plot_3d_data <- function(plot_data, label, group, color_options) {
  color <- label
  for (i in 1:length(group)) {
    color[color == group[i]] <- color_options[i]
  }
  
  G1 <- plot_data[, 1]
  G2 <- plot_data[, 2]
  G3 <- plot_data[, 3]
  plot3d(G1,
         G2,
         G3,
         col = color,
         type = "p",
         r = 0.2)
  # play3d(spin3d(axis = c(0, 0, 1), rpm = 20), duration = 10)

  # you have to create a subfolder under your working directory called 3dplot before running this
  # movie3d(
  #   movie="rotate",
  #   spin3d( axis = c(0, 0, 1), rpm = 7),
  #   duration = 10,
  #   type = "gif",
  #   clean = TRUE,
  #   dir = file.path(getwd(), "3dplot")
  # )
  
  dev.off()
}

# Example:
# ggtern is no longer availble in current version
# plot_ternary_data_deprecated(data_grnd, lbl2, c("blue", "green", "magenta", "orange", "red"))
plot_ternary_data_deprecated <-
  function(plot_data, label, color_options) {
    plot_data <-
      t(apply(plot_data, 1, function(x)
        x / sum(x)))
    
    p <- ggtern(data = as.data.frame(plot_data),
                aes(
                  x = plot_data[, 1],
                  y = plot_data[, 2],
                  z = plot_data[, 3],
                  color = label
                )) +
      geom_point() +
      labs(x = "Group 1",
           y = "Group 2",
           z = "Group 3") +
      scale_color_manual(values = color_options) +
      theme(
        title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank()
      )
    
    return(p)
    
  }


# Example:
# plot_ternary_data(data_grnd, lbl2, c("blue", "green", "magenta", "orange", "red"))
plot_ternary_data <-
  function(plot_data, label, group, color_options) {
    plot_data <-
      t(apply(plot_data, 1, function(x)
        x / sum(x)))
    
    color <- label
    for (i in 1:length(group)) {
      color[color == group[i]] <- color_options[i]
    }
    
    TernaryPlot()
    AddToTernary(points,
                 plot_data,
                 pch = 16,
                 cex = .8,
                 col = color)
    
    dev.off()
  }

########################## evaluation functions
# CEG dislocation: cosine_dislocation(data_grnd, CEG_grnd, norm_data_list)
# DEG dislocation (exclude sDEG): cosine_dislocation(data_grnd, DEG_grnd, norm_data_list)
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

# sDEG dislocation: groupwise_cosine_dislocation(data_grnd, norm_data_list, nGroup, pSMG)
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
# ROC (DEG (exclude SMG) vs non DEG (exclude SMG)) : ROC(norm_data_list, DEG_grnd, "DEG", c("red", "green", "blue", "orange", "magenta"))
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
      lst <-
        sort(cos_iDEG(norm_data_list[[i]]),
             index.return = TRUE,
             decreasing = TRUE)$ix
      lst <- lst[lst > sum(pSMG)]
      P <- sum(grnd_label)
      N <- nGene - sum(pSMG) - P
      limit <- N + P
      threshold <- 0.8
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
    legend = c("Total count", "Cosbin", "TMM", "DESeq2", "TCC"),
    col = color_list,
    lty = rep(2, length(norm_data_list)),
    cex = 0.8
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
