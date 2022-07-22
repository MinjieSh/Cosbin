#' Normalize data (in super sample) using Cosbin algorithm
#'
#' @param data data (super sample: mean of each group).
#' @return a list contains normalized data and normalization factor.
#' @examples
#' cosbin(data)
cosbin <- function(data) {
  result <- data
  row.names(result) <- 1:nrow(result)
  
  # Step 1: identify aDEG
  threshold <- 0.728284 #cos(arccos(0.5773503) - arccos(0.98))
  
  while (max(cos_iDEG(result)) >= threshold) {
    temp_sDEG_ind <- which.max(cos_iDEG(result))
    result <- result[-temp_sDEG_ind, ]
    result <- totalcount(result)$data
  }
  
  
  # STEP 2: identify iCEG, normalize based on iCEG
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
#' Low/high-expressed genes are filtered by their L2-norm ranks.
#' @param data Each row is a gene and each column is a sample.
#' @param thres.low The lower bound of percentage of genes to keep for Cosbin
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 0.05.
#' @param thres.high The higher bound of percentage of genes to keep for Cosbin
#'     with ranked norm. The value should be between 0 and 1.
#'     The default is 1.
data_cleaning <- function(data, thres.low = 0.05, thres.high = 1){
  if (thres.low >= thres.high) {
    stop("thres.low must be smaller than thres.high!")
  }
  
  sigNorm <- apply(data, 1, function(x) norm(matrix(x),"F") )
  Valid <- sigNorm >= quantile(sigNorm, thres.low) &
    sigNorm <= quantile(sigNorm, thres.high)
  data <- data[Valid,]
  return(data)
}
                   
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
