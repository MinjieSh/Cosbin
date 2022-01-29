########################## main function
cosbin <- function(data){
  result <- data
  row.names(result) <- 1:nrow(result)
  
  # Step 1 identify sDEG
  # mycosine(c(1, 0,0), c(1, 1, 1)) = 0.5773503
  # solve 0.5773503 = 0.99 * x - sin(arccos(x)) * sin(arccos(0.99))
  # cos(arccos(0.5773503) - arccos(0.99))
  
  #threshold <- 0.686758
  threshold <- 0.728284
  
  while(max(cos_iDEG(result)) >= threshold){
    temp_sDEG_ind <- which.max(cos_iDEG(result))
    result <- result[-temp_sDEG_ind, ]
    result <- apply(result, 2, function(x) x / sum(x))
  }
  print(dim(result))
  # View(result)
  
  # STEP 2 identify CEG, normalize based on CEG
  # converge: # not change
  while(length(which(cos_iCEG(result) >= 0.98)) != dim(result)[1]){
    temp_CEG_ind <- which(cos_iCEG(result) >= 0.98)
    result <- result[temp_CEG_ind,]
    result <- apply(result, 2, function(x) x / sum(x))
    # print(dim(result))
  }
  print(dim(result))
  # View(result)
  
  # Final step normalization
  ind <- as.numeric(row.names(result))
  scalar <- colMeans(data[ind, ] / result)
  print(scalar)
  for (i in 1:ncol(data)) {
    data[, i] <- data[, i] / scalar[i]
  }
  
  return(list(data = data, norm_factor = scalar))
}

########################## helper functions
mycosine <- function(A, B){
  return(sum(A*B)/sqrt(sum(A^2)*sum(B^2)))
}

cos_ref <- function(data, ref){
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  res <- array(0, dim = gene)
  for(i in 1:gene){
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}

cos_iCEG <- function(data){
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- rep(1, smp)
  
  res <- array(0, dim = gene)
  for(i in 1:gene){
    res[i] <- mycosine(ref, data[i, ])
  }
  return(res)
}


cos_iDEG <- function(data){
  gene <- dim(data)[1]
  smp <- dim(data)[2]
  
  ref <- vector("list", length = smp)
  for(i in 1:smp){
    ref[[i]] <- rep(0, smp)
    ref[[i]][i] <- 1
  }
  
  res <- array(0, dim = gene)
  
  for(i in 1:gene){
    temp <- NULL
    for(j in 1:smp){
      temp <- c(temp, mycosine(ref[[j]], data[i, ])) 
    }
    res[i] <- max(temp)
  }
  return(res)
}

totalcount <- function(data){
  scalar <- colSums(data)
  data <- apply(data, 2, function(x)
    x / sum(x)) * mean(colSums(data))
  
  return(list(data = data, norm_factor = scalar))
}
################################## super sample -> sample
cosbin_convert <- function(cosbin_out, data2, nGroup, nRep){
  group_factor <- cosbin_out$norm_factor / cosbin_out$norm_factor[1]
  
  within_group_factor <- NULL
  within_group_factor_mean <- NULL
  
  start <- 1
  for(i in 1:nGroup){
    temp <- colSums(data2[, start : (start + nRep[i] - 1)]) / sum(data2[, start])
    within_group_factor_mean <- c(within_group_factor_mean, mean(temp))
    within_group_factor <- c(within_group_factor, temp)
    start <- start + nRep[i]
  }
  
  f <- group_factor / within_group_factor_mean
  f <- f / f[1]
  
  f <- rep(f, nRep)
  f <- f * within_group_factor
  
  data_norm <- data2
  for(i in 1:sum(nRep)){
    data_norm[, i] <- data_norm[, i] / f[i]
  }
  data_norm <-  t(apply(data_norm, 1, function(x) x / sum(x))) * mean(nRep)
  
  return(data_norm)
}
