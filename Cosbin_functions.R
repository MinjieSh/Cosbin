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
    result <- result[-temp_sDEG_ind,]
    result <- apply(result, 2, function(x)
      x / sum(x))
  }
  print(dim(result))
  
  # STEP 2: identify CEG, normalize based on CEG
  # converge: # not change
  while (length(which(cos_iCEG(result) >= 0.98)) != dim(result)[1]) {
    temp_CEG_ind <- which(cos_iCEG(result) >= 0.98)
    result <- result[temp_CEG_ind, ]
    result <- apply(result, 2, function(x)
      x / sum(x))
    # print(dim(result))
  }
  print(dim(result))
  # View(result)
  
  # Final step: normalization
  ind <- as.numeric(row.names(result))
  scalar <- colMeans(data[ind,] / result)
  print(scalar)
  for (i in 1:ncol(data)) {
    data[, i] <- data[, i] / scalar[i]
  }
  
  return(list(data = data, norm_factor = scalar))
}

#' convert the normalized data (in super sample) to the correct scale (with replicates)
#'
#' @param cosbin_out normalized data, output of the cosbin function
#' @param data2 original data (with replicates, not in super sample)
#' @param nGroup number of group, e.g. 3
#' @param nRep number of replicates in each group, e.g. c(10,20,30)
#' @return final normalization results
#' @examples
#' cosbin_convert(cosbin_out, original_data, 3, c(10, 20, 30))
cosbin_convert <- function(cosbin_out, data2, nGroup, nRep) {
  group_factor <- cosbin_out$norm_factor / cosbin_out$norm_factor[1]
  
  within_group_factor <- NULL
  within_group_factor_mean <- NULL
  
  start <- 1
  for (i in 1:nGroup) {
    temp <-
      colSums(data2[, start:(start + nRep[i] - 1)]) / sum(data2[, start])
    within_group_factor_mean <-
      c(within_group_factor_mean, mean(temp))
    within_group_factor <- c(within_group_factor, temp)
    start <- start + nRep[i]
  }
  
  f <- group_factor / within_group_factor_mean
  f <- f / f[1]
  
  f <- rep(f, nRep)
  f <- f * within_group_factor
  
  data_norm <- data2
  for (i in 1:sum(nRep)) {
    data_norm[, i] <- data_norm[, i] / f[i]
  }
  data_norm <-
    t(apply(data_norm, 1, function(x)
      x / sum(x))) * mean(nRep)
  
  return(data_norm)
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
    res[i] <- mycosine(ref, data[i,])
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
    res[i] <- mycosine(ref, data[i,])
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
      temp <- c(temp, mycosine(ref[[j]], data[i,]))
    }
    res[i] <- max(temp)
  }
  return(res)
}

#' total count normalization
totalcount <- function(data) {
  scalar <- colSums(data)
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
  
  # movie3d(
  #   movie="rotate",
  #   spin3d( axis = c(0, 0, 1), rpm = 7),
  #   duration = 10,
  #   type = "gif",
  #   clean = TRUE,
  #   dir = file.path(getwd(), "3dplot")
  # )
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
  }
