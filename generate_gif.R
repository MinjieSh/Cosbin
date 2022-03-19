# those lines commented out are for ggtern (it's no longer availble in the newest ggplot2)
# dir_out <- file.path(getwd(), "cosbin_gif")
# dir.create(dir_out, recursive = TRUE)

result <- data3
plot_lbl <- lbl2

#########################################################

group <- c("sDEG G1", "sDEG G2", "sDEG G3", "iCEG (labeled)", "DEG (exclude sDEG) G1", "DEG (exclude sDEG) G2", "DEG (exclude sDEG) G3")
color_options <- c("magenta", "orange", "red", "green4", "steelblue1", "steelblue1", "steelblue1")

plot_ternary_data(result, plot_lbl, group, color_options)
plot_3d_data(result, plot_lbl, group, color_options)

row.names(result) <- 1:nrow(result)

#########################################################

threshold <- 0.728284 # cos(arccos(1/sqrt(3)) - arccos(0.98))

# Step 1: identify sDEG
count <- 1
while(max(cos_iDEG(result)) >= threshold){
  
  temp_SMG_ind <- which.max(cos_iDEG(result))
  result <- result[-temp_SMG_ind, ]
  result <- apply(result, 2, function(x) x / sum(x))
  
  plot_lbl <- plot_lbl[-temp_SMG_ind]
  
  if(count %% 10 == 1){
    plot_data <- result

    # p <- plot_ternay_data_deprecated(plot_data, plot_lbl, c("magenta", "orange", "red", "green4", "steelblue1"))
    # p <- plot_ternay_data_deprecated(plot_data, plot_lbl, c( "steelblue1","green4", "red", "red", "red"))
    
    
    # fp <- file.path(dir_out, paste0("0", (count-1) / 10, ".png"))
    # ggsave(plot = p,
    #        filename = fp,
    #        device = "png",
    #        width = 8,
    #        height = 5
    # )
    
    plot_ternary_data(plot_data, plot_lbl, group, color_options)
    
  }
  
  count <- count + 1
}

plot_3d_data(plot_data, plot_lbl, group, color_options)

# STEP 2 identify CEG, normalize based on CEG
# converge: # not change
count2 <- count + 1

while(length(which(cos_iCEG(result) >= 0.99)) != dim(result)[1]){
  
  temp_CEG_ind <- which(cos_iCEG(result) >= 0.99)
  result <- result[temp_CEG_ind,]
  result <- apply(result, 2, function(x) x / sum(x))
 
  plot_lbl <- plot_lbl[temp_CEG_ind]
  
  if (count2 %% 5 == 1) {
    plot_data <- result
    
    # p <- plot_ternay_data_deprecated(plot_data, plot_lbl, c("lightseagreen", "royalblue", "steelblue1", "green4", "red", "orange", "magenta"))
    # p <- plot_ternay_data_deprecated(plot_data, plot_lbl, c("magenta", "orange", "red", "green4", "steelblue1"))
    
    
    # fp <- file.path(dir_out, paste0(count2, ".png"))
    # ggsave(plot = p,
    #        filename = fp,
    #        device = "png",
    #        width = 8,
    #        height = 5
    # )
    
    plot_ternary_data(plot_data, plot_lbl, group, color_options)
    
  }
  
  count2 <- count2 + 1
  
}

plot_3d_data(plot_data, plot_lbl, group, color_options)

#########################################################

# step 3: normalization
data4 <- data3

ind <- as.numeric(row.names(result))
scalar <- colMeans(data4[ind, ] / result)
# print(scalar)

for (i in 1:ncol(data4)) {
  data4[, i] <- data4[, i] / scalar[i]
}

plot_data <- data4

# p <- plot_ternay_data_deprecated(plot_data, lbl2, c("lightseagreen", "royalblue", "steelblue1", "green4", "red", "orange", "magenta"))
# p <- plot_ternay_data_deprecated(plot_data, lbl2, c("magenta", "orange", "red", "green4", "steelblue1"))

# fp <- file.path(dir_out, paste0(count2 + 1, ".png"))
# ggsave(plot = p,
#        filename = fp,
#        device = "png",
#        width = 8,
#        height = 5
# )

plot_ternary_data(plot_data, lbl2, group, color_options)
plot_3d_data(plot_data, lbl2, group, color_options)

