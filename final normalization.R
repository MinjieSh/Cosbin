##########################################################
cosbin_out <- cosbin(data3)


# normalize the replicates
within_group_factor_grnd <- NULL
group_factor_grnd <- NULL

start <- 1
for(i in 1:nGroup){
  within_group_factor_grnd <- c(within_group_factor_grnd, factor[start : (start + nRep[i] - 1)] /  factor[start])
  group_factor_grnd <- c(group_factor_grnd,  mean(factor[start : (start + nRep[i] - 1)]))
  start <- start + nRep[i]
}

group_factor_grnd <- group_factor_grnd / group_factor_grnd[1]

f_grnd <- factor

############
# Given data2
# Find data_grnd

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

cor(f, factor)

###########################################

data_norm <- data2
for(i in 1:sum(nRep)){
  data_norm[, i] <- data_norm[, i] / f[i]
}
data_norm <-  t(apply(data_norm, 1, function(x) x / sum(x))) * mean(nRep)

data_norm_mean <- NULL
start <- 1
for(rep in nRep){
  data_norm_mean <- cbind(data_norm_mean, rowMeans(data_norm[, start : (start + rep - 1)]))
  start <- start + rep
}

# 100%
# data_grnd2 <- NULL
# start <- 1
# for(rep in nRep){
#   data_grnd2 <- cbind(data_grnd2, rowMeans(data_grnd[, start : (start + rep - 1)]))
#   start <- start + rep
# }

plot_data <- cosbin_out$data
plot_data <- data_norm_mean

plot_lbl <- lbl2
plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
  geom_point() +
  labs( x       = "Group 1",
        y       = "Group 2",
        z       = "Group 3") +
  scale_color_manual(values=c("blue", "green", "magenta", "orange", "red"))

length(which(cos_iCEG(plot_data) >= 0.99))
