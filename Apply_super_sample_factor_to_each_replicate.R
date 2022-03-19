# data2
toy_example_data2 <- matrix(c(1 ,1.02, 1, 1, 0.94, 1,
                        1, 1.02, 1.04, 1, 1, 1.03,
                        7, 1, 0.98, 1, 0.97, 1, 
                        1, 1, 1, 3.02, 3, 2.98),
                      nrow = 4, ncol = 6, byrow = TRUE)
nRep <- c(1, 2, 3)

# data_grnd
toy_example_data_grnd <- NULL
start <- 1
for(rep in nRep){
  if(rep == 1){
    toy_example_data_grnd <- cbind(toy_example_data_grnd, toy_example_data2[, start])
  }else {
    toy_example_data_grnd <- cbind(toy_example_data_grnd, rowMeans(toy_example_data2[, start : (start + rep - 1)]))
  }
  start <- start + rep
}

# toy_example_data_grnd <- matrix(c(1,	1.01,	0.98,
#                                   1,	1.03,	1.01,
#                                   7,	0.99,	0.99,
#                                   1,	1.00,	3.00),
#                                 nrow = 4, ncol = 3, byrow = TRUE)

# data2 -> between sample total count
toy_example_data2 <- totalcount(toy_example_data2)$data

# toy_example_data2 <- matrix(c(
# 0.6,	1.514851,	1.492537,	0.9966777,	0.9543147,	0.9983361,
# 0.6,	1.514851,	1.552239,	0.9966777,	1.0152284,	1.0282862,
# 4.2,	1.485149,	1.462687,	0.9966777,	0.9847716,	0.9983361,
# 0.6,	1.485149,	1.492537,	3.0099668,	3.0456853,	2.9750416),
# nrow = 4, ncol = 6, byrow = TRUE)


# data3
toy_example_data3 <- NULL
start <- 1
for(rep in nRep){
  if(rep == 1){
    toy_example_data3 <- cbind(toy_example_data3, toy_example_data2[, start])
  }else {
    toy_example_data3 <- cbind(toy_example_data3, rowMeans(toy_example_data2[, start : (start + rep - 1)]))
  }
  start <- start + rep
}

# toy_example_data3 <- matrix(c(
# 0.6,	1.503694,	0.9831095,
# 0.6,	1.533545,	1.0133975,
# 4.2,	1.473918,	0.9932618,
# 0.6,	1.488843,	3.0102312),
# nrow = 4, ncol = 3, byrow = TRUE)


which(cos_iCEG(toy_example_data_grnd) >= 0.99) # iCEG 1, 2
which(cos_iDEG(toy_example_data_grnd) >= 0.98) # iDEG 3

which(cos_iCEG(toy_example_data3) >= 0.99) # None detected
which(cos_iDEG(toy_example_data3) >= 0.98) # None detected


##########################################################

cosbin_out <- cosbin(toy_example_data3)

# cosbin_out$data = matrix(c(
# 1.093292,	1.082547,	1.076706,
# 1.093292,	1.104037,	1.109878,
# 7.653043,	1.061110,	1.087825,
# 1.093292,	1.071855,	3.296819),
# nrow = 4, ncol = 3, byrow = TRUE)

# cosbin_out$CEG_index: 1, 2
# cosbin_out$norm_factor: 0.5488013, 1.3890342, 0.9130713

##########################################################

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
