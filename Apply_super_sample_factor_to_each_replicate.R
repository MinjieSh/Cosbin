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

# First remove: 4.2,	1.473918,	0.9932618
# sample 1  sample 2  sample 3
# 1.259202 1.255029 0.7417606
# 1.259202 1.279944 0.7646130
# 1.259202 1.242634 2.2712332

# second remove: 1.259202 1.242634 2.2712332
# sample 1  sample 2  sample 3
# 1.093292 1.082547 1.076706
# 1.093292 1.104037 1.109878

# output
# cosbin_out$data = matrix(c(
# 1.093292,	1.082547,	1.076706,
# 1.093292,	1.104037,	1.109878,
# 7.653043,	1.061110,	1.087825,
# 1.093292,	1.071855,	3.296819),
# nrow = 4, ncol = 3, byrow = TRUE)

# Ground truth
# 1,	1.01,	0.98
# 1,	1.03,	1.01
# 7,	0.99,	0.99
# 1,	1.00,	3.00

# cosbin_out$CEG_index: 1, 2
# cosbin_out$norm_factor: 0.5488013, 1.3890342, 0.9130713

# (input - super sample)
# 0.6,	     1.503694,	   0.9831095 
# 0.6,	     1.533545,	   1.0133975 
# 4.2,	     1.473918,	   0.9932618 
# 0.6,	     1.488843,	   3.0102312
#/0.5488013  /1.3890342    /0.9130713 (normalization factor)

# (input - raw - with replicates)
# 0.6,	1.514851,	1.492537,	0.9966777,	0.9543147,	0.9983361
# 0.6,	1.514851,	1.552239,	0.9966777,	1.0152284,	1.0282862
# 4.2,	1.485149,	1.462687,	0.9966777,	0.9847716,	0.9983361
# 0.6,	1.485149,	1.492537,	3.0099668,	3.0456853,	2.9750416
# /0.5428227 /1.3704929 /1.3773113 /0.9016988 /0.8909272 /0.9167471 (normalization factor)

# (ground truth - with replicates)
# 1,  1.02, 1,    1,    0.94, 1
# 1,  1.02, 1.04, 1,    1,    1.03
# 7,  1,    0.98, 1,    0.97, 1 
# 1,  1,    1,    3.02, 3,    2.98

# Final result (with replicates)
#	1.105333	1.105333	1.083660	1.105333	1.071148	1.088998
#	1.105333	1.105333	1.127007	1.105333	1.139519	1.121668
#	7.737333	1.083660	1.061987	1.105333	1.105333	1.088998
#	1.105333	1.083660	1.083660	3.338107	3.418557	3.245215

##########################################################

toy_example_result <- cosbin_convert(cosbin_out, toy_example_data2)
