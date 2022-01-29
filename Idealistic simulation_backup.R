library(Ternary)
library(ggtern)
library(hrbrthemes)
library(ggplot2)
library(DirichletReg)
library(DescTools)
library(TCC)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)

source('Cosbin_functions.R')

set.seed(7)

# Generate ideal simulation data
nGroup <- 3
nGene <- 1000
#pCEG <- 0.1
#pNeither <- 0.3
psDEG <-  c(0.06, 0.18, 0.36)
ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}

# non_sDEG <- rdirichlet(20000, rep(3, nGroup))
# data_CEG <-
#  non_sDEG[sample(which(cos_iCEG(non_sDEG) >= 0.99), nGene * pCEG),]

# data_Neither <-
  # non_sDEG[sample(intersect(which(cos_iCEG(non_sDEG) < 0.99),
                            # which(cos_iDEG(non_sDEG) < 0.99)),
                 # nGene * pNeither),]

non_DEG <- rdirichlet(400, rep(3, nGroup))

# iDEG
data_sDEG <- NULL
for (i in 1:nGroup) {
  data_sDEG <-
    rbind(data_sDEG, matrix(rep(ref[[i]], nGene * psDEG[i]), ncol = nGroup, byrow = TRUE))
}
# sDEG (iDEG + noise)
noise <- rnorm(nGroup * nGene * sum(psDEG), mean = 0, sd = 0.04)
data_sDEG <- abs(data_sDEG + noise)


# noise cannot be too large
# length(which(cos_iDEG(data_sDEG) < 0.99)) = 0

# data <- rbind(data_sDEG, data_Neither, data_CEG)
data <- rbind(data_sDEG, non_DEG)

# shuffle
# data <- data[sample(nrow(data)), ]


####################################################



# Add replicates
nRep <- rep(10, nGroup)

data2 <- NULL
for (i in 1:nGroup) {
  data2 <-
    cbind(data2, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}
dim(data2)
noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data2 <- abs(data2 + noise)


# #### modify the rowSums
for(i in 1:nGene){
   data2[i,] <- data2[i,] * 1000
  # data2[i,] <- data2[i,] * tccrowsum[i]
}

data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

ind_CEG <- which(cos_iCEG(data_grnd) >= 0.95)
ind_sDEG <- which(cos_iDEG(data_grnd) >= 0.99)
ind_neither <- intersect(which(cos_iCEG(data_grnd) < 0.95),
                         which(cos_iDEG(data_grnd) < 0.99))

ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(data_grnd, ref[[i]]) >= 0.99)
}

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'
lbl[ind_neither] <- '/'
lbl[ind_sDEG] <- 'sDEG'

lbl2 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
}

# Ground truth
CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1
sum(CEG_grnd)

sDEG_grnd <- rep(0, nGene)
sDEG_grnd[ind_sDEG] <- 1
sum(sDEG_grnd)

#######################################################

# check result in the ternary plot
plot_data <- data_grnd
label <- lbl2
plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point() +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") +
  scale_color_manual(values=c("blue", "green", "magenta", "orange", "red"))

# length(which(cos_iCEG(data) >= 0.99))
####################################################

#############################################
test <- data_grnd
rs <- rowSums(test)
#test <- test[-which(rs > 1000), ]
#coltemp <- lbl2[-which(rs > 1000)]
coltemp <- lbl2

coltemp[coltemp == 'sDEG G1'] <- "magenta"
coltemp[coltemp == 'sDEG G2'] <- "orange"
coltemp[coltemp == 'sDEG G3'] <- "red"
coltemp[coltemp == 'iCEG (labeled)'] <- "green4"
coltemp[coltemp == '/'] <- "steelblue1"

# Static chart
G1 <- test[,1]
G2 <- test[,2]
G3 <- test[,3]
plot3d( G1, G2, G3, col = coltemp, type = "p", r = 0.2)

#################

# assume same total count
factor <- 1 / totalcount(data2)$norm_factor
data2 <- totalcount(data2)$data

#### within group mean
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

# test <- data3
# rs <- rowSums(test)
# test <- test[-which(rs > 1000), ]
# 
# #library(plotly)
# fig <- plot_ly(as.data.frame(test), x = test[,1], y = test[,2], z = test[,3], size = 1)
# fig


plot_data <- data_grnd
label <- lbl2
plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
  geom_point() +
  labs( x       = "Group 1",
        y       = "Group 2",
        z       = "Group 3") + 
  scale_color_manual(values=c("blue", "green", "magenta", "orange", "red"))


ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("green4", "steelblue1", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )
dev.off()
###########################################################



length(which(cos_iCEG(data3) >= 0.99))
###########################################################
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


###############



###########################################################
# super sample
cosbin_out <- cosbin(data3)

# normaliztaion
totalcount_out <- totalcount(data3)$data
cor(colSums(data2), factor)

# DESeq2 super sample
cond <- factor(1:nGroup)
data3_int <- round(data3)
dds <- DESeqDataSetFromMatrix(data3_int, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
DESeq2_out <- data3
for(i in 1:nGroup){
  DESeq2_out[, i] <- DESeq2_out[, i] / sizeFactors(dds)[i]
}

# # DESeq2 duplicate
# cond <- factor(rep(1:nGroup, nRep))
# data2_int <- round(data2)
# dds <- DESeqDataSetFromMatrix(data2_int, DataFrame(cond), ~ cond)
# dds <- estimateSizeFactors(dds)
# sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
# DESeq2_temp <- data2
# for(i in 1:sum(nRep)){
#   DESeq2_temp[, i] <- DESeq2_temp[, i] / sizeFactors(dds)[i]
# }
# DESeq2_out <- NULL
# start <- 1
# for(rep in nRep){
#   DESeq2_out <- cbind(DESeq2_out, rowMeans(DESeq2_temp[, start : (start + rep - 1)]))
#   start <- start + rep
# }

####
normFactors <- edgeR::calcNormFactors(data3)
normFactors <- normFactors/geometric.mean(normFactors)
tmm_out <- data3
for(i in 1:nGroup){
  tmm_out[, i] <- tmm_out[, i] / normFactors[i]
}
cor(normFactors, factor)

# edgeR normalization - TMM
# normFactors <- edgeR::calcNormFactors(data2)
# totCounts <- colSums(data2)
# normFactors <- normFactors * totCounts
# normFactors <- normFactors/geometric.mean(normFactors)
# tmm_temp <- data2
# for(i in 1:sum(nRep)){
#   tmm_temp[, i] <- tmm_temp[, i] / normFactors[i]
# }
# tmm_out <- NULL
# start <- 1
# for(rep in nRep){
#   tmm_out <- cbind(tmm_out, rowMeans(tmm_temp[, start : (start + rep - 1)]))
#   start <- start + rep
# }
# cor(normFactors, factor)

#
# tcc <- new("TCC", data3, 1:nGroup)
# tcc <-
#   TCC::calcNormFactors(
#     tcc,
#     norm.method = "tmm",
#     test.method = "wad",
#     iteration = 1,
#     FDR = 0.1,
#     floorPDEG = 0.05
#   )
# tccEsts <- tcc$norm.factors
# totCounts <- colSums(data3)
# tccEsts <- tccEsts * totCounts
# tccEsts <- tccEsts / geometric.mean(tccEsts)
# 
# tcc_out <- getNormalizedData(tcc)


# TCC with replicates
# tcc <- new("TCC", data2, rep(1:nGroup, each = 10))
# tcc <-
#   calcNormFactors(
#     tcc,
#     norm.method = "tmm",
#     test.method = "edger",
#     iteration = 1,
#     FDR = 0.1,
#     floorPDEG = 0.05
#   )
# tccEsts <- tcc$norm.factors
# totCounts <- colSums(data2)
# tccEsts <- tccEsts * totCounts
# tccEsts <- tccEsts / geometric.mean(tccEsts)
# 
# # tcc_temp <- getNormalizedData(tcc)
# 
# tcc_temp <- data2
# for(i in 1:sum(nRep)){
#   tcc_temp[, i] <- tcc_temp[, i] / tccEsts[i]
# }
# 
# tcc_out <- NULL
# start <- 1
# for(rep in nRep){
#   tcc_out <- cbind(tcc_out, rowMeans(tcc_temp[, start : (start + rep - 1)]))
#   start <- start + rep
# }
# 
# cor(tccEsts, factor)

# TCC with replicates
tcc <- new("TCC", data3, 1:nGroup)
tcc <-
  calcNormFactors(
    tcc,
    norm.method = "tmm",
    test.method = "edger",
    iteration = 1,
    FDR = 0.1,
    floorPDEG = 0.05
  )
tccEsts <- tcc$norm.factors
tccEsts <- tccEsts / geometric.mean(tccEsts)

# tcc_temp <- getNormalizedData(tcc)

tcc_out <- data3
for(i in 1:nGroup){
  tcc_out[, i] <- tcc_out[, i] / tccEsts[i]
}


###########################################################
# check if move the center back

###########################################################

# Evaluation
# correlation of normalization factor


###########################################################
dir_out <- file.path(getwd(), "secondtry")
dir.create(dir_out, recursive = TRUE)

result <- data3
row.names(result) <- 1:nrow(result)
#threshold <- 0.686758
# cos(arccos(1/sqrt(3)) - arccos(0.95))
threshold <- 0.803434
count <- 1
plot_lbl <- lbl2
while(max(cos_iDEG(result)) >= threshold){
  temp_sDEG_ind <- which.max(cos_iDEG(result))
  result <- result[-temp_sDEG_ind, ]
  result <- apply(result, 2, function(x) x / sum(x))
  
  plot_lbl <- plot_lbl[-temp_sDEG_ind]
  # 
  # 
  #   plot_data <- result
  #   
  #   plot_data <-
  #     t(apply(plot_data, 1, function(x) x / sum(x)))
  #   
  #   p <- ggtern(data = as.data.frame(plot_data),
  #          aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
  #     geom_point(alpha = 0.3) +
  #     labs( x       = "G1",
  #           y       = "G2",
  #           z       = "G3") + 
  #     scale_color_manual(values=c("steelblue1", "green4", "red", "red", "red")) +
  #     #scale_color_manual(values=c("magenta", "orange", "red", "green4", "steelblue1")) +
  #     theme(
  #       title = element_text(size = 12, face = "bold"),
  #       axis.title = element_text(size = 12, face = "bold"),
  #       legend.text = element_text(size = 12),
  #       strip.background = element_blank(),
  #       strip.text.x = element_blank(),
  #       strip.text.y = element_blank()
  #     )
  #   
  #   fp <- file.path(dir_out, paste0("00", count, ".png"))
  #   ggsave(plot = p, 
  #          filename = fp, 
  #          device = "png",
  #          width = 8,
  #          height = 5
  #   )
  # 
  # 
  count <- count + 1
  
}
print(dim(result))
print(count)
# View(result)

# STEP 2 identify CEG, normalize based on CEG
# converge: # not change
count2 <- count + 1
while(length(which(cos_iCEG(result) >= 0.95)) != dim(result)[1]){
  temp_CEG_ind <- which(cos_iCEG(result) >= 0.95)
  result <- result[temp_CEG_ind,]
  result <- apply(result, 2, function(x) x / sum(x))
  # print(dim(result))
  plot_lbl <- plot_lbl[temp_CEG_ind]
  
  plot_data <- result
  
  plot_data <-
    t(apply(plot_data, 1, function(x) x / sum(x)))
  
  p <-  ggtern(data = as.data.frame(plot_data),
               aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
    geom_point(alpha = 0.3) +
    labs( x       = "G1",
          y       = "G2",
          z       = "G3") + 
    #scale_color_manual(values=c("steelblue1", "green4", "red", "red", "red")) +
    scale_color_manual(values=c("green4", "steelblue1", "red", "red", "red")) +
    theme(
      title = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      strip.background = element_blank(),
      strip.text.x = element_blank(),
      strip.text.y = element_blank()
    )
  
  fp <- file.path(dir_out, paste0("00", count2, ".png"))
  ggsave(plot = p, 
         filename = fp, 
         device = "png",
         width = 8,
         height = 5
  )
  
  
  count2 <- count2 + 1
  
  
}
print(dim(result))
# View(result)
data4 <- data3

# Final step normalization
ind <- as.numeric(row.names(result))
scalar <- colMeans(data4[ind, ] / result)
print(scalar)
for (i in 1:ncol(data4)) {
  data4[, i] <- data4[, i] / scalar[i]
}


plot_data <- data4

plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

p <-  ggtern(data = as.data.frame(plot_data),
             aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = lbl2)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "red", "red", "red")) +
  #scale_color_manual(values=c("green4", "steelblue1", "red", "red", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )

fp <- file.path(dir_out, paste0(count2 + 1, ".png"))
ggsave(plot = p, 
       filename = fp, 
       device = "png",
       width = 8,
       height = 5
)



#####################################
library(magick)

imgs <- list.files(dir_out, full.names = TRUE)
img_list <- lapply(imgs, image_read)

## join the images together
img_joined <- image_join(img_list)

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 4)

## view animated image
img_animated

## save to disk
image_write(image = img_animated,
            path = "animation.gif")

system("magick.exe convert -delay 20 images521/*.png try521.gif")
###########################################################

outlis <- list(totalcount_out, cosbin_out$data, tmm_out, DESeq2_out, tcc_out)

i <- 1
plot_data <- outlis[[i]]
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

lst <- sort(cos_iCEG(plot_data), index.return = TRUE, decreasing = TRUE)
ind <- lst$ix[1 : (pCEG * nGene)]
#length(ind)
#length(ind_CEG)
length(intersect(ind, ind_CEG)) / length(ind_CEG)


##############################

# ROC (DEG vs iCEG)
outlis <- list(totalcount_out, cosbin_out$data, tmm_out, DESeq2_out, tcc_out)
TPRlis <- list()
FPRlis <- list()
AUClis <- list()
for(i in 1:length(outlis)){
  lst <- sort(cos_iDEG(outlis[[i]]), index.return = TRUE, decreasing = TRUE)
  res <- NULL
  P <- sum(psDEG) * nGene
  N <- nGene - P 
  TPR <- NULL
  FPR <- NULL
  res <- NULL
  for (j in 1:nGene) {
    res <- c(res, sum(sDEG_grnd[lst$ix[1:j]]))
    TP <- sum(sDEG_grnd[lst$ix[1:j]] == 1)
    FP <- sum(sDEG_grnd[lst$ix[1:j]] == 0)
    TPR <- c(TPR, TP / P)
    FPR <- c(FPR, FP / N)
  }
  TPRlis[[i]] <- TPR
  FPRlis[[i]] <- FPR
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR[[i]], TPR[[i]])) / 0.15
}

colorlis <- c("red", "green", "blue", "orange", "pink")

plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlim = c(0, 0.1),
  ylim = c(0, 0.1),
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

# save.image(file='492021.RData')

###########################################################
# average of pairwise LFC MSE
ind <- ind_CEG
for(i in 1:length(outlis)){
  plot_data <- outlis[[i]]
  res <- NULL
  nGroups <- 3
  for (i in 1:(nGroups - 1)) {
    for (j in (i + 1):nGroups) {
      normCounts1 <- plot_data[ind, i]
      normCounts2 <- plot_data[ind, j]
      res <- c(res, mean(log(normCounts1 / normCounts2) ^ 2))
    }
  }
  print(Mean(res))
}

###########################################################
# cosine iCEG displacement
for(i in 1:length(outlis)){
  plot_data <- outlis[[i]]
  res <- NULL
  nGroups <- 3
  for (i in ind) {
    ref <- rep(1, nGroups)
    res <-
      c(res, sum(ref * plot_data[i,]) / sqrt(sum(ref ^ 2) * sum(plot_data[i,] ^ 2)))
  }
  print(Mean(res))
}

# cosine sDEG displacement

# DEG testing
TPRlis <- list()
FPRlis <- list()
AUClis <- list()
for(i in 1:length(outlis)){
  lst <- sort(cos_iDEG(outlis[[i]]), index.return = TRUE, decreasing = TRUE)
  res <- NULL
  P <- nGene - length(ind_CEG)
  N <- nGene - P 
  TPR <- NULL
  FPR <- NULL
  res <- NULL
  for (j in 1:nGene) {
    res <- c(res, sum(sDEG_grnd[lst$ix[1:j]]))
    TP <- sum(CEG_grnd[lst$ix[1:j]] == 0)
    FP <- sum(CEG_grnd[lst$ix[1:j]] == 1)
    TPR <- c(TPR, TP / P)
    FPR <- c(FPR, FP / N)
  }
  TPRlis[[i]] <- TPR
  FPRlis[[i]] <- FPR
  print(DescTools::AUC(FPR, TPR))
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR, TPR)) / 0.95
}
