# generate simulation data
nGroup <- 3
nGene <- 1000
pCEG <- 0.1
pNeither <- 0.3
psDEG <- c(0.06, 0.18, 0.36)
nRep <- rep(10, nGroup)
data_TCC <- simulateReadCounts(Ngene = nGene, PDEG = sum(psDEG),
                          DEG.assign = psDEG / sum(psDEG),
                          DEG.foldchange = rep(15, nGroup),
                          replicates = nRep)

data <- data_TCC$count
sDEG_grnd <- data_TCC$simulation$trueDEG

lst <- sort(cos_iCEG(data), index.return = TRUE, decreasing = TRUE)
ind_CEG <- lst$ix[1 : (pCEG * nGene)]
CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1

lbl <- rep(0, nGene)
lbl[sDEG_grnd == 0] <- 'non-sDEG'
lbl[sDEG_grnd == 1] <- 'sDEG'
lbl[ind_CEG] <- 'iCEG'

#####################################

lbl2 <- lbl
start <- 1
for (i in 1:nGroup) {
  lbl2[start : (start + psDEG[i] * nGene - 1)] <- paste0('sDEG G', i)
  start <- start + (psDEG[i]* nGene)
}
#############################################3

# same total count assumption

data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data[, start : (start + rep - 1)]))
  start <- start + rep
}

norm_factor <- 1/ totalcount(data)$norm_factor
data2 <- totalcount(data)$data

# super sample
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

tccrowsum <- rowSums(data3)
# # check if there is a dislocation
# plot_data = data3
# plot_lbl = lbl
# 
# plot_data <-
#   t(apply(plot_data, 1, function(x) x / sum(x)))
# 
# ggtern(data = as.data.frame(plot_data),
#        aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = plot_lbl)) +
#   geom_point() +
#   labs( x       = "Group 1",
#         y       = "Group 2",
#         z       = "Group 3") + 
#   scale_color_manual(values=c("green", "blue", "red"))

# check if there is a dislocation
plot_data = outlis[[4]]
label = lbl2

plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

emf(file ="tcc1.emf", width = 6.4, height = 5)
ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point() +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("green", "blue", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )
dev.off()

emf(file ="tcc2.emf", width = 6.4, height = 5)
ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point() +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("green", "blue", "magenta", "orange", "red")) + theme_nolabels()+ theme(
    title = element_text(
      size = 12,
      face = "bold",
      colour = "white"
    ),
    axis.title = element_text(
      size = 12,
      face = "bold",
      colour = "white"
    ),
    axis.text.x = element_text(color = "white"),
    axis.text.y = element_text(color = "white"),
    legend.text = element_text(size = 12, colour = "white"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank())
dev.off()


###################################################
totalcount_out <- totalcount(data3)$data

# testMat <- matrix(c(1489, 906, 22, 13, 793, 410, 76, 42, 521, 1196), 
#                   nrow = 5, ncol = 2, byrow = TRUE)
# data2 <- testMat
# cond <- factor(c("A","B"))

# cosbin
cosbin_out <- cosbin(data3)
# cosbin_out_full <- cosbin_convert(cosbin_out, data2, nGroup, nRep)

# DESeq2 super sample
# only deal with house keeping genes by log(0) = -Inf (remove them)
# geometric mean removes outlier read counts
# median also downplay high read counts
# focus more on the moderately expressed genes
# cond <- factor(1:nGroup)
# data3_int <- round(data3) 
# dds <- DESeqDataSetFromMatrix(data3_int, DataFrame(cond), ~ cond)
# dds <- estimateSizeFactors(dds)
# sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
# DESeq2_out <- data3
# for(i in 1:nGroup){
#   DESeq2_out[, i] <- DESeq2_out[, i] / sizeFactors(dds)[i]
# }

# # DESeq2
cond <- factor(rep(1:nGroup, nRep))
data2_int <- round(data2)
dds <- DESeqDataSetFromMatrix(data2_int, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
DESeq2_temp <- data2
for(i in 1:sum(nRep)){
  DESeq2_temp[, i] <- DESeq2_temp[, i] / sizeFactors(dds)[i]
}
DESeq2_out <- NULL
start <- 1
for(rep in nRep){
  DESeq2_out <- cbind(DESeq2_out, rowMeans(DESeq2_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

# edgeR normalization - TMM
normFactors <- edgeR::calcNormFactors(data2)
totCounts <- colSums(data2)
normFactors <- normFactors * totCounts
normFactors <- normFactors/geometric.mean(normFactors)
tmm_temp <- data2
for(i in 1:sum(nRep)){
  tmm_temp[, i] <- tmm_temp[, i] / normFactors[i]
}
tmm_out <- NULL
start <- 1
for(rep in nRep){
  tmm_out <- cbind(tmm_out, rowMeans(tmm_temp[, start : (start + rep - 1)]))
  start <- start + rep
}


# TCC with replicates
tcc <- new("TCC", data2, rep(1:nGroup, each = 10))
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
totCounts <- colSums(data2)
tccEsts <- tccEsts * totCounts
tccEsts <- tccEsts / geometric.mean(tccEsts)

tcc_temp <- getNormalizedData(tcc)
tcc_out <- NULL
start <- 1
for(rep in nRep){
  tcc_out <- cbind(tcc_out, rowMeans(tcc_temp[, start : (start + rep - 1)]))
  start <- start + rep
}


###################################################
# check if fix the dislocation
outlis <- list(totalcount_out, cosbin_out$data, tmm_out, DESeq2_out, tcc_out)

i <- 5
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
  scale_color_manual(values=c("green", "blue", "magenta", "orange", "red"))

#######################################################


# overlap of CEG (total found vs ground truth)
# ind <- which(cos_iCEG(plot_data) >= 0.99)
lst <- sort(cos_iCEG(plot_data), index.return = TRUE, decreasing = TRUE)
ind <- lst$ix[1 : (pCEG * nGene)]
#length(ind)
#length(ind_CEG)
length(intersect(ind, ind_CEG)) / length(ind_CEG)

######################################################

# ROC (iCEG vs DEG)
#outlis <- list(totalcount_out, cosbin_out$data, tmm_out, DESeq2_out, tcc_out)
TPRlis <- list()
FPRlis <- list()
AUClis <- list()
for(i in 1:length(outlis)){
  lst <- sort(cos_iCEG(outlis[[i]]), index.return = TRUE, decreasing = TRUE)
  res <- NULL
  P <- pCEG * nGene
  N <- nGene - P
  TPR <- NULL
  FPR <- NULL
  res <- NULL
  for (j in 1:nGene) {
    res <- c(res, sum(CEG_grnd[lst$ix[1:j]]))
    TP <- sum(CEG_grnd[lst$ix[1:j]] == 1)
    FP <- sum(CEG_grnd[lst$ix[1:j]] == 0)
    TPR <- c(TPR, TP / P)
    FPR <- c(FPR, FP / N)
  }
  TPRlis[[i]] <- TPR
  FPRlis[[i]] <- FPR
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR, TPR)) / 0.15
}

colorlis <- c("red", "green", "blue", "orange", "magenta")

plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlim = c(0, 0.15),
  ylim = c(0, 1),
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

############################################


emf(file ="roc1.emf", width = 6.4, height = 5)
plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlim = c(0, 0.15),
  ylim = c(0, 1),
  type= "l",
  lty = 2,
  col = alpha(colorlis[1], 1))


for(i in 2:length(outlis)){
  lines(FPRlis[[i]], TPRlis[[i]], lty=2, col = alpha(colorlis[i], 1))
}


legend("bottomright",
       legend = c("Total count", "Cosbin", "TMM/edgeR", "DESeq2",
                  "DEGES/TCC"),
       col= colorlis,
       lty=rep(2, length(outlis)), cex=0.8)
dev.off()
########################################################

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
  print(DescTools::AUC(FPR, TPR))
  AUClis[[i]] <- 1 - (1 - DescTools::AUC(FPR, TPR)) / 0.95
}

colorlis <- c("red", "green", "blue", "orange", "magenta")


emf(file ="roc3.emf", width = 6.4, height = 5)

plot(
  FPRlis[[1]],
  TPRlis[[1]],
  xlim = c(0, 0.003),
  ylim = c(0, 1),
  type= "l",
  lty = 2,
  col = colorlis[1]
)


for(i in 2:length(outlis)){
  lines(FPRlis[[i]], TPRlis[[i]], lty=2, col = colorlis[i])
}


legend("bottomright",
       legend = c("Total count", "Cosbin", "TMM/edgeR", "DESeq2",
                  "DEGES/TCC"),
       col= colorlis,
       lty=rep(2, length(outlis)), cex=0.8)
dev.off()

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
ind <- ind_CEG
b <- c(20, 2, 10, 10, 3)
for(j in 1:length(outlis)){
  plot_data <- outlis[[j]]
  res <- NULL
  nGroups <- 3
  for (i in ind) {
    ref <- rep(1, nGroups)
    res <-
      c(res, sum(ref * plot_data[i,]) / sqrt(sum(ref ^ 2) * sum(plot_data[i,] ^ 2)))
  }
  # ggplot(as.data.frame(res), aes(x=res)) + 
  #   geom_histogram(binwidth=0.0005)
  hist(res, xlim = c(0.8, 1), breaks = seq(0.8, 1, 0.001))
  print(Mean(res))
}

#########################################################

# cosine sDEG displacement
for(k in 1:length(outlis)){
  plot_data <- outlis[[k]]
  for (j in 1:nGroup) {
    res <- NULL
    start <- 1
    for (i in start : (start + psDEG[j] * nGene - 1)) {
      res <-
        c(res, sum(data_grnd[i, ] * plot_data[i,]) / sqrt(sum(data_grnd[i, ]  ^ 2) * sum(plot_data[i,] ^ 2)))
    }
    start <- start + (psDEG[j]* nGene)
    hist(res, xlim = c(0.8, 1), breaks = seq(0.8, 1, 0.001))
    print(Mean(res))
  }
  print("*****")
}
