# super sample
start.time <- Sys.time()
cosbin_out <- cosbin(data3)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# normaliztaion
totalcount_out <- totalcount(data3)$data


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


#### TMM/edgeR
normFactors <- edgeR::calcNormFactors(data3)
normFactors <- normFactors/geometric.mean(normFactors)
tmm_out <- data3
for(i in 1:nGroup){
  tmm_out[, i] <- tmm_out[, i] / normFactors[i]
}


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

tcc_out <- data3
for(i in 1:nGroup){
  tcc_out[, i] <- tcc_out[, i] / tccEsts[i]
}

####################################################

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