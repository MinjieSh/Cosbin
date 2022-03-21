
#################################################### (group mean)

input_data <- data3
nGroup <- 3

#### Cosbin normaliztaion (group mean)
start.time <- Sys.time()
cosbin_out <- cosbin(input_data)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


#### Total count normaliztaion (group mean)
totalcount_out <- totalcount(input_data)$data


#### DESeq2 (group mean)
cond <- factor(1:nGroup)
input_data_int <- round(input_data)
dds <- DESeqDataSetFromMatrix(input_data_int, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
DESeq2_out <- input_data
for(i in 1:nGroup){
  DESeq2_out[, i] <- DESeq2_out[, i] / sizeFactors(dds)[i]
}


#### TMM/edgeR (group mean)
normFactors <- edgeR::calcNormFactors(input_data)
normFactors <- normFactors/geometric.mean(normFactors)
tmm_out <- input_data
for(i in 1:nGroup){
  tmm_out[, i] <- tmm_out[, i] / normFactors[i]
}


#### TCC (group mean)
tcc <- new("TCC", input_data, 1:nGroup)
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

tcc_out <- input_data
for(i in 1:nGroup){
  tcc_out[, i] <- tcc_out[, i] / tccEsts[i]
}

#################################################### (with replicates, then average to evaluate)

input_data2 <- data2
nGroup <- 3
nRep <- c(10, 10, 10)

#### DESeq2 (with replicates, then average to evaluate)
cond <- factor(rep(1:nGroup, nRep))
input_data2_int <- round(input_data2)
dds <- DESeqDataSetFromMatrix(input_data2_int, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
DESeq2_temp <- input_data2
for(i in 1:sum(nRep)){
  DESeq2_temp[, i] <- DESeq2_temp[, i] / sizeFactors(dds)[i]
}

## get average for evaluation
DESeq2_out <- NULL
start <- 1
for(rep in nRep){
  DESeq2_out <- cbind(DESeq2_out, rowMeans(DESeq2_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

#### edgeR normalization - TMM (with replicates, then average to evaluate)
normFactors <- edgeR::calcNormFactors(input_data2)
totCounts <- colSums(input_data2)
normFactors <- normFactors * totCounts
normFactors <- normFactors/geometric.mean(normFactors)
tmm_temp <- input_data2
for(i in 1:sum(nRep)){
  tmm_temp[, i] <- tmm_temp[, i] / normFactors[i]
}

## get average for evaluation
tmm_out <- NULL
start <- 1
for(rep in nRep){
  tmm_out <- cbind(tmm_out, rowMeans(tmm_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

#### TCC (group mean, different test.method)
tcc <- new("TCC", input_data, 1:nGroup)
tcc <-
  TCC::calcNormFactors(
    tcc,
    norm.method = "tmm",
    test.method = "wad",
    iteration = 1,
    FDR = 0.1,
    floorPDEG = 0.05
  )
tccEsts <- tcc$norm.factors
totCounts <- colSums(input_data)
tccEsts <- tccEsts * totCounts
tccEsts <- tccEsts / geometric.mean(tccEsts)

tcc_out <- getNormalizedData(tcc)


#### TCC (with replicates, then average to evaluate)
tcc <- new("TCC", input_data2, rep(1:nGroup, each = 10))
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
totCounts <- colSums(input_data2)
tccEsts <- tccEsts * totCounts
tccEsts <- tccEsts / geometric.mean(tccEsts)

tcc_temp <- input_data2
for(i in 1:sum(nRep)){
  tcc_temp[, i] <- tcc_temp[, i] / tccEsts[i]
}

## get average for evaluation
tcc_out <- NULL
start <- 1
for(rep in nRep){
  tcc_out <- cbind(tcc_out, rowMeans(tcc_temp[, start : (start + rep - 1)]))
  start <- start + rep
}