library(DirichletReg)
library(ggplot2)


test <- rdirichlet(1000, c(0.2, 0.2, 0.2))
pca_res <- prcomp(test, scale. = TRUE)
plot(pca_res$x[,1], pca_res$x[,2])


#############################################



res <- inCOT(data)

hist_res <- as.vector(inCOT(data))
hist(hist_res)

hist(hist_res[CEG_grnd], 50)
hist(hist_res[-CEG_grnd], 50)

thrld <- quantile(hist_res, k)

grnd_truth <- status2
grnd_truth[grnd_truth == 1] <- -1
grnd_truth[grnd_truth == 0] <- 1
grnd_truth[grnd_truth == -1] <- 0

hist_res[hist_res > thrld] <- 1
hist_res[hist_res <= thrld] <- 0
length(which(hist_res == 1))

CEG_grnd <- which(grnd_truth == 1)
CEG_cot <- which(hist_res == 1)
length(intersect(CEG_cot, CEG_grnd))









###############################################

#### Generate simulation data

# condition * n
nSmpCond <- 1
nCond <- 3
condList <- rep(nSmpCond, nCond)

# Generate baseline proportions for desired number of genes
ngenes <- 2000
# percentages of differential expression
percDE <- 0.6
ndiffExpress <- ngenes * percDE

# DEG index
# overall DEGs, some in some samples, some in other samples
indDEG <- sample(1:ngenes, ndiffExpress)

# Get distribution function of abundance proportions
# This distribution was generated from a real dataset
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))
baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
# baselineprop <- baselineprop/sum(baselineprop) ### 1000 gene baseline
# baseline proportions for each gene; DE genes get multiplied by the fold change
baselinepropList <- vector("list", length = nCond)
for (j in 1:nCond){
  baselinepropList[[j]] <- baselineprop
}
  
# will differential expression be symmetric?
symmetric <- FALSE
# if asymmetric, how asymmetric?
partiallyAsymm <- FALSE # most DE are up-reg under one cond.
completelyAsymm <- FALSE # all DE are up-reg under one cond.

fc <- 5 #fold change for DE genes
# we can either have asymmetric differential expression, with unequal proportions
# of up- and down-regulated genes, or symmetric differential expression
if (symmetric){

  nDEG <- floor(ndiffExpress / nCond)
  indDEGSub <- vector("list", length = nCond)
  x <- 1
  for (j in 1:nCond){
    if(j == nCond)
      indDEGSub[[j]] <- indDEG[x: ndiffExpress]
    else{
      indDEGSub[[j]] <- indDEG[x: (x + nDEG - 1)]
    }
    x <- x + nDEG
  }
  
  for (j in 1:nCond){
    baselinepropList[[j]][indDEGSub[[j]]] <- baselinepropList[[j]][indDEGSub[[j]]] * fc
  }
  
  
} else if (completelyAsymm){
  
  j <- 1
  baselinepropList[[j]][indDEGSub[[j]]] <- baselinepropList[[j]][indDEGSub[[j]]] * fc
  
} else {
  
  nDEG <- floor(c(0.6, 0.3, 0.1) * ndiffExpress)
  indDEGSub <- vector("list", length = nCond)
  x <- 1
  for (j in 1:nCond){
    if(j == nCond)
      indDEGSub[[j]] <- indDEG[x: ndiffExpress]
    else{
      indDEGSub[[j]] <- indDEG[x: (x + nDEG[j] - 1)]
    }
    x <- x + nDEG[j]
  }
  
  for (j in 1:nCond){
    baselinepropList[[j]][indDEGSub[[j]]] <- baselinepropList[[j]][indDEGSub[[j]]] * fc
  }
  
}

################

backup_groundtruth <- do.call("cbind", baselinepropList)


# make these actual proportions
baseSum <- vector("list", length = nCond)
for (j in 1:nCond){
  baseSum[[j]] <- sum(baselinepropList[[j]])
  baselinepropList[[j]] <- baselinepropList[[j]] / baseSum[[j]] 
}



nlibs <- sum(condList)
# Use equal or unequal library sizes
equal <- TRUE
# Library size 
if(equal){
  expected.lib.size <- rep(11e6, nlibs)
} else {
  expected.lib.size <- 20e6 * c(1,0.1,1,0.1,1,0.1)
}
# sequencing depths
sizeRatio <- rtruncnorm(nlibs, a = 0.6, b = 1.4, mean=1, sd=0.2)
# mean expression for each gene
meanGene <- vector("list", length = nCond)
start <- 1
end <- condList[1]
for(j in 1:nCond){
  meanGene[[j]] <- matrix(baselinepropList[[j]], ngenes, 1) %*% 
    matrix(expected.lib.size[start: end] * sizeRatio[start: end], 1, condList[j])
  start <- end + 1
  if (j == nCond){
    end <- ndiffExpress
  } else {
    end <- start + (condList[j] - 1)
  }
}
meanGeneCmb <- do.call("cbind", meanGene)

# status will keep track of which genes are deferentially expressed
status <- rep(0, ngenes)
status[indDEG] <- 1

# Use inverse chi-square or log-normal dispersion
invChisq <- TRUE
# Biological variation
BCV0 <- 0.2 + 1/sqrt(meanGeneCmb)
if(invChisq){
  df.BCV <- 40
  BCV <- BCV0 * sqrt(df.BCV/rchisq(ngenes, df=df.BCV))
} else {
  BCV <- BCV0 * exp( rnorm(ngenes, mean=0, sd=0.25)/2 )
}
if(NCOL(BCV)==1) BCV <- matrix(BCV,ngenes,nlibs)
shape <- 1/BCV^2
scale <- meanGeneCmb/shape
mu <- matrix(rgamma(ngenes*nlibs,shape=shape,scale=scale),ngenes,nlibs)

# Technical variation
counts <- matrix(rpois(ngenes*nlibs,lambda=mu),ngenes,nlibs)

# remove lowly expressed genes
#Filter: following limma/voom code, remove genes with
# rowsums < 10 in the read count matrix
keep <- rowSums(counts)>=10
nkeep <- sum(keep)
counts2 <- counts[keep,]
status2 <- status[keep]

############## ground truth
# cutoff for BH  adjust
cutoff <- 0.05

## naive-DESeq2
cond <- factor(rep(1:nCond, each = nSmpCond))
dds <- DESeqDataSetFromMatrix(counts2, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")
efdrDESeq2[j] <- 
  sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)

powerDESeq2[j] <- 
  sum(res$padj[status2==1] < cutoff)/sum(status2==1)

# want non-DE genes with no zero counts
countsToCompare <- as.data.frame(counts2[status2==0, ])
row_sub = apply(countsToCompare, 1, function(row) all(row !=0 ))
countsToCompare <- countsToCompare[row_sub, ]

# Oracle DESEq2
normFactors <- sizeRatio / unlist(rep(baseSum, each = nSmpCond))
sizeFactors(dds) <- normFactors / geometric.mean(normFactors)

dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
res <- results(dds)

res$padj <- p.adjust(res$pvalue, method="BH")
efdrOracle[j] <- 
  sum(res$padj[status2==0] < cutoff)/sum(res$padj < cutoff)

powerOracle[j] <- 
  sum(res$padj[status2==1] < cutoff)/sum(status2==1)

normCounts <- vector("list", length = nCond)
start <- 1
end <- condList[1]
for(j in 1:nCond){
  normCounts[[j]] <- countsToCompare[,start:end]/(sizeFactors(dds)[j])
  start <- end + 1
  if (j == nCond){
    end <- ndiffExpress
  } else {
    end <- start + (condList[j] - 1)
  }
}

normCounts <- do.call("cbind", normCounts)
  
#####  Evaluation (CEG among different conditions)
mseOracle <- mean(log(normCounts[,1:5]/normCounts[,6:10])^2)
