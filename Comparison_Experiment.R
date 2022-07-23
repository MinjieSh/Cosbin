library(Ternary)
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
library(ggtern)
source('Cosbin_functions.R')
source('Cosbin_functions_visualization.R')
source('Cosbin_functions_evaluation.R')
library(metagenomeSeq)
library(qsmooth)
####################################################

###### Generate ideal simulation data
###### 1000 genes, 600 SDEGs (corner, 60 + 180 + 360), iCEGs (center) determined by cosine
nGroup <- 3
nGene <- 1000
psDEG <-  c(0.06,0.18,0.36) ###### absolute percentage over all genes

ref <- vector("list", length = nGroup)
for (i in 1:nGroup) {
  ref[[i]] <- rep(0, nGroup)
  ref[[i]][i] <- 1
}

non_DEG <- rdirichlet(400, rep(3, nGroup))

###### iDEG
data_sDEG <- NULL
for (i in 1:nGroup) {
  data_sDEG <-
    rbind(data_sDEG, matrix(rep(ref[[i]], nGene * psDEG[i]), ncol = nGroup, byrow = TRUE))
}

###### sDEG (iDEG + noise)
noise <- rnorm(nGroup * nGene * sum(psDEG), mean = 0, sd = 0.04)
data_sDEG <- abs(data_sDEG + noise)


###### noise cannot be too large
# length(which(cos_iDEG(data_sDEG) < 0.99)) = 0


data <- rbind(data_sDEG, non_DEG)


####################################################

###### Add replicates, here: same number for each group, but you can do different: e.g. nRep <- c(50, 100, 1000)
nRep <- rep(10, nGroup)

data2 <- NULL
for (i in 1:nGroup) {
  data2 <-
    cbind(data2, matrix(rep(data[, i], nRep[i]), ncol = nRep[i], byrow = FALSE))
}

dim(data2)
noise <- rnorm(nRep * nGroup * nGene, mean = 0, sd = 0.01)
data2 <- abs(data2 + noise)


##### modify the rowSums, this constant 1000 can be any large integer
for(i in 1:nGene){
  data2[i,] <- data2[i,] * 1000
}

#### Ground truth data (super sample, to be compared with the normalized result)
data_grnd <- NULL
start <- 1
for(rep in nRep){
  data_grnd <- cbind(data_grnd, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

###### Add label
ind_CEG <- which(cos_iCEG(data_grnd) >= 0.95)
#ind_sDEG <- which(cos_iDEG(data_grnd) >= 0.99)
ind_sDEG <- seq(1,600,1)
#ind_neither <- intersect(which(cos_iCEG(data_grnd) < 0.95),
#                         which(cos_iDEG(data_grnd) < 0.99))
ind_neither <- setdiff(seq(1, 1000, 1), c(ind_CEG, ind_sDEG))


ind_sDEG2 <- list()
for (i in 1:nGroup) {
  ind_sDEG2[[i]] <- which(cos_ref(data_grnd, ref[[i]]) >= 0.98)
}

lbl <- rep(0, nGene)
lbl[ind_CEG] <- 'iCEG (labeled)'
lbl[ind_neither] <- 'DEG'
lbl[ind_sDEG] <- 'sDEG'

lbl2 <- lbl
for (i in 1:nGroup) {
  lbl2[ind_sDEG2[[i]]] <- paste0('sDEG G', i)
}

###### Ground truth label
CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1
sum(CEG_grnd)

sDEG_grnd <- rep(0, nGene)
sDEG_grnd[ind_sDEG] <- 1
sum(sDEG_grnd)

DEG_grnd <- rep(1, nGene)
DEG_grnd <- DEG_grnd-(sDEG_grnd + CEG_grnd)
sum(DEG_grnd)

DEG_and_sDEG_grnd <- rep(1, nGene)
DEG_and_sDEG_grnd <- DEG_and_sDEG_grnd-CEG_grnd
sum(DEG_and_sDEG_grnd)
###### check result in the ternary plot
plot_data <- data_grnd
label <- lbl2
plot_data <-
  t(apply(plot_data, 1, function(x) x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point() +
  labs( x = "G1",
        y = "G2",
        z = "G3") +
  scale_color_manual(values=c("blue", "green", "magenta", "orange", "red","pink"))

# use data2 as input
input_data <- data2
nGroup <- 3
nRep <- c(10, 10, 10)


input_data <- totalcount(input_data)$data # Input data for all methods are the total-count normalized one


### Total count
totalcount_out_temp <- totalcount(input_data)$data

# get the super-samples
totalcount_out <- NULL
start <- 1
for(rep in nRep){
  totalcount_out <- cbind(totalcount_out, rowMeans(totalcount_out_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

plot_data <- totalcount_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )

################################ Plot Show ####################################

# 
# # solid plot for ground truth
# plot_data <- data
# label <- lbl2
# plot_data <-
#   t(apply(plot_data, 1, function(x) x / sum(x)))
# 
# ggtern(data = as.data.frame(plot_data),
#        aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
#   geom_point() +
#   labs( x       = "Group 1",
#         y       = "Group 2",
#         z       = "Group 3") + 
#   scale_color_manual(values=c("blue", "green", "magenta", "orange", "red","pink"))
# 
# # hollow plot for ground truth
# ggtern(data = as.data.frame(plot_data),
#        aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
#   geom_point(alpha = 0.3) +
#   labs( x       = "G1",
#         y       = "G2",
#         z       = "G3") + 
#   scale_color_manual(values=c( "steelblue1", "green4","magenta", "orange", "red")) +
#   theme(
#     title = element_text(size = 12, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 12),
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank()
#   )





### CSS
OTU_read_count = as.data.frame(input_data)


metaSeqObject = newMRexperiment(OTU_read_count) 
metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )
OTU_read_count_CSS_temp = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=FALSE)) 
                                                        

# get the super-samples
OTU_read_count_CSS <- NULL
start <- 1
for(rep in nRep){
  OTU_read_count_CSS <- cbind(OTU_read_count_CSS, rowMeans(OTU_read_count_CSS_temp[, start : (start + rep - 1)]))
  start <- start + rep
}


# check the plot for CSS
plot_data <- OTU_read_count_CSS
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )





### Cosbin
cosbin_out <- cosbin(data)

plot_data <- cosbin_out$data
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )




# 
# # Ground truth
# plot_data <- data
# plot_data <-
#   t(apply(plot_data, 1, function(x)
#     x / sum(x)))
# 
# ggtern(data = as.data.frame(plot_data),
#        aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
#   geom_point(alpha = 0.3) +
#   labs( x       = "G1",
#         y       = "G2",
#         z       = "G3") + 
#   scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
#   theme(
#     title = element_text(size = 12, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     legend.text = element_text(size = 12),
#     strip.background = element_blank(),
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank()
#   )


### Total count
totalcount_out_temp <- totalcount(input_data)$data

# get the super-samples
totalcount_out <- NULL
start <- 1
for(rep in nRep){
  totalcount_out <- cbind(totalcount_out, rowMeans(totalcount_out_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

plot_data <- totalcount_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )



### Qsmooth
dat_qs <- qsmooth(object = input_data,
                  group_factor = c(rep(0, 10),rep(1, 10),rep(2, 10)))

# get the super-samples
qsmooth_out <- NULL
start <- 1
for(rep in nRep){
  qsmooth_out <- cbind(qsmooth_out, rowMeans(qsmoothData(dat_qs)[, start : (start + rep - 1)]))
  start <- start + rep
}


plot_data <- qsmooth_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )



#### DESeq2 (with replicates, then average to evaluate)
cond <- factor(rep(1:nGroup, nRep))
input_data_int <- round(input_data)
dds <- DESeqDataSetFromMatrix(input_data_int, DataFrame(cond), ~ cond)
dds <- estimateSizeFactors(dds)
sizeFactors(dds) <- sizeFactors(dds)/geometric.mean(sizeFactors(dds))
DESeq2_temp <- input_data
for(i in 1:sum(nRep)){
  DESeq2_temp[, i] <- DESeq2_temp[, i] / sizeFactors(dds)[i]
}

# get average for evaluation
DESeq2_out <- NULL
start <- 1
for(rep in nRep){
  DESeq2_out <- cbind(DESeq2_out, rowMeans(DESeq2_temp[, start : (start + rep - 1)]))
  start <- start + rep
}


plot_data <- DESeq2_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )





#### edgeR normalization - TMM (with replicates, then average to evaluate)
normFactors <- edgeR::calcNormFactors(input_data)
totCounts <- colSums(input_data)
normFactors <- normFactors * totCounts
normFactors <- normFactors/geometric.mean(normFactors)
tmm_temp <- input_data
for(i in 1:sum(nRep)){
  tmm_temp[, i] <- tmm_temp[, i] / normFactors[i]
}

# get average for evaluation
tmm_out <- NULL
start <- 1
for(rep in nRep){
  tmm_out <- cbind(tmm_out, rowMeans(tmm_temp[, start : (start + rep - 1)]))
  start <- start + rep
}

plot_data <- tmm_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )





### TCC/DEGES
tcc <- new("TCC", input_data, rep(1:nGroup, each = 10))
tcc <-
  TCC::calcNormFactors(
    tcc,
    norm.method = "tmm",
    test.method = "edger",
    iteration = 3,
    FDR = 0.1,
    floorPDEG = 0.05
  )
tccEsts <- tcc$norm.factors
totCounts <- colSums(input_data)
tccEsts <- tccEsts * totCounts
tccEsts <- tccEsts / geometric.mean(tccEsts)

tcc_temp <- input_data
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

plot_data <- tcc_out
plot_data <-
  t(apply(plot_data, 1, function(x)
    x / sum(x)))

ggtern(data = as.data.frame(plot_data),
       aes(x=plot_data[,1], y=plot_data[,2], z=plot_data[,3], color = label)) +
  geom_point(alpha = 0.3) +
  labs( x       = "G1",
        y       = "G2",
        z       = "G3") + 
  scale_color_manual(values=c("steelblue1", "green4", "magenta", "orange", "red")) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    strip.text.y = element_blank()
  )





#### Evaluation 
norm_data_list <- list(totalcount_out,
                       cosbin_out$data,
                       OTU_read_count_CSS,
                       qsmooth_out,
                       DESeq2_out,
                       tmm_out,
                       tcc_out)


# CEG dislocation: 
########################################################################################################
#|Name        |Cosine with the reference |Mean angle (arccos) with the reference |Max angle |Min angle |
#|Totalcount  |                          |                                       |          |          |
#|Cosbin      |                          |                                       |          |          |
#|CSS         |                          |                                       |          |          |
#|Qsmooth     |                          |                                       |          |          |
#|DESeq2      |                          |                                       |          |          |
#|TMM         |                          |                                       |          |          |
#|TCC         |                          |                                       |          |          |
########################################################################################################

CEG_dislocation <- cosine_dislocation(data_grnd, CEG_grnd, norm_data_list)




# ROC (DEG vs non DEG): 
ROC_DEG <- ROC(norm_data_list, DEG_grnd, "DEG", c("blue", "#00AA0066", "#0000AA66",  "pink", "orange","purple","red"))
ROC_DEG[3] # AUC scores
