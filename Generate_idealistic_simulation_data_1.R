setwd("C:/Users/Minjie Shen/Google Drive/inCOT simulation")

###### Prepare - install packages - R 3.6.3

# remove.packages("rlang")
# install.packages("rlang", INSTALL_opts = "--no-lock")
library(rlang)
library(Ternary)
library(ggtern)
library(hrbrthemes)
library(DirichletReg)
library(DescTools)
library(TCC)
library(psych)
library(scales)
library(devEMF)
library(rgl)
library(magick)
# devtools::install_version("ggplot2", version = "3.2.1", repos = "http://cran.us.r-project.org", type="source")

source('Cosbin_functions.R')

set.seed(7)

####################################################

###### Generate ideal simulation data
###### 1000 genes, 600 SDEGs (corner, 60 + 180 + 360), iCEGs (center) determined by cosine
nGroup <- 3
nGene <- 1000
psDEG <-  c(0.06, 0.18, 0.36) ###### absolute percentage over all genes

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

###### Ground truth label
CEG_grnd <- rep(0, nGene)
CEG_grnd[ind_CEG] <- 1
sum(CEG_grnd)

sDEG_grnd <- rep(0, nGene)
sDEG_grnd[ind_sDEG] <- 1
sum(sDEG_grnd)

#######################################################

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
  scale_color_manual(values=c("blue", "green", "magenta", "orange", "red"))

# length(which(cos_iCEG(data) >= 0.99))

#############################################

##### 3D plot
# test <- data_grnd
# rs <- rowSums(test)
# #test <- test[-which(rs > 1000), ]
# #coltemp <- lbl2[-which(rs > 1000)]
# coltemp <- lbl2
# 
# coltemp[coltemp == 'sDEG G1'] <- "magenta"
# coltemp[coltemp == 'sDEG G2'] <- "orange"
# coltemp[coltemp == 'sDEG G3'] <- "red"
# coltemp[coltemp == 'iCEG (labeled)'] <- "green4"
# coltemp[coltemp == '/'] <- "steelblue1"
# 
# # Static chart
# G1 <- test[,1]
# G2 <- test[,2]
# G3 <- test[,3]
# plot3d( G1, G2, G3, col = coltemp, type = "p", r = 0.2)

#################

#### assume same total count
factor <- 1 / totalcount(data2)$norm_factor
data2 <- totalcount(data2)$data

#### within group mean
data3 <- NULL
start <- 1
for(rep in nRep){
  data3 <- cbind(data3, rowMeans(data2[, start : (start + rep - 1)]))
  start <- start + rep
}

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
