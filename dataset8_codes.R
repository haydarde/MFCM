rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")


data <- read.csv("datasets/8_flower.csv")
load(file = "distMats/genDistdata8.RData")

contVars <- c(7,8)
nomVars <- c(4,5,6)
binVars <- c(1,2,3) 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)
weightStatus <- c(0.5, 0.2, 0,0.3)

weightStatus <- c(0.6, 0.2, 0,0.2)

data$V1 <- as.integer(data$V1)
data$V2 <- as.integer(data$V2)
data$V3 <- as.integer(data$V3)
data$V4 <- as.integer(data$V4)
data$V5 <- as.integer(data$V5)
data$V6 <- as.integer(data$V6)

new <- 0:1
old <- 1:2
data$V1[data$V1 %in% old] <- new[match(data$V1,old,nomatch = 0)]
data$V2[data$V2 %in% old] <- new[match(data$V2,old,nomatch = 0)]
data$V3[data$V3 %in% old] <- new[match(data$V3,old,nomatch = 0)]


ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
c <- 2
m <- 1.1
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes
