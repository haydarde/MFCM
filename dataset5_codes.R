rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/5_diagnosisUse.csv")
load(file = "distMats/genDistdata5.RData")

data <- data[,c(-4,-7,-8)] 

contVars <- c(1)
nomVars <-NULL
binVars <- c(2:5) 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)
weightStatus <- c(0.5, 0, 0, .5)

new <- 0:1
old <- 1:2
data$X2[data$X2 %in% old] <- new[match(data$X2,old,nomatch = 0)]
data$X3[data$X3 %in% old] <- new[match(data$X3,old,nomatch = 0)]
data$X4[data$X4 %in% old] <- new[match(data$X4,old,nomatch = 0)]
data$X5[data$X5 %in% old] <- new[match(data$X5,old,nomatch = 0)]
data$X6[data$X6 %in% old] <- new[match(data$X6,old,nomatch = 0)]

ranges <- apply(scale(data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

c <- 2
m <- 2
weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes
