rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/15_Blood1Use.csv")
load(file = "distMats/genDistdata15.RData")
load(file = "mixDist/data15MixDist.RData")

contVars <- c(1)
nomVars <- c(3)
binVars <- c(2) 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
c <- 3
m <- 1.5
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
