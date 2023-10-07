rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/17_beautyUse.csv")
load(file = "distMats/genDistdata17.RData")

contVars <- c(1,2,5,16)
nomVars <- c(6)
binVars <- c(7:15) 
ordVars <- c(17)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
c <- 2
m <- 1.25
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes
