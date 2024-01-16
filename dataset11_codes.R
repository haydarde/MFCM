rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")


data <- read.csv("datasets/11_FairUse.csv")
load(file = "distMats/genDistdata11.RData")
load(file = "mixDist/data11MixDist.RData")

contVars <- c(2,3,9)
nomVars <- c(5,7,8)
binVars <- c(1,4) # Binary variables need to be coded to 0 and 1. 
ordVars <- c(6)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$sex <- as.integer(as.factor(data$sex))
data$child <- as.integer(as.factor(data$child))

new <- 0:1
old <- 1:2
data$sex[data$sex %in% old] <- new[match(data$sex,old,nomatch = 0)]
data$child[data$child %in% old] <- new[match(data$child,old,nomatch = 0)]

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)


c <- 3
weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
m <- 1.25
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
