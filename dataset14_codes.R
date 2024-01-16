rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/14_HouseTasksUse.csv")
load(file = "distMats/genDistdata14.RData")
load(file = "mixDist/data14MixDist.RData")

contVars <- c(3)
nomVars <- c(1,2)
binVars <- NULL 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$Task <- as.integer(as.factor(data$Task))
data$Who <- as.integer(as.factor(data$Who))

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars)
c <- 4
useIndex <- "SIL.F"
m <- 1.25
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
