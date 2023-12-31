rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/14_HouseTasksUse.csv")
load(file = "distMats/genDistdata14.RData")

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
c <- 2
useIndex <- "SIL.F"
m <- 1.25
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes
