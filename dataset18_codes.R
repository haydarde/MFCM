rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/18_NCbirthsUse.csv")
load(file = "distMats/genDistdata18.RData")

contVars <- c(3,4,7,9,10)
nomVars <- c(1,6,13)
binVars <- c(2,5,8,11,12) 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)


data$MomRace <- as.integer(as.factor(data$MomRace))

new <- 0:1
old <- 1:2
data$Sex[data$Sex %in% old] <- new[match(data$Sex,old,nomatch = 0)]
data$Marital[data$Marital %in% old] <- new[match(data$Marital,old,nomatch = 0)]

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
