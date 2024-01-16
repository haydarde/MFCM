rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/10_HousePricesUse.csv")
load(file = "distMats/genDistdata10.RData")
load(file = "mixDist/data10MixDist.RData")

contVars <- c(1,2,3,4,5,11)
nomVars <- NULL
binVars <- c(6,7,8,9,10)
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$driveway <- as.integer(as.factor(data$driveway))
data$recreation <- as.integer(as.factor(data$recreation))
data$fullbase <- as.integer(as.factor(data$fullbase))
data$gasheat <- as.integer(as.factor(data$gasheat))
data$aircon <- as.integer(as.factor(data$aircon))

data$stories <- as.numeric(data$stories)

new <- 0:1
old <- 1:2
data$driveway[data$driveway %in% old] <- new[match(data$driveway,old,nomatch = 0)]
data$recreation[data$recreation %in% old] <- new[match(data$recreation,old,nomatch = 0)]
data$fullbase[data$fullbase %in% old] <- new[match(data$fullbase,old,nomatch = 0)]
data$gasheat[data$gasheat %in% old] <- new[match(data$gasheat,old,nomatch = 0)]
data$aircon[data$aircon %in% old] <- new[match(data$aircon,old,nomatch = 0)]

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")

optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
c <- 2
m <- 2.5
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
