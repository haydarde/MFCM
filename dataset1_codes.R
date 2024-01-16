rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/1_SouthGermanCreditUse.csv")
load(file = "distMats/genDistdata1.RData")
load(file = "mixDist/data1MixDist.RData")

contVars <- c(2,5,13,16)
nomVars <- c(1,3,4,6,9,10,14,15)
binVars <- c(18,19,20,21) 
ordVars <- c(7,8,11,12,17)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

new <- 0:1
old <- 1:2
data$pers[data$pers %in% old] <- new[match(data$pers,old,nomatch = 0)]
data$telef[data$telef %in% old] <- new[match(data$telef,old,nomatch = 0)]
data$gastarb[data$gastarb %in% old] <- new[match(data$gastarb,old,nomatch = 0)]

ranges <- apply(data,2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1,
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)
# 
# mixDist <- mixDistance(data  = data, H = data, varStatus = varStatus, weightStatus = weightStatus,
#                        contDist = "Euclidean", n = nrow(data), k = nrow(data), range = allRanges)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)
c <- 2
m <- 1.05

res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = T)
res$allRes

