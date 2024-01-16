rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/2_exasensUse.csv")
load(file = "distMats/genDistdata2.RData")
load(file = "mixDist/data2MixDist.RData")

data <- data[,c(-1,-7)] 
contVars <- c(1,2,4)
nomVars <- c(5)
binVars <- c(3) # Binary variables need to be coded to 0 and 1. kredit is already in os and 1s.
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

new <- 0:1
old <- 1:2
data$Gender[data$Gender %in% old] <- new[match(data$Gender,old,nomatch = 0)]

ranges <- apply(scale(data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars, ordVars)

# mixDist <- mixDistance(data  = data, H = data, varStatus = varStatus, weightStatus = weightStatus,
#                        contDist = "Euclidean", n = nrow(data), k = nrow(data), range = allRanges)
# save(mixDist,file = "data2MixDist.RData")

c <- 4
m <- 1.5

res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
