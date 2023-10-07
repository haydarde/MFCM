rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")


data <- read.csv("datasets/3_absenteeismUse.csv")
load(file = "distMats/genDistdata3.RData")

contVars <- c(5)
nomVars <- NULL
binVars <- c(1,2,4)
ordVars <- c(3)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$eth <- as.integer(as.factor(data$eth))
data$sex <- as.integer(as.factor(data$sex))
data$age <- as.integer(as.factor(data$age))
data$lrn <- as.integer(as.factor(data$lrn))

new <- 0:1
old <- 1:2
data$eth[data$eth %in% old] <- new[match(data$eth,old,nomatch = 0)]
data$sex[data$sex %in% old] <- new[match(data$sex,old,nomatch = 0)]
data$lrn[data$lrn %in% old] <- new[match(data$lrn,old,nomatch = 0)]

ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars)
c <- 2
m <- 1.5
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes
