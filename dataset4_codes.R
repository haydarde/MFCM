rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/4_tipsUse.csv")
load(file = "distMats/genDistdata4.RData")

contVars <- c(1,2)
nomVars <- c(5)
binVars <- c(3,4,6) 
ordVars <- c(7)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$sex <- as.integer(as.factor(data$sex))
data$smoker <- as.integer(as.factor(data$smoker))
data$day <- as.integer(as.factor(data$day))
data$time <- as.integer(as.factor(data$time))

new <- 0:1
old <- 1:2
data$sex[data$sex %in% old] <- new[match(data$sex,old,nomatch = 0)]
data$smoker[data$smoker %in% old] <- new[match(data$smoker,old,nomatch = 0)]
data$time[data$time %in% old] <- new[match(data$time,old,nomatch = 0)]


ranges <- apply((data),2,range)
allRanges <- ranges[2,]-ranges[1,]

dataCl <- data
rownames(dataCl) <- 1:nrow(dataCl)
clRes <- clValid(obj = dataCl, nClust = 2:6, maxitems = nrow(dataCl)+1, 
                 clMethods=c("hierarchical"), validation="internal")
optimalScores(clRes)

weightStatus <- calcEntropy(data, contVars, binVars, nomVars)
c <- 2
m <- 1.1
res <- runAllMethods(data = data, c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, seed = 1234)
res$allRes

