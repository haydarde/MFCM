rm(list=ls())
setwd("YOUR FOLDER COMSE HERE")
source("MFCM_Functions_V1.R")

data <- read.csv("datasets/12_CES11Use.csv")
load(file = "distMats/genDistdata12.RData")

contVars <- c(2)
nomVars <- c(1)
binVars <- c(3,4,7) 
ordVars <- c(5,6)
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)
weightStatus <- c(0.7, 0.1, 0.1, 0.1)


data$gender <- as.integer(as.factor(data$gender))
data$abortion <- as.integer(as.factor(data$abortion))
data$location <- as.integer(as.factor(data$location))
data$province <- as.integer(as.factor(data$province))
data$importance <- as.integer(as.factor(data$importance))
data$education <- as.integer(as.factor(data$education))

new <- 0:1
old <- 1:2
data$gender[data$gender %in% old] <- new[match(data$gender,old,nomatch = 0)]
data$abortion[data$abortion %in% old] <- new[match(data$abortion,old,nomatch = 0)]
data$location[data$location %in% old] <- new[match(data$location,old,nomatch = 0)]

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
res <- runAllMethods(data = (data), c = c, m = m, varStatus = varStatus, weightStatus = weightStatus, 
                     allRanges = allRanges, contVars = contVars, nomVars = nomVars, binVars = binVars, 
                     ordVars = ordVars, distanceMatrix = genDistdata, scale = TRUE, seed = 1234)
res$allRes
