# setwd("YOUR FOLDER COMSE HERE")
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/MFCM_R_Codes/")
source("MFCM_Functions_V2.R")

data <- read.csv("datasets/9_CPS1985Use.csv")
load(file = "distMats/genDistdata9.RData")
load(file = "mixDist/data9MixDist.RData")

contVars <- c(1,2,3,4)
nomVars <- c(5,8,9)
binVars <- c(6,7,10,11) 
ordVars <- NULL
varStatus <- list(cont = contVars, nom = nomVars, ord = ordVars, bin = binVars)

data$ethnicity <- as.integer(as.factor(data$ethnicity))
data$region <- as.integer(as.factor(data$region))
data$gender <- as.integer(as.factor(data$gender))
data$occupation <- as.integer(as.factor(data$occupation))
data$sector <- as.integer(as.factor(data$sector))
data$union <- as.integer(as.factor(data$union))
data$married <- as.integer(as.factor(data$married))

new <- 0:1
old <- 1:2
data$region[data$region %in% old] <- new[match(data$region,old,nomatch = 0)]
data$gender[data$gender %in% old] <- new[match(data$gender,old,nomatch = 0)]
data$union[data$union %in% old] <- new[match(data$union,old,nomatch = 0)]
data$married[data$married %in% old] <- new[match(data$married,old,nomatch = 0)]

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
                     ordVars = ordVars, distanceMatrix = genDistdata, mixDistanceMatrix = mixDist, 
                     seed = 1234, DoAll = TRUE)
res$allRes
