rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/10_HousePricesUse.csv")
load("distMats/data10EnsDist.RData")

contVars <- c(1,2,3,4,5,11)
binVars <- c(6,7,8,9,10)

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

# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data10EnsDist.RData")

c <- 2
m <- 2.5
tictoc::tic(quiet = TRUE)

ensembleFC <- ensemble_FC(data = data, contVars = contVars, binVars = binVars, c = c, m = m, distDat = distDat, seed = 12321)

a <- tictoc::toc(quiet = TRUE)
cat("Ensemble Fuzzy clustering done in ;", a$toc - a$tic, "; seconds ...", "\n")

clsStats <- list()
clsStats[[1]] <- cqcluster.stats(d = distDat, clustering = as.integer(ensembleFC$clust[,1]))
statsTake <- c(28,29,42,44)
df <- data.frame( FCMM_MD = clsStats[[1]][statsTake])
SILs <- mean(unlist(clsStats[[1]][20]))
df <- data.frame(SIL = SILs, df)
dfClusSizes <- data.frame(ensembleFC = toString(table(factor(ensembleFC$clust[,1]))))
df <- data.frame(df, ClustSize = t(dfClusSizes))
df
