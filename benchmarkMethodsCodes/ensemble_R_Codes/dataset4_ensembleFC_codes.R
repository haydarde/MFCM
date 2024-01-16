rm(list=ls())
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/4_tipsUse.csv")
load("distMats/data4EnsDist.RData")

contVars <- c(1,2)
binVars <- c(3,4,6) 


data$sex <- as.integer(as.factor(data$sex))
data$smoker <- as.integer(as.factor(data$smoker))
data$day <- as.integer(as.factor(data$day))
data$time <- as.integer(as.factor(data$time))

new <- 0:1
old <- 1:2
data$sex[data$sex %in% old] <- new[match(data$sex,old,nomatch = 0)]
data$smoker[data$smoker %in% old] <- new[match(data$smoker,old,nomatch = 0)]
data$time[data$time %in% old] <- new[match(data$time,old,nomatch = 0)]


# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data4EnsDist.RData")

c <- 2
m <- 1.1
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

