rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")


data <- read.csv("datasets/3_absenteeismUse.csv")
load("distMats/data3EnsDist.RData")

contVars <- c(5)
binVars <- c(1,2,4)


data$eth <- as.integer(as.factor(data$eth))
data$sex <- as.integer(as.factor(data$sex))
data$age <- as.integer(as.factor(data$age))
data$lrn <- as.integer(as.factor(data$lrn))

new <- 0:1
old <- 1:2
data$eth[data$eth %in% old] <- new[match(data$eth,old,nomatch = 0)]
data$sex[data$sex %in% old] <- new[match(data$sex,old,nomatch = 0)]
data$lrn[data$lrn %in% old] <- new[match(data$lrn,old,nomatch = 0)]

# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data3EnsDist.RData")

c <- 2
m <- 1.5

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
