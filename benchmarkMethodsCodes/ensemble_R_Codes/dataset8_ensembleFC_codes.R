rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")


data <- read.csv("datasets/8_flower.csv")
load("distMats/data8EnsDist.RData")

contVars <- c(7,8)
binVars <- c(1,2,3) 

data$V1 <- as.integer(data$V1)
data$V2 <- as.integer(data$V2)
data$V3 <- as.integer(data$V3)
data$V4 <- as.integer(data$V4)
data$V5 <- as.integer(data$V5)
data$V6 <- as.integer(data$V6)

new <- 0:1
old <- 1:2
data$V1[data$V1 %in% old] <- new[match(data$V1,old,nomatch = 0)]
data$V2[data$V2 %in% old] <- new[match(data$V2,old,nomatch = 0)]
data$V3[data$V3 %in% old] <- new[match(data$V3,old,nomatch = 0)]


# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data8EnsDist.RData")

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
