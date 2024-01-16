rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/5_diagnosisUse.csv")
load("distMats/data5EnsDist.RData")

data <- data[,c(-4,-7,-8)] 

contVars <- c(1)
binVars <- c(2:5) 

new <- 0:1
old <- 1:2
data$X2[data$X2 %in% old] <- new[match(data$X2,old,nomatch = 0)]
data$X3[data$X3 %in% old] <- new[match(data$X3,old,nomatch = 0)]
data$X4[data$X4 %in% old] <- new[match(data$X4,old,nomatch = 0)]
data$X5[data$X5 %in% old] <- new[match(data$X5,old,nomatch = 0)]
data$X6[data$X6 %in% old] <- new[match(data$X6,old,nomatch = 0)]

# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data5EnsDist.RData")

c <- 3
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
