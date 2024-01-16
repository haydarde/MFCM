rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/18_NCbirthsUse.csv")


contVars <- c(3,4,7,9,10)
binVars <- c(2,5,8,11,12) 

data$MomRace <- as.integer(as.factor(data$MomRace))
load("distMats/data18EnsDist.RData")

new <- 0:1
old <- 1:2
data$Sex[data$Sex %in% old] <- new[match(data$Sex,old,nomatch = 0)]
data$Marital[data$Marital %in% old] <- new[match(data$Marital,old,nomatch = 0)]


c <- 3
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
