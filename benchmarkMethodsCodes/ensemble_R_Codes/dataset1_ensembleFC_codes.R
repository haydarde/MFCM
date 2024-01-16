rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/1_SouthGermanCreditUse.csv")
load("distMats/data1EnsDist.RData")

contVars <- c(2,5,13,16)
binVars <- c(18,19,20,21) 

new <- 0:1
old <- 1:2
data$pers[data$pers %in% old] <- new[match(data$pers,old,nomatch = 0)]
data$telef[data$telef %in% old] <- new[match(data$telef,old,nomatch = 0)]
data$gastarb[data$gastarb %in% old] <- new[match(data$gastarb,old,nomatch = 0)]

# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data1EnsDist.RData")

c <- 2
m <- 1.075

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
