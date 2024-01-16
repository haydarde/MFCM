rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/12_CES11Use.csv")
load("distMats/data12EnsDist.RData")


contVars <- c(2)
binVars <- c(3,4,7) 

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
