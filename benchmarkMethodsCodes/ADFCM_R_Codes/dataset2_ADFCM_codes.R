rm(list=ls())
# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ADFCM_R_Codes/")
source("ADFCM_functions.R")

data <- read.csv("datasets/2_exasensUse.csv")

data <- data[,c(-1,-7)] 
contVars <- c(1,2,4)
binVars <- c(3)


c <- 2
m <- 1.5

tictoc::tic(quiet = TRUE)

catVars <- setdiff(c(1:ncol(data)) , c(contVars))
data <- cbind(data[,contVars], data[,catVars])
contVars <- 1:length(contVars)
catVars <- (length(contVars)+1):ncol(data)

adfcm <- ADFCM(X = data, contVars = contVars, catVars = catVars, k = c, m = m, UL = 4000, conv = 0.001, convThreshold = 10^-6, seed = 12321)

ccDistData <- catContDistanceData(data = data, H = data, U = adfcm$U, w = adfcm$w, contVars = contVars, 
                                  catVars = catVars, m = m, n = nrow(data), k = ncol(data))

a <- tictoc::toc(quiet = TRUE)
cat("ADFCM clustering done in ;", a$toc - a$tic, "; seconds ...", "\n")

clsStats <- list()
clsStats[[1]] <- cqcluster.stats(d = ccDistData, clustering = as.integer(adfcm$clus[,1]))
statsTake <- c(28,29,42,44)
df <- data.frame( ADFCM = clsStats[[1]][statsTake])
SILs <- mean(unlist(clsStats[[1]][20]))
df <- data.frame(SIL = SILs, df)
dfClusSizes <- data.frame(ADFCM = toString(table(factor(adfcm$clus[,1]))))
df <- data.frame(df, ClustSize = t(dfClusSizes))
df
