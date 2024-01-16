# setwd("YOUR FOLDER COMSE HERE")
setwd("~/Documents/makaleler/W_GND_FKMmix/analysis/ensemble_R_Codes/")
source("ensemble_Functions.R")

data <- read.csv("datasets/9_CPS1985Use.csv")
load("distMats/data9EnsDist.RData")

contVars <- c(1,2,3,4)
binVars <- c(6,7,10,11) 


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

# distDat <- modeDistanceData(data = data, H = data, contVars = contVars, n = nrow(data), k = nrow(data))
# save(distDat, file = "distMats/data9EnsDist.RData")

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
