# Implementation of FCM mix by D’Urso and Massari (2019): https://doi.org/10.1016/j.ins.2019.07.100

rm(list=ls())
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/FCMMd-MD_R_Codes/")
source("~/Documents/makaleler/W_GND_FKMmix/FCMMd-MD_functions.R")

dataOrg <- read.csv("datasets/18_NCbirthsUse.csv")

# Same types of variables must come next to each other. 
# So rearrange the columns of data based on the data types.

contVars <- c(3,4,7,9,10)
nomVars <- c(1,6,13)
binVars <- c(2,5,8,11,12) 
ordVars <- NULL

dataOrg$MomRace <- as.integer(as.factor(dataOrg$MomRace))

new <- 0:1
old <- 1:2
dataOrg$Sex[dataOrg$Sex %in% old] <- new[match(dataOrg$Sex,old,nomatch = 0)]
dataOrg$Marital[dataOrg$Marital %in% old] <- new[match(dataOrg$Marital,old,nomatch = 0)]

dat <- cbind(dataOrg[,contVars], dataOrg[,binVars], dataOrg[,nomVars], dataOrg[,ordVars])

X <- dat #dat is arranged to have the same types of variables next to each other.
p <- list(c(1:length(contVars)),c((length(contVars)+1):ncol(dat))) # Each element of p is a vector that shows the indices of 
                                                                   # each variable type among the columns of X
S <- 2 # There are S types of variables. D’Urso and Massari (2019) do not distinguish the types of categorical variables.
x <- list()
for ( s in 1:S){
  x[[s]] <- as.matrix(X[,p[[s]]])
}

C <- 2 # The number of clusters
m <- 1.1 # The degree of fuzziness
n <- nrow(X) # sample size
max.iter <- 1000
q <- matrix(NA, nrow = S, ncol = C)
init.medoids <- list() # init.medoids has one element for each s = 1,...,S, under each element, there is a Cxp[[s]] matrix
for (s in 1:S){
  q[s,] <- sample(1:n,C)
  init.medoids[[s]] <- x[[s]][q[s,],] # Randomly selected initial medoids for each variable type
}
init.w <- c(1,1) # Initial weights
w <- init.w

# distanceObs <- array(NA, dim = c(S,n,n))
# for (s in 1:S){
#   distanceObs[s,,] <- euclideanDistanceObs(data = x[[s]], n = n, square = TRUE)
#   distanceObs[s,,] <- (distanceObs[s,,] - min(distanceObs[s,,]))/(max(distanceObs[s,,])-min(distanceObs[s,,])) # Normalise to [0,1] range
# }
# save(distanceObs,file = "/distMats/distanceObs_data18.RData")
load("distMats/distanceObs_data18.RData")

FCMMd_MDres <- FCMMd_MD(x = x, n = n, C = C, S = S, init.medoids = init.medoids, max.iter = max.iter)

clusts <- assignClust(FCMMd_MDres$u)

clsStats <- list()
clsStats[[1]] <- cqcluster.stats(d = dist(dat), clustering = as.integer(clusts[,1]))
statsTake <- c(28,29,42,44)
df <- data.frame( FCMM_MD = clsStats[[1]][statsTake])
SILs <- mean(unlist(clsStats[[1]][20]))
df <- data.frame(SIL = SILs, df)
dfClusSizes <- data.frame(FCMM_MD = toString(table(factor(clusts[,1]))))
df <- data.frame(df, ClustSize = t(dfClusSizes))
df

