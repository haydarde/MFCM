# Implementation of FCM mix by D’Urso and Massari (2019): https://doi.org/10.1016/j.ins.2019.07.100

rm(list=ls())
setwd("/Users/haydardemirhan/Documents/makaleler/W_GND_FKMmix/analysis/FCMMd-MD_R_Codes/")
source("~/Documents/makaleler/W_GND_FKMmix/FCMMd-MD_functions.R")

dataOrg <- read.csv("datasets/1_SouthGermanCreditUse.csv")

# Same types of variables must come next to each other. 
# So rearrange the columns of data based on the data types.
contVars <- c(2,5,13,16)
nomVars <- c(1,3,4,6,9,10,14,15)
binVars <- c(18,19,20,21) 
ordVars <- c(7,8,11,12,17)

dat <- cbind(dataOrg[,contVars], dataOrg[,binVars], dataOrg[,nomVars], dataOrg[,ordVars])

tictoc::tic(quiet = TRUE)
X <- dat #dat is arranged to have the same types of variables next to each other.
p <- list(c(1:length(contVars)),c((length(contVars)+1):ncol(dat))) # Each element of p is a vector that shows the indices of 
                                                                   # each variable type among the columns of X
S <- 2 # There are S types of variables. D’Urso and Massari (2019) do not distinguish the types of categorical variables.
x <- list()
for ( s in 1:S){
  x[[s]] <- X[,p[[s]]]
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
# save(distanceObs,file = "/distMats/distanceObs_data1.RData")
load("distMats/distanceObs_data1.RData")

FCMMd_MDres <- FCMMd_MD(x = x, n = n, C = C, S = S, m = m, init.medoids = init.medoids, max.iter = max.iter)

clusts <- assignClust(FCMMd_MDres$u)
a <- tictoc::toc(quiet = TRUE)
cat("FCMMd-MD done in", a$toc - a$tic, "seconds ...", "\n")

clsStats <- list()
clsStats[[1]] <- cqcluster.stats(d = dist(dat), clustering = as.integer(clusts[,1]))
statsTake <- c(28,29,42,44)
df <- data.frame( FCMM_MD = clsStats[[1]][statsTake])
SILs <- mean(unlist(clsStats[[1]][20]))
df <- data.frame(SIL = SILs, df)
dfClusSizes <- data.frame(FCMM_MD = toString(table(factor(clusts[,1]))))
df <- data.frame(df, ClustSize = t(dfClusSizes))
df
