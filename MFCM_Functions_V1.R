#############################################################################
## This script includes the codes to run Mixed Fuzzy C-Means algorithm     ##
## developed by Haydar Demirhan and is a supplementary material to the     ##
## article Demirhan, H., Mixed fuzzy C-means clustering, under review.     ##
##                                                                         ##
## Functions related to the implementation of straightforward FCM          ##
## algorithm are adapted from the C++ codes provided by the 'fclust' R     ##
## package by "M.B. Ferraro, P. Giordani and A. Serafini (2019) fclust:    ##
## An R Package for Fuzzy Clustering, The R Journal, 11,                   ##
## https://journal.r-project.org/archive/2019/RJ-2019-017/RJ-2019-017.pdf" ##
#############################################################################

library(fclust); library(ggplot2); library(caret); library(class); library(mclust) 
library(cluster); library(tictoc); library(adana); library(e1071); library(clValid) 
library(fpc); library(ClustOfVar); library(CluMix); library(FD); library(kamila) 

centroidsFKM <- function(data, U, n, k, p, m){
  out <- matrix(0, nrow = k, ncol = p)
  Unew <- t(U) 
  for ( i in 1:k){
    out[i, ] <- ((Unew[i,]^m ) %*% data)/sum(U[,i]^m)
  }
  return(out)
}

# euclideanDistance <- function(data, H, n, k, square = FALSE){
#   out <- matrix(0, nrow = n, ncol = k)
#   for ( i in 1:n){
#     for ( j in 1:k){
#       out[i,j] <- sum((data[i,] - H[j,])^2)
#       if(square == TRUE){
#         out <- sqrt(out)
#       }
#     }
#   }
#   return(out)
# }

genMinkowskiDist <- function(x1, x2, range, p){  
  sum((abs(x1 - x2)/range)^p)^(1/p)  
}

membershipDegrees <- function(D, m, n, k, p){
  out <- matrix(0, nrow = n, ncol = k)
  d <- array(0,n)
  mm <- 0
  for(i in 1:n){
    d <- D[i,]
    if(min(d)== 0){
      mm <- which.min(d)
      out[i,mm] = 1
    }else{
      for(j in 1:k){
        out[i,j] = (1/D[i,j])^(1/(m-1)) / sum( (1/D[i,])^(1/(m-1)))
      }
    }
  }
  return(out)
}

mixDistance <- function(data, H, varStatus, weightStatus, contDist = c("Euclidean","Manhattan","Mahalanobis"), n, k, range = NULL){
  out <- matrix(0, nrow = n, ncol = k)
  for ( i in 1:n){
    for ( j in 1:k){
      
      if (!is.null(varStatus$cont)){
        if (contDist == "Euclidean"){
          out[i,j] <- weightStatus[1] * (sum((data[i,varStatus$cont] - H[j,varStatus$cont])^2))
        }else if (contDist == "Manhattan"){
          out[i,j] <- weightStatus[1] * sum(abs(data[i,varStatus$cont] - H[j,varStatus$cont]))
        }else if (contDist == "GenMinkowski"){
          out[i,j] <- weightStatus[1] * genMinkowskiDist(data[i,varStatus$cont], H[j,varStatus$cont], 
                                                         range = range[varStatus$cont], p = 2)
        }
          
      }
      
      if (!is.null(varStatus$bin)){
        sumDist <- 0
        for ( r in 1:length(varStatus$bin)){ 
          binaryCenter <- sum(H[j,varStatus$bin[r]] >= 0.5)
          sumDist <- sumDist + sum(xor(data[i,varStatus$bin[r]], binaryCenter)) 
        }
        out[i,j] <- out[i,j] + weightStatus[4] *  sumDist
      }
      
      if (!is.null(varStatus$nom)){
        sumDist <- 0
        for ( r in 1:length(varStatus$nom)){ 
          catCenter <- round(H[j,varStatus$nom[r]])
          sumDist <- sumDist + sum(hamming.distance(data[i,varStatus$nom[r]], catCenter))
        }
        out[i,j] <- out[i,j] + weightStatus[2] *  sumDist
      }
      
      if (!is.null(varStatus$ord)){
        out[i,j] <- out[i,j] + weightStatus[3] *  genMinkowskiDist(data[i,varStatus$ord], H[j,varStatus$ord],
                                                                   range = range[varStatus$ord], p = 2)
      }
      
    }
  }
  
  if ((contDist == "Mahalanobis") & (!is.null(varStatus$cont))){
    for ( j in 1:k){
      out[,j] <- weightStatus[1] * mahalanobis(x = data[ ,varStatus$cont], center = H[j,varStatus$cont],
                                              cov = cov(data[ ,varStatus$cont]) )
    }
  }

  return(out)
}


indices <- function(type, X, U, H, m, b, alpha, distance = FALSE){
  value <- 0
  if(type == "PC"){
    value = PC(U)
  } else if (type == "PE"){
    value <- PE(U,b)
  } else if (type == "MPC"){
    value <- MPC(U)
  } else if (type == "SIL"){
    value <- SIL(X = X, U = U,  distance = distance)$sil 
  } else if (type == "SIL.F"){
    value <- SIL.F(X = X, U = U, alpha = alpha, distance = distance) 
  } else if (type == "XB"){
    value = XB(X = X, U = U, H = H, m = m)
  } else {
    stop("No match names.")
  }
  return(value)
}

unifInit <- function(n, d){
  unif <- matrix(0,nrow = n, ncol = d)
  unif_temp <- array(0,n)
  for(i in 1:d) {
    unif_temp <- runif(n,0,1)
    unif[,i] <- unif_temp
  }
  for(j in 1:n) {
    unif[j,] <- unif[j,]/sum(unif[j,])
  }
  return(unif)
}

mainMFCM <- function(data, varStatus, weightStatus, contDist, m, n, p, k, rs, conv, maxit, index, alpha, range = NULL){
  value <- array(0, rs)
  it <- array(0, rs)
  H <- matrix(0, nrow = k, ncol = p)
  Hopt <- matrix(0, nrow = k, ncol = p)
  D <- matrix(0, nrow = n, ncol = k)
  U <- matrix(0, nrow = n, ncol = k)
  Uold <- matrix(0, nrow = n, ncol = k)
  func <- 0; funcOpt <- 0; ind <- 0; indMax <- 0
  Uopt <- matrix(0, nrow = n, ncol = k)
  for ( r in 1:rs){
    U <- unifInit(n, k)
    Uold <- U
    prova <- TRUE
    iter <- 0
    while ((prova == TRUE) & (iter < maxit)){
      iter <- iter + 1
      Uold <- U
      H <- centroidsFKM(data, U, n, k, p, m)
      D <- mixDistance(data, H, varStatus, weightStatus, contDist, n, k, range = range)
      U <- membershipDegrees(D, m, n, k, p)
      prova <- sum(abs(Uold - U)) > conv
    }
    func <- sum((U^m) * D)
    it[r] <- iter
    value[r] <- func
    if(!is.infinite(func)){
      if ((r == 1) | (func < funcOpt)){
        Uopt <- U
        Hopt <- H
        funcOpt <- func
      }
    }
  }
  ind <- indices(type = index, X = data, U = Uopt, H = Hopt, m = m, b = exp(1.0), alpha = alpha)
  if(index == "PE" || index == "XB"){
    indMax = -ind;
  }else{
    indMax = ind;
  }
  return(list( H = Hopt, U = Uopt, iter = it, value = value, index = ind, indexMax = indMax, k = k))
}

MFCM <- function (X, varStatus = NULL, weightStatus = NULL, contDist = "Euclidean",
                         k = c(2:6), m = 2, RS = 1, stand = 1, startU = NULL, 
                         index = "SIL.F", alpha = 1, conv = 1E-09, maxit = 1E+06, seed = NULL, range = NULL){
  X = data.matrix(X)
  n = nrow(X)
  p = ncol(X)
  nk <- length(k)
  
  set.seed(seed)
  
  Xraw = X
  rownames(Xraw) = rownames(X)
  colnames(Xraw) = colnames(X)
  if (stand == 1){ 
    X[,varStatus$cont] = scale(X[,varStatus$cont], center = TRUE, scale = TRUE)[, ]
  }

  crit.f <- array(NA, nk)
  crit = 0
  for (c in 1:nk){
    main.temp <- mainMFCM(data = X, varStatus = varStatus, weightStatus = weightStatus,
                                    contDist = contDist, m = m, index = index, alpha = alpha, 
                                    n = n, p = p, k = k[c], rs = RS, conv = conv, maxit = maxit, range = range)
    crit.temp = main.temp$indexMax
    crit.f[c] = main.temp$index
    if (c == 1 | crit < crit.temp) {
      main = main.temp
      crit = crit.temp
    }
  }
  value = as.vector(main$value)
  it = as.vector(main$iter)
  U.opt = main$U
  H.opt = main$H
  names(crit.f) = paste(index, " ", "k=", k, sep = "")
  k = main$k
  rownames(H.opt) = paste("Clus", 1:k, sep = " ")
  colnames(H.opt) = colnames(X)
  rownames(U.opt) = rownames(X)
  colnames(U.opt) = rownames(H.opt)
  names(value) = paste("Start", 1:RS, sep = " ")
  names(it) = names(value)
  names(k) = c("Number of clusters")
  names(m) = c("Parameter of fuzziness")
  if (stand != 1) 
    stand = 0
  names(stand) = c("Standardization (1=Yes, 0=No)")
  clus = cl.memb(U.opt)
  out = list(U = U.opt, H = H.opt, F = NULL, clus = clus, medoid = NULL, 
             value = value, criterion = crit.f, iter = it, k = k, 
             m = m, ent = NULL, b = NULL, vp = NULL, stand = stand, 
             Xca = X, X = Xraw, D = NULL, call = match.call())
  class(out) = c("fclust")
  return(out)
  
} 


nomOrdEntropy <- function(data, vars){
  n <- nrow(data)
  Ent <- array(NA, length(vars))
  for ( i in 1:length(vars)){
    distr <- table(data[,vars[i]])
    p <- distr/n
    k <- length(distr)
    Ent[i] <- (k-1)/2*log(2*pi*n*exp(1))-0.5*sum(log(p))
  }
  return(Ent)
}

calcEntropy <- function(data, contVars = NULL, binVars = NULL, nomVars = NULL, ordVars = NULL){
  n <- nrow(data)
  binEntAvr <- 0
  if (!is.null(binVars)){
    binEnt <- array(NA, length(binVars))
    for ( i in 1:length(binVars)){ 
      p <- sum(data[,binVars[i]])/n
      q <- 1-p
      binEnt[i] <- 0.5*log(2*pi*n*p*q*exp(1))
    }
    binEntAvr <- mean(binEnt)
  }
  
  nomEntAvr <- 0
  if (!is.null(nomVars)){
    nomEnt <- nomOrdEntropy(data, nomVars)
    nomEntAvr <- mean(nomEnt)
  }
  
  ordEntAvr <- 0
  if (!is.null(ordVars)){
    ordEnt <- array(NA, length(ordVars))
    ordEntX <- nomOrdEntropy(data, ordVars)
    for ( i in 1:length(ordVars)){
      N <- length(unique(data[,ordVars[i]])) 
      redOrdEnt <- log(N) + sum(log(choose(N-1,(1:N)-1)))/N - 0.5*(N-1) 
      ordEnt[i] <- ordEntX[i] - redOrdEnt
    }
    ordEntAvr <- mean(ordEnt)
  }
  
  contEntAvr <- 0
  if (!is.null(contVars)){
    contEnt <- array(NA, length(contVars))
    for ( i in 1:length(contVars)){
      nclass <- round(1+1.3*log(n))
      cont <- TRUE
      while (cont){ 
        distr <- hist(data[,contVars[i]], nclass = nclass, plot = FALSE)
        p <- distr$density
        nclass <- nclass - 1
        cont <- any(p == 0)
      }
      k <- nclass
      contEnt[i] <- (k-1)/2*log(2*pi*n*exp(1))-0.5*sum(log(p))
    }
    contEntAvr <- mean(contEnt)
  }

  totalEnt <- binEntAvr +  nomEntAvr +  ordEntAvr + contEntAvr
  
  AbinEnt <-  binEntAvr / totalEnt
  AnomEnt <-  nomEntAvr / totalEnt
  AordEnt <-  ordEntAvr / totalEnt
  AcontEnt <-  contEntAvr / totalEnt

  return(w = c(AcontEnt, AnomEnt, AordEnt, AbinEnt))
}

benchmarkClusteringMethods <- function(dat, k, index = "SIL.F", seed = 1234, distanceMatrix = NULL, 
                                       ranges = NULL, contVars = NULL, 
                                       nomVars = NULL, binVars = NULL, ordVars = NULL){
  
  distDat <- dist(scale(dat))
  tictoc::tic(quiet = TRUE)
  genMinkHierFit <- agnes(distanceMatrix, diss = TRUE,  method = "ward")
  genMinkHierClus <- cutree(genMinkHierFit, k = k)
  ss <- silhouette(genMinkHierClus, distanceMatrix)
  genMinkHier <- list(genMinkHierClus)
  genMinkHier$SIL <- mean(ss[, 3])
  genMinkHier$Stats <- cqcluster.stats(d = distanceMatrix, clustering = genMinkHierClus)
  genMinkHier$clusSizes <- table(factor(genMinkHierClus))
  
  a <- tictoc::toc(quiet = TRUE)
  cat("HwGM done in", a$toc - a$tic, "seconds....", "\n")
  
  tictoc::tic(quiet = TRUE)
  distDat2 <- gowdis((dat))
  mixedHclust.out <- hclust(distDat2, method = "ward.D2") 
  mixedClustHier <- cutree(mixedHclust.out, k = k) 
  ss <- silhouette(mixedClustHier, distDat2)
  GowMixedHier <- list(clustHier = mixedClustHier)
  GowMixedHier$SIL <- mean(ss[, 3])
  GowMixedHier$Stats <- cqcluster.stats(d = distDat2, clustering = mixedClustHier)
  GowMixedHier$clusSizes <- table(factor(mixedClustHier))
  
  a <- tictoc::toc(quiet = TRUE)
  cat("MHwG done in", a$toc - a$tic, "seconds ...", "\n")
  
  tictoc::tic(quiet = TRUE)
  dat2 <- dat
  dat2[,c(nomVars,binVars,ordVars)] <- as.character(dat2[,c(nomVars,binVars,ordVars)])
  mixedHierFit <-  hclustvar(X.quanti = dat2[,c(contVars)],X.quali =  as.matrix(dat2[,c(nomVars,binVars,ordVars)]))
  clstr <- cutreevar(mixedHierFit, k = k)
  clstrScores <- clstr$scores
  mixedHierClus <- apply(clstrScores,1,which.max)
  ss <- silhouette(mixedHierClus, distDat)
  mixedHier <- list(mixedHierClus)
  mixedHier$SIL <- mean(ss[, 3])
  mixedHier$Stats <- cqcluster.stats(d = distDat, clustering = mixedHierClus)
  mixedHier$clusSizes <- table(factor(mixedHierClus))
  
  a <- tictoc::toc(quiet = TRUE)
  cat("MixedHier done in", a$toc - a$tic, "seconds....", "\n")
  
  tictoc::tic(quiet = TRUE)
  clustPam <- cluster::pam(distDat2, k = k)
  ss <- silhouette(clustPam$clustering, distDat2)
  PAM <- clustPam
  PAM$SIL <- mean(ss[, 3])
  PAM$Stats <- cqcluster.stats(d = distDat2, clustering = PAM$clustering)
  PAM$clusSizes <- table(factor(clustPam$clustering))
  
  a <- tictoc::toc(quiet = TRUE)
  cat("PAMwG done in", a$toc - a$tic, "seconds ...", "\n")
  
  tictoc::tic(quiet = TRUE)
  catVars <- c(nomVars,ordVars,binVars)
  catDf <- data.frame(apply(as.matrix(dat[,catVars]), 2, factor), stringsAsFactors = TRUE)
  conDf <- data.frame(scale(dat[,contVars]), stringsAsFactors = TRUE)
  kamRes <- kamila(conVar = conDf, catFactor = catDf, numClust = k, numInit = 10)
  ss <- silhouette(kamRes$finalMemb, distDat)
  kamilaClust <- list(kamRes$finalMemb)
  kamilaClust$SIL <- mean(ss[, 3])
  kamilaClust$Stats <- cqcluster.stats(d = distDat, clustering = kamRes$finalMemb)
  kamilaClust$clusSizes <- table(factor(kamRes$finalMemb))

  a <- tictoc::toc(quiet = TRUE)
  cat("kamila done in", a$toc - a$tic, "seconds ...", "\n")
  
  return(list(PAMwG = PAM, HwGM = genMinkHier, MixedHier = mixedHier, MHwG = GowMixedHier, kamila = kamilaClust))
}

# Set DoAll = FALSE to only run MFCM.
runAllMethods <- function(data, c, m, varStatus, weightStatus, allRanges, contVars, nomVars, binVars, ordVars, 
                          distanceMatrix, scale = FALSE, scaleCont = FALSE, seed = 1234, DoAll = TRUE){
  tictoc::tic(quiet = TRUE)
  if (scale == TRUE){
    MFCMdata <- scale(data)
    if (scaleCont == TRUE){
      MFCMdata <- cbind(scale(data[,contVars]),data[,c(nomVars,binVars,ordVars)])
    }
  } else{
    MFCMdata <- data
  }
  
  clustMFCM <- MFCM(MFCMdata, varStatus = varStatus, weightStatus = weightStatus,
                    contDist = "Euclidean", k = c , m = m, RS = 1, index = "SIL.F", 
                    stand = 1, seed = seed, range = allRanges)
  a <- tictoc::toc(quiet = TRUE)
  cat("MFCM done in", a$toc - a$tic, "seconds ...", "\n")
  if (DoAll == FALSE){
    clsStats <- list()
    clsStats[[1]] <- cqcluster.stats(d = dist(data), clustering = as.integer(clustMFCM$clus[,1]))
    statsTake <- c(28,29,42,44)
    df <- data.frame( MFCM = clsStats[[1]][statsTake])
    SILs <- mean(unlist(clsStats[[1]][20]))
    df <- data.frame(SIL = SILs, df)
    dfClusSizes <- data.frame(MFCM = toString(table(factor(clustMFCM$clus[,1]))))
    df <- data.frame(df, ClustSize = t(dfClusSizes))
    return(list(allRes = df, clustMFCM = clustMFCM, clustStats = clsStats))
  }
  
  if (DoAll == TRUE){
    others <- benchmarkClusteringMethods(data, k = c, index = "SIL", seed = seed, 
                                       distanceMatrix = distanceMatrix, ranges = allRanges,
                                       contVars = contVars, nomVars = nomVars, binVars = binVars, 
                                       ordVars = ordVars)
    clsStats <- list()
    clsStats[[1]] <- cqcluster.stats(d = dist(data), clustering = as.integer(clustMFCM$clus[,1]))
    clsStats[[2]] <- others$HwGM$Stats
    clsStats[[3]] <- others$MixedHier$Stats
    clsStats[[4]] <- others$MHwG$Stats
    clsStats[[5]] <- others$kamila$Stats
    clsStats[[6]] <- others$PAMwG$Stats 
    
    df <- list(data.frame())
    statsTake <- c(28,29,42,44)
    rn <-  c("MFCM", "HwGM", "MixedHier", "MHwG", "kamila", "PAMwG")
    dfrn <- array()
    for (j in 1:length(statsTake)){
      df.temp <- data.frame( clsStats[[1]][statsTake[j]])
      dfrn <- rn[1]
      for (i in 2:6){
        if (!is.null(clsStats[[i]])){
          df.temp <- data.frame(df.temp, clsStats[[i]][statsTake[j]])
          dfrn <- cbind(dfrn,rn[i])
        }
      }
      df[[j]] <- df.temp
      df.temp <- NULL
    }
    df <- rbindlist(df,use.names=FALSE)
    df <- t(df)
    rownames(df) <- dfrn
    colnames(df) <- c("wb",	"ch", "dindex", "highdgap")
    SILs <- data.frame(MFCM = mean(unlist(clsStats[[1]][20])), HwGM= mean(unlist(clsStats[[2]][20])), 
                       MixedHier= mean(unlist(clsStats[[3]][20])), MHwG= mean(unlist(clsStats[[4]][20])), 
                       kamila = mean(unlist(clsStats[[5]][20])), PAMwG = mean(unlist(clsStats[[6]][20])))
    
    df <- data.frame(SIL = t(SILs), df)
    
    dfClusSizes <- data.frame(MFCM = toString(table(factor(clustMFCM$clus[,1]))), HwGM= toString(others$HwGM$clusSizes), 
                              MixedHier= toString(others$MixedHier$clusSizes), MHwG= toString(others$MHwG$clusSizes), 
                              kamila = toString(others$kamila$clusSizes), PAMwG = toString(others$PAMwG$clusSizes))
    df <- data.frame(df, ClustSize = t(dfClusSizes))
    return(list(allRes = df, clustMFCM = clustMFCM, benchmark = others, clustStats = clsStats))
  }
}
