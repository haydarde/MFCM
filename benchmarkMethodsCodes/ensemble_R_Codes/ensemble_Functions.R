library(fclust)
library(fastDummies)

# This script implements the ensamble method of "Suguna, J. and Selve, A. (2012). Ensemble fuzzy clustering for mixed 
# numeric and categorical data. International Journal of Computer Applications, 42(3), 19-23."


modefunc <- function(x){
  tabresult <- tabulate(x)
  mode <- which(tabresult == max(tabresult))
  if(sum(tabresult == max(tabresult))>1){
    mode <- mean(mode)
  } 
  return(mode)
}

centroids <- function(data, U, n, k, p, m){
  out <- matrix(0, nrow = k, ncol = p)
  temp.cluster <- cl.memb(U)[,1]
  for ( i in 1:k){
    out[i, ] <- apply(data[which(temp.cluster == i),], 2, modefunc)
  }
  return(out)
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

modeDistance <- function(data, H, n, k){
  out <- matrix(0, nrow = n, ncol = k)
  for ( i in 1:n){
    for ( j in 1:k){
      out[i,j] <- sum(data[i,] != H[j,])
    }
  }
  return(out)
}

modeDistanceData <- function(data, H, contVars, n, k){
  out <- matrix(0, nrow = n, ncol = k)
  catVars <- setdiff(c(1:ncol(data)) , contVars)
  for ( i in 1:n){
    for ( j in 1:k){
      sum1 <- 0
      for ( r in contVars){ 
        sum1 <- sum1 + sum((data[i,r] - H[j,r])^2)
      }
      sum2 <- 0 
      for ( r in catVars){ 
        sum2 <- sum2 + sum(data[i,r] == H[j,r])
      }
      out[i,j] <- sqrt(sum1) + sum2
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

mainFcModes <- function(data, m, n, p, k, rs, conv, maxit, index, alpha){
  # ###############################################
  # # Run Armadillo accu() function
  # # https://stackoverflow.com/questions/38589705/difference-between-rs-sum-and-armadillos-accu
  # require(RcppArmadillo)
  # suppressMessages(require(inline))
  # code <- 'arma::vec ax = Rcpp::as<arma::vec>(x);
  #          return Rcpp::wrap(arma::accu(ax));'
  # ## create the compiled function
  # armasum <- cxxfunction(signature(x="numeric"),
  #                        code,plugin="RcppArmadillo")
  # ###############################################
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
      H <- centroids(data, U, n, k, p, m)
      D <- modeDistance(data, H, n, k)
      U <- membershipDegrees(D, m, n, k, p)
      print(sum(abs(Uold - U)))
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


FcModes <- function (X, k = c(2:6), m = 2,  RS = 1, startU = NULL, 
                     index = "SIL.F", alpha = 1, conv = 1E-09, maxit = 1E+06, seed = NULL){
  X = data.matrix(X)
  n = nrow(X)
  p = ncol(X)
  # k = 2:6
  nk <- length(k)
  
  set.seed(seed)
  
  Xraw = X
  rownames(Xraw) = rownames(X)
  colnames(Xraw) = colnames(X)

  crit.f <- array(NA, nk)
  # crit.f <- rep(NA, nk)
  crit = 0
  for (c in 1:nk) {
      main.temp <- mainFcModes(data = X, m = m, index = index, alpha = alpha, 
                              n = n, p = p, k = k[c], rs = RS, conv = conv, maxit = maxit)
 
    
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
  clus = cl.memb(U.opt)
  out = list(U = U.opt, H = H.opt, F = NULL, clus = clus,
             value = value, criterion = crit.f, iter = it, k = k, 
             m = m, ent = NULL, b = NULL, vp = NULL, Xca = X, X = Xraw, 
             D = main$D, call = match.call())
  class(out) = c("fclust")
  return(out)
  
} 


ensemble_FC <- function(data, contVars, binVars, c = 2, m = 2, distDat, index = "SIL.F", stand = 1, seed = 12321){
  
  catVars <- setdiff(c(1:ncol(data)) , c(contVars))
  
  catData <- data[,-contVars]
  
  # Fuzzy C-Means clustering algorithm is applied for numeric data
  clustFCM <- FKM(data[,contVars], k=c, m=m, RS=1, index = index, stand = stand, seed = seed)
  clustFcModes <- FcModes(X = catData, k = c, m = m, seed = seed)
  ensembleU <- (clustFCM$U + clustFcModes$U)/2
  ensembleClust <- cl.memb(ensembleU)
  
  return(list(u = ensembleU, clust = ensembleClust))
}
  