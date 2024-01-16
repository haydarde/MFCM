library(fclust)
library(fpc)
# This script implements the method of "Ahmad, A., & Dey, L. (2005, December). Algorithm 
# for fuzzy clustering of mixed data with numeric and categorical attributes. In International 
# Conference on Distributed Computing and Internet Technology (pp. 561-572). Berlin, 
# Heidelberg: Springer Berlin Heidelberg.."


centroids <- function(data, U, n, k, p, m){
  out <- matrix(0, nrow = k, ncol = p)
  Unew <- t(U) 
  for ( i in 1:k){
    out[i, ] <- ((Unew[i,]^m ) %*% data)/sum(U[,i]^m)
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

weightedDistance <- function(data, w){
  n <- nrow(data)
  k <- ncol(data)
  if (any(abs(w)>10^10)){
    w <- (w-min(w))/(max(w)-min(w))
    data <- scale(data)
  }
  wDist <- matrix(0, nrow = n , ncol = n)
  for ( i in 1:(n-1)){
    for ( j in (i+1):n){
      wDist[i,j] <- sqrt(sum((w^2)*(data[i,]-data[j,])^2))
    }
  }
  wDist[lower.tri(wDist)] <- wDist[lower.tri(wDist)]
  wDist <- (wDist - min(wDist))/(max(wDist)-min(wDist))
  return(wDist)
}

calcEw <- function(data, distData, w, diff, lambda, beta, w_new = NULL){
  N <- nrow(data)
  if (is.null(N)){
    N <- length(N)
  }
  if (is.null(w_new)){
    w_new <- w - lambda * diff
  }
  wDist_new <- weightedDistance(data, w_new)
  Ew <- 0
  for (p in 1:(N-1)){
    for (q in (p+1):N){
      rho_pqw <- (1/(1+beta*wDist_new[p,q]))
      rho_pq1 <- (1/(1+beta*distData[p,q]))
      Ew <- Ew + 0.5*( rho_pqw * (1-rho_pq1) + rho_pq1 * (1-rho_pqw) )
    }
  }
  Ew <- 2 * Ew / (N*(N-1))
  return(Ew)
}


featureWeights <- function(data, k, step = 0.001, UL = 1000, convThreshold = 10^-9){
  # Implements the feature weights of "Yeung, D. S., & Wang, X. Z. (2002). Improving 
  # performance of similarity-based clustering by feature weight learning. IEEE transactions 
  # on pattern analysis and machine intelligence, 24(4), 556-561."
  
  # data: only continuous features.
  # k: number of continuous features.
  # step: step size to search for beta
  
  # Step 1:
  N <- nrow(data)
  if (is.null(N)){
    N <- length(data)
  }
  distData <- as.matrix(dist(data))
  distData <- (distData - min(distData))/(max(distData)-min(distData))
  difference <- 1
  beta <- 0
  iter <- 0
  differenceOld <- 1
  while ((difference > 0.00001) & (iter < 10000)){
    iter <- iter + 1
    
    sum1 <- 0
    for (p in 1:(N-1)){
      for (q in (p+1):N){
        sum1 <- sum1 + (1/(1+beta*distData[p,q]))
      }
    }
    difference <- abs(0.5 - ((2*sum1)/(N*(N-1))))
    # print(difference)
    if (differenceOld > difference){
      differenceOld <- difference
      beta <- beta + step
    } else {
      iter <- 10001
    }
    
    
  } # When we are out of this loop, beta will be specified
  
  w <- runif(k,0,10) # Step 2
  iter <- 0
  cont <- TRUE
  EwOld <- 0
  while(cont){
    iter <- iter + 1
    wDist <- weightedDistance(data, w)
    # Step 3:
    diff <- array(0, k)
    for ( j in 1:k){
      for (p in 1:(N-1)){
        for (q in (p+1):N){
          diff1 <- -beta / ((1+beta*wDist[p,q])^2)
          if (wDist[p,q] == 0){
            diff2 <- 0
          } else {
            diff2 <- (w[j]*(data[p,j]-data[q,j])^2)/wDist[p,q]
          }
          diff[j] <- diff[j] + ((1-(2/(1+beta*distData[p,q])))*diff1*diff2)/(N*(N-1))
          if (is.nan(diff[j])){
            diff[j]<-diff[j]
            }
        }
      }
    }
    
    # Step 4:
    cont <- TRUE

    #https://indrag49.github.io/Numerical-Optimization/solving-one-dimensional-optimization-problems.html#golden-section-search-method
    
    xl <- 0
    xr <- UL
    epsilon <- 0.01
    L <- xr-xl
    phi <- (sqrt(5)-1)/2
    x1 <- xl+(phi^2)*L
    x2 <- xl+phi*L
    while (L>epsilon){
      f1 <- calcEw(data = data, distData = distData, w = w, diff = diff, lambda = x1, beta = beta)
      f2 <- calcEw(data = data, distData = distData, w = w, diff = diff, lambda = x2, beta = beta)
      if (is.nan(f1) ){
        f1 <- f1
        print(diff)
        print("burada")
      } 
      if  (is.nan(f2)){
        f2 <- f2
      }
  
      if (f1 > f2){
        xl <- x1
        x1 <- x2
        L <- xr-x1
        x2 <- xl+phi*L
      } else {
        xr <- x2
        x2 <- x1
        L <- xr - xl
        x1 <- xl+(phi^2)*L
      }
    }
    lambda <- (xl+xr)/2
    
    
    # Step 5:
    candidate <- (w - lambda * diff)
    # print(candidate)
    w[which(candidate > 0)] <- candidate[which(candidate > 0)]
    # w <- (w-min(w))/(max(w)-min(w))
    Ew <- calcEw(data, distData, w, diff, lambda, beta, w_new = w)
    # print(abs(EwOld - Ew))
    if ((abs(EwOld - Ew) <= convThreshold) | (iter > 1000)) {
      cont <- FALSE
      cat("Final Ew:", abs(EwOld - Ew),"\n")
    } else {
      EwOld <- Ew
    }
    
  }
  return(w)
}

omega <- function(X, C, U, m){ 
  # Only one cluster will come in this function at a time. So, C and U are for cluster j.
  # U is a vector for cluster j.
  # C is a vector for cluster j.
  # X is a vector for a particular feature.
  Nc <- sum(U^m)
  tableX <- table(X)
  catNum <- length(tableX)
  n <- length(X)
  Nac <- array(NA, catNum)
  for ( p in 1:catNum){
    Nac[p] <- sum(U[which(X == as.numeric(names(tableX))[p])]^m) # N_{a,Aa,p,c}
  }
  
  delta <- matrix(NA, nrow = n , ncol = catNum)
  for ( p in 1:catNum){
    for (i in 1:n){
      if (X[i] == p){
        delta[i,p] <- 0 
      } else {
        delta[i,p] <- abs((Nac[which(as.numeric(names(tableX)) == X[i])] - X[i])/Nc) # C can be used instead of X in this calculation,.  
      }
    }
  }
  res <- array(0,n)
  for ( i in 1:n){
    res[i] <- sum(Nac*delta[i,])/Nc
    # for ( p in 1:catNum){
    #   res[i] <- res[i] + Nac[p]*delta[i,p]/Nc
    # }
  }
  return(res)
}

catContDistance <- function(data, H, U, w, contVars = NULL, catVars = NULL, m, n, k){
  out <- matrix(0, nrow = n, ncol = k)
  clusters <- cl.memb(U)
  for ( i in 1:n){
    for ( j in 1:k){
      if (!is.null(contVars)){
        sum1 <- 0
        for ( t in 1:length(contVars)){
          sum1 <- sum1 + (w[t]*(data[i,contVars[t]] - H[j,t]))^2
        }
      }
      out[i,j] <- sum1
    }
  }
  for ( j in 1:k){     
    if (!is.null(catVars)){
      sum2 <- array(0,length(which(clusters[,1] == j)))
      for ( t in catVars){
        sum2 <- sum2 + (omega(X = data[which(clusters[,1] == j),t], U = U[which(clusters[,1] == j),j], m = m))^2
      }
      out[which(clusters[,1] == j),j] <- sqrt(out[which(clusters[,1] == j),j] + sum2)
    }
    
  }
  return(out)
}

catContDistanceData <- function(data, H, U, w, contVars = NULL, catVars = NULL, m, n, k){
  out <- matrix(0, nrow = n, ncol = n)

  for ( i in 1:n){
    for ( j in 1:n){
      if (!is.null(contVars)){
        sum1 <- 0
        for ( t in 1:length(contVars)){
          sum1 <- sum1 + (w[t]*(data[i,contVars[t]] - H[j,contVars[t]]))^2
        }
      }
      out[i,j] <- sum1
    }
  }
  if (!is.null(catVars)){
    Nc <- sum(U^m)
    outCat <- list()
    countCat <- 0
    for ( j in 1:(length(catVars)-1)){
      for (r in (j+1):length(catVars)){
        countCat <- countCat + 1
        
        X <- data[,catVars[j]]
        Y <- data[,catVars[r]]
        
        tableX <- table(X)
        catNumX <- length(tableX)
        n <- length(X)
        NacX <- array(NA, catNumX)
        for ( p in 1:catNumX){
          NacX[p] <- sum(U[which(X == as.numeric(names(tableX))[p])]^m) # N_{a,Aa,p,c}
        }
        
        tableY <- table(Y)
        catNumY <- length(tableY)
        n <- length(X)
        NacY <- array(NA, catNumY)
        for ( p in 1:catNumY){
          NacY[p] <- sum(U[which(Y == as.numeric(names(tableY))[p])]^m) # N_{a,Aa,p,c}
        }
        outCat[[countCat]] <- matrix(NA, nrow = n, ncol = n)
        for (t in 1:n){
          for (s in 1:n){
            outCat[[countCat]][t,s] <- abs(NacX[which(as.numeric(names(tableX)) == data[t,catVars[j]])]-
                                             NacY[which(as.numeric(names(tableY)) == data[s,catVars[r]])])/Nc
          }
        }
      }
    }
    out2 <-  matrix(0, nrow = n, ncol = n)
    for ( i in 1:countCat){
      out2 <- out2 + outCat[[countCat]]
    }
  }
  out <- sqrt(out + out2)
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

mainADFCM <- function(data, contVars, catVars, m, n, p, k, rs, conv, maxit, index, alpha, UL, convThreshold){
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
      H <- centroids(data[,contVars], U, n, k, p = length(contVars), m)
      weights <- featureWeights(data = as.matrix(data[,contVars]), k =length(contVars), step = 0.001, UL = UL, convThreshold = convThreshold)
      weights <- weights/max(weights)
      D <- catContDistance(data, H, U, weights, contVars, catVars, m, n, k)
      U <- membershipDegrees(D, m, n, k, p)
      cat("ADFCM convergence trace:", sum(abs(Uold - U)),"\n")
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
  return(list( H = Hopt, U = Uopt, iter = it, value = value, index = ind, indexMax = indMax, k = k, w = weights))
}


ADFCM <- function (X, contVars, catVars, k = c(2:6), m = 2,  RS = 1, startU = NULL, 
                     index = "SIL.F", alpha = 1, conv = 1E-09, maxit = 1E+06, UL = 1000, convThreshold = 10^-9, seed = NULL){
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
    main.temp <- mainADFCM(data = X, contVars = contVars, catVars = catVars,
                             m = m, index = index, alpha = alpha, 
                             n = n, p = p, k = k[c], rs = RS, conv = conv, maxit = maxit,
                             UL = UL, convThreshold =convThreshold)
    
    
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
  colnames(H.opt) = colnames(X[,contVars])
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
             D = main$D, w = main$w, call = match.call())
  class(out) = c("fclust")
  return(out)
  
} 


