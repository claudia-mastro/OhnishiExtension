calcSigma <- function(N, phi, w) {
  Sigma <- matrix(0, nrow=N, ncol=N)
  for (k in 1:ncol(w)) {
    diff <- abs(outer(w[, k], w[, k], "-"))
    Sigma <- Sigma - phi[k] * diff
  }
  return(exp(Sigma))
}

calcRa<-function(N, R, h4, l5, h6, l6) {
  R_a <- rep(NA, N)
  R_a[R==1 | (R==5 & a_long >= l5) | (R==6 & a_long >= l6)] <- 1
  R_a[R==2 | (R==4 & a_long < h4) | (R==6 & a_long < h6)] <- 2
  R_a[R==3 | (R==4 & a_long >= h4) | (R==5 & a_long < l5) | (R==6 & (a_long >= h6) & (a_long < l6))] <- 3
  return(R_a)
}

CADE.A <- function(alpha, R, h4, l5, h6, l6, beta, sig2, mc_iter) {
  ## Calc R(alpha)
  
  R_a_long <- rep(0, sum(N))
  R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
  R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
  R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | ((R == 6 & alpha >= h6 & alpha < l6))]<-3 
  R_a <- list(0)
  R_a[[1]]<-R_a_long[1:N[1]]
  for(j in 2:J){
    R_a[[j]]<-R_a_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
  }
  eff <- rep(NA, mc_iter)
  for (s in 1:mc_iter) {
    ## Sample Z
    set.seed(i)
    Z<-list(0)
    for(j in 1:J){
      trt <- sample(1:N[j], size=alpha*N[j], replace=FALSE)
      Z[[j]] <- rep(0, N[j])
      Z[[j]][trt] <- 1
      T[j]<-sum(Z[[j]])
      a[j]<-T[j]/N[j]
    }
    ## Calc S
    D<-list(0)
    S<-rep(NA,
           times = J)
    for(j in 1:J){
      D[[j]]<-rep(NA,
                  times = N[j])
      D[[j]][R_a[[j]] == 1]<-1
      D[[j]][R_a[[j]] == 2]<-0
      D[[j]][R_a[[j]] == 3]<-Z[[j]][R_a[[j]] == 3]
      S[j]<-sum(D[[j]])
    }
    S_long <- rep(S, each=sum(N)/J)
    Y1 <- rep(NA, sum(N))
    Y0 <- rep(NA, sum(N))
    ij <- 0
    
    for (j in 1:J) {
      for (i in 1:N[J]) {
        ij <- ij + 1
        R.eff.a <- R_a[[j]][i]
        W0 <- W1 <- W[ij,]
        W0[3] <- W1[3] <- S_long[ij]
        W0[4] <- W1[4] <- alpha
        W0[5] <- 0
        W1[5] <- 1
        mu0<-W0%*%beta[R.eff.a,]
        var0<-sig2[R.eff.a]
        mu1<-W1%*%beta[R.eff.a,]
        var1<-sig2[R.eff.a]
        
        Y0[ij] <- rnorm(n = 1, mean = mu0, sd = sqrt(var0))
        
        Y1[ij]<-rnorm(n = 1, mean = mu1, sd = sqrt(var1))     
      }
    }
    eff[s] <- sum((Y1 - Y0)*(R_a_long==3))/sum(R_a_long==3)
  }
  return(mean(eff))
}

CADE.S <- function(s, R, h4, l5, h6, l6, beta, sig2) {
  
  eff <- rep(NA, length(unique(a)))
  i <- 0
  for (alpha in unique(a)) {
    ## Calc R(alpha)
    R_a_long <- rep(0, sum(N))
    R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
    R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
    R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | (R == 6 & (alpha >= h6) & (alpha < l6))]<-3 
    
    Y0 <- rep(NA, sum(N))
    Y1 <- rep(NA, sum(N))
    
    for (n in 1:sum(N)) {
      W0 <- W1 <- W[n,]
      W0[3] <- W1[3] <- s
      W0[4] <- W1[4] <- alpha
      W0[5] <- 0
      W1[5] <- 1
      
      mu0<-W0%*%beta[R_a_long[n],]
      var0<-sig2[R_a_long[n]]
      mu1<-W1%*%beta[R_a_long[n],]
      var1<-sig2[R_a_long[n]]
      
      Y0[n] <- rnorm(n = 1, mean = mu0, sd = sqrt(var0))
      
      Y1[n]<-rnorm(n = 1, mean = mu1, sd = sqrt(var1))   
    }
    i <- i + 1
    eff[i] <- (sum((Y1 - Y0)*(R_a_long==3)))/sum(R_a_long==3)
  }
  return(mean(eff))
}

CASE.S <- function(s, Sp, Z, g, R, h4, l5, h6, l6, beta, sig2) {
  
  eff <- rep(NA, length(a))
  i <- 0
  for (alpha in a) {
    ## Calc R(alpha)
    R_a_long <- rep(0, sum(N))
    R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
    R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
    R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | (R == 6 & (alpha >= h6) & (alpha < l6))]<-3 
    
    Y0 <- rep(NA, sum(N))
    Y1 <- rep(NA, sum(N))
    
    for (n in 1:sum(N)) {
      W0 <- W1 <- W[n,]
      W0[3] <- Sp
      W1[3] <- s
      W0[4] <- W1[4] <- alpha
      W0[5] <- Z
      W1[5] <- Z
      
      mu0<-W0%*%beta[R_a_long[n],]
      var0<-sig2[R_a_long[n]]
      mu1<-W1%*%beta[R_a_long[n],]
      var1<-sig2[R_a_long[n]]
      
      Y0[n] <- rnorm(n = 1, mean = mu0, 
                     sd = sqrt(var0))
      Y1[n] <- rnorm(n = 1, mean = mu1, 
                     sd = sqrt(var1))      
    }
    i <- i + 1
    eff[i] <- (sum((Y1 - Y0)*(R_a_long==g)))/sum(R_a_long==g)
  }
  return(mean(eff))
}

CADE.CASE <- function(alpha, s, Sp, Z, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2) {
  Y1 <- rep(NA, sum(N))
  Y0 <- rep(NA, sum(N))
  Y0p <- rep(NA, sum(N))
  C <- rep(NA, sum(N))
  
  Whl <- list(W_true[[1]], W_true[[1]], W_true[[1]], 
              cbind(W_true[[1]], h4), 
              cbind(W_true[[1]], l5),
              cbind(W_true[[1]], h6, l6))
  Sigma <- list(calcSigma(sum(N), phi[,1], Whl[[1]]),
                calcSigma(sum(N), phi[,2], Whl[[2]]),
                calcSigma(sum(N), phi[,3], Whl[[3]]),
                calcSigma(sum(N), phi[,4], Whl[[4]]),
                calcSigma(sum(N), phi[,5], Whl[[5]]),
                calcSigma(sum(N), phi[,6], Whl[[6]]))
  
  eff.a <- alpha
  eff.s <- s
  eff.sp <- Sp
  ij <- 0
  for (j in 1:J) {
    for (i in 1:N[J]) {
      ij <- ij + 1
      ## Need to figure out G(eff.a)
      if (R[ij] %in% 1:3) {
        G.eff.a <- R[ij]
      }  else if (R[ij]==4) {
        if (eff.a < h4[ij]) {
          G.eff.a <- 2
        } else {
          G.eff.a <- 3
        }
      } else if (R[ij]==5) {
        if (eff.a < l5[ij]) {
          G.eff.a <- 3
        } else {
          G.eff.a <- 1
        }
      } else if (R[ij]==6) {
        if (eff.a < h6[ij]) {
          G.eff.a <- 2
        } else if (eff.a < l6[ij]) {
          G.eff.a <- 3
        } else {
          G.eff.a <- 1
        }
      }
      
      if (G.eff.a == 3) {
        C[ij] <- 1
      } else {
        C[ij] <- 0
      }
      
      if (C[ij]==1) {
        W0p <- W0 <- W1 <- W[ij,]
        W0[2] <- W1[2] <- eff.s
        W0p[2] <- eff.sp
        W0p[3] <- W0[3] <- W1[3] <- eff.a
        W0p[4] <- W0[4] <- 0
        W1[4] <- 1
        
        Sig0 <- updateSigma(sum(N), phi[,R[ij]], W0, Sigma[[R[ij]]], ij)
        theta0 <- rmnorm(1, mean = 0, varcov = psi2[R[ij]]*Sig0)
        Y0[ij] <- rnorm(n = 1,
                        mean = mu[R[ij]] + theta0[ij],
                        sd = sqrt(sig2[R[ij]]))
        Sig1 <- updateSigma(sum(N), phi[,R[ij]], W1, Sigma[[R[ij]]], ij)
        theta1 <- rmnorm(1, mean = 0, varcov = psi2[R[ij]]*Sig1)
        Y1[ij] <- rnorm(n = 1,
                        mean = mu[R[ij]] + theta1[ij],
                        sd = sqrt(sig2[R[ij]]))
        
        Sig0p <- updateSigma(sum(N), phi[,R[ij]], W0p, Sigma[[R[ij]]], ij)
        theta0p <- rmnorm(1, mean = 0, varcov = psi2[R[ij]]*Sig0p)
        Y0p[ij] <- rnorm(n = 1,
                         mean = mu[R[ij]] + theta0p[ij],
                         sd = sqrt(sig2[R[ij]]))
      }
    }
  }
  CADE <- sum(C*(Y1-Y0), na.rm=TRUE)/sum(C)
  CASE <- sum(C*(Y0-Y0p), na.rm=TRUE)/sum(C)
  
  return(list(CADE, CASE))
}

