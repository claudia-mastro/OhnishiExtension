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

CADE.A <- function(alpha, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2, mc_iter) {
  ## Calc R(alpha)
  R_a_long <- rep(0, sum(N))
  R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
  R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
  R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | (R == 6 & (alpha >= h6) & (alpha < l6))]<-3 
  R_a <- list(0)
  R_a[[1]]<-R_a_long[1:N[1]]
  for(j in 2:J){
    R_a[[j]]<-R_a_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
  }
  eff <- rep(NA, mc_iter)
  for (i in 1:mc_iter) {
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
    Y0 <- rep(NA, sum(N))
    Y1 <- rep(NA, sum(N))
    Whl <- list(W_true[[1]], W_true[[1]], W_true[[1]], 
                cbind(W_true[[1]], h4), 
                cbind(W_true[[1]], l5),
                cbind(W_true[[1]], h6, l6))
    for (n in 1:sum(N)) {
      W0 <- W1 <- Whl[[R[n]]]
      W0[n,2] <- W1[n,2] <- S_long[n]
      W0[n,3] <- W1[n,3] <- alpha
      W0[n,4] <- 0
      W1[n,4] <- 1
      Sigma0 <- calcSigma(sum(N), phi[,R[n]], W0)
      Sigma1 <- calcSigma(sum(N), phi[,R[n]], W1)
      theta0 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig0)
      theta1 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig1)
      Y0[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                      sd = sqrt(sig2[R[n]]))
      Y1[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                     sd = sqrt(sig2[R[n]]))      
    }
    eff[i] <- sum((Y1 - Y0)*(R_a_long==3))/sum(R_a_long==3)
  }
  return(mean(eff))
}

CADE.S <- function(S, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2) {
  
  eff <- rep(NA, length(a))
  
  for (alpha in a) {
    ## Calc R(alpha)
    R_a_long <- rep(0, sum(N))
    R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
    R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
    R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | (R == 6 & (alpha >= h6) & (alpha < l6))]<-3 

    Y0 <- rep(NA, sum(N))
    Y1 <- rep(NA, sum(N))
    Whl <- list(W_true[[1]], W_true[[1]], W_true[[1]], 
                cbind(W_true[[1]], h4), 
                cbind(W_true[[1]], l5),
                cbind(W_true[[1]], h6, l6))
    for (n in 1:sum(N)) {
      W0 <- W1 <- Whl[[R[n]]]
      W0[n,2] <- W1[n,2] <- S
      W0[n,3] <- W1[n,3] <- alpha
      W0[n,4] <- 0
      W1[n,4] <- 1
      Sigma0 <- calcSigma(sum(N), phi[,R[n]], W0)
      Sigma1 <- calcSigma(sum(N), phi[,R[n]], W1)
      theta0 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig0)
      theta1 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig1)
      Y0[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                     sd = sqrt(sig2[R[n]]))
      Y1[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                     sd = sqrt(sig2[R[n]]))      
    }
    eff[i] <- (sum(Y1 - Y0)*R_a_long==3)/sum(R_a_long==3)
  }
  return(mean(eff))
}

CASE.S <- function(S, Sp, Z, G, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2) {
  
  eff <- rep(NA, length(a))
  
  for (alpha in a) {
    ## Calc R(alpha)
    R_a_long <- rep(0, sum(N))
    R_a_long[(R == 1) | ((R == 5) & (alpha >= l5)) | ((R == 6) & (alpha >= l6))]<-1
    R_a_long[(R == 2) | ((R == 4) & (alpha < h4)) | ((R == 6) & (alpha < h6))]<-2
    R_a_long[(R == 3) | ((R == 4) & (alpha >= h4)) | ((R == 5) & (alpha < l5)) | (R == 6 & (alpha >= h6) & (alpha < l6))]<-3 
    
    Y0 <- rep(NA, sum(N))
    Y1 <- rep(NA, sum(N))
    Whl <- list(W_true[[1]], W_true[[1]], W_true[[1]], 
                cbind(W_true[[1]], h4), 
                cbind(W_true[[1]], l5),
                cbind(W_true[[1]], h6, l6))
    for (n in 1:sum(N)) {
      W0 <- W1 <- Whl[[R[n]]]
      W0[n,2] <- S
      W1[n,2] <- Sp
      W0[n,3] <- W1[n,3] <- alpha
      W0[n,4] <- Z
      W1[n,4] <- Z
      Sigma0 <- calcSigma(sum(N), phi[,R[n]], W0)
      Sigma1 <- calcSigma(sum(N), phi[,R[n]], W1)
      theta0 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig0)
      theta1 <- rmnorm(1, mean = 0, varcov = psi2[R[n]]*Sig1)
      Y0[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                     sd = sqrt(sig2[R[n]]))
      Y1[n] <- rnorm(n = 1, mean = mu[R[n]] + theta0[n], 
                     sd = sqrt(sig2[R[n]]))      
    }
    eff[i] <- (sum(Y1 - Y0)*R_a_long==G)/sum(R_a_long==G)
  }
  return(mean(eff))
}

CADE.CASE <- function(alpha, S, Sp, Z, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2) {
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
  eff.s <- S
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
        W0p <- W0 <- W1 <- Whl[[R[ij]]]
        W0[ij,2] <- W1[ij,2] <- eff.s
        W0p[ij,2] <- eff.sp
        W0p[ij,3] <- W0[ij,3] <- W1[ij,3] <- eff.a
        W0p[ij,4] <- W0[ij,4] <- 0
        W1[ij,4] <- 1
  
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

CADE.CASE.ERT <- function(alpha, S, Sp, Z, R, h4, l5, h6, l6, phi, theta, mu, sig2, psi2) {
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
  eff.s <- S
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
        W0p <- W0 <- W1 <- Whl[[R[ij]]]
        W0[ij,2] <- W1[ij,2] <- eff.s
        W0p[ij,2] <- eff.sp
        W0p[ij,3] <- W0[ij,3] <- 0
        W1[ij,3] <- 1
        
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


