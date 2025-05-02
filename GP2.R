###########
# Packages
###########
library(mnormt)
library(matrixStats)
library(LaplacesDemon)
library(truncnorm)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
print(id)
source("~/project/OhnishiExtension/JWCode/Data_Simulation_GP2.R")
p <- 4
saveRDS(list(theta_true, sigma2_true, mu_true, R_long_true, phi_true, psi2_true, 
             h4_true, l5_true, h6_true, l6_true, gamma_true, 
             delta_h4_true, delta_l5_true, delta_h6_true, delta_l6_true, 
             tau2_h4_true, tau2_l5_true, tau2_h6_true, tau2_l6_true),
        paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/true_vals", 
               id, ".rds"))

#######
# Seed
#######
set.seed(id)

############
# Functions
############
logSumExp<-function(x){
  m<-max(x)
  val<-m + 
    log(sum(exp(x - m)))
  return(val)
}

# calcSigma<-function(N, phi, w) {
#   Sigma <- matrix(0, nrow=N, ncol=N)
#   for (i in 1:N) {
#     for (j in 1:N){
#       for (k in 1:ncol(w)) {
#           Sigma[i,j] = Sigma[i,j] - phi[k] * abs(w[i,k] - w[j,k])
#       }
#       Sigma[i,j] = exp(Sigma[i,j])
#     }
#   }
#   return(Sigma)
# }

# calcSigma <- function(N, phi, w) {
#   Sigma <- matrix(0, nrow=N, ncol=N)
#   diff <- matrix(0, nrow=N, ncol=N)
#   diff <- list(diff)[rep(1,ncol(w))]
#   for (k in 1:ncol(w)) {
#     for (i in 1:N) {
#       for (j in 1:N) {
#         if (i <= j) {
#           diff[[k]][i,j] = diff[[k]][j,i] = abs(w[i,k] - w[j,k])
#         }
#       }
#     }
#     Sigma <- Sigma - phi[k]*diff[[k]]
#   }
#   return(exp(Sigma))
# }

# calcSigma <- function(N, phi, w) {
#   Sigma <- matrix(0, nrow=N, ncol=N)
#   diff <- matrix(0, nrow=N, ncol=N)
#   diff <- list(diff)[rep(1,ncol(w))]
#   for (k in 1:ncol(w)) {
#     for (i in 1:N) {
#       for (j in i:N) {
#           diff[[k]][i,j] = diff[[k]][j,i] = abs(w[i,k] - w[j,k])
#       }
#     }
#     Sigma <- Sigma - phi[k]*diff[[k]]
#   }
#   return(exp(Sigma))
# }

calcSigma <- function(N, phi, w) {
  Sigma <- matrix(0, nrow=N, ncol=N)
  for (k in 1:ncol(w)) {
    diff <- abs(outer(w[, k], w[, k], "-"))
    Sigma <- Sigma - phi[k] * diff
  }
  return(exp(Sigma))
}

updateSigma<-function(N, phi, w, Sigma, j) {
  for (i in 1:N) {
    Sigma[i,j] = Sigma[j,i] = 0
    for (k in 1:ncol(w)) {
     Sigma[i,j] = Sigma[i,j] - phi[k] * abs(w[i,k] - w[j,k])
    }
    Sigma[i,j] = Sigma[j,i] = exp(Sigma[i,j])
  }
  return(Sigma)
}

# calcPi<-function(N, gamma, v) {
#   Pi <- matrix(0, nrow=N, ncol=6)
#     for (k in 1:6) {
#       Pi[,k] <- exp(v%*%gamma[k,])
#     }
#     Pi <- Pi/rowSums(Pi)
#   return(Pi)
# }

# calcPi<-function(N, gamma, v) {
#   Pi <- matrix(0, nrow=N, ncol=6)
#   for (k in 1:6) {
#     Pi[,k] <- exp(crossprod(t(v), gamma[k,]))
#   }
#   Pi <- Pi/rowSums(Pi)
#   return(Pi)
# }

calcPi<-function(N, gamma, v) {
  Pi <- exp(crossprod(t(v_long), t(gamma)))
  Pi <- Pi/rowSums(Pi)
  return(Pi)
}

calcRa<-function(N, R, h4, l5, h6, l6) {
  R_a <- rep(NA, N)
  R_a[R==1 | (R==5 & a_long >= l5) | (R==6 & a_long >= l6)] <- 1
  R_a[R==2 | (R==4 & a_long < h4) | (R==6 & a_long < h6)] <- 2
  R_a[R==3 | (R==4 & a_long >= h4) | (R==5 & a_long < l5) | (R==6 & (a_long >= h6) & (a_long < l6))] <- 3
  return(R_a)
}

##################
# Global Settings
##################
mcmc_samples <- 10000

mu_priormean <- 0
mu_priorvar <- 100^2
sig_priora <- 3
sig_priorb <- 2
psi_priora <- 3
psi_priorb <- 2
phi_priormean <- 0
phi_priorvar <- 1
gam_priormean <- 0
gam_priorvar <- 1
del_priormean <- 0
del_priorvar <- 1
tau_priora <- 3
tau_priorb <- 2

#############
# Parameters
#############

mu <- rep(NA, 6)
mu <- list(mu)[rep(1,mcmc_samples)]

theta <- matrix(data=NA, nrow=sum(N), ncol=6)
theta <- list(theta)[rep(1,mcmc_samples)]

sig2 <- rep(NA, 6)
sig2 <- list(sig2)[rep(1,mcmc_samples)]

psi2 <- rep(NA, 6)
psi2 <- list(psi2)[rep(1,mcmc_samples)]

phi <- matrix(data=NA, nrow=p+2, ncol=6)
phi <- list(phi)[rep(1,mcmc_samples)]

logit_h4 <- rep(NA, sum(N))
logit_h4 <- list(logit_h4)[rep(1,mcmc_samples)]

logit_l5 <- rep(NA, sum(N))
logit_l5 <- list(logit_l5)[rep(1,mcmc_samples)]

logit_h6 <- rep(NA, sum(N))
logit_h6 <- list(logit_h6)[rep(1,mcmc_samples)]

logit_l6 <- rep(NA, sum(N))
logit_l6 <- list(logit_l6)[rep(1,mcmc_samples)]

delta <- matrix(data=NA, nrow=4, ncol=2)
delta <- list(delta)[rep(1,mcmc_samples)]

tau2 <- rep(NA, 4)
tau2 <- list(tau2)[rep(1,mcmc_samples)]

R <- rep(NA, sum(N))
R <- list(R)[rep(1,mcmc_samples)]

gamma <- matrix(data=NA, nrow=6, ncol=2)
gamma <- list(gamma)[rep(1,mcmc_samples)]

CADE<-rep(NA, mcmc_samples)
CASE<-rep(NA, mcmc_samples)

#################
# Initial Values
#################

mu[[1]] <- rep(mean(Y_long), 6)
psi2[[1]] <- rinvgamma(6, shape=psi_priora, scale=psi_priorb)
sig2[[1]] <- rinvgamma(6, shape=sig_priora, scale=sig_priorb)
phi[[1]] <- matrix(runif(6*(p+2), 0, 1), nrow=p+2, ncol=6)
phi[[1]][5:6,1:3] <- 0
phi[[1]][6,4:5] <- 0

delta[[1]] <- matrix(rnorm(4*2, del_priormean, sqrt(del_priorvar)), nrow=4)
tau2[[1]] <- rinvgamma(4, shape=tau_priora, scale=tau_priorb)

logit_h4[[1]] <- rnorm(sum(N), crossprod(t(v_long), delta[[1]][1,]), sqrt(tau2[[1]][1]))
logit_l5[[1]] <- rnorm(sum(N), crossprod(t(v_long), delta[[1]][2,]), sqrt(tau2[[1]][2]))
logit_h6[[1]] <- rnorm(sum(N), crossprod(t(v_long), delta[[1]][3,]), sqrt(tau2[[1]][3]))
logit_l6[[1]] <- rnorm(sum(N), crossprod(t(v_long), delta[[1]][4,]), sqrt(tau2[[1]][4]))

h4 <- invlogit(logit_h4[[1]])
l5 <- invlogit(logit_l5[[1]])
h6 <- invlogit(logit_h6[[1]])
l6 <- (h6 + exp(logit_l6[[1]]))/(1+exp(logit_l6[[1]]))

W <- W_true[[1]]
Whl <- list(W, W, W, 
            cbind(W, h4), 
            cbind(W, l5),
            cbind(W, h6, l6))
Sigma <- list(calcSigma(sum(N), phi[[1]][,1], Whl[[1]]),
              calcSigma(sum(N), phi[[1]][,2], Whl[[2]]),
              calcSigma(sum(N), phi[[1]][,3], Whl[[3]]),
              calcSigma(sum(N), phi[[1]][,4], Whl[[4]]),
              calcSigma(sum(N), phi[[1]][,5], Whl[[5]]),
              calcSigma(sum(N), phi[[1]][,6], Whl[[6]]))

for (k in 1:6) {
  theta[[1]][,k] <- rmnorm(1, 0, psi2[[1]][k]*Sigma[[k]])
}

gamma[[1]] <- rbind(matrix(rnorm(5*2, gam_priormean, sqrt(gam_priorvar)), 
                           nrow=5, ncol=2), c(0,0))
Pi <- calcPi(sum(N), gamma[[1]], v_long)
R[[1]] <- sapply(c(1:sum(N)), 
                 function(idx){
                   sample(c(1:6), 
                          size = 1, 
                          prob = Pi[idx,],
                          replace = TRUE)
                 })

######################
# Metropolis Settings
######################

metrop_sd_phi <- rep(1, 6)
acctot_phi <- list(matrix(0, nrow=6, ncol=6))[rep(1,mcmc_samples)]

metrop_sd_gamma <- matrix(0.2, nrow = 5, ncol = 2)
acctot_gamma <- list(matrix(0, nrow=6, ncol=2))[rep(1,mcmc_samples)]

metrop_sd_h4 <- rep(0.5, times = sum(N))
acctot_h4 <- list(rep(0, 6))[rep(1,mcmc_samples)]

metrop_sd_l5 <- rep(0.5, times = sum(N))
acctot_l5 <- list(rep(0, 6))[rep(1,mcmc_samples)]

metrop_sd_h6 <- rep(0.5, times = sum(N))
acctot_h6 <- list(rep(0, 6))[rep(1,mcmc_samples)]

metrop_sd_l6 <- rep(0.5, times = sum(N))
acctot_l6 <- list(rep(0, 6))[rep(1,mcmc_samples)]

#####################
# Main Sampling Loop
#####################

for(s in 2:mcmc_samples){
  
  print(paste('iteration', s))
  
  ##########################
  # theta, mu, sigma2, psi2
  ##########################
  sig2R <- diag(sig2[[s-1]][R[[s-1]]])
  sig2R_inv <- diag(1.00/sig2[[s-1]][R[[s-1]]])
  theta[[s]] <- theta[[s-1]]
  mu[[s]] <- mu[[s-1]]
  sig2[[s]] <- sig2[[s-1]]
  gamma[[s]] <- gamma[[s-1]]
  for (k in 1:6) {
    ## THETA
    IR1 <- diag(R[[s-1]]==k)
    Sig_inv <- chol2inv(chol(Sigma[[k]]))
    Y_star1 <- Y_long - mu[[s]][R[[s-1]]] - theta[[s]][cbind(seq_along(R[[s-1]]), R[[s-1]])]*(R[[s-1]]!=k)
    
    theta_cov <- chol2inv(chol(crossprod(IR1, (sig2R_inv%*%IR1)) + Sig_inv/psi2[[s-1]][k]))
    theta_mean <- theta_cov%*%crossprod(IR1, (sig2R_inv%*%Y_star1))
    
    theta[[s]][,k] <- rmnorm(1, theta_mean, theta_cov)
    
    ## MU
    IR2 <- as.matrix(diag(IR1))
    Y_star2 <- Y_long - mu[[s]][R[[s-1]]]*(R[[s-1]]!=k) - theta[[s]][cbind(seq_along(R[[s-1]]), R[[s-1]])]
    
    mu_cov <- chol2inv(chol(crossprod(IR2, (sig2R_inv%*%IR2)) + 1/mu_priorvar))
    mu_mean <- mu_cov%*%crossprod(IR2, (sig2R_inv%*%Y_star2))

    mu[[s]][k] <- rnorm(1, mu_mean, sqrt(mu_cov))
    
    ## SIGMA2
    Y_k <- Y_long[R[[s-1]]==k]
    theta_k <- theta[[s]][R[[s-1]]==k,k]
    Y_star3 <- Y_k - mu[[s]][k] - theta_k
    
    sig_a <- sig_priora + sum(R[[s-1]]==k)/2
    sig_b <- sig_priorb + 1/2*crossprod(Y_star3)
    
    sig2[[s]][k] <- rinvgamma(1, sig_a, sig_b)
    sig2R <- diag(sig2[[s]][R[[s-1]]])
    sig2R_inv <- diag(1.00/sig2[[s]][R[[s-1]]])
    
    ## PSI2
    psi_a <- psi_priora + sum(N)/2
    psi_b <- psi_priorb + 1/2*crossprod(theta[[s]][,k], (Sig_inv%*%theta[[s]][,k]))
    
    psi2[[s]][k] <- rinvgamma(1, psi_a, psi_b)

    ## PHI
    phi[[s]][,k] <- phi[[s-1]][,k]

    for (h in 1:6) {
      if (h < 5 | (h==5 & k > 3) | (h==6 & k == 6)) {
        lphi <- log(phi[[s]][h,k])
        lphi_new <- rnorm(1, mean=lphi, sd=metrop_sd_phi)
        
        phi[[s]][h,k] <- exp(lphi_new)
  
        Sig_new <- calcSigma(sum(N), phi[[s]][,k], Whl[[k]])
  
        denom <- dmnorm(theta[[s]][,k], mean=0, varcov=psi2[[s]][k]*Sigma[[k]], log=TRUE) +
          dnorm(lphi, mean=phi_priormean, sd=sqrt(phi_priorvar), log=TRUE)
        numer <- dmnorm(theta[[s]][,k], mean=0, varcov=psi2[[s]][k]*Sig_new, log=TRUE) +
          dnorm(lphi_new, mean=phi_priormean, sd=sqrt(phi_priorvar), log=TRUE)
        
        
        ratio <- exp(numer - denom)
        uni_draw <- runif(n = 1, min = 0.00, max = 1.00)

        if(ratio > uni_draw){
          Sigma[[k]] <- Sig_new
          acctot_phi[[s]][h,k] <- acctot_phi[[s]][h,k] + 1
        } else {
          phi[[s]][h,k] <- phi[[s-1]][h,k]
        }
      }
    }
    
    ## GAMMA
    
    if (k!=6) {
      for (l in 1:2){ 
        gamma[[s]][k,l] <- rnorm(1, mean=gamma[[(s-1)]][k,l], sd=metrop_sd_gamma)
        Pi_new <- calcPi(sum(N), gamma[[s]], v_long)
        
        denom <- sum(log(Pi[cbind(seq_along(R[[s-1]]), R[[s-1]])])) + 
          dnorm(gamma[[s-1]][k,l], mean=gam_priormean, sd=sqrt(gam_priorvar), log=T)
        
        numer <- sum(log(Pi_new[cbind(seq_along(R[[s-1]]), R[[s-1]])])) + 
          dnorm(gamma[[s]][k,l], mean=gam_priormean, sd=sqrt(gam_priorvar), log=T)
        
        ratio<-exp(numer - denom)
        uni_draw<-runif(n = 1,
                        min = 0.00,
                        max = 1.00)
 
        if(ratio > uni_draw){
          Pi <- Pi_new
          acctot_gamma[[s]][k,l] <- acctot_gamma[[s]][k,l] + 1 
        } else {
          gamma[[s]][k,l] <- gamma[[s-1]][k,l]
        }
      }
    }
  }
 
  ## DELTA
  v_4 <- v_long[R[[s-1]]==4,,drop=FALSE]
  v_5 <- v_long[R[[s-1]]==5,,drop=FALSE]
  v_6 <- v_long[R[[s-1]]==6,,drop=FALSE]
  
  delta_cov_h4 <- chol2inv(chol(crossprod(v_4)/tau2[[s-1]][1] +
                                  diag(2)/del_priorvar))
  delta_mean_h4 <- crossprod(delta_cov_h4, 
                             crossprod(v_4, logit_h4[[s-1]][R[[s-1]]==4])/tau2[[s-1]][1])
  delta[[s]][1,] <- rmnorm(n = 1, mean = delta_mean_h4, varcov = delta_cov_h4)
  
  delta_cov_l5 <- chol2inv(chol(crossprod(v_5)/tau2[[s-1]][2] +
                                  diag(2)/del_priorvar))
  delta_mean_l5 <- crossprod(delta_cov_l5, 
                             crossprod(v_5, logit_l5[[s-1]][R[[s-1]]==5])/tau2[[s-1]][2])
  delta[[s]][2,] <- rmnorm(n = 1, mean = delta_mean_l5, varcov = delta_cov_l5)
  
  delta_cov_h6 <- chol2inv(chol(crossprod(v_6)/tau2[[s-1]][3] +
                                  diag(2)/del_priorvar))
  delta_mean_h6 <- crossprod(delta_cov_h6, 
                             crossprod(v_6, logit_h6[[s-1]][R[[s-1]]==6])/tau2[[s-1]][3])
  delta[[s]][3,] <- rmnorm(n = 1, mean = delta_mean_h6, varcov = delta_cov_h6)   
  
  delta_cov_l6 <- chol2inv(chol(crossprod(v_6)/tau2[[s-1]][4] +
                                  diag(2)/del_priorvar))
  delta_mean_l6 <- crossprod(delta_cov_l6, 
                             crossprod(v_6, logit_l6[[s-1]][R[[s-1]]==6])/tau2[[s-1]][4])
  delta[[s]][4,] <- rmnorm(n = 1, mean = delta_mean_l6, varcov = delta_cov_l6)   
  
  ## TAU
  n_4 <- sum(R[[s-1]]==4)
  n_5 <- sum(R[[s-1]]==5)
  n_6 <- sum(R[[s-1]]==6)
  
  tau_b_h4 <- crossprod(logit_h4[[s-1]]-crossprod(t(v_long), delta[[s]][1,]))/2 + tau_priorb
  tau2[[s]][1] <- rinvgamma(1, tau_priora + sum(N)/2, tau_b_h4)
  
  tau_b_l5 <- crossprod(logit_l5[[s-1]]-crossprod(t(v_long), delta[[s]][2,]))/2 + tau_priorb
  tau2[[s]][2] <- rinvgamma(1, tau_priora + sum(N)/2, tau_b_l5)
  
  tau_b_h6 <- crossprod(logit_h6[[s-1]]-crossprod(t(v_long), delta[[s]][3,]))/2 + tau_priorb
  tau2[[s]][3] <- rinvgamma(1, tau_priora + sum(N)/2, tau_b_h6)
  
  tau_b_l6 <- crossprod(logit_l6[[s-1]]-crossprod(t(v_long), delta[[s]][4,]))/2 + tau_priorb
  tau2[[s]][4] <- rinvgamma(1, tau_priora + sum(N)/2, tau_b_l6)
  
  ## H/L
  for (k in 1:6) {
    logit_h4_new <- rnorm(sum(R[[s-1]]==k), logit_h4[[s-1]][R[[s-1]]==k], metrop_sd_h4)
    h4_new <- invlogit(logit_h4_new)
    
    W_h4_new <- Whl[[4]]
    W_h4_new[R[[s-1]]==k,5] <- h4_new
    
    Sig_new <- calcSigma(sum(N), phi[[s]][,4], W_h4_new)
    
    Sig_inv_new <- chol2inv(chol(Sigma[[4]]))
    
    denom <- sum(dmnorm(theta[[s]][,4], mean=0, varcov=psi2[[s]][4]*Sigma[[4]], log=TRUE) +
      dnorm(logit_h4[[s-1]][R[[s-1]]==k], mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][1,]), sd=sqrt(tau2[[s]][1]), log=TRUE))
    numer <- sum(dmnorm(theta[[s]][,4], mean=0, varcov=psi2[[s]][4]*Sig_new, log=TRUE) +
      dnorm(logit_h4_new, mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][1,]), sd=sqrt(tau2[[s]][1]), log=TRUE))

    ratio <- exp(numer - denom)
    uni_draw <- runif(n = 1,
                    min = 0.00,
                    max = 1.00)
    
    if(ratio > uni_draw){
      logit_h4[[s]][R[[s-1]]==k] <- logit_h4_new
      h4[R[[s-1]]==k] <- h4_new
      Sigma[[4]] <- Sig_new
      Sig_inv <- Sig_inv_new
      Whl[[4]] <- W_h4_new
      acctot_h4[[s]][k] <- acctot_h4[[s]][k] + 1
    } else {
      logit_h4[[s]][R[[s-1]]==k] <- logit_h4[[s-1]][R[[s-1]]==k]
    }
    
    logit_l5_new <- rnorm(sum(R[[s-1]]==k), logit_l5[[s-1]][R[[s-1]]==k], metrop_sd_l5)
    l5_new <- invlogit(logit_l5_new)
    
    W_l5_new <- Whl[[5]]
    W_l5_new[R[[s-1]]==k,5] <- l5_new
    
    Sig_new <- calcSigma(sum(N), phi[[s]][,5], W_l5_new)
    
    Sig_inv_new <- chol2inv(chol(Sigma[[5]]))

    denom <- sum(dmnorm(theta[[s]][,5], mean=0, varcov=psi2[[s]][5]*Sigma[[5]], log=TRUE) +
      dnorm(logit_l5[[s-1]][R[[s-1]]==k], mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][2,]), sd=sqrt(tau2[[s]][2]), log=TRUE))
    numer <- sum(dmnorm(theta[[s]][,5], mean=0, varcov=psi2[[s]][5]*Sig_new, log=TRUE) +
      dnorm(logit_l5_new, mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][2,]), sd=sqrt(tau2[[s]][2]), log=TRUE))
    
    ratio <- exp(numer - denom)
    uni_draw <- runif(n = 1,
                      min = 0.00,
                      max = 1.00)
    
    if(ratio > uni_draw){
      logit_l5[[s]][R[[s-1]]==k] <- logit_l5_new
      l5[R[[s-1]]==k] <- l5_new
      Sigma[[5]] <- Sig_new
      Sig_inv <- Sig_inv_new
      Whl[[5]] <- W_l5_new
      acctot_l5[[s]][k] <- acctot_l5[[s]][k] + 1
    } else {
      logit_l5[[s]][R[[s-1]]==k] <- logit_l5[[s-1]][R[[s-1]]==k]
    }
    
    logit_h6_new <- rnorm(sum(R[[s-1]]==k), logit_h6[[s-1]][R[[s-1]]==k], metrop_sd_h6)
    h6_new <- invlogit(logit_h6_new)
    
    W_h6_new <- Whl[[6]]
    W_h6_new[R[[s-1]]==k,5] <- h6_new
    
    Sig_new <- calcSigma(sum(N), phi[[s]][,6], W_h6_new)
    
    Sig_inv_new <- chol2inv(chol(Sigma[[6]]))

    denom <- sum(dmnorm(theta[[s]][,6], mean=0, varcov=psi2[[s]][6]*Sigma[[6]], log=TRUE) +
      dnorm(logit_h6[[s-1]][R[[s-1]]==k], mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][3,]), sd=sqrt(tau2[[s]][3]), log=TRUE))
    numer <- sum(dmnorm(theta[[s]][,6], mean=0, varcov=psi2[[s]][6]*Sig_new, log=TRUE) +
      dnorm(logit_h6_new, mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][3,]), sd=sqrt(tau2[[s]][3]), log=TRUE))
    
    ratio <- exp(numer - denom)
    uni_draw <- runif(n = 1,
                      min = 0.00,
                      max = 1.00)

    if(ratio > uni_draw){
      logit_h6[[s]][R[[s-1]]==k] <- logit_h6_new
      h6[R[[s-1]]==k] <- h6_new
      Sigma[[6]] <- Sig_new
      Sig_inv <- Sig_inv_new
      Whl[[6]] <- W_h6_new
      acctot_h6[[s]][k] <- acctot_h6[[s]][k] + 1
    } else {
      logit_h6[[s]][R[[s-1]]==k] <- logit_h6[[s-1]][R[[s-1]]==k]
    }
    
    logit_l6_new <- rnorm(sum(R[[s-1]]==k), logit_l6[[s-1]][R[[s-1]]==k], metrop_sd_l6)
    l6_new <- (logit_h6[[s]][R[[s-1]]==k] + exp(logit_l6_new))/(1.00 + exp(logit_l6_new))
    
    W_l6_new <- Whl[[6]]
    W_l6_new[R[[s-1]]==k,6] <- l6_new
    
    Sig_new <- calcSigma(sum(N), phi[[s]][,6], W_l6_new)
    
    Sig_inv_new <- chol2inv(chol(Sigma[[6]]))

    denom <- sum(dmnorm(theta[[s]][,6], mean=0, varcov=psi2[[s]][6]*Sigma[[6]], log=TRUE) +
      dnorm(logit_l6[[s-1]][R[[s-1]]==k], mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][4,]), sd=sqrt(tau2[[s]][4]), log=TRUE))
    numer <- sum(dmnorm(theta[[s]][,6], mean=0, varcov=psi2[[s]][6]*Sig_new, log=TRUE) +
      dnorm(logit_l6_new, mean=crossprod(t(v_long[R[[s-1]]==k,]), delta[[s]][4,]), sd=sqrt(tau2[[s]][4]), log=TRUE))
    
    ratio <- exp(numer - denom)
    uni_draw <- runif(n = 1,
                      min = 0.00,
                      max = 1.00)
    
    if(ratio > uni_draw){
      logit_l6[[s]][R[[s-1]]==k] <- logit_l6_new
      l6[R[[s-1]]==k] <- l6_new
      Sigma[[6]] <- Sig_new
      Sig_inv <- Sig_inv_new
      Whl[[6]] <- W_l6_new
      acctot_l6[[s]][k] <- acctot_l6[[s]][k] + 1
    } else {
      logit_l6[[s]][R[[s-1]]==k] <- logit_l6[[s-1]][R[[s-1]]==k]
    }
  }
  
  ## R
  D_mat <- cbind(D_long==1, D_long==0, D_long==Z_long,
                 (calcRa(sum(N), rep(4, sum(N)), h4, l5, h6, l6)==2 & D_long==0) |
                  (calcRa(sum(N), rep(4, sum(N)), h4, l5, h6, l6)==3 & D_long==Z_long),
                 (calcRa(sum(N), rep(5, sum(N)), h4, l5, h6, l6)==1 & D_long==1) |
                   (calcRa(sum(N), rep(5, sum(N)), h4, l5, h6, l6)==3 & D_long==Z_long),
                 (calcRa(sum(N), rep(6, sum(N)), h4, l5, h6, l6)==1 & D_long==1) |
                   (calcRa(sum(N), rep(6, sum(N)), h4, l5, h6, l6)==2 & D_long==0) |
                   (calcRa(sum(N), rep(6, sum(N)), h4, l5, h6, l6)==3 & D_long==Z_long))
  
  numer <- matrix(NA, nrow=sum(N), ncol=6)
  for (k in 1:6) {
    numer[,k] <- dnorm(Y_long, mean=theta[[s]][,k] + mu[[s]][k], sd=sqrt(sig2[[s]][k]), log=T)
  }
  numer <- numer + log(D_mat) + log(Pi)
  probs <- sapply(c(1:6), function(k) 1.00/rowSums(exp(numer - numer[,k])))
  probs[is.na(probs) == 1]<-0.00
  
  R[[s]]<-sapply(c(1:sum(N)),
                 function(idx){
                   sample(c(1:6),
                          size = 1,
                          prob = probs[idx,],
                          replace = TRUE)
                 })

  ##########
  #Estimands
  ##########
  
  # Y1 <- rep(NA, sum(N))
  # Y0 <- rep(NA, sum(N))
  # Y0p <- rep(NA, sum(N))
  # C <- rep(NA, sum(N))
  # 
  # eff.a <- 0.8
  # eff.s <- 0.4
  # eff.sp <- 0.8
  # ij <- 0
  # for (j in 1:J) {
  #   for (i in 1:N[J]) {
  #     ij <- ij + 1
  #     ## Need to figure out G(eff.a)
  #     if (R[[s]][ij] %in% 1:3) {
  #       G.eff.a <- R[[s]][ij]
  #     }  else if (R[[s]][ij]==4) {
  #       if (eff.a < h4[ij]) {
  #         G.eff.a <- 2
  #       } else {
  #         G.eff.a <- 3
  #       }
  #     } else if (R[[s]][ij]==5) {
  #       if (eff.a < l5[ij]) {
  #         G.eff.a <- 3
  #       } else {
  #         G.eff.a <- 1
  #       }
  #     } else if (R[[s]][ij]==6) {
  #       if (eff.a < h6[ij]) {
  #         G.eff.a <- 2
  #       } else if (eff.a < l6[ij]) {
  #         G.eff.a <- 3
  #       } else {
  #         G.eff.a <- 1
  #       }
  #     }
  # 
  #     if (G.eff.a == 3) {
  #       C[ij] <- 1
  #     } else {
  #       C[ij] <- 0
  #     }
  # 
  #     if (C[ij]==1) {
  #       W0p <- W0 <- W1 <- Whl[[R[[s]][ij]]]
  #       W0[ij,2] <- W1[ij,2] <- eff.s
  #       W0p[ij,2] <- eff.sp
  #       W0p[ij,3] <- W0[ij,3] <- W1[ij,3] <- eff.a
  #       W0p[ij,4] <- W0[ij,4] <- 0
  #       W1[ij,4] <- 1
  # 
  #       Sig0 <- updateSigma(sum(N), phi[[s]][,R[[s]][ij]], W0, Sigma[[R[[s]][ij]]], ij)
  #       theta0 <- rmnorm(1, mean = 0, varcov = psi2[[s]][R[[s]][ij]]*Sig0)
  #       Y0[ij] <- rnorm(n = 1,
  #                     mean = mu[[s]][R[[s]][ij]] + theta0[ij],
  #                     sd = sqrt(sig2[[s]][R[[s]][ij]]))
  #       Sig1 <- updateSigma(sum(N), phi[[s]][,R[[s]][ij]], W1, Sigma[[R[[s]][ij]]], ij)
  #       theta1 <- rmnorm(1, mean = 0, varcov = psi2[[s]][R[[s]][ij]]*Sig1)
  #       Y1[ij] <- rnorm(n = 1,
  #                     mean = mu[[s]][R[[s]][ij]] + theta1[ij],
  #                     sd = sqrt(sig2[[s]][R[[s]][ij]]))
  # 
  #       Sig0p <- updateSigma(sum(N), phi[[s]][,R[[s]][ij]], W0p, Sigma[[R[[s]][ij]]], ij)
  #       theta0p <- rmnorm(1, mean = 0, varcov = psi2[[s]][R[[s]][ij]]*Sig0p)
  #       Y0p[ij] <- rnorm(n = 1,
  #                      mean = mu[[s]][R[[s]][ij]] + theta0p[ij],
  #                      sd = sqrt(sig2[[s]][R[[s]][ij]]))
  #     }
  #   }
  # }
  # CADE[s] <- sum(C*(Y1-Y0), na.rm=TRUE)/sum(C)
  # CASE[s] <- sum(C*(Y0-Y0p), na.rm=TRUE)/sum(C)
  
  if (s %in% c(1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)){
    print("SAVING")
    #write.csv(cbind(CADE, CASE), "/home/cim24/project/OhnishiExtension/Results/GP2_4.10/CADE_CASE.csv")
    saveRDS(theta, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/theta",
                          id, ".rds"))
    saveRDS(sig2, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/sig2",
                         id, ".rds"))
    saveRDS(mu, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/mu",
                       id, ".rds"))
    saveRDS(R, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/R",
                      id, ".rds"))
    saveRDS(phi, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/phi",
                        id, ".rds"))
    saveRDS(psi2, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/psi2",
                         id, ".rds"))
    saveRDS(logit_h4, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/h4",
                             id, ".rds"))
    saveRDS(logit_l5, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/l5",
                             id, ".rds"))
    saveRDS(logit_h6, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/h6",
                             id, ".rds"))
    saveRDS(logit_l6,paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/l6",
                            id, ".rds"))
    saveRDS(gamma, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/gamma",
                          id, ".rds"))
    saveRDS(Pi, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/Pi",
                       id, ".rds"))
    saveRDS(delta, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/delta",
                       id, ".rds"))
    saveRDS(tau2, paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/tau2",
                       id, ".rds"))
    saveRDS(list(acctot_phi, acctot_gamma, acctot_h4, acctot_l5, acctot_h6, acctot_l6),
            paste0("/home/cim24/project/OhnishiExtension/Results/GP2_4.10/acctots", 
                   id, ".rds"))
  }
}




