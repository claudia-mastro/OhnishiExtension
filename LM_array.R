###############
#Packages
###############
library(mnormt)
library(matrixStats)
args <- commandArgs(trailingOnly = TRUE)
id <- as.integer(args[1])
print(id)
J <- as.integer(args[2])
print(J)
Nj <- as.integer(args[3])
print(Nj)
v <- paste0("LM_4.28_J", J, "Nj", Nj)
print(v)
source("~/project/OhnishiExtension/JWCode/Data_Simulation_LM.R")

#########
#Function
#########
logSumExp<-function(x){
  m<-max(x)
  val<-m + 
    log(sum(exp(x - m)))
  return(val)
  
}

##############
#Seed
##############
set.seed(id)

################
#Global Settings
################
mcmc_samples<-100000
burnin <- 50000
thin <- 50
iters <- burnin:mcmc_samples
iters <- iters[seq(1, mcmc_samples-burnin + thin, thin)]

sigma2_beta<-100.00^2
a_sigma2<-0.01
b_sigma2<-0.01
sigma2_alpha<-1.00^2
sigma2_delta<-1.00^2
a_tau2<-0.01
b_tau2<-0.01

shape_tau2_update<-sum(N)/2.00 +
  a_tau2

###########
#Parameters
###########
beta <- matrix(NA, nrow=3, ncol=5)
beta <- list(beta)[rep(1,mcmc_samples)]

sigma2 <- rep(NA, 3)
sigma2 <- list(sigma2)[rep(1,mcmc_samples)]

G <- matrix(NA, nrow=J, ncol=Nj)
G <- list(G)[rep(1,mcmc_samples)]

alpha <- matrix(NA, nrow=2, ncol=5)
alpha <- list(alpha)[rep(1,mcmc_samples)]

h0 <- l0 <- h1 <- l1 <- rep(NA, sum(N))
h0 <- list(h0)[rep(1,mcmc_samples)]
l0 <- list(l0)[rep(1,mcmc_samples)]
h1 <- list(h1)[rep(1,mcmc_samples)]
l1 <- list(l1)[rep(1,mcmc_samples)]

delta_h0 <- delta_l0 <- delta_h1 <- delta_l1 <- rep(NA, 2)
delta_h0 <- list(delta_h0)[rep(1,mcmc_samples)]
delta_l0 <- list(delta_l0)[rep(1,mcmc_samples)]
delta_h1 <- list(delta_h1)[rep(1,mcmc_samples)]
delta_l1 <- list(delta_l1)[rep(1,mcmc_samples)]

tau2_h0 <- tau2_l0 <- tau2_h1 <- tau2_l1 <- NA
tau2_h0 <- list(tau2_h0)[rep(1,mcmc_samples)]
tau2_h1 <- list(tau2_h1)[rep(1,mcmc_samples)]
tau2_l0 <- list(tau2_l0)[rep(1,mcmc_samples)]
tau2_l1 <- list(tau2_l1)[rep(1,mcmc_samples)]

CADE<-rep(NA, mcmc_samples)
CASE<-rep(NA, mcmc_samples)

###############
#Initial Values
###############
for(k in 1:3){
  beta[[1]][k,]<- rep(0.00, 5)
}

sigma2[[1]] <- rep(1.00, 3)

log_pi_mat_temp <- matrix(NA, nrow = sum(N), ncol = 6)
log_pi_mat_temp[,6] <- rep(0.00, sum(N))

for(k in 1:5){
  alpha[[1]][,k] <- rep(0.00, 2)
  log_pi_mat_temp[,k] <- q_long%*%alpha[[1]][,k]
}
log_pi_mat <- pi_mat<-matrix(NA, nrow = sum(N), ncol = 6)

for(k in 1:6){
  log_pi_mat[,k] <- log_pi_mat_temp[,k] -
    apply(log_pi_mat_temp, 1, logSumExp) 
  pi_mat[,k] <- 1.00/rowSums(exp(log_pi_mat_temp - log_pi_mat_temp[,k]))
}

for(j in 1:J){
  for(i in 1:N[j]){
    
    if(Z[[j]][i] == 1 & D[[j]][i] == 1){
      G[[1]][j, i] <- sample(c(1, 3, 4, 5, 6), size = 1, replace = TRUE)
    }
    
    if(Z[[j]][i] == 1 & D[[j]][i] == 0){
      G[[1]][j, i] <- sample(c(2, 4, 6), size = 1, replace = TRUE)
    }
    
    if(Z[[j]][i] == 0 & D[[j]][i] == 1){
      G[[1]][j, i]<-sample(c(1, 5, 6), size = 1, replace = TRUE)
    }
    
    if(Z[[j]][i] == 0 & D[[j]][i] == 0){
      G[[1]][j, i]<-sample(c(2, 3, 5, 6), size = 1, replace = TRUE)
    }
    
  }
}

G_long<-G[[1]][1,]
for(j in 2:J){
  G_long <- c(G_long,
              G[[1]][j,])
}

delta_h0[[1]] <- rep(0.00, 2)
tau2_h0[[1]] <- 0.01
logit_h0 <- rnorm(n = sum(N),
                  mean = (v_long%*%delta_h0[[1]]),
                  sd = sqrt(tau2_h0[[1]]))
h0[[1]] <- 1.00/(1.00 + exp(-logit_h0))

delta_l0[[1]] <- rep(0.00, 2)
tau2_l0[[1]] <- 0.01
logit_l0 <- rnorm(n = sum(N),
                  mean = (v_long%*%delta_l0[[1]]),
                  sd = sqrt(tau2_l0[[1]]))
l0[[1]] <- 1.00/(1.00 + exp(-logit_l0))

delta_h1[[1]] <- rep(0.00, 2)
tau2_h1[[1]] <- 0.01
logit_h1 <- rnorm(n = sum(N),
                  mean = (v_long%*%delta_h1[[1]]),
                  sd = sqrt(tau2_h1[[1]]))
h1[[1]] <- 1.00/(1.00 + exp(-logit_h1))

delta_l1[[1]] <- rep(0.00, 2)
tau2_l1[[1]] <- 0.01
logit_l1 <- rnorm(n = sum(N),
                  mean = (v_long%*%delta_l1[[1]]),
                  sd = sqrt(tau2_l1[[1]]))
l1[[1]] <- (h1[[1]] + exp(logit_l1))/(1.00 + exp(logit_l1))

G_a_long<-rep(NA, times = sum(N))
G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[1]])) | ((G_long == 6) & (a_long >= l1[[1]]))] <- 1
G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[1]])) | ((G_long == 6) & (a_long < h1[[1]]))] <- 2
G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[1]])) | ((G_long == 5) & (a_long < l0[[1]])) | ((G_long == 6) & (a_long >= h1[[1]]) & (a_long < l1[[1]]))] <- 3

####################
#Metropolis Settings
####################
metrop_sd_alpha <- matrix(0.50, nrow = 5, ncol = 2)
acctot_alpha <- matrix(1, nrow = 5, ncol = 2)

metrop_sd_h0 <- rep(1.00, times = sum(N))
acctot_h0 <- rep(1, times = sum(N))

metrop_sd_l0 <- rep(1.00, times = sum(N))
acctot_l0 <- rep(1, times = sum(N))

metrop_sd_h1 <- rep(1.00, times = sum(N))
acctot_h1 <- rep(1, times = sum(N))

metrop_sd_l1 <- rep(1.00, times = sum(N))
acctot_l1 <- rep(1, times = sum(N))

###################
#Main Sampling Loop
###################
for(s in 2:mcmc_samples){ 
  
  #############
  #beta, sigma2
  #############
  for(k in 1:3){
    Y_k <- Y_long[G_a_long == k]
    W_k <- W[(G_a_long == k),]
    cov_beta <- chol2inv(chol(crossprod(W_k)/sigma2[[s-1]][k] + diag(5)/sigma2_beta))
    mu_beta <- cov_beta%*%(crossprod(W_k, Y_k))/sigma2[[s-1]][k]
    beta[[s]][k,] <- rmnorm(n = 1, mean = mu_beta, varcov = cov_beta)
    
    rate <- crossprod(Y_k - W_k%*%beta[[s]][k,])/2.00 + b_sigma2
    shape <- length(Y_k)/2.00 + a_sigma2
    sigma2[[s]][k] <- 1.00/rgamma(n = 1, shape = shape, rate = rate)
  }
  
  #######
  #G, G_a
  #######
  num <- matrix(NA, nrow = sum(N), ncol = 6)
  
  for(k in 1:6) {
    G_long_temp<-rep(k, times = sum(N))
    G_a_long_temp<-rep(NA, times = sum(N))
    G_a_long_temp[(G_long_temp == 1) | ((G_long_temp == 5) & (a_long >= l0[[s-1]])) | ((G_long_temp == 6) & (a_long >= l1[[s-1]]))] <- 1
    G_a_long_temp[(G_long_temp == 2) | ((G_long_temp == 4) & (a_long < h0[[s-1]])) | ((G_long_temp == 6) & (a_long < h1[[s-1]]))] <- 2
    G_a_long_temp[(G_long_temp == 3) | ((G_long_temp == 4) & (a_long >= h0[[s-1]])) | ((G_long_temp == 5) & (a_long < l0[[s-1]])) | ((G_long_temp == 6) & (a_long >= h1[[s-1]]) & (a_long < l1[[s-1]]))] <- 3
    
    mu <- rep(NA, times = sum(N))
    var <- rep(NA,times = sum(N))
    
    for(j in 1:sum(N)) {
      mu[j] <- W[j,]%*%beta[[s]][G_a_long_temp[j],]
      var[j] <- sigma2[[s]][G_a_long_temp[j]]
    }
    
    num[,k] <- dnorm(x = Y_long, mean = mu, sd = sqrt(var), log = TRUE)
    if(k < 6) {
      num[,k] <- num[,k] + q_long%*%alpha[[s-1]][,k]
    }
    
    num[((D_long == 0) & (G_a_long_temp == 1)), k] <- -Inf
    num[((D_long == 1) & (G_a_long_temp == 2)), k] <- -Inf
    num[((D_long != Z_long) & (G_a_long_temp == 3)), k] <- -Inf
  }
  
  probs<-matrix(NA, nrow = sum(N), ncol = 6)
  for(k in 1:6) {
    probs[,k]<-1.00/rowSums(exp(num - num[,k]))
  }
  probs[is.na(probs) == 1] <- 0.00
  
  G_long <- sapply(c(1:sum(N)), 
                   function(idx){
                     sample(c(1:6), 
                            size = 1, 
                            prob = probs[idx,],
                            replace = TRUE)
                   })
  
  G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[s-1]])) | ((G_long == 6) & (a_long >= l1[[s-1]]))] <- 1
  G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[s-1]])) | ((G_long == 6) & (a_long < h1[[s-1]]))] <- 2
  G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[s-1]])) | ((G_long == 5) & (a_long < l0[[s-1]])) | ((G_long == 6) & (a_long >= h1[[s-1]]) & (a_long < l1[[s-1]]))] <- 3
  
  G[[s]][1,]<-G_long[1:N[1]]
  for(j in 2:J){
    G[[s]][j,]<-G_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
  }
  
  ######
  #alpha
  ######
  for(k in 1:5) {
    alpha[[s]][,k] <- alpha[[s-1]][,k]
    for(l in 1:2) {
      log_pi_mat_temp_old <- log_pi_mat_temp
      log_pi_mat_old <- log_pi_mat
      pi_mat_old <- pi_mat
      denom <- sum(log_pi_mat_old[cbind(seq_along(G_long), G_long)]) +
        dnorm(x = alpha[[s-1]][l, k],
              mean = 0.00,
              sd = sqrt(sigma2_alpha),
              log = TRUE)
      
      alpha[[s]][l,k] <- rnorm(n = 1,
                             mean = alpha[[s-1]][l, k],
                             sd = metrop_sd_alpha[k,l])
      
      log_pi_mat_temp[,k] <- q_long%*%alpha[[s]][,k]
      
      for(m in 1:6) {
        log_pi_mat[,m] <- log_pi_mat_temp[,m] -
          apply(log_pi_mat_temp, 1, logSumExp)
        pi_mat[,m] <- 1.00/rowSums(exp(log_pi_mat_temp - log_pi_mat_temp[,m]))
      }
      
      numer <- sum(log_pi_mat[cbind(seq_along(G_long), G_long)]) +
        dnorm(x = alpha[[s]][l, k],
              mean = 0.00,
              sd = sqrt(sigma2_alpha),
              log = TRUE)
      
      accept <- 1
      ratio <- exp(numer - denom)
      uni_draw <- runif(n = 1, min = 0.00, max = 1.00)
      
      if(ratio < uni_draw){
        log_pi_mat_temp<-log_pi_mat_temp_old
        log_pi_mat<-log_pi_mat_old
        pi_mat<-pi_mat_old
        alpha[[s]][l, k] <- alpha[[s-1]][l, k]
        accept<-0
      }
      
      acctot_alpha[k,l] <- acctot_alpha[k,l] + accept
    }
  }
  
  ###
  #h0
  ###
  logit_h0_old <- logit_h0
  G_a_long_old <- G_a_long
  mu_old <- sapply(1:sum(N), 
                 function(j){
                   W[j,]%*%beta[[s]][G_a_long_old[j],]
                 })
  
  var_old <- sigma2[[s]][G_a_long_old]
  denom <- dnorm(x = Y_long,
               mean = mu_old,
               sd = sqrt(var_old),
               log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
    dnorm(x = logit_h0_old,
          mean = (v_long%*%delta_h0[[s-1]]),
          sd = sqrt(tau2_h0[[s-1]]),
          log = TRUE)
  
  logit_h0 <- rnorm(n = sum(N), mean = logit_h0_old, sd = metrop_sd_h0)
  h0[[s]] <- 1.00/(1.00 + exp(-logit_h0))
  G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[s-1]])) | ((G_long == 6) & (a_long >= l1[[s-1]]))] <- 1
  G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[s]])) | ((G_long == 6) & (a_long < h1[[s-1]]))] <- 2
  G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[s]])) | ((G_long == 5) & (a_long < l0[[s-1]])) | ((G_long == 6) & (a_long >= h1[[s-1]]) & (a_long < l1[[s-1]]))] <- 3
  mu <- sapply(1:sum(N), 
             function(j){
               W[j,]%*%beta[[s]][G_a_long[j],]
             })
  var <- sigma2[[s]][G_a_long]
  numer <- dnorm(x = Y_long, mean = mu, sd = sqrt(var), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
    dnorm(x = logit_h0, mean = (v_long%*%delta_h0[[s-1]]),
          sd = sqrt(tau2_h0[[s-1]]),
          log = TRUE)       
  
  accept <- rep(1,
              times = sum(N))
  ratio <- exp(numer - denom)
  uni_draw <- runif(n = sum(N), min = 0.00, max = 1.00)
  
  logit_h0[ratio < uni_draw] <- logit_h0_old[ratio < uni_draw]
  h0[[s]][(ratio < uni_draw)] <- h0[[(s-1)]][(ratio < uni_draw)]
  G_a_long[ratio < uni_draw] <- G_a_long_old[ratio < uni_draw]
  accept[ratio < uni_draw] <- 0
  
  acctot_h0 <- acctot_h0 + accept
  
  ###
  #l0
  ###
  logit_l0_old <- logit_l0
  G_a_long_old <- G_a_long
  mu_old <- sapply(1:sum(N), 
                 function(j){
                   W[j,]%*%beta[[s]][G_a_long_old[j],]
                 })
  
  var_old <- sigma2[[s]][G_a_long_old]
  denom <- dnorm(x = Y_long, mean = mu_old, sd = sqrt(var_old), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
    dnorm(x = logit_l0_old,
          mean = (v_long%*%delta_l0[[s-1]]),
          sd = sqrt(tau2_l0[[s-1]]),
          log = TRUE)
  
  logit_l0 <- rnorm(n = sum(N), mean = logit_l0_old, sd = metrop_sd_l0)
  l0[[s]] <- 1.00/(1.00 + exp(-logit_l0))
  G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[s]])) | ((G_long == 6) & (a_long >= l1[[s-1]]))] <- 1
  G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[s]])) | ((G_long == 6) & (a_long < h1[[s-1]]))] <- 2
  G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[s]])) | ((G_long == 5) & (a_long < l0[[s]])) | ((G_long == 6) & (a_long >= h1[[s-1]]) & (a_long < l1[[s-1]]))] <- 3
  mu<-sapply(1:sum(N), 
             function(j){
               W[j,]%*%beta[[s]][G_a_long[j],]
             })
  var<-sigma2[[s]][G_a_long]
  numer<-dnorm(x = Y_long, mean = mu, sd = sqrt(var), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
    dnorm(x = logit_l0,
          mean = (v_long%*%delta_l0[[s-1]]),
          sd = sqrt(tau2_l0[[s-1]]),
          log = TRUE)       
  
  accept <- rep(1, times = sum(N))
  ratio <- exp(numer - denom)
  uni_draw <- runif(n = sum(N), min = 0.00, max = 1.00)
  
  logit_l0[ratio < uni_draw] <- logit_l0_old[ratio < uni_draw]
  l0[[s]][(ratio < uni_draw)] <- l0[[s-1]][(ratio < uni_draw)]
  G_a_long[ratio < uni_draw] <- G_a_long_old[ratio < uni_draw]
  accept[ratio < uni_draw] <- 0
  
  acctot_l0<-acctot_l0 + accept
  
  ###
  #h1
  ###
  logit_h1_old <- logit_h1
  G_a_long_old <- G_a_long
  mu_old <- sapply(1:sum(N), 
                 function(j){
                   W[j,]%*%beta[[s]][G_a_long_old[j],]
                 })
  
  var_old <- sigma2[[s]][G_a_long_old]
  denom <- dnorm(x = Y_long, mean = mu_old, sd = sqrt(var_old), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
    dnorm(x = logit_h1_old,
          mean = (v_long%*%delta_h1[[s-1]]),
          sd = sqrt(tau2_h1[[s-1]]),
          log = TRUE)
  
  logit_h1 <- rnorm(n = sum(N),mean = logit_h1_old, sd = metrop_sd_h1)
  h1[[s]] <- 1.00/(1.00 + exp(-logit_h1))
  G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[s]])) | ((G_long == 6) & (a_long >= l1[[s-1]]))] <- 1
  G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[s]])) | ((G_long == 6) & (a_long < h1[[s]]))] <- 2
  G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[s]])) | ((G_long == 5) & (a_long < l0[[s]])) | ((G_long == 6) & (a_long >= h1[[s]]) & (a_long < l1[[s-1]]))] <- 3
  mu <- sapply(1:sum(N), 
             function(j){
               W[j,]%*%beta[[s]][G_a_long[j],]
             })
  var <- sigma2[[s]][G_a_long]
  numer <- dnorm(x = Y_long, mean = mu, sd = sqrt(var), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
    dnorm(x = logit_h1,
          mean = (v_long%*%delta_h1[[s-1]]),
          sd = sqrt(tau2_h1[[s-1]]),
          log = TRUE)       
  
  accept <- rep(1, times = sum(N))
  ratio <- exp(numer - denom)
  uni_draw <- runif(n = sum(N), min = 0.00, max = 1.00)
  
  logit_h1[ratio < uni_draw] <- logit_h1_old[ratio < uni_draw]
  h1[[s]][(ratio < uni_draw)] <- h1[[s-1]][(ratio < uni_draw)]
  G_a_long[ratio < uni_draw] <- G_a_long_old[ratio < uni_draw]
  accept[ratio < uni_draw] <- 0
  
  acctot_h1 <- acctot_h1 + accept
  
  ###
  #l1
  ###
  logit_l1_old <- logit_l1
  G_a_long_old <- G_a_long
  mu_old <- sapply(1:sum(N), 
                 function(j){
                   W[j,]%*%beta[[s]][G_a_long_old[j],]
                 })
  
  var_old <- sigma2[[s]][G_a_long_old]
  denom <- dnorm(x = Y_long, mean = mu_old, sd = sqrt(var_old), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long_old == 1) | (D_long == 0 & G_a_long_old == 2) | (D_long == Z_long & G_a_long_old == 3))) +
    dnorm(x = logit_l1_old,
          mean = (v_long%*%delta_l1[[s-1]]),
          sd = sqrt(tau2_l1[[s-1]]),
          log = TRUE)
  
  logit_l1 <- rnorm(n = sum(N), mean = logit_l1_old, sd = metrop_sd_l1)
  l1[[s]] <- (h1[[s]] + exp(logit_l1))/(1.00 + exp(logit_l1))
  G_a_long[(G_long == 1) | ((G_long == 5) & (a_long >= l0[[s]])) | ((G_long == 6) & (a_long >= l1[[s]]))] <- 1
  G_a_long[(G_long == 2) | ((G_long == 4) & (a_long < h0[[s]])) | ((G_long == 6) & (a_long < h1[[s]]))] <- 2
  G_a_long[(G_long == 3) | ((G_long == 4) & (a_long >= h0[[s]])) | ((G_long == 5) & (a_long < l0[[s]])) | ((G_long == 6) & (a_long >= h1[[s]]) & (a_long < l1[[s]]))] <- 3
  mu <- sapply(1:sum(N), 
             function(j){
               W[j,]%*%beta[[s]][G_a_long[j],]
             })
  var <- sigma2[[s]][G_a_long]
  numer <- dnorm(x = Y_long, mean = mu, sd = sqrt(var), log = TRUE) +
    log(as.numeric((D_long == 1 & G_a_long == 1) | (D_long == 0 & G_a_long == 2) | (D_long == Z_long & G_a_long == 3))) +
    dnorm(x = logit_l1,
          mean = (v_long%*%delta_l1[[s-1]]),
          sd = sqrt(tau2_l1[[s-1]]),
          log = TRUE)       
  
  accept <- rep(1, times = sum(N))
  ratio <- exp(numer - denom)
  uni_draw <- runif(n = sum(N), min = 0.00, max = 1.00)
  
  logit_l1[ratio < uni_draw] <- logit_l1_old[ratio < uni_draw]
  l1[[s]][(ratio < uni_draw)] <- l1[[s-1]][(ratio < uni_draw)]
  G_a_long[ratio < uni_draw] <- G_a_long_old[ratio < uni_draw]
  accept[ratio < uni_draw] <- 0
  
  acctot_l1 <- acctot_l1 + accept
  
  #########
  #delta_h0
  #########
  cov_delta_h0 <- chol2inv(chol(crossprod(v_long)/tau2_h0[[s-1]] + diag(2)/sigma2_delta))
  mu_delta_h0 <- cov_delta_h0%*%(crossprod(v_long, logit_h0))/tau2_h0[[s-1]]
  delta_h0[[s]]<-rmnorm(n = 1, mean = mu_delta_h0, varcov = cov_delta_h0)
  
  ########
  #tau2_h0
  ########
  rate <- crossprod(logit_h0 - v_long%*%delta_h0[[s]])/2.00 + b_tau2
  tau2_h0[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = rate)
  
  #########
  #delta_l0
  #########
  cov_delta_l0 <- chol2inv(chol(crossprod(v_long)/tau2_l0[[s-1]] + diag(2)/sigma2_delta))
  mu_delta_l0 <- cov_delta_l0%*%(crossprod(v_long, logit_l0))/tau2_l0[[s-1]]
  delta_l0[[s]] <- rmnorm(n = 1, mean = mu_delta_l0, varcov = cov_delta_l0)

  ########
  #tau2_l0
  ########
  rate <- crossprod(logit_l0 - v_long%*%delta_l0[[s]])/2.00 + b_tau2
  tau2_l0[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = rate)
  
  #########
  #delta_h1
  #########
  cov_delta_h1 <- chol2inv(chol(crossprod(v_long)/tau2_h1[[s-1]] + diag(2)/sigma2_delta))
  mu_delta_h1 <- cov_delta_h0%*%(crossprod(v_long, logit_h1))/tau2_h1[[s-1]]
  delta_h1[[s]] <- rmnorm(n = 1, mean = mu_delta_h1, varcov = cov_delta_h1)
  
  ########
  #tau2_h1
  ########
  rate <- crossprod(logit_h1 - v_long%*%delta_h1[[s]])/2.00 + b_tau2
  tau2_h1[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = rate)

  #########
  #delta_l1
  #########
  cov_delta_l1 <- chol2inv(chol(crossprod(v_long)/tau2_l1[[s-1]] + diag(2)/sigma2_delta))
  mu_delta_l1 <- cov_delta_l1%*%(crossprod(v_long, logit_l1))/tau2_l1[[s-1]]
  delta_l1[[s]] <- rmnorm(n = 1, mean = mu_delta_l1, varcov = cov_delta_l1)
  
  ########
  #tau2_l1
  ########
  rate <- crossprod(logit_l1 - v_long%*%delta_l1[[s]])/2.00 + b_tau2
  tau2_l1[[s]] <- 1.00/rgamma(n = 1, shape = shape_tau2_update, rate = rate)
  
  print("#####")
  print(round(100*(s/mcmc_samples), 2))
  
  ##########
  #Estimands
  ##########
  
  if (s %in% iters) { 
    Y1 <- rep(NA, sum(N))
    Y0 <- rep(NA, sum(N))
    Y0p <- rep(NA, sum(N))
    C <- rep(NA, sum(N))
    
    eff.a <- 0.8
    eff.s <- 0.4
    eff.sp <- 0.8
    ij <- 0
    for (j in 1:J) {
      for (i in 1:N[J]) {
        ij <- ij + 1
        ## Need to figure out G(eff.a)
        if (G[[s]][j,i] %in% 1:3) {
          G.eff.a <- G[[s]][j,i]
        }  else if (G[[s]][j,i]==4) {
          if (eff.a < h0[[s]][ij]) {
            G.eff.a <- 2
          } else {
            G.eff.a <- 3
          }
        } else if (G[[s]][j,i]==5) {
          if (eff.a < l0[[s]][ij]) {
            G.eff.a <- 3
          } else {
            G.eff.a <- 1
          }
        } else if (G[[s]][j,i]==6) {
          if (eff.a < h1[[s]][ij]) {
            G.eff.a <- 2
          } else if (eff.a < l1[[s]][ij]) {
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
        if (a_long[j] == eff.a & S_long[j] == eff.s & G_a_long_true[j] == 3) {
          if (Z_long[j] == 0) {
            Y0 <- Y_long[j]
            W0p <- W1 <- W[ij,]
            W1[3] <- eff.s
            W0p[3] <- eff.sp
            W0p[4] <- W1[4] <- eff.a
            W0p[5] <- 0
            W1[5] <- 1

            mu1<-W1%*%beta[[s]][G.eff.a,]
            var1<-sigma2[[s]][G.eff.a]
            Y1[ij]<-rnorm(n = 1,
                          mean = mu1,
                          sd = sqrt(var1))
            
            mu0p<-W0p%*%beta[[s]][G.eff.a,]
            var0p<-sigma2[[s]][G.eff.a]
            Y0p[ij]<-rnorm(n = 1,
                           mean = mu0p,
                           sd = sqrt(var0p))  
          } else if (Z_long[j]==1) {
          Y1 <- Y_long[j]
          W0p <- W0 <- W[ij,]
          W0[3] <- eff.s
          W0p[3] <- eff.sp
          W0p[4] <- W0[4] <- eff.a
          W0p[5] <- W0[5] <- 0
          
          mu0<-W0%*%beta[[s]][G.eff.a,]
          var0<-sigma2[[s]][G.eff.a]
          Y0[ij]<-rnorm(n = 1,
                        mean = mu0,
                        sd = sqrt(var0)) 
          
          mu0p<-W0p%*%beta[[s]][G.eff.a,]
          var0p<-sigma2[[s]][G.eff.a]
          Y0p[ij]<-rnorm(n = 1,
                         mean = mu0p,
                         sd = sqrt(var0p)) 
          }
        } else {
        W0p <- W0 <- W1 <- W[ij,]
        W0[3] <- W1[3] <- eff.s
        W0p[3] <- eff.sp
        W0p[4] <- W0[4] <- W1[4] <- eff.a
        W0p[5] <- W0[5] <- 0
        W1[5] <- 1
        
        mu0<-W0%*%beta[[s]][G.eff.a,]
        var0<-sigma2[[s]][G.eff.a]
        Y0[ij]<-rnorm(n = 1,
                      mean = mu0,
                      sd = sqrt(var0)) 
        
        mu1<-W1%*%beta[[s]][G.eff.a,]
        var1<-sigma2[[s]][G.eff.a]
        Y1[ij]<-rnorm(n = 1,
                      mean = mu1,
                      sd = sqrt(var1))
        
        mu0p<-W0p%*%beta[[s]][G.eff.a,]
        var0p<-sigma2[[s]][G.eff.a]
        Y0p[ij]<-rnorm(n = 1,
                       mean = mu0p,
                       sd = sqrt(var0p)) 
        }
      }
    }
    CADE[s] <- sum(C*(Y1-Y0))/sum(C)
    CASE[s] <- sum(C*(Y0-Y0p))/sum(C)
  }
}
print("saving")
print(paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/CADE",
             id, ".rds"))
saveRDS(CADE, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/CADE",
                     id, ".rds"))
saveRDS(CASE, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/CASE",
                     id, ".rds"))
saveRDS(beta, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/beta",
                     id, ".rds"))
saveRDS(sigma2, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/sigma2",
                       id, ".rds"))
saveRDS(G, paste0("/home/cim24/project/OhnishiExtension/Results/", v, "/G",
                  id, ".rds"))
saveRDS(alpha, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/alpha",
                      id, ".rds"))
saveRDS(h0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/h0",
                   id, ".rds"))
saveRDS(l0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/l0",
                   id, ".rds"))
saveRDS(h1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/h1",
                   id, ".rds"))
saveRDS(l1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/l1",
                   id, ".rds"))
saveRDS(delta_h0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/deltah0",
                         id, ".rds"))
saveRDS(delta_l0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/deltal0",
                         id, ".rds"))
saveRDS(delta_h1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/deltah1",
                         id, ".rds"))
saveRDS(delta_l1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/deltal1",
                         id, ".rds"))
saveRDS(tau2_h0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/tau2h0",
                        id, ".rds"))
saveRDS(tau2_l0, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/tau2l0",
                        id, ".rds"))
saveRDS(tau2_h1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/tau2h1",
                        id, ".rds"))
saveRDS(tau2_l1, paste0("/home/cim24/project/OhnishiExtension/Results/", v,"/tau2l1",
                        id, ".rds"))
