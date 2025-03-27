##############
#Seed
##############
set.seed(id)

################
#Global Settings
################
J<-20
N<-rep(25,
       times = J)

M<-3
pz1=0.0
pz2=0.2
pz3=0.4
pz4=0.6
pz5=0.8
pz6=1.0
pz=c(pz1,pz2,pz3,pz4,pz5,pz6)
a<-sample(pz,
          size = J,
          replace = TRUE)
Z<-list(0)
T<-rep(NA,
       times = J)
for(j in 1:J){
   trt <- sample(1:N[j], size=a[j]*N[j], replace=FALSE)
   Z[[j]] <- rep(0, N[j])
   Z[[j]][trt] <- 1
   T[j]<-sum(Z[[j]])
   a[j]<-T[j]/N[j]

   }
Z_long<-Z[[1]]
a_long<-rep(a[1], 
            times = N[1])
N_long<-rep(N[1],
            times = N[1])
for(j in 2:J){

   Z_long<-c(Z_long,
             Z[[j]])
   a_long<-c(a_long,
             rep(a[j],
                 times = N[j]))
   N_long<-c(N_long,
             rep(N[j],
                 times = N[j]))

}

# Function
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

###########
#Predictors
###########
x_long<-rnorm(n = sum(N))
q_long<-cbind(1, 
              x_long)
v_long<-q_long

################
#True Parameters
################
gamma_true<-matrix(NA,
                   nrow = 6,
                   ncol = 2)
for(k in 1:5){
  gamma_true[k,]<-rnorm(n = 2)
}
gamma_true[6,]<-0.00

pi_mat_temp<-matrix(NA,
                    nrow = sum(N),
                    ncol = 6)
for(k in 1:6){
  pi_mat_temp[,k]<-exp(q_long%*%gamma_true[k,])
}
pi_mat<-matrix(NA,
               nrow = sum(N),
               ncol = 6)
for(k in 1:6){
  pi_mat[,k]<-pi_mat_temp[,k]/rowSums(pi_mat_temp)
}

R_long_true<-rep(NA,
                 times = sum(N))
for(j in 1:sum(N)){
  R_long_true[j]<-sample(c(1:6),
                         size = 1,
                         prob = pi_mat[j,])
}
R_true<-list(0)
R_true[[1]]<-R_long_true[1:N[1]]
for(j in 2:J){
  R_true[[j]]<-R_long_true[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
}

delta_h4_true <- rnorm(n = 2)
delta_l5_true <- rnorm(n = 2)
delta_h6_true <- rnorm(n = 2)
delta_l6_true <- rnorm(n = 2)
tau2_h4_true<-0.01
tau2_l5_true<-0.01
tau2_h6_true<-0.01
tau2_l6_true<-0.01

logit_h4_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_h4_true),
                     sd = sqrt(tau2_h4_true))
h4_true<-1.00/(1.00 + exp(-logit_h4_true))

logit_l5_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_l5_true),
                     sd = sqrt(tau2_l5_true))
l5_true<-1.00/(1.00 + exp(-logit_l5_true))

logit_h6_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_h6_true),
                     sd = sqrt(tau2_h6_true))
h6_true<-1.00/(1.00 + exp(-logit_h6_true))

logit_l6_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_l6_true),
                     sd = sqrt(tau2_l6_true))
l6_true<-(h6_true + exp(logit_l6_true))/(1.00 + exp(logit_l6_true))


R_a_long_true<-rep(NA,
                   times = sum(N))
R_a_long_true[(R_long_true == 1) | ((R_long_true == 5) & (a_long >= l5_true)) | ((R_long_true == 6) & (a_long >= l6_true))]<-1
R_a_long_true[(R_long_true == 2) | ((R_long_true == 4) & (a_long < h4_true)) | ((R_long_true == 6) & (a_long < h6_true))]<-2
R_a_long_true[(R_long_true == 3) | ((R_long_true == 4) & (a_long >= h4_true)) | ((R_long_true == 5) & (a_long < l5_true)) | (R_long_true == 6 & (a_long >= h6_true) & (a_long < l6_true))]<-3
R_a_true<-list(0)
R_a_true[[1]]<-R_a_long_true[1:N[1]]
for(j in 2:J){
  R_a_true[[j]]<-R_a_long_true[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
}

D<-list(0)
S<-rep(NA,
       times = J)
for(j in 1:J){
  
  D[[j]]<-rep(NA,
              times = N[j])
  D[[j]][R_a_true[[j]] == 1]<-1
  D[[j]][R_a_true[[j]] == 2]<-0
  D[[j]][R_a_true[[j]] == 3]<-Z[[j]][R_a_true[[j]] == 3]
  S[j]<-sum(D[[j]])
  
}
D_long<-D[[1]]
S_long<-rep(S[1],
            times = N[1])
for(j in 2:J){
  
  D_long<-c(D_long,
            D[[j]])
  S_long<-c(S_long,
            rep(S[j],
                times = N[j]))
  
}

W_true<-list(cbind(x_long, S_long/N_long, a_long, Z_long),
        cbind(x_long, S_long/N_long, a_long, Z_long),
        cbind(x_long, S_long/N_long, a_long, Z_long),
        cbind(x_long, S_long/N_long, a_long, Z_long, h4_true),
        cbind(x_long, S_long/N_long, a_long, Z_long, l5_true),
        cbind(x_long, S_long/N_long, a_long, Z_long, h6_true, l6_true))


psi2_true <- runif(6)
phi_true <- list(runif(4), runif(4), runif(4), runif(5), runif(5), runif(6))
Sigma_true <- list(calcSigma(sum(N), phi_true[[1]], W_true[[1]]),
                   calcSigma(sum(N), phi_true[[2]], W_true[[2]]),
                   calcSigma(sum(N), phi_true[[3]], W_true[[3]]),
                   calcSigma(sum(N), phi_true[[4]], W_true[[4]]),
                   calcSigma(sum(N), phi_true[[5]], W_true[[5]]),
                   calcSigma(sum(N), phi_true[[6]], W_true[[6]]))

theta_true <- cbind(rmnorm(1, 0, varcov=Sigma_true[[1]]*psi2_true[1]),
                   rmnorm(1, 0, varcov=Sigma_true[[2]]*psi2_true[2]),
                   rmnorm(1, 0, varcov=Sigma_true[[3]]*psi2_true[3]),
                   rmnorm(1, 0, varcov=Sigma_true[[4]]*psi2_true[4]),
                   rmnorm(1, 0, varcov=Sigma_true[[5]]*psi2_true[5]),
                   rmnorm(1, 0, varcov=Sigma_true[[6]]*psi2_true[6]))

mu_true <- rnorm(6, 0, sd=4)

sigma2_true <- runif(6)

Y_long <- rmnorm(1,
                 mu_true[R_long_true] + theta_true[cbind(seq_along(R_long_true), R_long_true)],
                 varcov=diag(sigma2_true[R_long_true]))

Y<-list(0)
Y[[1]]<-Y_long[1:N[1]]
for(j in 2:J){
   Y[[j]]<-Y_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
   }

Y0_long <- rep(NA, sum(N))
Y0p_long <- rep(NA, sum(N))
Y1_long <- rep(NA, sum(N))
C <- rep(NA, sum(N))
CADE.G <- rep(NA, sum(N))
eff.a <- 0.8
eff.s <- 0.4
eff.sp <- 0.8

for(j in 1:sum(N)) {
  if (a_long[j] == eff.a & S_long[j] == eff.s & R_a_long_true[j] == 3) {
    CADE.G[j] <- 3
    if (Z_long[j] == 0) {
      Y0_long[j] <- Y_long[j]
      W1 <- W_true[[R_long_true[j]]]
      W1[j,2] <- eff.s
      W1[j,3] <- eff.a
      W1[j,4] <- 1
      Sig1 <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W1)
      theta1 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig1)
      Y1_long[j] <- rnorm(n = 1,
                      mean = mu_true[R_long_true[j]] + theta1[j],
                      sd = sigma2_true[R_long_true[j]]) 

      W0p <- W_true[[R_long_true[j]]]
      W0p[j,2] <- eff.sp
      W0p[j,3] <- eff.a
      W0p[j,4] <- 0
      Sig0p <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0p) 
      theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
      Y0p_long[j] <- rnorm(n = 1,
                       mean = mu_true[R_long_true[j]] + theta0p[j],
                       sd = sigma2_true[R_long_true[j]])     
    } else {
      Y1_long[j] <- Y_long[j]
      W0 <- W_true[[R_long_true[j]]]
      W0[j,2] <- eff.s
      W0[j,3] <- eff.a
      W0[j,4] <- 0
      Sig0 <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0)
      theta0 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0)
      Y0_long[j] <- rnorm(n = 1,
                      mean = mu_true[R_long_true[j]] + theta0[j],
                      sd = sigma2_true[R_long_true[j]]) 
      W0p <- W_true[[R_long_true[j]]]
      W0p[j,2] <- eff.s
      W0p[j,3] <- eff.a
      W0p[j,4] <- 0
      Sig0p <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0p) 
      theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
      Y0p_long[j] <- rnorm(n = 1,
                       mean = mu_true[R_long_true[j]] + theta0p[j],
                       sd = sigma2_true[R_long_true[j]])   
    }
  } else {
    if (R_a_long_true[j] %in% 1:3) {
      CADE.G[j] <- R_a_long_true[j]
    }  else if (R_a_long_true[j] == 4) {
      if (eff.a < h4_true[j]) {
        CADE.G[j] <- 2
      } else {
        CADE.G[j] <- 3
      }
    } else if (R_a_long_true[j] == 5) {
      if (eff.a < l5_true[j]) {
        CADE.G[j] <- 3
      } else {
        CADE.G[j] <- 1
      }
    } else if (R_a_long_true[j] == 6) {
      if (eff.a < h6_true[j]) {
        CADE.G[j] <- 2
      } else if (eff.a < l6_true[j]) {
        CADE.G[j] <- 3
      } else {
        CADE.G[j] <- 1
      }
    }
    
    if (CADE.G[j] == 3) {
      C[j] <- 1
    } else {
      C[j] <- 0
    }
  }
  W0 <- W0p <- W1 <- W_true[[R_long_true[j]]]
  W0[j,2] <- W1[j,2] <- eff.s
  W0p[j,2] <- eff.sp
  W0[j,3] <- W0p[j,3] <- W1[j,3] <- eff.a
  W0[j,4] <- W0p[j,4] <- 0
  W1[j,4] <- 1
  Sig0 <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W0, Sigma_true[[R_long_true[j]]], j)
  theta0 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0)
  Y0_long[j] <- rnorm(n = 1,
                      mean = mu_true[R_long_true[j]] + theta0[j],
                      sd = sigma2_true[R_long_true[j]]) 
  Sig1 <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W1, Sigma_true[[R_long_true[j]]], j)
  theta1 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig1)
  Y1_long[j] <- rnorm(n = 1,
                      mean = mu_true[R_long_true[j]] + theta1[j],
                      sd = sigma2_true[R_long_true[j]]) 
  Sig0p <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W0p, Sigma_true[[R_long_true[j]]], j) 
  theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
  Y0p_long[j] <- rnorm(n = 1,
                       mean = mu_true[R_long_true[j]] + theta0p[j],
                       sd = sigma2_true[R_long_true[j]]) 
}

CADE.true <- sum((CADE.G==3)*(Y1_long - Y0_long))/(sum(CADE.G==3))
CASE.true <- sum((CADE.G==3)*(Y0_long - Y0p_long))/(sum(CADE.G==3))


# Y0_long <- matrix(NA, nrow=sum(N), ncol=length(pz))
# Y1_long <- matrix(NA, nrow=sum(N), ncol=length(pz))
# C <- matrix(NA, nrow=sum(N), ncol=length(pz))
# CADE.G <- matrix(NA, nrow=sum(N), ncol=length(pz))
# 
# for(j in 1:sum(N)) {
#   for(t in pz) {
#     if (a_long[j] == eff.a & S_long[j] == eff.s & R_a_long_true[j] == 3) {
#       CADE.G[j,t] <- 3
#       if (Z_long[j] == 0) {
#         Y0_long[j,t] <- Y_long[j]
#         W1 <- W_true[[R_long_true[j]]]
#         W1[j,2] <- eff.s
#         W1[j,3] <- t
#         W1[j,4] <- 1
#         Sig1 <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W1)
#         theta1 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig1)
#         Y1_long[j,t] <- rnorm(n = 1,
#                             mean = mu_true[R_long_true[j]] + theta1[j],
#                             sd = sigma2_true[R_long_true[j]]) 
#         
#         W0p <- W_true[[R_long_true[j]]]
#         W0p[j,2] <- eff.sp
#         W0p[j,3] <- t
#         W0p[j,4] <- 0
#         Sig0p <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0p) 
#         theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
#         Y0p_long[j,t] <- rnorm(n = 1,
#                              mean = mu_true[R_long_true[j]] + theta0p[j],
#                              sd = sigma2_true[R_long_true[j]])     
#       } else {
#         Y1_long[j,t] <- Y_long[j]
#         W0 <- W_true[[R_long_true[j]]]
#         W0[j,2] <- eff.s
#         W0[j,3] <- t
#         W0[j,4] <- 0
#         Sig0 <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0)
#         theta0 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0)
#         Y0_long[j,t] <- rnorm(n = 1,
#                             mean = mu_true[R_long_true[j]] + theta0[j],
#                             sd = sigma2_true[R_long_true[j]]) 
#         W0p <- W_true[[R_long_true[j]]]
#         W0p[j,2] <- eff.s
#         W0p[j,3] <- t
#         W0p[j,4] <- 0
#         Sig0p <- calcSigma(sum(N), phi_true[[R_long_true[j]]], W0p) 
#         theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
#         Y0p_long[j,t] <- rnorm(n = 1,
#                              mean = mu_true[R_long_true[j]] + theta0p[j],
#                              sd = sigma2_true[R_long_true[j]])   
#       }
#     } else {
#       if (R_a_long_true[j] %in% 1:3) {
#         CADE.G[j,t] <- R_a_long_true[j]
#       }  else if (R_a_long_true[j] == 4) {
#         if (eff.a < h4_true[j]) {
#           CADE.G[j,t] <- 2
#         } else {
#           CADE.G[j,t] <- 3
#         }
#       } else if (R_a_long_true[j] == 5) {
#         if (eff.a < l5_true[j]) {
#           CADE.G[j,t] <- 3
#         } else {
#           CADE.G[j,t] <- 1
#         }
#       } else if (R_a_long_true[j] == 6) {
#         if (eff.a < h6_true[j]) {
#           CADE.G[j,t] <- 2
#         } else if (eff.a < l6_true[j]) {
#           CADE.G[j,t] <- 3
#         } else {
#           CADE.G[j,t] <- 1
#         }
#       }
#       
#       if (CADE.G[j,t] == 3) {
#         C[j,t] <- 1
#       } else {
#         C[j,t] <- 0
#       }
#     }
#     W0 <- W0p <- W1 <- W_true[[R_long_true[j]]]
#     W0[j,2] <- W1[j,2] <- eff.s
#     W0p[j,2] <- eff.sp
#     W0[j,3] <- W0p[j,3] <- W1[j,3] <- t
#     W0[j,4] <- W0p[j,4] <- 0
#     W1[j,4] <- 1
#     Sig0 <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W0, Sigma_true[[R_long_true[j]]], j)
#     theta0 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0)
#     Y0_long[j,t] <- rnorm(n = 1,
#                         mean = mu_true[R_long_true[j]] + theta0[j],
#                         sd = sigma2_true[R_long_true[j]]) 
#     Sig1 <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W1, Sigma_true[[R_long_true[j]]], j)
#     theta1 <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig1)
#     Y1_long[j,t] <- rnorm(n = 1,
#                         mean = mu_true[R_long_true[j]] + theta1[j],
#                         sd = sigma2_true[R_long_true[j]]) 
#     Sig0p <- updateSigma(sum(N), phi_true[[R_long_true[j]]], W0p, Sigma_true[[R_long_true[j]]], j) 
#     theta0p <- rmnorm(1, mean = 0, varcov = psi2_true[R_long_true[j]]*Sig0p)
#     Y0p_long[j,t] <- rnorm(n = 1,
#                          mean = mu_true[R_long_true[j]] + theta0p[j],
#                          sd = sigma2_true[R_long_true[j]]) 
#   }
# }












                    