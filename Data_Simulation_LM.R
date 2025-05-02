##############
#Seed
##############
set.seed(id)

################
#Global Settings
################

N<-rep(Nj,
       times = J)

M<-J
a<-sample(seq(0, 1, 1/(nalpha-1)),
          size = J,
          replace = TRUE)
Z<-list(0)
T<-rep(NA,
       times = J)
for(j in 1:J){

   Z[[j]]<-rbinom(n = N[j],
                  size = 1,
                  prob = a[j])
   T[j]<-sum(Z[[j]])
   a[j]<-T[j]/N[j]

   }
Z_long<-Z[[1]]
a_long<-rep(a[1], 
            times = N[1])
for(j in 2:J){

   Z_long<-c(Z_long,
             Z[[j]])
   a_long<-c(a_long,
             rep(a[j],
                 times = N[j]))

   }

###########
#Predictors
###########
set.seed(22)
x_long<-rnorm(n = sum(N))
q_long<-cbind(1, x_long)
v_long<-q_long

################
#True Parameters
################
beta_true<-matrix(NA,
                  nrow = 3,
                  ncol = 5)
sigma2_true<-rep(NA,
                 times = 3)
for(k in 1:3){

   beta_true[k,]<-c(rnorm(n = 2), rnorm(n=1, mean=2, sd=0.5),
                    rnorm(n=1, mean=1, sd=0.5), rnorm(n=1, mean=2, sd=0.5))
   sigma2_true[k]<-runif(n = 1,
                         min = 0.00,
                         max = 0.10)

   }

alpha_true<-matrix(NA,
                   nrow = 6,
                   ncol = 2)
for(k in 1:5){
   alpha_true[k,]<-rnorm(n = 2)
   }
alpha_true[6,]<-0.00

pi_mat_temp<-matrix(NA,
                    nrow = sum(N),
                    ncol = 6)
for(k in 1:6){
   pi_mat_temp[,k]<-exp(q_long%*%alpha_true[k,])
   }
pi_mat<-matrix(NA,
               nrow = sum(N),
               ncol = 6)
for(k in 1:6){
   pi_mat[,k]<-pi_mat_temp[,k]/rowSums(pi_mat_temp)
   }

G_long_true<-rep(NA,
                 times = sum(N))
for(j in 1:sum(N)){
   G_long_true[j]<-sample(c(1:6),
                          size = 1,
                          prob = pi_mat[j,])
   }
G_true<-list(0)
G_true[[1]]<-G_long_true[1:N[1]]
for(j in 2:J){
   G_true[[j]]<-G_long_true[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
   }

delta_h0_true<-rnorm(n = 2)
tau2_h0_true<-0.01
logit_h0_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_h0_true),
                     sd = sqrt(tau2_h0_true))
h0_true<-1.00/(1.00 + exp(-logit_h0_true))

delta_l0_true<-rnorm(n = 2)
tau2_l0_true<-0.01
logit_l0_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_l0_true),
                     sd = sqrt(tau2_l0_true))
l0_true<-1.00/(1.00 + exp(-logit_l0_true))

delta_h1_true<-rnorm(n = 2)
tau2_h1_true<-0.01
logit_h1_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_h1_true),
                     sd = sqrt(tau2_h1_true))
h1_true<-1.00/(1.00 + exp(-logit_h1_true))

delta_l1_true<-rnorm(n = 2)
tau2_l1_true<-0.01
logit_l1_true<-rnorm(n = sum(N),
                     mean = (v_long%*%delta_l1_true),
                     sd = sqrt(tau2_l1_true))
l1_true<-(h1_true + exp(logit_l1_true))/(1.00 + exp(logit_l1_true))

G_a_long_true<-rep(NA,
                   times = sum(N))
G_a_long_true[(G_long_true == 1) | ((G_long_true == 5) & (a_long >= l0_true)) | ((G_long_true == 6) & (a_long >= l1_true))]<-1
G_a_long_true[(G_long_true == 2) | ((G_long_true == 4) & (a_long < h0_true)) | ((G_long_true == 6) & (a_long < h1_true))]<-2
G_a_long_true[(G_long_true == 3) | ((G_long_true == 4) & (a_long >= h0_true)) | ((G_long_true == 5) & (a_long < l0_true)) | (G_long_true == 6 & (a_long >= h1_true) & (a_long < l1_true))]<-3
G_a_true<-list(0)
G_a_true[[1]]<-G_a_long_true[1:N[1]]
for(j in 2:J){
   G_a_true[[j]]<-G_a_long_true[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
   }

#####
#Data
#####
D<-list(0)
S<-rep(NA,
       times = J)
for(j in 1:J){

   D[[j]]<-rep(NA,
               times = N[j])
   D[[j]][G_a_true[[j]] == 1]<-1
   D[[j]][G_a_true[[j]] == 2]<-0
   D[[j]][G_a_true[[j]] == 3]<-Z[[j]][G_a_true[[j]] == 3]
   S[j]<-sum(D[[j]])/N[j]

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

W<-cbind(1, x_long, S_long, a_long, Z_long)
   
Y0_long <- rep(NA, sum(N))
Y1_long <- rep(NA, sum(N))
Y0p_long <- rep(NA, sum(N))
C <- rep(NA, sum(N))
CADE.G <- rep(NA, sum(N))
eff.a <- 0.8
eff.s <- 0.4
eff.sp <- 0.8

for(j in 1:sum(N)) {
    if (G_long_true[j] %in% 1:3) {
      CADE.G[j] <- G_long_true[j]
    }  else if (G_long_true[j] == 4) {
      if (eff.a < h0_true[j]) {
        CADE.G[j] <- 2
      } else {
        CADE.G[j] <- 3
      }
    } else if (G_long_true[j] == 5) {
      if (eff.a < l0_true[j]) {
        CADE.G[j] <- 3
      } else {
        CADE.G[j] <- 1
      }
    } else if (G_long_true[j] == 6) {
      if (eff.a < h1_true[j]) {
        CADE.G[j] <- 2
      } else if (eff.a < l1_true[j]) {
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
    W0 <- W0p <- W1 <- W[j, ]
    W0[3] <- W1[3] <- eff.s
    W0p[3] <- eff.sp
    W0[4] <- W0p[4] <- W1[4] <- eff.a
    W0[5] <- W0p[5] <- 0
    W1[5] <- 1
    mu0 <- W0 %*% beta_true[CADE.G[j], ]
    var0 <- sigma2_true[CADE.G[j]]
    Y0_long[j] <- rnorm(n = 1,
                        mean = mu0,
                        sd = sqrt(var0))
    mu1 <- W1 %*% beta_true[CADE.G[j], ]
    var1 <- sigma2_true[CADE.G[j]]
    Y1_long[j] <- rnorm(n = 1,
                        mean = mu1,
                        sd = sqrt(var1))
    mu0p <- W0p %*% beta_true[CADE.G[j], ]
    var0p <- sigma2_true[CADE.G[j]]
    Y0p_long[j] <- rnorm(n = 1,
                         mean = mu0p,
                         sd = sqrt(var0p))
}

Y_long<-rep(NA,
            time = sum(N))

for(j in 1:sum(N)){
  if (a_long[j] == eff.a & S_long[j] == eff.s & G_a_long_true[j] == 3) {
    if (Z_long[j]==0) {
      Y_long[j] <- Y0_long[j]
    } else if (Z_long[j]==1) {
      Y_long[j] <- Y1_long[j]
    }
  } else if (a_long[j] == eff.a & S_long[j] == eff.sp & G_a_long_true[j] == 3 & Z_long[j]==0) {
    Y_long[j] <- Y0p_long[j]
  } else {
    mu<-W[j,]%*%beta_true[G_a_long_true[j],]
    var<-sigma2_true[G_a_long_true[j]]
    Y_long[j]<-rnorm(n = 1,
                     mean = mu,
                     sd = sqrt(var))
  }
}
Y<-list(0)
Y[[1]]<-Y_long[1:N[1]]
for(j in 2:J){
  Y[[j]]<-Y_long[(1 + sum(N[1:(j-1)])):sum(N[1:j])]
}

CADE.true <- sum((CADE.G==3)*(Y1_long - Y0_long))/(sum(CADE.G==3))
CASE.true <- sum((CADE.G==3)*(Y0_long - Y0p_long))/(sum(CADE.G==3))










                    