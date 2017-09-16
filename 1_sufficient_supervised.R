# clears workspace: 
rm(list=ls()) 
setwd("/git/stan_hmm_outbreaks/")
library(rstan)

K <- 2;
N <- 500;

lambdas <- c(0.01,0.1) # per-person prob of sickness
consults <- rpois(N,400)
z <- c(rep(1,N/4),rep(2,N/2),rep(1,N/4))
betas <- c(0.01,0.11)

outbreak <- z-1

lambda <- betas[1] + betas[2]*outbreak
y <- rpois(N,lambda=lambda*consults)

alpha <- rep(1,K);

data <- list(N=N, K=K, z=z, alpha=alpha, consults=consults, outbreak=outbreak,y=y)  # To be passed on to Stan

a1 <- Sys.time()
fit1 <- stan('1_sufficient_supervised.stan',   
             data=data,iter=1000, chains=1, init=0)
b1 <- Sys.time()

