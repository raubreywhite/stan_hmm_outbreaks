# clears workspace: 
rm(list=ls()) 
setwd("/git/example-models/misc/hmm")
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

a0 <- Sys.time()
#fit0 <- stan('continuous.stan',   
#             data=data,iter=1000, chains=1, init=0)
b0 <- Sys.time()

a1 <- Sys.time()
fit1 <- stan('continuous_sufficient.stan',   
             data=data,iter=1000, chains=1, init=0)
b1 <- Sys.time()


x <- rstanarm::stan_glm(y ~ outbreak,family=poisson, offset=log(consults))

# Now the values for the monitored parameters are in the "samples" object, 
# ready for inspection.

print(samples, digits=3)
