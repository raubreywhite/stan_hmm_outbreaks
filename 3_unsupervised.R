# clears workspace: 
rm(list=ls()) 
setwd("/git/stan_hmm_outbreaks/")
library(rstan)
library(ggplot2)
library(data.table)

rt <- stanc(file = '3_unsupervised.stan',model_name = "foo")
sm <- stan_model(stanc_ret = rt) 

K <- 2;
N <- 800;

lambdas <- c(0.01,0.1) # per-person prob of sickness
consults <- rpois(N,400)
z <- c(rep(1,N/4),rep(2,N/2),rep(1,N/4))
outbreak <- z-1
betas <- c(0.01,0.11)

outbreak <- z-1

lambda <- betas[1] + betas[2]*outbreak
y <- rpois(N,lambda=lambda*consults)

alpha <- rep(1,K);

data <- list(N=N, K=K, z=z, alpha=alpha, consults=consults, outbreak=outbreak,y=y)  # To be passed on to Stan

retval <- vector("list",100)
for(i in 1:length(retval)){
  T <- 28+7*i
  stan.data<-list(K=2,T=T,y=round(y[1:T]),consults=consults[1:T])
  fit2<-sampling(sm,data=stan.data,iter=2000,chains=1)
  tmp2<-extract(fit2)
  index <- (T-6):T
  retval[[i]] <- data.frame(p=apply(tmp2$prob_grp,2,mean)[index],outbreak=outbreak[index],day=index)
}

plotData <- rbindlist(retval)
q <- ggplot(plotData,aes(x=day))
q <- q + geom_line(aes(y=p))
q <- q + geom_ribbon(data=plotData[plotData$outbreak==1,],aes(ymax=Inf,ymin=-Inf),fill="red",alpha=0.3)
q

plot(tmp2[[1]])
plot(tmp2[[2]])
plot(tmp2[[3]])
plot(tmp2[[2]][,1])
plot(tmp2[[2]][,2])
plot(tmp2[[3]][,1])
plot(tmp2[[3]][,2])
plot(tmp2[[4]][,1])
plot(tmp2[[4]][,2])
plot(tmp2[[5]][,1,1],main=paste("X:",X[1,1],"Y:",Y[1,1],sep=" "))
plot(tmp2[[5]][,1,10],main=paste("X:",X[1,10],"Y:",Y[1,10],sep=" "))

# Convergance and estimation is perfect.
# Don't be concerned about label switching.





