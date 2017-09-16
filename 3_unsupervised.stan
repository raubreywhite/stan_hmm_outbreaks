data {
  int<lower=1> K; // number of groups
  int<lower=1> T; // length of timeseries
  int y[T]; // observations
  vector[T] consults;
}
transformed data {
  vector[T] log_expo;
  log_expo = log(consults);
}
parameters {
  real<lower=0,upper=1> p0 ;     //initial prob grp 1
  real<lower=0,upper=1> TP[K] ;  //transition probs of staying in group
  real mu[K]; // locations of mixture components
  real<lower=0> sigma[K]; // scales of mixture components
}
transformed parameters {
  real<lower=0,upper=1> prob_grp[T];  //smoother estimate probability of group membership
  real<lower=0,upper=1> pred[T];   //one-step filter prediction of probabililty of group membership
  {
    real F[T];   //filter forwards group membership prob
    real B1[T];  //backwards information filter from grp 1
    real B2[T];  //backwards information filter from grp 2
    real Z1[T];  //intermediate data
    real Z2[T];  //intermediate data
    real like1;    
    real like2;
    real p1;
    real p2;
    real k;
    int i;
    
    //Forwards algorithm
    
      F[1]=p0; 
      pred[1]=F[1];
    
    for (t in 1:T){
        //update prior using data
        like1=exp(poisson_lpmf(y[t] | log_expo[t]+mu[1]));
        like2=exp(poisson_lpmf(y[t] | log_expo[t]+mu[2]));
        p1=F[t]*like1;
        p2=(1-F[t])*like2;
        F[t]=p1/(p1+p2);
        
        //predict forward one timestep
        if (t != T) {
          p1=F[t]*TP[1]+(1-F[t])*(1-TP[2]);
          p2=F[t]*(1-TP[1])+(1-F[t])*TP[2];
          F[t+1]=p1/(p1+p2);
          pred[t+1]=F[t+1];
        }
    }  
    //backwards algorithm
    B1[T]=1; 
    B2[T]=1; 

    for (t in 1:(T-1)){
      i=t*(-1)+T;      // transform t to get a backwards loop
      like1=exp(poisson_lpmf(y[i+1] | log_expo[i+1]+mu[1]));
      like2=exp(poisson_lpmf(y[i+1] | log_expo[i+1]+mu[2]));  
      
      B1[i]=TP[1]*like1*B1[(i+1)]+(1-TP[2])*like2*B2[(i+1)];
      B2[i]=(1-TP[1])*like1*B1[(i+1)]+TP[2]*like2*B2[(i+1)];
      
      k=B1[i]+B2[i];
      B1[i]=B1[i]/k;
      B2[i]=B2[i]/k;
    }
    // put it all together
    for (t in 1:T){
        Z1[t]=F[t]*B1[t];
        Z2[t]=(1-F[t])*B2[t];
        prob_grp[t]=Z1[t]/(Z1[t]+Z2[t]);
    }
    
  }
}
model {
  real ps; // temp for log component densities
  sigma ~ cauchy(0,2.5);
  mu[1] ~ normal(0.01,0.001);
  mu[2] ~ normal(0.11,0.001);
  
  for (t in 1:T){
      ps = pred[t]*exp(poisson_lpmf(y[t] | log_expo[t]+mu[1]))+
        (1-pred[t])*exp(poisson_lpmf(y[t] | log_expo[t]+mu[2]));
      increment_log_prob(log(ps));
  }
}

