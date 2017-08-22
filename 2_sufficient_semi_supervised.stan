data {
  int<lower=0> N; 
  int<lower=1> K;  // num categories
  int<lower=1,upper=K> z[N]; // categories
  vector<lower=0>[K] alpha;  // transit prior
  vector[N] consults;
  vector[N] outbreak;
  int y[N];
}
transformed data {
  vector[N] log_expo;
  int<lower=0> trans[K,K];
  
  log_expo = log(consults);
  
  for (k1 in 1:K) 
    for (k2 in 1:K)
      trans[k1,k2] = 0;
    for (i in 2:N)
      trans[z[i - 1], z[i]] = 1 + trans[z[i - 1], z[i]];
    
}
parameters {
  simplex[K] theta[K];  // transit probs
  vector[K] phi;    // category effects
} 
model {
  for(k in 1:K)
    theta[k] ~ dirichlet(alpha);
  phi ~ normal(0,1);
  
  //z[1] ~ categorical(theta[1]);
  for (k in 1:K)
    trans[k] ~ multinomial(theta[k]);
  
  y ~ poisson_log(log_expo + phi[z]);
  //for(i in 1:N){
  //  y[i] ~ poisson_log(log_expo + phi[z[i]]);
  //}
}

