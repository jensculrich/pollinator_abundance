data {
  int<lower=0> M; // Size of augumented data set
  int<lower=0> nind; // Number of observed individuals
  int<lower=0> n_visits; // Number of sampling occasions
  int<lower=0> n_sites;
  int sites[M];
  int<lower=0> y[M, n_visits]; // capture histories matrix
  vector[n_sites] X; // site restored
}

transformed data {
  
  // s == totals in each row // i.e. was a particular individual ever detected?
  int<lower=0> s[M]; 
  
  for (i in 1:M) {
    s[i] = sum(y[i]);
  }
  
}

parameters {
  
  real mu_alpha0;
  real mu_alpha1;
  real mu_beta0;
  real mu_beta1;
  
}

transformed parameters {
  
  vector[M] logit_p;
  vector[n_sites] lambda;
  vector[n_sites] probs;
  
  real<lower=0, upper=1> omega; // Inclusion probability
  
  for(i in 1:n_sites){
    lambda[i] = exp(mu_alpha0 + mu_alpha1 * X[i]);
    probs[i] = lambda[i] / sum(lambda[]);
  }
  
  for(i in 1:M) {
    logit_p[i] = mu_beta0 + mu_beta1 * X[sites[i]]; 
  }
  
  omega = sum(lambda[]) / M;
  
}

model {
  
  // Priors
  // omega ~ uniform(0, 1); // implicitly defined by lower/upper bounds
  mu_alpha0 ~ normal(0, 2);
  mu_alpha1 ~ uniform(0, 2);
  mu_beta0 ~ normal(0, 2);
  mu_beta1 ~ normal(0, 2);
  
  for (i in 1:M) {
    
    sites[i] ~ categorical(probs[]); 
    
    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega)
                + bernoulli_logit_lpmf(y[i] | logit_p[i]);
    } else // s[i] == 0
    {
      target += log_sum_exp(bernoulli_lpmf(1 | omega)
                            // z[i] == 1
                            + bernoulli_logit_lpmf(y[i] | logit_p[i]),
                            bernoulli_lpmf(0 | omega) // z[i] == 0
                            );
    }
  }
}


