// Multinomial N-mixture model
// for estimating effects of covariate on abundance and detection rate
// given a set of capture/encounter histories (y) from N (value unknown) individuals
// across a spatially stratified system consisting of n_sites units.

// there are 3 repeat surveys per site -> 
// 8 possible detection histories (including 000 never detected)

data {
  
  int<lower=0> n_sites; // number of sites
  int<lower=0> y[n_sites, 7]; // capture histories matrix (7 observable histories possible)
  int<lower=0> nobs[n_sites]; // total number of individuals observed at each site 
                                // (with any history) // i.e. nobs = sum(y[,1:7])
  int<lower=0> K; // some arbitrarilly large maximum site-specific abundance to search up until
  vector[n_sites] X; // site covariate (categorical 0 or 1)

}

parameters {
  
  real mu_alpha0; // abundance intercept
  real mu_alpha1; // effect of covariate on abundance 
  real mu_beta0; // detection intercept
  real mu_beta1; // effect of covariate on detection
  
}

transformed parameters {
  
  vector<lower=0, upper=1>[n_sites] p; // site-specific detection, ranges between 0 and 1 
  vector[n_sites] lambda; // site-specific abundance
  
  real cellprobs[n_sites,8]; // cell probabilities (probability of landing in a cell)
  real cellprobs_conditional[n_sites,7]; // cell probs conditional on the sum of all other observable probabilities
  real p_det[n_sites]; // probability of detecting an individual >= once
  
  for(i in 1:n_sites){
    
    lambda[i] = mu_alpha0 + mu_alpha1 * X[i]; // linear model for abundance
    p[i] = inv_logit(mu_beta0 + mu_beta1 * X[i]); // linear model for detection
    
    // define cell probabilities (capture histories based on site-specific detection prob.)
    cellprobs[i, 1] = p[i]*p[i]*p[i]; // 111 (observed, observed, observed)
    cellprobs[i, 2] = p[i]*p[i]*(1-p[i]); // 110 (observed, observed, not observed)
    cellprobs[i, 3] = p[i]*(1-p[i])*p[i]; // 101
    cellprobs[i, 4] = p[i]*(1-p[i])*(1-p[i]); // 100
    cellprobs[i, 5] = (1-p[i])*p[i]*p[i]; // 011
    cellprobs[i, 6] = (1-p[i])*p[i]*(1-p[i]); // 010
    cellprobs[i, 7] = (1-p[i])*(1-p[i])*p[i]; // 001
    cellprobs[i, 8] = 1-sum(cellprobs[i,1:7]); // 000 // INDIVIDUAL NEVER OBSERVED
    // define conditional cell probabilities (each cell prob relative contribution to total prob.)
    cellprobs_conditional[i,1] = cellprobs[i, 1]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,2] = cellprobs[i, 2]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,3] = cellprobs[i, 3]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,4] = cellprobs[i, 4]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,5] = cellprobs[i, 5]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,6] = cellprobs[i, 6]/sum(cellprobs[i, 1:7]);
    cellprobs_conditional[i,7] = cellprobs[i, 7]/sum(cellprobs[i, 1:7]);
    // probability of detecting an individual >= 1 times
    p_det[i] = sum(cellprobs[i,1:7]); 
    
  }
  
}

model {
  
  // Priors
  mu_alpha0 ~ normal(0, 2); // abundance intercept
  mu_alpha1 ~ normal(0, 2); // effect of covariate on abundance
  mu_beta0 ~ normal(0, 2); // detection intercept
  mu_beta1 ~ normal(0, 2); // effect of covariate on detection
  
  // Likelihood
  for (i in 1:n_sites) {
    
    // Roughly how I would write this if I were using Jags:
    //y[i,1:7] ~ multinomial(cellprobs_conditional[i,1:7]);
    //nobs[i] ~ binomial(N[i], p_det[i]);
    //N[i] ~ poisson(lambda[i]);
    
    // marginalized across possible total abundances nobs:K
    vector[K - nobs[i] + 1] lp; // lp vector of length of possible abundances 
      
      for (j in 1:(K - nobs[i] + 1)) { // for each possible abundance:
      
        // lp of capture histories conditional on
        // abundance model and
        // observation model
        lp[j] = 
          // prob of N individuals per site given lambda per site
          poisson_log_lpmf(nobs[i] + j - 1 | lambda[i]) +
          // prob of observing nobs individuals given N per site and p_det odds of detecting individuals >= 1 times
          binomial_lpmf(nobs[i] | nobs[i] + j - 1, p_det[i]) +
          // prob of the site-specific detection histories, given the site-specific conditional cell probabilities
          multinomial_lpmf(y[i,1:7] | to_vector(cellprobs_conditional[i,1:7]))
        ; // end lp[j]
          
        target += log_sum_exp(lp);
      }
      
  }
  
}

generated quantities{
  
  // predict total abundance in the system 
  // not really an important result (for me) but it's a nice initial check of
  // whether the model agrees with the simulation
  // develop posterior predictive check next
  
  vector[n_sites] N; // abundance per site (to be simulated)
  real totalN; // total abundance across all sites
  
  for(i in 1:n_sites){
    N[i] = poisson_rng(lambda[i]); // simulate site-specific abundance
  }

  totalN = sum(N); // sum abundance across sites

}
