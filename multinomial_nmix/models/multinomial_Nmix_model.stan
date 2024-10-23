// Multinomial N-mixture model
// for estimating effects of covariate on abundance and detection rate
// given a set of capture/encounter histories (y) from N (value unknown) individuals
// across a spatially stratified system consisting of n_sites units.
// There are 3 repeat surveys per site -> 
// 7 observable detection histories plus an unobservable 8th history (never detected)
data {
  int<lower=0> n_sites; // number of sites
  int<lower=0> y[n_sites, 7]; // capture histories matrix (7 observable histories possible)
  int<lower=0> nobs[n_sites]; // total number of individuals observed at each site, i.e. nobs = sum(y[,1:7])
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
  vector<lower=0, upper=1>[n_sites] p = inv_logit(mu_beta0 + mu_beta1 * X); // site-specific detection, ranges between 0 and 1 
  vector[n_sites] lambda = mu_alpha0 + mu_alpha1 * X;; // site-specific abundance
  real cellprobs[n_sites, 8]; // cell probabilities (probability of landing in a cell)
  real cellprobs_conditional[n_sites, 7]; // cell probs conditional on the sum of all other observable probabilities
  real p_det[n_sites]; // probability of detecting an individual >= once
  
  for(i in 1:n_sites){ // define site-specific cell probabilities (prob. of capture histories)
    real one_minus_p = 1 - p[i];
    cellprobs[i, 1] = p[i]^3; // 111 (observed, observed, observed)
    cellprobs[i, 2] = p[i]^2*one_minus_p; // 110 (observed, observed, not observed)
    cellprobs[i, 3] = cellprobs[i, 2]; // 101
    cellprobs[i, 4] = p[i]*one_minus_p^2; // 100
    cellprobs[i, 5] = cellprobs[i, 2]; // 011
    cellprobs[i, 6] = cellprobs[i, 4]; // 010
    cellprobs[i, 7] = cellprobs[i, 4]; // 001
    cellprobs[i, 8] = 1-sum(cellprobs[i, 1:7]); // 000 // INDIVIDUAL NEVER OBSERVED
    // define conditional cell probabilities (each cell prob relative contribution to total prob.)
    vector[7] qs = to_vector(cellprobs[i, 1:7]);
    cellprobs_conditional[i, 1:7] = to_array_1d(qs / sum(qs));
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
    int length_lp = K - nobs[i] + 1;
    vector[length_lp] lp; // lp vector of length of possible abundances at site i
    vector[7] site_cellprobs = to_vector(cellprobs_conditional[i,1:7]);
    for (j in 1:length_lp) { // for each possible abundance:
      int latent_abundance = nobs[i] + j - 1;
      lp[j] = poisson_log_lpmf(latent_abundance | lambda[i]) + // abundance
            binomial_lpmf(nobs[i] | latent_abundance, p_det[i]) + // detection of n ind
            multinomial_lpmf(y[i,1:7] | site_cellprobs); // observed detection histories
    }
    target += log_sum_exp(lp); // sum across possible latent states
  }
}
generated quantities {
  array[n_sites] real N = poisson_log_rng(lambda);
  real totalN = sum(N);
}
