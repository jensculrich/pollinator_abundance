// Multinomial N-mixture model
// for estimating effects of covariate on abundance and detection rate
// given a set of capture/encounter histories (y) from N (value unknown) individuals
// across a spatially stratified system consisting of n_sites units.
// There are 3 repeat surveys per site -> 
// 7 observable detection histories plus an unobservable 8th history (never detected)
data {
  int<lower=0> R; // length of the data set (unique species*site*year combinations)
  int<lower=0> n_sites; // number of sites
  int<lower=1> sites[R]; // vector of sites
  int<lower=0> n_species; // number of species
  int<lower=1> species[R]; // vector of species
  int<lower=0> y[R, 7]; // capture histories matrix (7 observable histories possible)
  int<lower=0> nobs[R]; // total number of individuals observed at each site, i.e. nobs = sum(y[,1:7])
  int<lower=0> K[R]; // some arbitrarilly large maximum site-specific abundance to search up until
  vector[R] X; // site covariate (categorical 0 or 1)
  vector[R] year; // year covariate (categorical 0 or 1)
}
parameters {
  // ABUNDANCE
  real mu_alpha0; // global interecept for abundance
  vector[n_species] alpha0_species_raw; // species-specific intercepts
  real<lower=0> sigma_alpha0_species; // among species variation in intercepts
  real mu_alpha1; // community mean effect of restoration
  vector[n_species] alpha1_species_raw; // vector of species specific restoration effects
  real<lower=0> sigma_alpha1_species; // among species variation in restoration effects
  real alpha2; // community mean effect of year
  vector[n_sites] alpha3_site_raw; // vector of site specific intercept effects
  real<lower=0> sigma_alpha3_site; // among site variation in intercepts
  real<lower=0> scale_param;
  // DETECTION
  real mu_beta0; // global interecept for abundance
  vector[n_species] beta0_species_raw; // species-specific intercepts
  real<lower=0> sigma_beta0_species; // among species variation in intercepts
  real mu_beta1; // community mean effect of restoration
  vector[n_species] beta1_species_raw; // vector of species specific restoration effects
  real<lower=0> sigma_beta1_species; // among species variation in restoration effects
  real beta2; //   effect of year
  vector[n_sites] beta3_site_raw; // vector of site specific intercept effects
  real<lower=0> sigma_beta3_site; // among site variation in intercepts
}
transformed parameters {
  vector<lower=0, upper=1>[R] p; // site-specific detection, ranges between 0 and 1 
  vector[R] log_lambda; // site-specific abundance
  real cellprobs[R, 8]; // cell probabilities (probability of landing in a cell)
  real cellprobs_conditional[R, 7]; // cell probs conditional on the sum of all other observable probabilities
  real p_det[R]; // probability of detecting an individual >= once
  
  // non-centered species-specific effects
  vector[n_species] alpha0_species;
  vector[n_species] alpha1_species;
  vector[n_sites] alpha3_site;
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  alpha0_species = mu_alpha0 + sigma_alpha0_species * alpha0_species_raw;
  alpha1_species = mu_alpha1 + sigma_alpha1_species * alpha1_species_raw;
  alpha3_site = 0 + sigma_alpha3_site * alpha3_site_raw; // centered on 0 (not adding another intercept)
  // non-centered species-specific effects
  vector[n_species] beta0_species;
  vector[n_species] beta1_species;
  vector[n_sites] beta3_site;
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  beta0_species = mu_beta0 + sigma_beta0_species * beta0_species_raw;
  beta1_species = mu_beta1 + sigma_beta1_species * beta1_species_raw;
  beta3_site = 0 + sigma_beta3_site * beta3_site_raw; // centered on 0 (not adding another intercept)
  
  for (i in 1:R) { // for each site*species combination
    log_lambda[i] = // log abundance is equal to..
      alpha0_species[species[i]] + // a species specific intercept  
      (alpha1_species[species[i]] * X[i]) + // a species specific effect (restoration)
      (alpha2 * year[i]) + // a species specific effect (year)
      alpha3_site[sites[i]]; // site random effect
    p[i] = inv_logit( // inv logit (detection probability) is equal to..
        beta0_species[species[i]] + // a species specific intercept  
        (beta1_species[species[i]] * X[i]) + // a species specific effect (restoration)
        (beta2 * year[i]) + // an effect (year)
        beta3_site[sites[i]]); // site random effect 
  }
  
  for(i in 1:R){ // define site-specific cell probabilities (prob. of capture histories)
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
  // ABUNDANCE
  // intercepts
  mu_alpha0 ~ normal(0, 2); // persistence intercept
  alpha0_species_raw ~ std_normal();
  sigma_alpha0_species ~ normal(0, 2);
  // habitat effects
  mu_alpha1 ~ normal(0, 2); // persistence intercept
  alpha1_species_raw ~ std_normal();
  sigma_alpha1_species ~ normal(0, 2);
  // year effects
  alpha2 ~ normal(0, 1); // effect of year
  // site random effect on intercept
  alpha3_site_raw ~ std_normal();
  sigma_alpha3_site ~ normal(0, 1);
  // neg bin scale param
  scale_param ~ normal(0, 1);
  // DETECTION
  // intercepts
  mu_beta0 ~ normal(0, 2); // persistence intercept
  beta0_species_raw ~ std_normal();
  sigma_beta0_species ~ normal(0, 2);
  // habitat effects
  mu_beta1 ~ normal(0, 2); // effect of restoration
  beta1_species_raw ~ std_normal();
  sigma_beta1_species ~ normal(0, 2);
  // year effects
  beta2 ~ normal(0, 1); // effect of year
  // site random effect on intercept
  beta3_site_raw ~ std_normal();
  sigma_beta3_site ~ normal(0, 1);
  
  // Likelihood
  for (i in 1:R) {
    int length_lp = K[i] - nobs[i] + 1;
    vector[length_lp] lp; // lp vector of length of possible abundances at site i
    vector[7] site_cellprobs = to_vector(cellprobs_conditional[i,1:7]);
    for (j in 1:length_lp) { // for each possible abundance:
      int latent_abundance = nobs[i] + j - 1;
      lp[j] = poisson_log_lpmf(latent_abundance | log_lambda[i]) + // abundance
            binomial_lpmf(nobs[i] | latent_abundance, p_det[i]) + // detection of n ind
            multinomial_lpmf(y[i,1:7] | site_cellprobs); // observed detection histories
    }
    target += log_sum_exp(lp); // sum across possible latent states
  }
}

generated quantities {
  
  int<lower=0> N[R]; // predicted abundance at each site for each species\
  vector[R] eval; // Expected values
  int y_new[R]; // new data for counts generated from eval
  vector[R] E; // squared scaled distance of real data from expected value
  vector[R] E_new; // squared scaled distance of new data from expected value
  vector[n_species] totalN; // total pop size PER SPECIES

  real fit = 0; // sum squared distances of real data across all observation intervals
  real fit_new = 0; // sum squared distances of new data across all observation intervals
 
  // predict abundance given log_lambda
  for (i in 1:R) {
    N[i] = poisson_log_rng(log_lambda[i]);
  }
  
  // Bayesian p-value fit. 
    
  // Initialize E and E_new
  for (i in 1:1) {
      E[i] = 0;
      E_new[i] = 0;
  }
  
  for (i in 2:R) {
    E[i] = E[i - 1];
    E_new[i] = E_new[i - 1];
  }
  
  for (i in 1:R) {

      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
      eval[i] = p_det[i] * poisson_log_rng(log_lambda[i]); // expected value at observation i for visit j 
        // (probabilty across visits is fixed) is = expected detection prob * expected abundance
      // Compute fit statistic E_new for real data (y)
      E[i] = square(nobs[i] - eval[i]) / (eval[i] + 0.5);
      // Generate new replicate data and
      y_new[i] = binomial_rng(N[i], p_det[i]); // always detect if there
      // Compute fit statistic E_new for replicate data
      E_new[i] = square(y_new[i] - eval[i]) / (eval[i] + 0.5);
    
    fit = fit + E[i]; // descrepancies for each siteXspecies combo 
    fit_new = fit_new + E_new[i]; // descrepancies for generated data for each 
                                      // siteXspecies combos 
  }

  for (i in 1:n_species){
    vector[(n_sites*2)] totalSpecies; // *2 because there are 2 years
    for (j in 1:(n_sites*2)){ // *2 because there are 2 years
      totalSpecies[j] = N[(i+((n_species*j)-n_species))];
    }
    totalN[i] = sum(totalSpecies);
  }
}

