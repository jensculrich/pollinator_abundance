// Poisson GLM
// with species-specific intercepts and slopes

data {
  
  int<lower=0> R; // Number of observation units (species*sites*years)
  int<lower=0> n_visits; // number of rep visits per site
  int<lower=0> y[R, n_visits]; // Counts observed on observation 1:R
  vector[R] X; // Habitat covariate (1 = restored or 0 = control)
  vector[R] year; // Year covariate (1 = 2023 or 0 = 2022)
  int<lower=0> n_sites; // Number of sites
  int<lower=1> sites[R]; // vector of sites
  int<lower=0> n_species; // Number of species
  int<lower=1> species[R]; // vector of species
  int<lower=0> K[R];       // Upper bound of population size
}
transformed data {
  int<lower=0> max_y[R];
  for (i in 1:R) {
    max_y[i] = max(y[i]);
  }
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
  // OTHER
  real<lower=0> scale_param;
}
transformed parameters {
  vector[R] log_mu; // Log population size (mu as the centering parameter
    // for the negative binomial count distribution).
  vector[R] logit_p; // logit-scaled detection probability
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
    log_mu[i] = // log abundance is equal to..
      alpha0_species[species[i]] + // a species specific intercept  
      (alpha1_species[species[i]] * X[i]) + // a species specific effect (restoration)
      (alpha2 * year[i]) + // a species specific effect (year)
      alpha3_site[sites[i]]; // a site specific intercept  
    
    logit_p[i] = // logit (detection probability) is equal to..
        beta0_species[species[i]] + // a species specific intercept  
        (beta1_species[species[i]] * X[i]) + // a species specific effect (restoration)
        (beta2 * year[i]) + // an effect (year)
        beta3_site[sites[i]]; // a site specific intercept  
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

  // OTHER
  // neg bin scale param
  scale_param ~ normal(0, 1);
  
  // Likelihood
  for (i in 1:R) { // for each siteXspecies
    vector[K[i] - max_y[i] + 1] lp; // lp vector of length of possible abundances 
    // (from max observed to K) 
    
    for (j in 1:(K[i] - max_y[i] + 1)) // for each possible abundance:
    // lp of abundance given ecological model and observational model
    lp[j] = neg_binomial_2_log_lpmf(max_y[i] + j - 1 | log_mu[i], scale_param)
      + binomial_logit_lpmf(y[i] | max_y[i] + j - 1, logit_p[i]); // vectorized over T
    target += log_sum_exp(lp);
  }
  
}

generated quantities {
  
  int<lower=0> N[R]; // predicted abundance at each site for each species
  vector[R] p; // detection probability of speciesXsite combos

  matrix[R, n_visits] eval; // Expected values
  
  int y_new[R, n_visits]; // new data for counts generated from eval
    
  matrix[R, n_visits] E; // squared scaled distance of real data from expected value
  matrix[R, n_visits] E_new; // squared scaled distance of new data from expected value
  
  real fit = 0; // sum squared distances of real data across all observation intervals
  real fit_new = 0; // sum squared distances of new data across all observation intervals

  vector[n_species] totalN; // total pop size PER SPECIES
 
  // predict abundance given log_mu
  for (i in 1:R) {
    N[i] = neg_binomial_2_log_rng(log_mu[i], scale_param);
    p[i] = inv_logit(logit_p[i]);
  }
  
  // Bayesian p-value fit. 
    
  // Initialize E and E_new
  for (i in 1:1) {
    for(j in 1:n_visits) {
      E[i,j] = 0;
      E_new[i,j] = 0;
    }
  }
  
  for (i in 2:R) {
    E[i] = E[i - 1];
    E_new[i] = E_new[i - 1];
  }
  
  for (i in 1:R) {
    for (j in 1:n_visits) {
      
      // Assess model fit using Chi-squared discrepancy
      // Compute fit statistic E for observed data
      eval[i,j] = p[i] * neg_binomial_2_log_rng(log_mu[i], scale_param); // expected value at observation i for visit j 
        // (probabilty across visits is fixed) is = expected detection prob * expected abundance
      // Compute fit statistic E_new for real data (y)
      E[i,j] = square(y[i,j] - eval[i,j]) / (eval[i,j] + 0.5);
      // Generate new replicate data and
      y_new[i,j] = binomial_rng(N[i], p[i]); 
      // Compute fit statistic E_new for replicate data
      E_new[i,j] = square(y_new[i,j] - eval[i,j]) / (eval[i,j] + 0.5);
      
    }
    
    fit = fit + sum(E[i]); // descrepancies for each siteXspecies combo 
    fit_new = fit_new + sum(E_new[i]); // descrepancies for generated data for each 
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
