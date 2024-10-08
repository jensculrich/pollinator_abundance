// Poisson GLM
// with species-specific intercepts and slopes

data {
  
  int<lower=0> R; // Number of observation units (species*sites)
  
  int<lower=0> y[R]; // Counts observed on observation 1:R
  
  vector[R] X; // Habitat covariate (1 = restored or 0 = control)
  vector[R] years; // Year covariate (1 = 2023 or 0 = 2022)
  
  int<lower=0> n_sites; // Number of sites
  int<lower=1> sites[R]; // vector of sites
  
  int<lower=0> n_species; // Number of species
  int<lower=1> species[R]; // vector of species
  
}

parameters {
  
  real mu_alpha0; // global interecept for abundance
  vector[n_species] alpha0_species_raw; // species-specific intercepts
  real<lower=0> sigma_alpha0_species; // among species variation in intercepts

  real mu_alpha1; // community mean effect of restoration
  vector[n_species] alpha1_species_raw; // vector of species specific restoration effects
  real<lower=0> sigma_alpha1_species; // among species variation in restoration effects
  
  real mu_alpha2; // community mean effect of year
  vector[n_species] alpha2_species_raw; // vector of species specific year effects
  real<lower=0> sigma_alpha2_species; // among species variation in year effects

}

transformed parameters {
  
  vector[R] log_mu; // Log population size (mu as the centering parameter
    // for the negative binomial count distribution).
  
  // non-centered species-specific effects
  vector[n_species] alpha0_species;
  vector[n_species] alpha1_species;
  vector[n_species] alpha2_species;
  // implies: xprocess_species ~ normal(mu_xprocess_species, sigma_xprocess_species)
  alpha0_species = mu_alpha0 + sigma_alpha0_species * alpha0_species_raw;
  alpha1_species = mu_alpha1 + sigma_alpha1_species * alpha1_species_raw;
  alpha2_species = mu_alpha2 + sigma_alpha2_species * alpha2_species_raw;
  
  for (i in 1:R) { // for each site*species combination
    log_mu[i] = // log abundance is equal to..
      alpha0_species[species[i]] + // a species specific intercept plus 
      (alpha1_species[species[i]] * X[i]) + // a species specific slope
      (alpha2_species[species[i]] * years[i]); // a species specific slope
  }
}

model {
  
  // Priors

  // intercepts
  mu_alpha0 ~ normal(0, 2); // persistence intercept
  alpha0_species_raw ~ std_normal();
  sigma_alpha0_species ~ normal(0, 2);
  
  // habitat effects
  mu_alpha1 ~ normal(0, 2); // persistence intercept
  alpha1_species_raw ~ std_normal();
  sigma_alpha1_species ~ normal(0, 2);
  
  // year effects
  mu_alpha2 ~ normal(0, 2); // persistence intercept
  alpha2_species_raw ~ std_normal();
  sigma_alpha2_species ~ normal(0, 2);
  
  // Likelihood
  for (i in 1:R) { // for each siteXspecies
    // (from max observed to K) 
    // lp of abundance given ecological model and observational model
    target += poisson_log_lpmf(y[i] | log_mu[i]);
  }
  
}

generated quantities {
  
  int<lower=0> N[R]; // predicted abundance at each site for each species

  vector[R] eval; // Expected values
  
  int y_new[R]; // new data for counts generated from eval
    
  vector[R] E; // squared scaled distance of real data from expected value
  vector[R] E_new; // squared scaled distance of new data from expected value
  
  real fit = 0; // sum squared distances of real data across all observation intervals
  real fit_new = 0; // sum squared distances of new data across all observation intervals

  vector[n_species] totalN; // total pop size PER SPECIES
 
  // predict abundance given log_mu
  for (i in 1:R) {
    N[i] = poisson_log_rng(log_mu[i]);
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
      eval[i] = poisson_log_rng(log_mu[i]); // expected value at observation i for visit j 
        // (probabilty across visits is fixed) is = expected detection prob * expected abundance
      // Compute fit statistic E_new for real data (y)
      E[i] = square(y[i] - eval[i]) / (eval[i] + 0.5);
      // Generate new replicate data and
      y_new[i] = binomial_rng(N[i], 1); // always detect if there
      // Compute fit statistic E_new for replicate data
      E_new[i] = square(y_new[i] - eval[i]) / (eval[i] + 0.5);
    
    fit = fit + E[i]; // descrepancies for each siteXspecies combo 
    fit_new = fit_new + E_new[i]; // descrepancies for generated data for each 
                                      // siteXspecies combos 
  }
  
  // to find total abundance PER SPECIES:
  // sum abundance by rows for each species i
  //for (i in 1:n_species){
    //vector[n_sites] totalSpecies;
    //for (j in 1:n_sites){
     //totalSpecies[j] = N[(i+((n_species*j)-n_species))];
    //}
    //totalN[i] = sum(totalSpecies);
  //}
  // do it with year
  for (i in 1:n_species){
    vector[(n_sites*2)] totalSpecies; // *2 because there are 2 years
    for (j in 1:(n_sites*2)){ // *2 because there are 2 years
      totalSpecies[j] = N[(i+((n_species*j)-n_species))];
    }
    totalN[i] = sum(totalSpecies);
  }
}
