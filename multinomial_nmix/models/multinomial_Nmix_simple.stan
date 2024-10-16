// multinomial nmix model for negbinomial distributed abundance
// with species-specific intercepts and slopes

// https://github.com/stan-dev/example-models/blob/master/BPA/Ch.06/M0.stan

data {
  
  int<lower=0> M; // size of augmented data set
  int<lower=0> n_visits; // number of rep visits per site
  
  int<lower=0> y[M, n_visits]; // Counts observed on observation 1:R
  
}
transformed data {
  array[M] int<lower=0> s; // Totals in each row
  int<lower=0> C; // Size of observed data set
  
  C = 0;
  for (i in 1 : M) {
    s[i] = sum(y[i]);
    if (s[i] > 0) {
      C = C + 1;
    }
  }
}

parameters {
  real<lower=0, upper=1> omega; // Inclusion probability
  real<lower=0, upper=1> p; // Detection probability
}


model {
  // Priors are imlicitly defined;
  //  omega ~ uniform(0, 1);
  //  p ~ uniform(0, 1);
  
  // Likelihood
  for (i in 1 : M) {
    if (s[i] > 0) {
      // z[i] == 1
      target += bernoulli_lpmf(1 | omega) + binomial_lpmf(s[i] | n_visits, p);
    } else // s[i] == 0
    {
      target += log_sum_exp(bernoulli_lpmf(1 | omega)
                            // z[i] == 1
                            + binomial_lpmf(0 | n_visits, p),
                            bernoulli_lpmf(0 | omega) // z[i] == 0
                            );
    }
  }
}

generated quantities {
  // prob present given never detected
  real omega_nd = (omega * (1 - p) ^ n_visits) / (omega * (1 - p) ^ n_visits + (1 - omega));
  int<lower=C, upper=M> totalN = C + binomial_rng(M - C, omega_nd);
}
