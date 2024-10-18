data {
  int<lower=0> M; // Size of augumented data set
  int<lower=0> n_visits; // Number of sampling occasions
  int<lower=0> C; // Size of observed data set
  array[M, n_visits] int<lower=0, upper=1> y; // Capture-history matrix
  array[C] int<lower=0, upper=1> X; // site restored
}

transformed data {
  array[M] int<lower=0> s; // Totals in each row
  
  for (i in 1 : M) {
    s[i] = sum(y[i]);
  }
}
parameters {
  real<lower=0, upper=1> omega; // Inclusion probability
  vector<lower=0, upper=1>[n_visits] mean_p; // Mean detection probability
  real alpha;
  real beta;
  real mu_type;
  real<lower=0> sd_type;
  
  array[M - C] real<lower=0, upper=1> X_mis; // Missing data

}
transformed parameters {
  vector[M] logit_p;
  
  for (i in 1 : C) {
    logit_p[i] = alpha + beta * X[i];
  }
  for (i in (C + 1) : M) {
    logit_p[i] = alpha + beta * X_mis[i - C];
  }
}
model {
  // Priors
  //  omega ~ uniform(0, 1);
  //  mean_p ~ uniform(0, 1);
  alpha ~ normal(0, 2);
  beta ~ normal(0, 2);
  mu_type ~ normal(0, 1);
  sd_type ~ uniform(0, 1);
  
  // Likelihood
  for (i in 1 : C) {
    X[i] ~ normal(mu_type, sd_type) T[0, 1];
  }
  for (i in (C + 1) : M) {
    X_mis[i - C] ~ normal(mu_type, sd_type) T[0, 1];
  }
  
  for (i in 1 : M) {
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
generated quantities {
  vector<lower=0, upper=1>[M] p = inv_logit(logit_p);
  array[M] int<lower=0, upper=1> z;
  int<lower=C> totalN;
  
  for (i in 1 : M) {
    if (s[i] > 0) {
      // species was detected at least once
      z[i] = 1;
    } else {
      // species was never detected
      // prob never detected given present
      real pr = prod(rep_vector(1, n_visits) - p[i]);
      z[i] = bernoulli_rng(omega * pr / (omega * pr + (1 - omega)));
    }
  }
  totalN = sum(z);
}
