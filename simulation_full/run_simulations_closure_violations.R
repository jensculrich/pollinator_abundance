## simulate data for abundance from a Poisson or Negative Binomial distribution
## with multinomial detection error, to fit with a multinomial (capture history) nmix model
## eventually I want to expand this into a multi-species multi-year model,
## hence the options to add species and year effects. Setting the species and year effects
## will simulate data for a single species for a singe year across n_sites spatial units.

library(dplyr) # required for data carpentry
library(rstan) # required to fit stan model

##########################
### Simulate data ########
##########################

# type: choose 1 of 3 options
type <- "multimix"
type <- "binmix"

n_sims <- 20
n_draws_per_sim <- 250
mu_alpha1 <- matrix(nrow=n_sims, ncol=n_draws_per_sim)
precision_mu_alpha1 <- vector(length=n_sims)
mu_alpha0 <- matrix(nrow=n_sims, ncol=n_draws_per_sim)
mu_beta0 <- matrix(nrow=n_sims, ncol=n_draws_per_sim)

# define study dimensions and some predictor variable values
# consider the effect of site covariates on abundance
n_sites = 20 # number of sites # must be an even number
n_species = 16 # number of species
n_visits = 3 # number of repeat visits (= number of temporal reps)
n_years = 2 # number of years

mu_alpha0 = 3 # abundance intercept
sigma_alpha0_species = 0.5 # community variation in abundance intercept
mu_alpha1 = 1 # community mean abundance response to management
sigma_alpha1_species = 0.5 # community variation in abundance response to management
mu_alpha2 = 0 # community mean in abundance response to year (only two years so treat as a binary, not categorcial effect)
sigma_alpha2_species = 0 # community variation in abundance response to year
sigma_alpha3_site = 0.5 # among site variation in abundance random effect

# violate abundance closure assumptions?
# set theta == 1 to perfectly satisfy closure assumptions.
# lower values mean that individuals in some true population have some probability < 1
# of being available for detection in the plot during the survey date.
# do these on the logit scale
theta0 = -1 # in the range of...  -1, 0, 1, 2
theta_habitat_effect = 1 # in the range of... -2,-1, 0, 1, 2

mu_beta0 = -2 # detection intercept
sigma_beta0_species = 0.5 # community variation in detection intercept
mu_beta1 = 0 # community mean detection response to management
sigma_beta1_species = 0.5 # community variation in detection response to management
beta2 = 0 # effect of year on detection rate (i.e., maybe we get better over time)

poisson = TRUE # if false, simulate data from a negative binomial distribution
phi = 0.7 #if using a negative binomial distribution, how much dispersion to use?

# set seed makes the random simmed data be identical whether using the GLM, binmix or multimix
seed <- seq(1:n_sims)

# now start the simulation
source("./simulation_full/simulate_function_closure_violations.R")

for(i in 1:n_sims){
  
  set.seed(seed[i])
  my_simulated_data <- simulate_data(n_sites, 
                                     n_species,
                                     n_visits,
                                     
                                     mu_alpha0,
                                     sigma_alpha0_species, 
                                     mu_alpha1, 
                                     sigma_alpha1_species, 
                                     mu_alpha2, 
                                     sigma_alpha2_species, 
                                     sigma_alpha3_site, 
                                     
                                     theta0,
                                     theta_habitat_effect,
                                     
                                     mu_beta0,
                                     sigma_beta0_species, 
                                     mu_beta1, 
                                     sigma_beta1_species, 
                                     beta2, 
                                     
                                     poisson,
                                     phi
  )
  
  R <- my_simulated_data$R
  y <- my_simulated_data$y
  y_names <- my_simulated_data$y_w_names
  nobs <- my_simulated_data$nobs
  X <- my_simulated_data$X
  K <- my_simulated_data$K
  totalN <- my_simulated_data$totalN
  sites <- my_simulated_data$sites
  species <- my_simulated_data$species
  year <- my_simulated_data$year
  
  ##########################
  ### Run model ############
  ##########################
  
  # MCMC settings
  n_iterations <- 300
  n_thin <- 1
  n_burnin <- n_iterations / 2
  n_chains <- 6
  n_cores <- n_chains
  
  if(type == "multimix"){ # prep stan for multinomial nmix
    
    # stan data
    stan_data <- c("R", "K", 
                   "nobs",
                   "X",
                   "n_sites",
                   "sites",
                   "species",
                   "n_species",
                   "year",
                   "y")
    
    # Parameters monitored
    params <- c("mu_alpha0",
                "sigma_alpha0_species",
                "mu_alpha1",
                "sigma_alpha1_species",
                "alpha2",
                "sigma_alpha3_site",
                "mu_beta0",
                "sigma_beta0_species",
                "mu_beta1",
                "sigma_beta1_species",
                "beta2",
                "sigma_beta3_site",
                "fit",
                "fit_new",
                "totalN"
    )
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      list(mu_alpha0 = runif(1, -1, 1),
           sigma_alpha0_species = runif(1, 0, 1),
           mu_alpha1 = runif(1, -1, 1),
           sigma_alpha1_species = runif(1, 0, 1),
           alpha2 = runif(1, -1, 1),
           sigma_alpha3_site = runif(1, 0.25, 0.5),
           mu_beta0 = runif(1, -1, 1),
           sigma_beta0_species = runif(1, 0, 1),
           mu_beta1 = runif(1, -1, 1),
           sigma_beta1_species = runif(1, 0, 1),
           beta2 = runif(1, -0.5, 0.5),
           sigma_beta3_site = runif(1, 0, 1)
      )
    )
    
    # Call STAN model from R 
    stan_model <- "./multimix/models/multimix_poisson.stan"
    
  } else if(type == "binmix"){ # prep stan for binomial nmix
    
    # stan data
    stan_data <- c("R", "K", 
                   "X",
                   "n_sites",
                   "sites",
                   "species",
                   "n_species",
                   "year",
                   "y",
                   "n_visits")
    
    # Parameters monitored
    params <- c("mu_alpha0",
                "sigma_alpha0_species",
                "mu_alpha1",
                "sigma_alpha1_species",
                "alpha2",
                "sigma_alpha3_site",
                #"scale_param",
                "mu_beta0",
                "sigma_beta0_species",
                "mu_beta1",
                "sigma_beta1_species",
                "beta2",
                "sigma_beta3_site",
                "fit",
                "fit_new",
                "totalN"
    )
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      list(mu_alpha0 = runif(1, -1, 1),
           sigma_alpha0_species = runif(1, 0, 1),
           mu_alpha1 = runif(1, -1, 1),
           sigma_alpha1_species = runif(1, 0, 1),
           alpha2 = runif(1, -1, 1),
           sigma_alpha3_site = runif(1, 0.25, 0.5),
           mu_beta0 = runif(1, -1, 1),
           sigma_beta0_species = runif(1, 0, 1),
           mu_beta1 = runif(1, -1, 1),
           sigma_beta1_species = runif(1, 0, 1),
           beta2 = runif(1, -0.5, 0.5),
           sigma_beta3_site = runif(1, 0, 1)
      )
    )
    
    # Call STAN model from R 
    stan_model <- "./binmix/models/binmix_poisson.stan"
    
  } else { # prep stan for glm
    
    # stan data
    stan_data <- c("R",
                   "X",
                   "n_sites",
                   "sites",
                   "species",
                   "n_species",
                   "year",
                   "y",
                   "n_visits")
    
    # Parameters monitored
    params <- c("mu_alpha0",
                "sigma_alpha0_species",
                "mu_alpha1",
                "sigma_alpha1_species",
                "sigma_alpha3_site",
                "scale_param",
                "fit",
                "fit_new",
                "totalN"
    )
    
    ## Initial values
    # given the number of parameters, the chains need some decent initial values
    # otherwise sometimes they have a hard time starting to sample
    inits <- lapply(1:n_chains, function(i)
      list(mu_alpha0 = runif(1, -1, 1),
           sigma_alpha0_species = runif(1, 0, 1),
           mu_alpha1 = runif(1, -1, 1),
           sigma_alpha1_species = runif(1, 0, 1),
           mu_alpha2 = runif(1, -1, 1),
           sigma_alpha2_species = runif(1, 0, 1),
           sigma_alpha3_site = runif(1, 0.25, 0.5)
      )
    )
    
    # Call STAN model from R 
    stan_model <- "./GLM/models/glm_poisson.stan"
    
  }
  
  ## Call Stan from R
  stan_out_sim <- stan(stan_model,
                       data = stan_data, 
                       init = inits, 
                       pars = params,
                       chains = n_chains, iter = n_iterations, 
                       warmup = n_burnin, thin = n_thin,
                       seed = 1,
                       open_progress = FALSE,
                       cores = n_cores)
  
  print(stan_out_sim, digits = 3)
  
  # capture the mean estimate
  #fit_summary <- rstan::summary(stan_out_sim)
  list_of_draws_mu_alpha1 <- as.data.frame(stan_out_sim)[,3]
  mu_alpha1[i, 1:n_draws_per_sim] <- sample(list_of_draws_mu_alpha1, size = n_draws_per_sim)
  
  quantiles <- quantile(list_of_draws_mu_alpha1, c(0.05, 0.95))
  precision_mu_alpha1[i] <- quantiles[2] - quantiles[1]
  
  list_of_draws_mu_alpha0 <- as.data.frame(stan_out_sim)[,1]
  mu_alpha0[i, 1:n_draws_per_sim] <- sample(list_of_draws, size = n_draws_per_sim)
  
  list_of_draws_mu_beta0 <- as.data.frame(stan_out_sim)[,8]
  mu_beta0[i, 1:n_draws_per_sim] <- sample(list_of_draws, size = n_draws_per_sim)
  
  print(i)
  
}

temp_list <- list(mu_alpha1, precision_mu_alpha1, mu_alpha0, mu_beta0)

# save outputs
saveRDS(temp_list, paste0(
  "./simulation_full/simulation_outputs/closure_violations/",
  type, 
  "_alpha1=", mu_alpha1,
  "_theta0=", theta0, 
  "_theta1=", theta_habitat_effect, 
  ".rds"))
