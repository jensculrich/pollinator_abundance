## simulate data for abundance from a Poisson or Negative Binomial distribution
## with multinomial detection error, to fit with a multinomial (capture history) nmix model

library(tidyverse) # for data carpentry

##########################
### Simulate data ########
##########################

# define study dimensions and some predictor variable values
# consider the effect of site covariates on abundance
n_sites = 32 # number of sites # must be an even number
n_species = 1 # number of species
n_visits = 3 # number of repeat visits (= number of temporal reps)
n_years = 1 # number of years

mu_alpha0 = 1.5 # abundance intercept
sigma_alpha0_species = 0 # community variation in abundance intercept
mu_alpha1 = 1 # community mean abundance response to management
sigma_alpha1_species = 0 # community variation in abundance response to management
mu_alpha2 = 0 # community mean in abundance response to year (only two years so treat as a binary, not categorcial effect)
sigma_alpha2_species = 0 # community variation in abundance response to year
sigma_alpha3_site = 0 # among site variation in abundance random effect

mu_beta0 = -0.5 # detection intercept
sigma_beta0_species = 0 # community variation in detection intercept
mu_beta1 = 0.5 # community mean detection response to management
sigma_beta1_species = 0 # community variation in detection response to management
beta2 = 0 # effect of year on detection rate (i.e., maybe we get better over time)

poisson = TRUE # if false, simulate data from a negative binomial distribution
phi = 0.7 #if using a negative binomial distribution, how much dispersion to use?

M_multiplier = 2.5 # multiplier for data augmentation (positive encounters * "" = max pop size)

# Define function for generating binom-mix model data
simulate_data <- function(
    n_sites, 
    n_species,
    n_visits,
    
    mu_alpha0,
    sigma_alpha0_species, 
    mu_alpha1, 
    sigma_alpha1_species, 
    mu_alpha2, 
    sigma_alpha2_species, 
    sigma_alpha3_site, 
    
    mu_beta0,
    sigma_beta0_species, 
    mu_beta1, 
    sigma_beta1_species, 
    beta2,
    
    poisson,
    phi,
    
    M_multiplier
){
 
  # number of unique data observation points (observations of unique species at unique site)
  R = n_sites*n_species*n_years
  
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = (R/2))) # half sites in each category
  #X <- rep(X, times = 2)
  
  sites <- rep(1:n_sites, each = n_species, times = n_years)
  #year <- rep(c(0,1), each = (R/2))
  species <- rep(1:n_species, times = (R / n_species))
   
  ## Ecological process
  
  ## add a species-level random effect for abundance
  # both on the intercept (some species more abundant than others) centered on 0
  species_intercept <- rnorm(n_species, mu_alpha0, sigma_alpha0_species)
  (mean(species_intercept)) # should be near mu_alpha0
  # and on the slope (change in abundance of different species responds differently to management)
  # centered on community mean
  species_slope <- rnorm(n_species, mu_alpha1, sigma_alpha1_species)
  (mean(species_slope)) # should be near mu_alpha1
  # and on the slope (change in abundance of different species responds differently to year)
  # centered on community mean
  species_year_effect <- rnorm(n_species, mu_alpha2, sigma_alpha2_species)
  (mean(species_year_effect)) # should be near mu_alpha2
  # add a site-level random effect for abundance
  # centered on 0 (because we already have a global intercept mu_alpha0)
  site_intercept <- rnorm(n_sites, 0, sigma_alpha3_site)
  (mean(site_intercept)) # should be near mu_alpha2
  
  species_intercept_data <- as.numeric(vector(length=R))
  species_intercept_data <- rep(species_intercept[1:n_species], 
                                times=R/n_species)
  
  species_slope_data <- as.numeric(vector(length=R))
  species_slope_data <- rep(species_slope[1:n_species], 
                            times=R/n_species)
  
  species_year_data <- as.numeric(vector(length=R))
  species_year_data <- rep(species_year_effect[1:n_species], 
                            times=R/n_species)
  
  site_intercept_data <- as.numeric(vector(length=R))
  site_intercept_data <- rep(site_intercept[1:n_sites], 
                                each=n_species, times=n_years)
  
  df <- as.data.frame(cbind(sites, X, species, 
                            species_intercept_data, species_slope_data, species_year_data,
                            site_intercept_data))
  
  # Relationship expected abundance – covariate + variation in species
  # determined in part by species-level effect on abundance
  lambda <- vector(length=R)
  
  for(i in 1:length(lambda)){ 
    lambda[i] <- exp(
                  species_intercept_data[i] +
                  species_slope_data[i] * X[i] +
                  site_intercept_data[i]
    )
  }
  
  abundance <- vector(length=R)
  if(poisson == TRUE){
    
    for(i in 1:length(lambda)){
      abundance[i] <- rpois(n = 1, lambda = lambda[i]) 
    }
    
  } else {
    
    # Add dispersion: draw N from neg_bin(mu, phi)
    for(i in 1:length(lambda)){
      abundance[i] <- rnbinom(n = 1, mu = lambda[i], size = phi)
    }

  }
  
  df <- as.data.frame(cbind(df, abundance))
  
  totalN = sum(abundance)
  
  ##############################################################################
  ## add imperfect detection process
  
  ## add a species-level random effect for abundance
  # both on the intercept (some species more abundant than others) centered on 0
  species_intercept_detection <- rnorm(n_species, mu_beta0, sigma_beta0_species)
  (mean(species_intercept_detection)) # should be near mu_beta0
  # and on the slope (change in abundance of different species responds differently to management)
  # centered on community mean
  species_slope_detection <- rnorm(n_species, mu_beta1, sigma_beta1_species)
  (mean(species_slope_detection)) # should be near mu_beta1
  
  species_intercept_detection_data <- as.numeric(vector(length=R))
  species_intercept_detection_data <- rep(species_intercept_detection[1:n_species], 
                                times=R/n_species)
  
  species_slope_detection_data <- as.numeric(vector(length=R))
  species_slope_detection_data <- rep(species_slope_detection[1:n_species], 
                            times=R/n_species)
  
  df <- as.data.frame(cbind(df, 
                            species_intercept_detection_data, 
                            species_slope_detection_data))
  
  
  ilogit <- function(x) exp(x)/(1+exp(x))
  # now get a detection rate per species*site*year
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- ilogit(species_intercept_detection_data[i] +
                     species_slope_detection_data[i] * X[i] #+
                   #beta2 * year[i]
                   )
  }
  
  df <- as.data.frame(cbind(df, p))
  
  # construct multinomial cell probabilities for three potential detection events
  #cellprobs <- c(
    #p*p*p, # 111
    #p*p*(1-p), # 110
    #p*(1-p)*p, # 101
    #p*(1-p)*(1-p), # 100
    #(1-p)*p*p, # 011
    #(1-p)*p*(1-p), # 010
    #(1-p)*(1-p)*p, # 001
    #(1-p)*(1-p)*(1-p) # 000
  #)
  
  df <- df %>%
    mutate("p111" = p*p*p, # 111
           "p110" = p*p*(1-p), # 110
           "p101" = p*(1-p)*p, # 101
           "p100" = p*(1-p)*(1-p), # 100
           "p011" = (1-p)*p*p, # 011
           "p010" = (1-p)*p*(1-p), # 010
           "p001" = (1-p)*(1-p)*p, # 001
           "p000" = (1-p)*(1-p)*(1-p) # 000)
    )
  
  # Make a 'census' (i.e., go out and count things)
  y <- matrix(NA, nrow = length(abundance), ncol = 8) # 8 cols = 8 different possible detection histories 
  for(i in 1:nrow(y)){
    
    # rmultinom(n, size, prob)
    y[i, 1:8] <- rmultinom(1, size = abundance[i], prob = as.numeric(df[i, 12:19]))
      
  }
  
  # we will use a data augmentation approach
  # we will need to stretch out the encounter frequencies into a row per individual
  # First, remove the 8th column (not detected individuals, i.e., the data we don't see)
  y <- y[,-8]
  
  y_cell_matrix <- y
  
  nobs <- apply(y_cell_matrix, 1, sum)
  
  K = max(nobs) * 2
  
  colnames(y) <- c("111", "110", "101","100", "011", "010", "001")
  
  # Return stuff
  return(list(
    # simulated data outcomes 
    y_w_names = y,
    y = y_cell_matrix,
    nobs = nobs,
    K = K,
    X = X, # simulated covariate values,
    sites = sites,
    species = species,
    year = year,
    df = df, # combined year, site, species, covariate data
    abundance = abundance, # simulated abundance
    totalN = totalN # simulated abundance per species across sites
  ))
     
}

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
                                   
                                   mu_beta0,
                                   sigma_beta0_species, 
                                   mu_beta1, 
                                   sigma_beta1_species, 
                                   beta2, 
                                   
                                   poisson,
                                   phi,
                                   
                                   M_multiplier)

y <- my_simulated_data$y
y_names <- my_simulated_data$y_w_names
nobs <- my_simulated_data$nobs
X <- my_simulated_data$X
K <- my_simulated_data$K
totalN <- my_simulated_data$totalN
sites <- my_simulated_data$sites

##########################
### Run model ############
##########################

stan_data <- c("K", 
               "nobs",
               "X",
               "n_visits",
               "n_sites",
               "sites",
               "y")

# Parameters monitored
params <- c("mu_alpha0",
            "mu_alpha1"
)


# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- 150
n_chains <- 4
n_cores <- 4

## Initial values
# given the number of parameters, the chains need some decent initial values
# otherwise sometimes they have a hard time starting to sample
inits <- lapply(1:n_chains, function(i)
  list(mu_alpha0 = runif(1, -1, 1),
       mu_alpha1 = runif(1, -1, 1),
       mu_beta0 = runif(1, -1, 1),
       mu_beta1 = runif(1, -1, 1)
  )
)

# Call STAN model from R 
stan_model <- "./multinomial_nmix/models/multinomial_Nmix_simple4.stan"

## Call Stan from R
library(rstan)
stan_out_sim <- stan(stan_model,
                   data = stan_data, 
                   init = inits, 
                   pars = params,
                   chains = n_chains, iter = n_iterations, 
                   warmup = n_burnin, thin = n_thin,
                   seed = 1,
                   open_progress = FALSE,
                   cores = n_cores)

# Output:
#"Chain 1: Rejecting initial value:
#Chain 1:   Log probability evaluates to log(0), i.e. negative infinity.
#Chain 1:   Stan can't start sampling from this initial value."

print(stan_out_sim, digits = 3)