## simulate data for abundance from a Poisson or Negative Binomial distribution
## with detection error, to fit with a Nmix model

##########################
### Simulate data ########
##########################

# define study dimensions and some predictor variable values
# consider the effect of site covariates on abundance
n_sites = 28 # number of sites # must be an even number
n_species = 15 # number of species
n_visits = 4 # number of repeat visits (= number of temporal reps)
n_years = 2 # number of years

mu_alpha0 = 0.5 # abundance intercept
sigma_alpha0_species = 1 # community variation in abundance intercept
mu_alpha1 = 1.5 # community mean abundance response to management
sigma_alpha1_species = 0.8 # community variation in abundance response to management
mu_alpha2 = -0.2 # community mean in abundance response to year (only two years so treat as a binary, not categorcial effect)
sigma_alpha2_species = 0.4 # community variation in abundance response to year

mu_beta0 = 0 # detection intercept
sigma_beta0_species = 0.5 # community variation in detection intercept
mu_beta1 = -0.75 # community mean detection response to management
sigma_beta1_species = 0.25 # community variation in detection response to management
beta2 = 0.2 # effect of year on detection rate (i.e., maybe we get better over time)

poisson = FALSE # if false, simulate data from a negative binomial distribution
phi = 0.7 #if using a negative binomial distribution, how much dispersion to use?

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
    
    mu_beta0,
    sigma_beta0_species, 
    mu_beta1, 
    sigma_beta1_species, 
    beta2,
    
    poisson,
    phi                
){
 
  # number of unique data observation points (observations of unique species at unique site)
  R = n_sites*n_species*n_years
  
  # Covariate values: sort for ease of presentation
  X <- sort(rep(c(0, 1), times = (R/4))) # half sites in each category
  X <- rep(X, times = 2)
  
  sites <- rep(1:n_sites, each = n_species, times = n_years)
  year <- rep(c(0,1), each = (R/2))
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
  
  species_intercept_data <- as.numeric(vector(length=R))
  species_intercept_data <- rep(species_intercept[1:n_species], 
                                times=R/n_species)
  
  species_slope_data <- as.numeric(vector(length=R))
  species_slope_data <- rep(species_slope[1:n_species], 
                            times=R/n_species)
  
  species_year_data <- as.numeric(vector(length=R))
  species_year_data <- rep(species_year_effect[1:n_species], 
                            times=R/n_species)
  
  df <- as.data.frame(cbind(year, sites, X, species, 
                            species_intercept_data, species_slope_data, species_year_data))
  
  # Relationship expected abundance – covariate + variation in species
  # determined in part by species-level effect on abundance
  lambda <- vector(length=R)
  for(i in 1:length(lambda)){ 
    lambda[i] <- exp(
                  species_intercept_data[i] +
                  species_slope_data[i] * X[i] +
                  species_year_data[i] * year[i]
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
  
  totalN1 <- vector(length=n_species)
  totalN1 <- rowsum(abundance[1:(R/2)], rep(1:n_species, times = n_sites))
  totalN2 <- vector(length=n_species)
  totalN2 <- rowsum(abundance[(R/2 + 1):R], rep(1:n_species, times = n_sites))
  totalN = totalN1 + totalN2
  
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
                     species_slope_detection_data[i] * X[i] +
                   beta2 * year[i])
  }
  
  df <- as.data.frame(cbind(df, p))
  
  par(mfrow = c(2, 1))
  hist(p[1:(R/4)], xlim = c(0,1))
  hist(p[(R/4):(R/2)], xlim = c(0,1))
  
  # Make a 'census' (i.e., go out and count things)
  y <- matrix(NA, nrow = length(p), ncol = n_visits) # Array for counts
  for(i in 1:nrow(y)){
    for(j in 1:n_visits){
      y[i,j] <- rbinom(n = 1, size = abundance[i], prob = p[i])
    }
  }
  
  # Return stuff
  return(list(
    # simulated data outcomes 
    R = R,
    X = X, # simulated covariate values,
    sites = sites,
    species = species,
    year = year,
    df = df, # combined year, site, species, covariate data
    abundance = abundance, # simulated abundance
    totalN = totalN, # simulated abundance per species across sites
    y = y # observed abundance
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
                                   
                                   mu_beta0,
                                   sigma_beta0_species, 
                                   mu_beta1, 
                                   sigma_beta1_species, 
                                   beta2, 
                                   
                                   poisson,
                                   phi)

R <- my_simulated_data$R
X <- my_simulated_data$X
species_data <- my_simulated_data$species_data
abundance <- my_simulated_data$abundance
totalN <- my_simulated_data$totalN
y <- my_simulated_data$y
sites <- my_simulated_data$sites
species <- my_simulated_data$species
years <- my_simulated_data$year
df <- my_simulated_data$df

y_max <- apply(y, 1, max)
K <- (y_max + 3 * 5)

##########################
### Run model ############
##########################

stan_data <- c("R", "K",
               "sites",
               "n_sites",
               "species",
               "years",
               "n_species", 
               "n_visits", 
               "y", "X")

# Parameters monitored
params <- c(
            "mu_alpha0",
            "sigma_alpha0_species",
            "mu_alpha1",
            "sigma_alpha1_species",
            "mu_alpha2",
            "sigma_alpha2_species",
            "mu_beta0",
            "sigma_beta0_species",
            "mu_beta1",
            "sigma_beta1_species",
            "beta2",
            "scale_param",
            #"alpha0_species",
            #"alpha1_species",
            #"alpha2_species",
            "fit",
            "fit_new",
            "totalN"
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
       sigma_alpha0_species = runif(1, 0, 1),
       mu_alpha1 = runif(1, -1, 1),
       sigma_alpha1_species = runif(1, 0, 1),
       mu_alpha2 = runif(1, -1, 1),
       sigma_alpha2_species = runif(1, 0, 1),
       mu_beta0 = runif(1, -1, 1),
       sigma_beta0_species = runif(1, 0, 1),
       mu_beta1 = runif(1, -1, 1),
       sigma_beta1_species = runif(1, 0, 1),
       beta2 = runif(1, -1, 1)
  )
)

# Call STAN model from R 
stan_model <- "./N_mixture/models/Nmix1_negbin.stan"

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

print(stan_out_sim, digits = 3)

traceplot(stan_out_sim, pars = c("mu_alpha0", "sigma_alpha0_species",
                                 "mu_alpha1", "sigma_alpha1_species", 
                                 "mu_alpha2", "sigma_alpha2_species",
                                 "scale_param"))


library(bayesplot)
library(tidyverse)

# p: how good is our model at estimating the parameters?
# plot posterior distribution
p <- mcmc_hist(stan_out_sim, pars = c("mu_alpha0"))
p <- p + labs(x = "mu_alpha0",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha0, linetype = "solid", size = 1)
p

p <- mcmc_hist(stan_out_sim, pars = c("mu_alpha1"))
p <- p + labs(x = "mu_alpha1",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(stan_out_sim, pars = c("mu_alpha2"))
p <- p + labs(x = "mu_alpha2",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha2, linetype = "solid", size = 1)
p

# q: how good is our model at estimating the true abundance?
(q <- mcmc_hist(stan_out_sim, pars = c("totalN[1]", "totalN[2]", 
                                      "totalN[3]", "totalN[4]",
                                      "totalN[5]", "totalN[6]"))
)

# plot posterior distribution
q <- mcmc_hist(stan_out_sim, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = my_simulated_data$totalN[1], linetype = "solid", size = 1)
q

q <- mcmc_hist(stan_out_sim, pars = c("totalN[12]"))
q <- q + labs(x = "totalN[12]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = my_simulated_data$totalN[12], linetype = "solid", size = 1)
q

# Evaluation of fit
list_of_draws <- as.data.frame(stan_out_sim)

par(mfrow = c(1, 1))

plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(20000, 150000),
     xlim = c(20000, 150000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)
