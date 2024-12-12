## simulate data for abundance from a Poisson or Negative Binomial distribution
## with multinomial detection error, to fit with a multinomial (capture history) nmix model
## eventually I want to expand this into a multi-species multi-year model,
## hence the options to add species and year effects. Setting the species and year effects
## will simulate data for a single species for a singe year across n_sites spatial units.

library(tidyverse) # required for data carpentry
library(rstan) # required to fit stan model

##########################
### Simulate data ########
##########################

# define study dimensions and some predictor variable values
# consider the effect of site covariates on abundance
n_sites = 20 # number of sites # must be an even number
n_species = 16 # number of species
n_visits = 3 # number of repeat visits (= number of temporal reps)
n_years = 2 # number of years

mu_alpha0 = 1.5 # abundance intercept
sigma_alpha0_species = 0.5 # community variation in abundance intercept
mu_alpha1 = 1 # community mean abundance response to management
sigma_alpha1_species = 1 # community variation in abundance response to management
mu_alpha2 = -0.3 # community mean in abundance response to year (only two years so treat as a binary, not categorcial effect)
sigma_alpha2_species = 0.5 # community variation in abundance response to year
sigma_alpha3_site = 0.5 # among site variation in abundance random effect

mu_beta0 = -1 # detection intercept
sigma_beta0_species = 0.5 # community variation in detection intercept
mu_beta1 = 0.5 # community mean detection response to management
sigma_beta1_species = 0.5 # community variation in detection response to management
beta2 = 0.5 # effect of year on detection rate (i.e., maybe we get better over time)

poisson = FALSE # if false, simulate data from a negative binomial distribution
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
                  species_year_data[i] * year[i] +
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
                   beta2 * year[i]
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
  
  cell_probs = matrix(nrow=nrow(df), ncol=8)
  for(i in 1:nrow(df)){
    cell_probs[i,1] = df$p[i]*df$p[i]*df$p[i] # 111
    cell_probs[i,2] = df$p[i]*df$p[i]*(1-df$p[i]) # 110
    cell_probs[i,3] = df$p[i]*(1-df$p[i])*df$p[i] # 101
    cell_probs[i,4] = df$p[i]*(1-df$p[i])*(1-df$p[i]) # 100
    cell_probs[i,5] = (1-df$p[i])*df$p[i]*df$p[i] # 011
    cell_probs[i,6] = (1-df$p[i])*df$p[i]*(1-df$p[i]) # 010
    cell_probs[i,7] = (1-df$p[i])*(1-df$p[i])*df$p[i] # 001
    cell_probs[i,8] = (1-df$p[i])*(1-df$p[i])*(1-df$p[i]) # 000
  }
  
  # Make a 'census' (i.e., go out and count things)
  y <- matrix(NA, nrow = length(abundance), ncol = 8) # 8 cols = 8 different possible detection histories 
  for(i in 1:nrow(y)){
    
    # rmultinom(n, size, prob)
    y[i, 1:8] <- rmultinom(1, size = abundance[i], prob = as.numeric(cell_probs[i, 1:8]))
      
  }
  
  # we will use a data augmentation approach
  # we will need to stretch out the encounter frequencies into a row per individual
  # First, remove the 8th column (not detected individuals, i.e., the data we don't see)
  y_cell_matrix <- y
  
  y_cell_matrix <- y_cell_matrix[,-8]
  
  # get nobs by species*site*year
  nobs <- apply(y_cell_matrix, 1, sum)
  
  K <- (nobs + 3 * 5)
  
  colnames(y) <- c("111", "110", "101","100", "011", "010", "001", "000")
  
  # Return stuff
  return(list(
    # simulated data outcomes 
    R = R,
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
            "mu_alpha2",
            "sigma_alpha2_species",
            "sigma_alpha3_site",
            "scale_param",
            "mu_beta0",
            "sigma_beta0_species",
            "mu_beta1",
            "sigma_beta1_species",
            "beta2",
            "fit",
            "fit_new",
            "totalN"
)


# MCMC settings
n_iterations <- 300
n_thin <- 1
n_burnin <- n_iterations / 2
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
       sigma_alpha3_site = runif(1, 0.25, 0.5),
       mu_beta0 = runif(1, -1, 1),
       sigma_beta0_species = runif(1, 0, 1),
       mu_beta1 = runif(1, -1, 1),
       sigma_beta1_species = runif(1, 0, 1),
       beta2 = runif(1, -0.5, 0.5)
  )
)

# Call STAN model from R 
stan_model <- "./multinomial_nmix/models/multinomial_Nmix_model_negbin.stan"

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

traceplot(stan_out_sim)
pairs(stan_out_sim, pars = c("mu_alpha0", "sigma_alpha0_species",
                      "mu_alpha1", "sigma_alpha1_species",
                      "sigma_alpha3_site",
                      "mu_beta0", "sigma_beta0_species",
                      "mu_beta1", "sigma_beta1_species"))

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

p <- mcmc_hist(stan_out_sim, pars = c("mu_beta1"))
p <- p + labs(x = "mu_beta1",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_beta1, linetype = "solid", size = 1)
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
     ylim = c(400, 1500),
     xlim = c(400, 1500))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)

