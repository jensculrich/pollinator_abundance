## simulate data for abundance from a Poisson or Negative Binomial distribution
## with multinomial detection error, to fit with a multinomial (capture history) nmix model
## eventually I want to expand this into a multi-species multi-year model,
## hence the options to add species and year effects. Setting the species and year effects
## will simulate data for a single species for a singe year across n_sites spatial units.

library(tidyverse) # required for data carpentry

##########################
### Simulate data ########
##########################

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
    
    theta0,
    theta_habitat_effect,
    
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
  
  theta <- vector(length=R)
  for(i in 1:length(theta)){
    theta[i] <- theta0 +
      theta_habitat_effect * X[i]
  }
  
  abundance <- vector(length=R)
  abundance_avail <- matrix(nrow=R, ncol=n_visits)
  if(poisson == TRUE){
    
    for(i in 1:length(lambda)){
      # simulate true abundance
      abundance[i] <- rpois(n = 1, lambda = lambda[i]) 
      # violate closure assumptions
      for(j in 1:n_visits){
        abundance_avail[i, j] <- rbinom(n = 1, size = abundance[i], prob = theta[i])
      }
    }
    
  } else {
    
    # Add dispersion: draw N from neg_bin(mu, phi)
    for(i in 1:length(lambda)){
      # simulate true abundance
      abundance[i] <- rnbinom(n = 1, mu = lambda[i], size = phi) 
      # violate closure assumptions
      for(j in 1:n_visits){
        abundance_avail[i, j] <- rbinom(n = 1, size = abundance[i], prob = theta[i])
      }
    }
    
  }
  
  df <- as.data.frame(cbind(df, abundance, abundance_avail))
  
  # splitting in half allows to add across the two years of data
  totalN1 <- vector(length=n_species) # year 1 abundances
  totalN1 <- rowsum(abundance[1:(R/2)], rep(1:n_species, times = n_sites))
  totalN2 <- vector(length=n_species) # year 2 abundances
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
  
  df <- as.data.frame(cbind(df, theta, p))
  
  if(type == "multimix"){ # prep data for multinomial Nmix
    
    # construct cell probs
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
    
    # prob of encounter = probability that the individual is in the survey area (theta)
    # times the probability of detecting the individual given that it is in the survey area (p)
    cell_probs = matrix(nrow=nrow(df), ncol=8)
    for(i in 1:nrow(df)){
      cell_probs[i,1] = (theta[i]*df$p[i])*(theta[i]*df$p[i])*(theta[i]*df$p[i]) # 111
      cell_probs[i,2] = (theta[i]*df$p[i])*(theta[i]*df$p[i])*(1-(theta[i]*df$p[i])) # 110
      cell_probs[i,3] = (theta[i]*df$p[i])*(1-(theta[i]*df$p[i]))*(theta[i]*df$p[i]) # 101
      cell_probs[i,4] = (theta[i]*df$p[i])*(1-(theta[i]*df$p[i]))*(1-(theta[i]*df$p[i])) # 100
      cell_probs[i,5] = (1-(theta[i]*df$p[i]))*(theta[i]*df$p[i])*(theta[i]*df$p[i]) # 011
      cell_probs[i,6] = (1-(theta[i]*df$p[i]))*(theta[i]*df$p[i])*(1-(theta[i]*df$p[i])) # 010
      cell_probs[i,7] = (1-(theta[i]*df$p[i]))*(1-(theta[i]*df$p[i]))*(theta[i]*df$p[i]) # 001
      cell_probs[i,8] = (1-(theta[i]*df$p[i]))*(1-(theta[i]*df$p[i]))*(1-(theta[i]*df$p[i])) # 000
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
    y_w_names <- y # save one that has the column names (for readability in the output)
    
    y <- y[,-8]
    
    # get nobs by species*site*year
    nobs <- apply(y, 1, sum)
    
    # set a ceiling on likelihood search
    if(mu_beta0 < -1.1){
      K <- ((nobs + 3) * 10)
    } else{
      K <- ((nobs + 1) * 10)
    }
    
    colnames(y_w_names) <- c("111", "110", "101","100", "011", "010", "001", "000")
    
  } else if(type == "binmix") { # prep data for binomial Nmix
    
    # Make a 'census' (i.e., go out and count things)
    y <- matrix(NA, nrow = length(p), ncol = n_visits) # Array for counts
    for(i in 1:nrow(y)){
      for(j in 1:n_visits){
        y[i,j] <- rbinom(n = 1, size = abundance_avail[i], prob = p[i])
      }
    }
    
    # set a ceiling on likelihood search
    y_max <- apply(y, 1, max)
    
    if(mu_beta0 < -1.1){
      K <- ((y_max + 3) * 10)
    } else{
      K <- ((y_max + 1) * 10)
    }
    
    nobs <- NULL
    y_w_names <- NULL
    
  } else { # prep data for GLM
    
    # Make a 'census' (i.e., go out and count things)
    y <- matrix(NA, nrow = length(p), ncol = n_visits) # Array for counts
    for(i in 1:nrow(y)){
      for(j in 1:n_visits){
        y[i,j] <- rbinom(n = 1, size = abundance_avail[i], prob = p[i])
      }
    }
    
    K <- NULL
    nobs <- NULL
    y_w_names <- NULL
    
  }
  
  
  # Return stuff
  return(list(
    # simulated data outcomes 
    R = R,
    y_w_names = y_w_names,
    y = y,
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