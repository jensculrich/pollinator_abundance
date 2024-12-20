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
  
  ilogit <- function(x) exp(x)/(1+exp(x))
  
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
    theta[i] <- ilogit(
      theta0 +
      theta_habitat_effect * X[i]
    )
  }
  
  abundance <- vector(length=R)
  if(poisson == TRUE){
    
    for(i in 1:length(lambda)){
      # simulate true abundance
      abundance[i] <- rpois(n = 1, lambda = lambda[i]) 
    }
    
  } else {
    
    # Add dispersion: draw N from neg_bin(mu, phi)
    for(i in 1:length(lambda)){
      # simulate true abundance
      abundance[i] <- rnbinom(n = 1, mu = lambda[i], size = phi) 
    }
    
  }
  
  
  
  df <- as.data.frame(cbind(df, abundance))
  
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
  
  
  
  # now get a detection rate per species*site*year
  p <- vector(length=R)
  for(i in 1:length(p)){
    p[i] <- ilogit(species_intercept_detection_data[i] +
                     species_slope_detection_data[i] * X[i] +
                     beta2 * year[i]
    )
  }
  
  df <- as.data.frame(cbind(df, theta, p))
  
  # simulate some detection histories
  det_histories <- matrix(nrow = R, ncol = 8)
  possible_histories <- as.data.frame(rbind(
    c('1', '1', '1'),
    c('1', '1', '0'),
    c('1', '0', '1'),
    c('1', '0', '0'),
    c('0', '1', '1'),
    c('0', '1', '0'),
    c('0', '0', '1'),
    c('0', '0', '0')
  )) %>%
    mutate(V1 = as.integer(V1),
           V2 = as.integer(V2),
           V3 = as.integer(V3))
  
  # matrix of counts (for binmix or glm)
  y_counts <- matrix(NA, nrow = length(p), ncol = n_visits) # counts
  # matrix of histories (for multimix)
  y_histories <- matrix(NA, nrow = length(abundance), ncol = 8) # 8 cols = 8 different possible detection histories 
  
  for(r in 1:R){
    available <- matrix(nrow = df$abundance[r], ncol = n_visits)
    detected <- matrix(nrow = df$abundance[r], ncol = n_visits)
    if(length(available) > 0){
      for(i in 1:nrow(available)){
        for(j in 1:n_visits){
          available[i,j] <- rbinom(1, 1, theta[r]) 
          detected[i,j] = available[i,j] * rbinom(1, 1, p[r]) 
        }
      }
      temp <- as.data.frame(detected) %>%
        count(V1, V2, V3) %>%
        full_join(., possible_histories) %>%
        mutate(n = replace_na(n, 0)) %>%
        mutate(history = paste0(V1,V2,V3)) %>%
        select(history, n) %>%
        arrange(history, levels = c('000', '001', '010', '011', '100', '101', '110', '111')) %>%
        arrange(desc(history))
      
      y_counts[r, 1:n_visits] = colSums(detected)
      y_histories[r, 1:8] = temp$n
    }
    
    # when there were zero individuals present at the site, 
    # the above protocol generates NAs for those y_counts and histories. 
    # we can just replace na's with zeros since by definition we detect zero of something that's not present
    y_counts <- replace_na(y_counts, 0)
    y_histories <- replace_na(y_histories, 0)
      
  }
  
  # now prep data for model run
  if(type == "multimix"){ # prep data for multinomial Nmix
    
    y_w_names <- y_histories
    y <- y_histories[,-8]
    
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
    
    # use the count data
    y <- y_counts
      
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
    y <- y_counts
    
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