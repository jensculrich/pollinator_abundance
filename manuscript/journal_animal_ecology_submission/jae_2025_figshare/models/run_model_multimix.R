# fit the Nmix model to my real data
# load required libraries
library(tidyverse)
library(rstan)

################################################################################
# load and format data

####################################
# analysis of real data ############
####################################

# read abundance data
df <- read.csv("./data/abundance_data.csv")
str(df)

df <- df %>%
  # calculate number of individuals with each observable encounter history (by species*site*year group)
  group_by(year, site, species) %>%
  mutate(
    "111" = obs_round_1_and_2_double_recaptures[3],
    "110" = obs_round_1_recaptures[2] - obs_round_1_and_2_double_recaptures[3],
    "101" = obs_round_1_recaptures[3],
    "100" = count[1] - obs_round_1_recaptures[3] - obs_round_1_and_2_double_recaptures[3] 
      - (obs_round_1_and_2_double_recaptures[3] - obs_round_1_recaptures[2]), 
    "011" = obs_round_2_recaptures[3],
    "010" = count[2] - obs_round_1_recaptures[2] - obs_round_2_recaptures[3],
    "001" = count[3] - obs_round_1_recaptures[3] - obs_round_2_recaptures[3] - obs_round_1_and_2_double_recaptures[3] 
    ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(missing_data = ifelse(is.na(.[13]), 0, 1)) # add a missing data variable

n_species <- nrow(distinct(df, species)) # number of species
n_sites <- nrow(distinct(df, site)) # number of sites

R <- nrow(df) # number of observation events
y <- as.matrix(df[,13:19]) %>%
  replace_na(-99) # give stan some numeric value (it will be treated as an NA)
y_names <- df[,13:19] # encounter history names
nobs <- apply(y, 1, sum) # total number of observed individuals at each site/species/year
X <- df$treatment # 0 = control; 1 = restored
K <- ((nobs + 5) * 12) %>% # search possible abundance up to K = ...
  replace_na(-99) # give stan some numeric value (it will be treated as an NA)
nobs <- nobs %>% replace_na(-99)# give stan some numeric value (it will be treated as an NA)
sites <- as.integer(as.factor(df$site)) # unique integer values to represent site names
site_names <- df$site # real site names (character version)
species <- as.integer(as.factor(df$species)) # unique integer values to represent species names
species_names <- df$species # real species names (character version)
year <- df$year # year 0 = 2022; 1 = 2023
missing_data <- as.vector(df$missing_data) # vector indicating whether the row should be treated as NA 

# keep track of integer / character name conversions
species_names_table <- as.data.frame(cbind(cbind(1:n_species), unique(species_names)))
site_names_unique <- unique(site_names)
site_names_table <- as.data.frame(cbind(cbind(1:n_sites), unique(site_names)))


##########################
### Run model ############
##########################

# pack some data for Stan
stan_data <- c("R", "K", 
               "nobs",
               "X",
               "n_sites",
               "sites",
               "species",
               "n_species",
               "year",
               "y",
               "missing_data")

# Parameters monitored
params <- c("mu_alpha0",
            "sigma_alpha0_species",
            "mu_alpha1",
            "sigma_alpha1_species",
            "alpha2",
            "sigma_alpha3_site",
            "scale_param",
            "mu_beta0",
            "sigma_beta0_species",
            "mu_beta1",
            "sigma_beta1_species",
            "beta2",
            "sigma_beta3_site",
            "alpha0_species",
            "alpha1_species",
            "alpha3_site",
            "beta0_species",
            "beta1_species",
            "fit",
            "fit_new",
            "totalN"
)


# MCMC settings
n_iterations <- 4000
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
stan_model <- "./multimix/models/multimix_negbin_NAs.stan"

## Call Stan from R
stan_out <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     control=list(adapt_delta=0.95),
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

# Save and view model outputs
saveRDS(stan_out, "./model_outputs/real_data/multinomial_Nmix.rds")
print(stan_out, digits = 3)

traceplot(stan_out, pars = c("mu_alpha0",
                             "sigma_alpha0_species",
                             "mu_alpha1",
                             "sigma_alpha1_species",
                             "alpha2",
                             "sigma_alpha3_site",
                             "scale_param",
                             "mu_beta0",
                             "sigma_beta0_species",
                             "mu_beta1",
                             "sigma_beta1_species",
                             "beta2",
                             "sigma_beta3_site"))

pairs(stan_out, pars = c("mu_alpha0",
                             "sigma_alpha0_species",
                             "mu_alpha1",
                             "sigma_alpha1_species",
                             "alpha2",
                             "sigma_alpha3_site",
                             "scale_param"))

# Evaluation of fit
list_of_draws <- as.data.frame(stan_out)
plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(000, 150000),
     xlim = c(000, 150000))
abline(0, 1, lwd = 2, col = "black")

Bp <- signif(mean(list_of_draws$fit_new > list_of_draws$fit), 4)
title(paste0("Freeman Tukey P = ", Bp)) 
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)

# some summary views of the param estimates
library(bayesplot)
library(tidyverse)

# p: parameter distributions
# plot posterior distribution
p <- mcmc_hist(stan_out, pars = c("mu_alpha0"))
p <- p + labs(x = "mu_alpha0 (abundance intercept)",
              y = "Frequency in Posterior Distribution") +
  geom_vline(xintercept = mean(list_of_draws$mu_alpha0), linetype = "solid", size = 1)
p

p <- mcmc_hist(stan_out, pars = c("mu_alpha1"))
p <- p + labs(x = "mu_alpha1",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha1, linetype = "solid", size = 1)
p

p <- mcmc_hist(stan_out, pars = c("mu_alpha2"))
p <- p + labs(x = "mu_alpha2",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha2, linetype = "solid", size = 1)
p

# q: how good is our model at estimating the true abundance?
(q <- mcmc_hist(stan_out, pars = c("totalN[1]", "totalN[2]", 
                                       "totalN[3]", "totalN[4]"))
)

# plot posterior distribution
q <- mcmc_hist(stan_out, pars = c("totalN[5]"))
q <- q + labs(x = "totalN[5]",
              y = "Frequency in 1000 Draws")
  # xlim(nobs_in_sim, 150) +
q

q <- mcmc_hist(stan_out, pars = c("totalN[2]"))
q <- q + labs(x = "totalN[2]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = my_simulated_data$totalN[2], linetype = "solid", size = 1)
q

par(mfrow = c(1, 1))


