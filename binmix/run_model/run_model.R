# fit the Nmix model to my real data
library(tidyverse)

################################################################################
# load and format data

####################################
# analysis of real data ############
####################################

df <- read.csv("./data/abundance_data.csv")
str(df)

## Clean and prep data for model fitting
# select needed columns
df <- df %>%
  dplyr::select(site, treatment, species, year, count)

# add visit number
df <- df %>%
  group_by(site, species, year) %>%
  mutate(visit = row_number()) %>%
  ungroup()

# for now we will filter out the species that we only looked at in one year (A. prunorum)
#df <- filter(df, species != "Andrena prunorum")
df1 <- df %>%
  filter(year == 1) %>% # add a missing data variable
  mutate(missing_data = ifelse(species == "Andrena prunorum", 0, 1))
  
df2 <- df %>%
  filter(year == 2) %>% # add a missing data variable
  mutate(missing_data = 1)

df <- rbind(df1, df2)

n_species <- nrow(distinct(df, species)) # number of species
n_sites <- nrow(distinct(df, site)) # number of sites

df <- pivot_wider(df, names_from = visit, values_from = count)

y <- as.matrix(df[,5:7])
X <- as.vector(df$treatment)
R <- nrow(y)

species_names <- as.vector(df$species) # vector of species
site_names <- as.vector(df$site) # vector of sites
species <- as.integer(as.factor(species_names))
sites <- as.integer(as.factor(site_names))
year <- as.vector(df$year)

y_max <- apply(y, 1, max)
K <- ((y_max + 5) * 12)
View(as.data.frame(cbind(y_max, K)))

n_visits <- ncol(y)

missing_data <- as.vector(df$missing_data)

names <- rbind("Agapostemon texanus",
               "Andrena prunorum",
               "Anthidium oblongatum",
               "Bombus flavifrons", 
               "Bombus mixtus", 
               "Halictus rubicundus",
               "Megachile montivaga",
               "Melissodes microsticus")

species_names_table <- as.data.frame(cbind(cbind(1:8), names))
site_names_unique <- unique(site_names)
factor_sites <- as.integer(as.factor(site_names_unique))
site_names_table <- as.data.frame(factor_sites, site_names_unique)

##########################
### Run model ############
##########################

stan_data <- c("R", "K",
               "missing_data",
               "sites",
               "n_sites",
               "species",
               "year",
               "n_species", 
               "n_visits", 
               "y", "X")

# Parameters monitored
params <- c(
  "mu_alpha0",
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
  "scale_param",
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
n_burnin <- 0.5*n_iterations
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
       mu_beta0 = runif(1, -1, 1),
       sigma_beta0_species = runif(1, 0, 1),
       mu_beta1 = runif(1, -1, 1),
       sigma_beta1_species = runif(1, 0, 1),
       beta2 = runif(1, -1, 1)
  )
)

# Call STAN model from R 
stan_model <- "./binmix/models/binmix_negbin_NAs.stan"

## Call Stan from R
library(rstan)
stan_out <- stan(stan_model,
                     data = stan_data, 
                     init = inits, 
                     pars = params,
                     chains = n_chains, iter = n_iterations, 
                     warmup = n_burnin, thin = n_thin,
                     seed = 1,
                     open_progress = FALSE,
                     cores = n_cores)

saveRDS(stan_out, "./model_outputs/real_data/binmix.rds")

stan_out1 <- stan_out
stan_out2 <- readRDS("./model_outputs/real_data/binomial_Nmix.rds")
print(stan_out, digits = 3)

traceplot(stan_out, pars = c("mu_alpha0", "sigma_alpha0_species",
                             "mu_alpha1", "sigma_alpha1_species", 
                             #"mu_alpha2", "sigma_alpha2_species",
                             "sigma_alpha3_site"))


library(bayesplot)
library(tidyverse)

# p: how good is our model at estimating the parameters?
# plot posterior distribution
p <- mcmc_hist(stan_out, pars = c("mu_alpha0"))
p <- p + labs(x = "mu_alpha0",
              y = "Frequency in 1000 Draws") +
  geom_vline(xintercept = mu_alpha0, linetype = "solid", size = 1)
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
q <- mcmc_hist(stan_out, pars = c("totalN[1]"))
q <- q + labs(x = "totalN[1]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = my_simulated_data$totalN[1], linetype = "solid", size = 1)
q

q <- mcmc_hist(stan_out, pars = c("totalN[2]"))
q <- q + labs(x = "totalN[2]",
              y = "Frequency in 1000 Draws") +
  # xlim(nobs_in_sim, 150) +
  geom_vline(xintercept = my_simulated_data$totalN[2], linetype = "solid", size = 1)
q

# Evaluation of fit
list_of_draws <- as.data.frame(stan_out)

par(mfrow = c(1, 1))

plot(list_of_draws$fit, list_of_draws$fit_new, main = "", xlab =
       "Discrepancy actual data", ylab = "Discrepancy replicate data",
     frame.plot = FALSE,
     ylim = c(000, 50000),
     xlim = c(000, 50000))
abline(0, 1, lwd = 2, col = "black")

mean(list_of_draws$fit_new > list_of_draws$fit)
mean(list_of_draws$fit) / mean(list_of_draws$fit_new)
