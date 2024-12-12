# plot intercepts and effects of habitat on abundance and detection

library(tidyverse)
library(rstan)

# get species names
n_species <- 7

df <- read.csv("./data/abundance_data.csv") %>%
  filter(species != "Andrena prunorum")
species_names <- df$species
species_names_table <- as.data.frame(cbind(cbind(1:n_species), unique(species_names))) %>%
  rename("X" = "V1",
         "species_name" = "V2") %>%
  mutate(X = as.integer(X))

species_names <- as.vector(species_names_table$species_name)

## --------------------------------------------------
# get multi mix species estimates

multinomial <- readRDS("./model_outputs/real_data/multinomial_Nmix.rds")

fit_summary1 <- rstan::summary(multinomial)
#View(cbind(1:nrow(fit_summary1$summary), fit_summary1$summary)) # View to see which row corresponds to the parameter of interest

X <- as.integer(seq(1:n_species)) 
Y <- vector(length = length(X))
Y_low95 <- vector(length = length(X))
Y_up95 <- vector(length = length(X))
Y_low50 <- vector(length = length(X))
Y_up50 <- vector(length = length(X))

# parameter means
for(i in 1:n_species){
  Y[i] <- fit_summary1$summary[44+i,6] # beta0 species starts at 45
  Y_low95[i] <- fit_summary1$summary[44+i,4] # beta0 species starts at 45
  Y_up95[i] <- fit_summary1$summary[44+i,8] # beta0 species starts at 45
  Y_low50[i] <- fit_summary1$summary[44+i,5] # beta0 species starts at 45
  Y_up50[i] <- fit_summary1$summary[44+i,7] # beta0 species starts at 45
}

out <- as.data.frame(cbind(X, Y, Y_low95, Y_up95, Y_low50, Y_up50)) %>%
  left_join(species_names_table) %>%
  mutate(X = as.factor(X))

## --------------------------------------------------
## Draw species detection parameter plot

(p <- ggplot(out) +
   theme_bw() +
   geom_point(aes(x=X, y=Y),
              size = 5, alpha = 0.8) +
   geom_errorbar(aes(x=X_eco, ymin=Y_low95, ymax=Y_up95),
                 color="black",width=0.1,size=1,alpha=0.5) +
   geom_errorbar(aes(x=X_eco, ymin=Y_low50, ymax=Y_up50),
                 color="black",width=0,size=3,alpha=0.8) +
   scale_x_discrete(name="", breaks = seq(1:n_species),
                    labels=species_names) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3.5, 0)) +
   #guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)



## --------------------------------------------------
# get binom mix species estimates

binomial <- readRDS("./model_outputs/real_data/binomial_Nmix.rds")

fit_summary2 <- rstan::summary(binomial)
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest

X <- as.integer(seq(1:n_species)) 
Y <- vector(length = length(X))
Y_low95 <- vector(length = length(X))
Y_up95 <- vector(length = length(X))
Y_low50 <- vector(length = length(X))
Y_up50 <- vector(length = length(X))

# parameter means
for(i in 1:n_species){
  Y[i] <- fit_summary2$summary[44+i,6] # beta0 species starts at 45
  Y_low95[i] <- fit_summary2$summary[44+i,4] # beta0 species starts at 45
  Y_up95[i] <- fit_summary2$summary[44+i,8] # beta0 species starts at 45
  Y_low50[i] <- fit_summary2$summary[44+i,5] # beta0 species starts at 45
  Y_up50[i] <- fit_summary2$summary[44+i,7] # beta0 species starts at 45
}

out <- as.data.frame(cbind(X, Y, Y_low95, Y_up95, Y_low50, Y_up50)) %>%
  left_join(species_names_table) %>%
  mutate(X = as.factor(X))

## --------------------------------------------------
## Draw species detection parameter plot

(q <- ggplot(out) +
   theme_bw() +
   geom_point(aes(x=X, y=Y),
              size = 5, alpha = 0.8) +
   geom_errorbar(aes(x=X_eco, ymin=Y_low95, ymax=Y_up95),
                 color="black",width=0.1,size=1,alpha=0.5) +
   geom_errorbar(aes(x=X_eco, ymin=Y_low50, ymax=Y_up50),
                 color="black",width=0,size=3,alpha=0.8) +
   scale_x_discrete(name="", breaks = seq(1:n_species),
                    labels=species_names) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-3.5, 0)) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)


