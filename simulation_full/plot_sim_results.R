library(tidyverse)
library(gridExtra)

# plot mean estimates for mu_alpha1, with mu_beta1 = 0 versus mu_beta1 = 1 
# (for each of the three different models)

mu_alpha1 <- 1

# read data
estimates_glm <- readRDS("./simulation_full/simulation_outputs/estimates/glm_alpha1=1_beta0=0_beta1=-0.5.rds")
estimates_bmix <- readRDS("./simulation_full/simulation_outputs/estimates/binomial_nmix_alpha1=1_beta0=-3_beta1=0.rds")
estimates_mmix <- readRDS("./simulation_full/simulation_outputs/estimates/multimix_alpha1=1_beta0=-1_beta1=0.rds")

#quantile(estimates_mmix, c(0.05, 0.95))

df <- as.data.frame(as.vector(estimates_glm)) %>%
  rename("estimates_glm" = "as.vector(estimates_glm)") %>%
  mutate(dist = estimates_glm - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95)) 
df2 <- as.data.frame(as.vector(estimates_bmix)) %>%
  rename("estimates_bmix" = "as.vector(estimates_bmix)") %>%
  mutate(dist = estimates_bmix - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95))  
df3 <- as.data.frame(as.vector(estimates_mmix)) %>%
  rename("estimates_mmix" = "as.vector(estimates_mmix)") %>%
  mutate(dist = estimates_mmix - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95)) 

# Basic histogram
#text <- as.character(signif(df[1,3], 3))
(p <- ggplot(df, aes(x=estimates_glm)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

(q <- ggplot(df2, aes(x=estimates_bmix)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df2[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df2[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df2[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

(r <- ggplot(df3, aes(x=estimates_mmix)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df3[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df3[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df3[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

# read data
estimates_glm <- readRDS("./simulation_full/simulation_outputs/glm_alpha1=1_beta0=0_beta1=1.rds")
estimates_bmix <- readRDS("./simulation_full/simulation_outputs/binomial_nmix_alpha1=1_beta0=0_beta1=1.rds")
estimates_mmix <- readRDS("./simulation_full/simulation_outputs/multinomial_nmix_alpha1=1_beta0=0_beta1=1.rds")

df4 <- as.data.frame(as.vector(estimates_glm)) %>%
  rename("estimates_glm" = "as.vector(estimates_glm)") %>%
  mutate(dist = estimates_glm - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95)) 
df5 <- as.data.frame(as.vector(estimates_bmix))%>%
  rename("estimates_bmix" = "as.vector(estimates_bmix)") %>%
  mutate(dist = estimates_bmix - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95)) 
df6 <- as.data.frame(as.vector(estimates_mmix)) %>%
  rename("estimates_mmix" = "as.vector(estimates_mmix)") %>%
  mutate(dist = estimates_mmix - 1,
         mean_dist = mean(dist), 
         lower_90 = quantile(dist, 0.05),
         upper_90 = quantile(dist, 0.95)) 

# Basic histogram
(s <- ggplot(df4, aes(x=estimates_glm)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df4[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df4[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df4[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

(t <- ggplot(df5, aes(x=estimates_bmix)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df5[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df5[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df5[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

(u <- ggplot(df6, aes(x=estimates_mmix)) + geom_histogram(binwidth=0.1) +
    geom_vline(aes(xintercept=mean(mu_alpha1)),
               color="skyblue2", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df6[1,3])),
               color="black", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df6[1,4])),
               color="#F8766D", linetype="dashed", size=2) +
    geom_vline(aes(xintercept=mean(mu_alpha1 + df6[1,5])),
               color="#F8766D", linetype="dashed", size=2) +
    xlim(-0.5, 3) +
    theme_classic())

(u <- ggplot(df6, aes(x=estimates_mmix)) +
    xlim(-0.5, 3) +
    ylim(0, 500) +
    theme_classic())

grid.arrange(p, s, q, t, r, u, ncol = 2)

#-------------------------------------------------------------
