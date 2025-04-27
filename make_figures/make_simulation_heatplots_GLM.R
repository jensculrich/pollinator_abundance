library(tidyverse)
library(gridExtra)

# plot mean estimates for mu_alpha1, with mu_beta1 = 0 versus mu_beta1 = 1 
# (for each of the three different models)

#-------------------------------------------------------------------------------
## changing baseline detection

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
list1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=0_beta1=-0.5.rds")
list2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=0_beta1=0.rds")
list3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=0_beta1=0.5.rds")
list4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=0_beta1=1.rds")
list5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-1_beta1=-0.5.rds")
list6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-1_beta1=0.rds")
list7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-1_beta1=0.5.rds")
list8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-1_beta1=1.rds")
list9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=-0.5.rds")
list10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=0.rds")
list11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
list12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=1.rds")
list13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-3_beta1=-0.5.rds")
list14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-3_beta1=0.rds")
list15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-3_beta1=0.5.rds")
list16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-3_beta1=1.rds")

# read data
n_datasets <- 16
list_dim <- 1 # which part of the list to access?
for(i in 1:n_datasets){
  name = paste0('list', as.character(i))
  x <- get(name)
  y <- x[[list_dim]]
  assign(  paste0("df", i), y)
  rm(name, x, y); gc()
}

# now calculate the means and BCI's across the range of sim situations tested
means <- vector(length=n_datasets)
lower90 <- vector(length=n_datasets)
upper90 <- vector(length=n_datasets)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}


df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

p <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,-0.5,0,0.5,1), limits = c(-1, 1)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab("Baseline detection rate (logit-scaled)") +
  ylab("Effect of habitat on\ndetection (logit-scaled)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
p


# precision

list_dim <- 2 # which part of the list to access?
for(i in 1:n_datasets){
  name = paste0('list', as.character(i))
  x <- get(name)
  y <- x[[list_dim]]
  assign(  paste0("df", i), y)
  rm(name, x, y); gc()
}

# now calculate the means and BCI's across the range of sim situations tested
means <- vector(length=n_datasets)
lower90 <- vector(length=n_datasets)
upper90 <- vector(length=n_datasets)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:n_datasets){
  precision <- get(paste0("df", i))
  means[i] <- mean(precision)
  lower90[i] <- quantile(precision, 0.05)
  upper90[i] <- quantile(precision, 0.95)
}


df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

p2 <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "ivory", high = "ivory3",, na.value = NA,
                      breaks=c(0,1), limits = c(0.5, 1.25)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab("Baseline detection rate (logit-scaled)") +
  ylab("Effect of habitat on\ndetection (logit-scaled)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"))
p2

#-------------------------------------------------------------------------------
## changing baseline detection

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
list1 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=-0.5_nsites=10.rds")
list2 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=0_nsites=10.rds")
list3 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=0.5_nsites=10.rds")
list4 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=1_nsites=10.rds")
list5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=-0.5.rds")
list6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=0.rds")
list7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
list8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/glm_alpha1=1_beta0=-2_beta1=1.rds")
list9 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=-0.5_nsites=30.rds")
list10 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=0_nsites=30.rds")
list11 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=0.5_nsites=30.rds")
list12 <- readRDS("./simulation_full/simulation_outputs/n_sites/glm_alpha1=1_beta0=-2_beta1=1_nsites=30.rds")

# read data
n_datasets <- 12
list_dim <- 1 # which part of the list to access?
for(i in 1:n_datasets){
  name = paste0('list', as.character(i))
  x <- get(name)
  y <- x[[list_dim]]
  assign(  paste0("df", i), y)
  rm(name, x, y); gc()
}

# now calculate the means and BCI's across the range of sim situations tested
means <- vector(length=n_datasets)
lower90 <- vector(length=n_datasets)
upper90 <- vector(length=n_datasets)
n_sites <- rep(c(10, 20, 30), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=3)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}


df <- as.data.frame(cbind(n_sites, beta1, means, lower90, upper90))

q <- ggplot(data=df, aes(x=n_sites, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,-0.5,0,0.5,1), limits = c(-1, 1)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = n_sites, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = n_sites, y = beta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab("n sites") +
  ylab("Effect of habitat on\ndetection (logit-scaled)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
q


# precision

list_dim <- 2 # which part of the list to access?
for(i in 1:n_datasets){
  name = paste0('list', as.character(i))
  x <- get(name)
  y <- x[[list_dim]]
  assign(  paste0("df", i), y)
  rm(name, x, y); gc()
}

# now calculate the means and BCI's across the range of sim situations tested
means <- vector(length=n_datasets)
lower90 <- vector(length=n_datasets)
upper90 <- vector(length=n_datasets)
n_sites <- rep(c(10, 20, 30), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=3)

for(i in 1:n_datasets){
  precision <- get(paste0("df", i))
  means[i] <- mean(precision)
  lower90[i] <- quantile(precision, 0.05)
  upper90[i] <- quantile(precision, 0.95)
}


df <- as.data.frame(cbind(n_sites, beta1, means, lower90, upper90))

q2 <- ggplot(data=df, aes(x=n_sites, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "ivory", high = "ivory3", na.value = NA,
                      breaks=c(0.5,1), limits = c(0.5, 1.25)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = n_sites, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = n_sites, y = beta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab("n sites") +
  ylab("Effect of habitat on\ndetection (logit-scaled)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"))
q2


#-------------------------------------------------------------------------------
# plot grid

cowplot::plot_grid(p, p2, q, q2, ncol=2,
                   labels = c("a)", "b)", "c)","d)"),
                   label_size = 20)

cowplot::plot_grid(p, q, ncol=2,
                   labels = c("a)", "b)"),
                   label_size = 20)

cowplot::plot_grid(p2, q2, ncol=2,
                   labels = c("a)", "b)"),
                   label_size = 20)
