library(tidyverse)
library(gridExtra)

# plot mean estimates for mu_alpha1, with mu_beta1 = 0 versus mu_beta1 = 1 
# (for each of the three different models)

n_datasets <- 12

#-------------------------------------------------------------------------------
## binmix

# read data
list1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=-1_theta1=-1.rds")
list2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=-1_theta1=0.rds")
list3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=-1_theta1=1.rds")
list4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=0_theta1=-1.rds")
list5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=0_theta1=0.rds")
list6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=0_theta1=1.rds")
list7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=1_theta1=-1.rds")
list8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=1_theta1=0.rds")
list9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=1_theta1=1.rds")
list10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=2_theta1=-1.rds")
list11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=2_theta1=0.rds") 
list12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/binmix_alpha1=1_theta0=2_theta1=1.rds") 

#-------------------------------------------------------------------------------
# bias abundance covariate

# get real value of effect
mu_alpha1 <- 1

# read data
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

# plot the results
q <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  ggtitle("abundance covariate bias (binmix)") +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=18)) 
q

#-------------------------------------------------------------------------------
# total abundance

# get real value of effect
mu_alpha0 <- 3

# read data
list_dim <- 3 # which part of the list to access?
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha0
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))


q2 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  ggtitle("abundance intercept bias (binmix)") +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=18)) 
q2

#-------------------------------------------------------------------------------
# detection rate

# get real value of effect
mu_beta0 <- -2

# read data
list_dim <- 4 # which part of the list to access?
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_beta0
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))


q3 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
q3

grid.arrange(q, q2, q3, ncol = 2)


#-------------------------------------------------------------------------------
## multimix

# read data
list1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=-1_theta1=-1.rds")
list2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=-1_theta1=0.rds")
list3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=-1_theta1=1.rds")
list4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=0_theta1=-1.rds")
list5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=0_theta1=0.rds")
list6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=0_theta1=1.rds")
list7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=1_theta1=-1.rds")
list8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=1_theta1=0.rds")
list9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=1_theta1=1.rds")
list10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=2_theta1=-1.rds")
list11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=2_theta1=0.rds") 
list12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/multimix_alpha1=1_theta0=2_theta1=1.rds") 

#-------------------------------------------------------------------------------
# bias abundance covariate

# get real value of effect
mu_alpha1 <- 1

# read data
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

#df[4:9, 3:5] <- NA # add NAs if we don't have all the data yet

# plot the results
q <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
q

#-------------------------------------------------------------------------------
# total abundance

# get real value of effect
mu_alpha0 <- 3

# read data
list_dim <- 3 # which part of the list to access?
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha0
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))


q2 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
q2

#-------------------------------------------------------------------------------
# detection rate

# get real value of effect
mu_beta0 <- -2

# read data
list_dim <- 4 # which part of the list to access?
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
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:n_datasets){
  temp <- get(paste0("df", i))
  bias <- temp - mu_beta0
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))


q3 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,0,1,2), limits = c(-1.5, 1.5)) +
  labs(fill="bias") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,2)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"),
        plot.title = element_text(face="bold", size=20)) 
q3

grid.arrange(q, q2, q3, ncol = 2)
#-------------------------------------------------------------------------------
## plot all

cowplot::plot_grid(q2, q, ncol=2,
             labels = c("a)", "b)"),
             label_size = 20)

