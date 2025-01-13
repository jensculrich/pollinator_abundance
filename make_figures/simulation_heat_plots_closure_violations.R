library(tidyverse)
library(gridExtra)

# plot mean estimates for mu_alpha1, with mu_beta1 = 0 versus mu_beta1 = 1 
# (for each of the three different models)

#-------------------------------------------------------------------------------
## binmix

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=-1_theta1=-1.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=-1_theta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=-1_theta1=1.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=0_theta1=-1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=0_theta1=0.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=0_theta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=1_theta1=-1.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=1_theta1=0.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=1_theta1=1.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=2_theta1=-1.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=2_theta1=-1.rds") # readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=2_theta1=0.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=2_theta1=-1.rds") # readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/binmix_alpha1=1_theta0=2_theta1=1.rds")

means <- vector(length=12)
lower90 <- vector(length=12)
upper90 <- vector(length=12)
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:12){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

df[11:12, 3:5] <- NA

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

# precision

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=-1_theta1=-1.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=-1_theta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=-1_theta1=1.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=0_theta1=-1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=0_theta1=0.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=0_theta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=1_theta1=-1.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=1_theta1=0.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=1_theta1=1.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=2_theta1=-1.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=2_theta1=0.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/binmix_alpha1=1_theta0=2_theta1=1.rds")

means <- vector(length=12)
lower90 <- vector(length=12)
upper90 <- vector(length=12)
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:12){
  temp <- get(paste0("df", i))
  means[i] <- median(temp)
  lower90[i] <- quantile(temp, 0.05)
  upper90[i] <- quantile(temp, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

q2 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "#fbf3ff", high = "#b595c5", na.value = NA,
                      breaks=c(0,1,2), limits = c(0, 2)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,4)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"))
q2

grid.arrange(q, q2, ncol = 2)


#-------------------------------------------------------------------------------
## multimix

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=-1_theta1=-1.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=-1_theta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=-1_theta1=1.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=0_theta1=-1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=0_theta1=0.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=0_theta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=1_theta1=-1.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=1_theta1=0.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=1_theta1=1.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=2_theta1=-1.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=2_theta1=0.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/estimates/multimix_alpha1=1_theta0=2_theta1=1.rds")

means <- vector(length=12)
lower90 <- vector(length=12)
upper90 <- vector(length=12)
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:12){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

r <- ggplot(data=df, aes(x=theta0, y=theta1)) +
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
r

# precision

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=-1_theta1=-1.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=-1_theta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=-1_theta1=1.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=0_theta1=-1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=0_theta1=0.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=0_theta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=1_theta1=-1.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=1_theta1=0.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=1_theta1=1.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=2_theta1=-1.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=2_theta1=0.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/closure_violations/precision/multimix_alpha1=1_theta0=2_theta1=1.rds")

means <- vector(length=12)
lower90 <- vector(length=12)
upper90 <- vector(length=12)
theta0 <- rep(c(-1, 0, 1, 2), each=3)
theta1 <- rep(c(-1, 0, 1), times=4)

for(i in 1:12){
  temp <- get(paste0("df", i))
  means[i] <- median(temp)
  lower90[i] <- quantile(temp, 0.05)
  upper90[i] <- quantile(temp, 0.95)
}

df <- as.data.frame(cbind(theta0, theta1, means, lower90, upper90))

r2 <- ggplot(data=df, aes(x=theta0, y=theta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "#fbf3ff", high = "#b595c5", na.value = NA,
                      breaks=c(0,1,2), limits = c(0, 2)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = theta0, y = theta1, label = signif(means,4)),
            vjust = -1.5, size = 4.5) +
  xlab(bquote(theta[0])) +
  ylab(bquote(theta[1])) +
  theme_classic() +
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14, face="bold"))
r2

grid.arrange(r, r2, ncol = 2)

#-------------------------------------------------------------------------------
## plot all

cowplot::plot_grid(q, r, ncol=2,
             labels = c("a)", "b)"),
             label_size = 20)

