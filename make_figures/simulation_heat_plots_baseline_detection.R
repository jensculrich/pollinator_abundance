library(tidyverse)
library(gridExtra)

# plot mean estimates for mu_alpha1, with mu_beta1 = 0 versus mu_beta1 = 1 
# (for each of the three different models)

#-------------------------------------------------------------------------------
## glm

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
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

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/glm_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
  temp <- get(paste0("df", i))
  means[i] <- median(temp)
  lower90[i] <- quantile(temp, 0.05)
  upper90[i] <- quantile(temp, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

p2 <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "#fbf3ff", high = "#b595c5", na.value = NA,
                      breaks=c(0,1,2), limits = c(0, 2)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = signif(means,4)),
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

grid.arrange(p, p2, ncol = 2)

#-------------------------------------------------------------------------------
## binmix

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/binmix_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

q <- ggplot(data=df, aes(x=beta0, y=beta1)) +
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
q

# precision

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/binmix_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
  temp <- get(paste0("df", i))
  means[i] <- median(temp)
  lower90[i] <- quantile(temp, 0.05)
  upper90[i] <- quantile(temp, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

q2 <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "#fbf3ff", high = "#b595c5", na.value = NA,
                      breaks=c(0,1,2), limits = c(0, 2)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = signif(means,4)),
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
q2

grid.arrange(q, q2, ncol = 2)


#-------------------------------------------------------------------------------
## multimix

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/multimix_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

r <- ggplot(data=df, aes(x=beta0, y=beta1)) +
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
r

# precision

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=0_beta1=-0.5.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=0_beta1=0.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=0_beta1=0.5.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=0_beta1=1.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-1_beta1=-0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-1_beta1=0.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-1_beta1=0.5.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-1_beta1=1.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-2_beta1=-0.5.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-2_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-2_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-2_beta1=1.rds")
df13 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-3_beta1=-0.5.rds")
df14 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-3_beta1=0.rds")
df15 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-3_beta1=0.5.rds")
df16 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/precision/multimix_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=16)
lower90 <- vector(length=16)
upper90 <- vector(length=16)
beta0 <- rep(c(0, -1, -2, -3), each=4)
beta1 <- rep(c(-0.5, 0, 0.5, 1), times=4)

for(i in 1:16){
  temp <- get(paste0("df", i))
  means[i] <- median(temp)
  lower90[i] <- quantile(temp, 0.05)
  upper90[i] <- quantile(temp, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

r2 <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient(low = "#fbf3ff", high = "#b595c5", na.value = NA,
                      breaks=c(0,1,2), limits = c(0, 2)) +
  labs(fill="precision") +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = paste0(
              #"90% BCI: [", signif(lower90,2), ", ", signif(upper90,2), "]")),
              "[", signif(lower90,2), ", ", signif(upper90,2), "]")),
            vjust = 1.5, size = 4.5) +
  geom_text(data = df, 
            aes(x = beta0, y = beta1, label = signif(means,4)),
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
r2

grid.arrange(r, r2, ncol = 2)

#-------------------------------------------------------------------------------
## plot all

cowplot::plot_grid(p, p2, q, q2, r, r2, ncol=2,
             labels = c("a)", "b)", "c)","d)", "e)", "f)"),
             label_size = 20)







#-------------------------------------------------------------------------------
## glm without labels (for presentation)

# bias

# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=0.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=0.5.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=0_beta1=1.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=0.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-1_beta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=1.rds")
df10 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=0.rds")
df11 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=0.5.rds")
df12 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-3_beta1=1.rds")

means <- vector(length=12)
lower90 <- vector(length=12)
upper90 <- vector(length=12)
beta0 <- rep(c(0, -1, -2, -3), each=3)
beta1 <- rep(c(0, 0.5, 1), times=4)

for(i in 1:12){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(beta0, beta1, means, lower90, upper90))

p_intercept <- ggplot(data=df, aes(x=beta0, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,-0.5,0,0.5,1), limits = c(0, 1)) +
  labs(fill="bias") +
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
p_intercept

#-------------------------------------------------------------------------------
## glm


# get real value of effect
mu_alpha1 <- 1

# read data
df1 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=0_nsites=10.rds")
df2 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=0.5_nsites=10.rds")
df3 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=1_nsites=10.rds")
df4 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.rds")
df5 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=0.5.rds")
df6 <- readRDS("./simulation_full/simulation_outputs/baseline_detection/estimates/glm_alpha1=1_beta0=-2_beta1=1.rds")
df7 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=0_nsites=30.rds")
df8 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=0.5_nsites=30.rds")
df9 <- readRDS("./simulation_full/simulation_outputs/n_sites/estimates/glm_alpha1=1_beta0=-2_beta1=1_nsites=30.rds")

means <- vector(length=9)
lower90 <- vector(length=9)
upper90 <- vector(length=9)
n_sites <- rep(c(10, 20, 30), each=3)
beta1 <- rep(c(0, 0.5, 1), times=3)

for(i in 1:9){
  temp <- get(paste0("df", i))
  bias <- temp - mu_alpha1
  means[i] <- mean(bias)
  lower90[i] <- quantile(bias, 0.05)
  upper90[i] <- quantile(bias, 0.95)
}

df <- as.data.frame(cbind(n_sites, beta1, means, lower90, upper90))

p_sites <- ggplot(data=df, aes(x=n_sites, y=beta1)) +
  geom_tile(aes(fill=means)) +
  scale_fill_gradient2(low = "#6bb4ff", mid = "white", high = "#ff6b6b", na.value = NA,
                       breaks=c(-1,-0.5,0,0.5,1), limits = c(0, 1)) +
  labs(fill="bias") +
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
p_sites

cowplot::plot_grid(p_intercept, p_sites, ncol=2)
