# plot intercepts and effects of habitat on abundance and detection

library(tidyverse)
library(rstan)

multimix <- as.data.frame(readRDS("./model_outputs/real_data/multimix.rds"))
binmix <- as.data.frame(readRDS("./model_outputs/real_data/binmix.rds"))
glm <- as.data.frame(readRDS("./model_outputs/real_data/GLM.rds"))

n_species <- length(8)

## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
# mu_alpha0, mu_alpha1, mu_beta0, mu_beta1
X_eco <- seq(1:10) # 10 ecological params of interest
# mean of eco params
Y_eco <- c(quantile(multimix$mu_alpha0, 0.5), # mu_alpha0
           quantile(binmix$mu_alpha0, 0.5), # mu_alpha0
           quantile(glm$mu_alpha0, 0.5), # mu_alpha0
           quantile(multimix$mu_alpha1, 0.5), # mu_alpha1
           quantile(binmix$mu_alpha1, 0.5), # mu_alpha1
           quantile(glm$mu_alpha1, 0.5), # mu_alpha1
           quantile(multimix$mu_beta0, 0.5), # mu_beta0
           quantile(binmix$mu_beta0, 0.5), # mu_beta0
           quantile(multimix$mu_beta1, 0.5), # mu_beta1
           quantile(binmix$mu_beta1, 0.5) # mu_beta1
)

# confidence intervals
lower_90_eco <- c(quantile(multimix$mu_alpha0, 0.05), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.05), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.05), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.05), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.05), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.05), # mu_alpha1
                  quantile(multimix$mu_beta0, 0.05), # mu_beta0
                  quantile(binmix$mu_beta0, 0.05), # mu_beta0
                  quantile(multimix$mu_beta1, 0.05), # mu_beta1
                  quantile(binmix$mu_beta1, 0.05) # mu_beta1
)

upper_90_eco <- c(quantile(multimix$mu_alpha0, 0.95), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.95), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.95), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.95), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.95), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.95), # mu_alpha1
                  quantile(multimix$mu_beta0, 0.95), # mu_beta0
                  quantile(binmix$mu_beta0, 0.95), # mu_beta0
                  quantile(multimix$mu_beta1, 0.95), # mu_beta1
                  quantile(binmix$mu_beta1, 0.95) # mu_beta1
)

# confidence intervals
lower_50_eco <- c(quantile(multimix$mu_alpha0, 0.25), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.25), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.25), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.25), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.25), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.25), # mu_alpha1
                  quantile(multimix$mu_beta0, 0.25), # mu_beta0
                  quantile(binmix$mu_beta0, 0.25), # mu_beta0
                  quantile(multimix$mu_beta1, 0.25), # mu_beta1
                  quantile(binmix$mu_beta1, 0.25) # mu_beta1
)

upper_50_eco <- c(quantile(multimix$mu_alpha0, 0.75), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.75), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.75), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.75), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.75), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.75), # mu_alpha1
                  quantile(multimix$mu_beta0, 0.75), # mu_beta0
                  quantile(binmix$mu_beta0, 0.75), # mu_beta0
                  quantile(multimix$mu_beta1, 0.75), # mu_beta1
                  quantile(binmix$mu_beta1, 0.75) # mu_beta1
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_90_eco, upper_90_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Draw ecological parameter plot

(ggplot <- ggplot(df_estimates_eco) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = seq(1:10),
                    labels=c(bquote(alpha[0]),
                             bquote(alpha[0]),
                             bquote(alpha[0]),
                             bquote(alpha[1]),
                             bquote(alpha[1]),
                             bquote(alpha[1]),
                             bquote(beta[0]),
                             bquote(beta[0]),
                             bquote(beta[1]),
                             bquote(beta[1])
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                      limits = c(-5, 5)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
   theme(legend.text=element_text(size=10),
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 24, angle=45, vjust=-0.5),
         axis.title.x = element_text(size = 18),
         axis.title.y = element_text(size = 18),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
   coord_flip()
)

ggplot <- ggplot +
  geom_errorbar(aes(x=X_eco, ymin=lower_90_eco, ymax=upper_90_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 


ggplot

## --------------------------------------------------
## Draw ecological parameter plot by param

# alpha0
temp <- df_estimates_eco[1:3,]

(p <- ggplot(temp) +
   theme_bw() +
   # scale_color_viridis(discrete=TRUE) +
   scale_x_discrete(name="", breaks = seq(1:3),
                    labels=c("multimix",
                             "binmix",
                             "GLM"
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (log-scaled)", width = 30),
                      limits = c(-2, 5)) +
   guides(color = guide_legend(title = "")) +
   geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(alpha[0])) +
    theme(plot.title = element_text(size = 32, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip() 
)

p <- p +
  geom_errorbar(aes(x=X_eco, ymin=lower_90_eco, ymax=upper_90_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 
p

# alpha1
temp <- df_estimates_eco[4:6,]

(q <- ggplot(temp) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = (seq(1:3)+3),
                     labels=c("multimix",
                              "binmix",
                              "GLM"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (log-scaled)", width = 30),
                       limits = c(-1, 3)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(mu[alpha[4]])) +
    theme(plot.title = element_text(size = 32, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip() 
)

q <- q +
  geom_errorbar(aes(x=X_eco, ymin=lower_90_eco, ymax=upper_90_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 
q

# beta0
temp <- df_estimates_eco[7:8,]

(r <- ggplot(temp) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = (seq(1:2)+6),
                     labels=c("multimix",
                              "binmix"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-4, 1)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(beta[0])) +
    theme(plot.title = element_text(size = 32, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip() 
)

r <- r +
  geom_errorbar(aes(x=X_eco, ymin=lower_90_eco, ymax=upper_90_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 
r

# beta1
temp <- df_estimates_eco[9:10,]

(s <- ggplot(temp) +
    theme_bw() +
    # scale_color_viridis(discrete=TRUE) +
    scale_x_discrete(name="", breaks = (seq(1:2)+8),
                     labels=c("multimix",
                              "binmix"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-1, 2)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(mu[beta[4]])) +
    theme(plot.title = element_text(size = 32, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    coord_flip() 
)

s <- s +
  geom_errorbar(aes(x=X_eco, ymin=lower_90_eco, ymax=upper_90_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 
s

cowplot::plot_grid(p,q,r,s, ncol = 2, 
                   labels = c('a)', 'b)', 'c)', 'd)'), 
                   label_size = 20)

## --------------------------------------------------
## compare overlap between distributions

glm <- as.data.frame(glm) 
binmix <- as.data.frame(binmix)  
multimix <- as.data.frame(multimix) 

quantile(binmix$mu_beta1, c(0.25, 0.5, 0.75))

SS## is estimate larger?
mean(binmix$mu_alpha1 > glm$mu_alpha1)
mean(multimix$mu_alpha1 > glm$mu_alpha1)

ilogit <- function(x) exp(x)/(1+exp(x))
quantile(multimix$mu_beta1, c(0.05, 0.5, 0.95))
quantile(multimix$mu_beta0, c(0.05, 0.5, 0.95))
detection_restored <- multimix$mu_beta0 + multimix$mu_beta1
ilogit(quantile(detection_restored, c(0.05, 0.5, 0.95)))
ilogit(quantile(multimix$mu_beta0, c(0.05, 0.5, 0.95)))

mean(detection_restored > multimix$mu_beta0)
