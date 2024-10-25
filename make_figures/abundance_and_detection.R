# plot intercepts and effects of habitat on abundance and detection

library(tidyverse)
library(rstan)

multinomial <- readRDS("./model_outputs/real_data/multinomial_Nmix.rds")
binomial <- readRDS("./model_outputs/real_data/binomial_Nmix.rds")
glm <- readRDS("./model_outputs/real_data/GLM.rds")

fit_summary1 <- rstan::summary(multinomial)
fit_summary2 <- rstan::summary(binomial)
fit_summary3 <- rstan::summary(glm)

View(cbind(1:nrow(fit_summary1$summary), fit_summary1$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary2$summary), fit_summary2$summary)) # View to see which row corresponds to the parameter of interest
View(cbind(1:nrow(fit_summary3$summary), fit_summary3$summary)) # View to see which row corresponds to the parameter of interest

n_species <- length(7)

## --------------------------------------------------
## Plot ecological paramter means and variation

# parameter means
# mu_alpha0, mu_alpha1, mu_beta0, mu_beta1
X_eco <- seq(1:10) # 10 ecological params of interest
# mean of eco params
Y_eco <- c(fit_summary1$summary[1,1], # mu_alpha0
           fit_summary2$summary[1,1], # mu_alpha0
           fit_summary3$summary[1,1], # mu_alpha0
           fit_summary1$summary[3,1], # mu_alpha1
           fit_summary2$summary[3,1], # mu_alpha1
           fit_summary3$summary[3,1], # mu_alpha1
           fit_summary1$summary[9,1], # mu_beta0
           fit_summary2$summary[8,1], # mu_beta0
           fit_summary1$summary[11,1], # mu_beta1
           fit_summary2$summary[10,1] # mu_beta1
)

# confidence intervals
lower_95_eco <- c(fit_summary1$summary[1,4], # mu_alpha0
                  fit_summary2$summary[1,4], # mu_alpha0
                  fit_summary3$summary[1,4], # mu_alpha0
                  fit_summary1$summary[3,4], # mu_alpha1
                  fit_summary2$summary[3,4], # mu_alpha1
                  fit_summary3$summary[3,4], # mu_alpha1
                  fit_summary1$summary[9,4], # mu_beta0
                  fit_summary2$summary[8,4], # mu_beta0
                  fit_summary1$summary[11,4], # mu_beta1
                  fit_summary2$summary[10,4] # mu_beta1
)

upper_95_eco <- c(fit_summary1$summary[1,8], # mu_alpha0
                  fit_summary2$summary[1,8], # mu_alpha0
                  fit_summary3$summary[1,8], # mu_alpha0
                  fit_summary1$summary[3,8], # mu_alpha1
                  fit_summary2$summary[3,8], # mu_alpha1
                  fit_summary3$summary[3,8], # mu_alpha1
                  fit_summary1$summary[9,8], # mu_beta0
                  fit_summary2$summary[8,8], # mu_beta0
                  fit_summary1$summary[11,8], # mu_beta1
                  fit_summary2$summary[10,8] # mu_beta1
)

# confidence intervals
lower_50_eco <- c(fit_summary1$summary[1,5], # mu_alpha0
                  fit_summary2$summary[1,5], # mu_alpha0
                  fit_summary3$summary[1,5], # mu_alpha0
                  fit_summary1$summary[3,5], # mu_alpha1
                  fit_summary2$summary[3,5], # mu_alpha1
                  fit_summary3$summary[3,5], # mu_alpha1
                  fit_summary1$summary[9,5], # mu_beta0
                  fit_summary2$summary[8,5], # mu_beta0
                  fit_summary1$summary[11,5], # mu_beta1
                  fit_summary2$summary[10,5] # mu_beta1
)

upper_50_eco <- c(fit_summary1$summary[1,7], # mu_alpha0
                  fit_summary2$summary[1,7], # mu_alpha0
                  fit_summary3$summary[1,7], # mu_alpha0
                  fit_summary1$summary[3,7], # mu_alpha1
                  fit_summary2$summary[3,7], # mu_alpha1
                  fit_summary3$summary[3,7], # mu_alpha1
                  fit_summary1$summary[9,7], # mu_beta0
                  fit_summary2$summary[8,7], # mu_beta0
                  fit_summary1$summary[11,7], # mu_beta1
                  fit_summary2$summary[10,7] # mu_beta1
)


df_estimates_eco <- as.data.frame(cbind(X_eco, Y_eco, 
                                        lower_95_eco, upper_95_eco,
                                        lower_50_eco, upper_50_eco))

df_estimates_eco$X_eco <- as.factor(df_estimates_eco$X_eco)

## --------------------------------------------------
## Get species specific estimates

species_estimates <- data.frame()

first_species_nat_green = 121

for(i in 1:n_species){
  
  # row is one before the row of the first species estimate
  #species_estimates[1,i] <- NA # psi species
  #species_estimates[1,i] <- fit_summary$summary[23+i,1] # psi species
  
  species_estimates[1,i] <- fit_summary$summary[(first_species_nat_green-1)+i,1] # herb shrub forest
  species_estimates[2,i] <- NA # dev green
  species_estimates[3,i] <- NA # income
  species_estimates[4,i] <- NA # race
  species_estimates[5,i] <- NA # site area
}

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
                      limits = c(-4, 4)) +
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

df_estimates_eco_species <- cbind(df_estimates_eco, species_estimates)

for(i in 1:n_species){
  
  test <- as.data.frame(cbind(X_eco, rev(df_estimates_eco_species[,6+i])))
  #test[1,2] <- NA
  #test[2,2] <- NA
  #test[4,2] <- NA
  colnames(test) <- c("X_eco", "Y_eco")
  
  s <- s + geom_point(data = test, aes(x=X_eco, y=Y_eco), 
                      col = "skyblue", size = 6, shape = "|", alpha = 0.75)
  
}

ggplot <- ggplot +
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
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
                    labels=c("mrc Nmix",
                             "binomial Nmix",
                             "GLM"
                    )) +
   scale_y_continuous(str_wrap("Posterior model estimate (log-scaled)", width = 30),
                      limits = c(-1, 4)) +
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
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
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
                     labels=c("mrc Nmix",
                              "binomial Nmix",
                              "GLM"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (log-scaled)", width = 30),
                       limits = c(-1, 4)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(alpha[1])) +
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
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
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
                     labels=c("mrc Nmix",
                              "binomial Nmix"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-2.5, 1)) +
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
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
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
                     labels=c("mrc Nmix",
                              "binomial Nmix"
                     )) +
    scale_y_continuous(str_wrap("Posterior model estimate (logit-scaled)", width = 30),
                       limits = c(-1, 2.5)) +
    guides(color = guide_legend(title = "")) +
    geom_hline(yintercept = 0, lty = "dashed") +
    ggtitle(bquote(beta[1])) +
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
  geom_errorbar(aes(x=X_eco, ymin=lower_95_eco, ymax=upper_95_eco),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=X_eco, ymin=lower_50_eco, ymax=upper_50_eco),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=X_eco, y=Y_eco),
             size = 5, alpha = 0.8) 
s

gridExtra::grid.arrange(p,q,r,s)
