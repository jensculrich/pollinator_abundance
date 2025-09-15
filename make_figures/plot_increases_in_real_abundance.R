abundance_intercept <- c(2,2,2)

abundance_effect <- c(1,1,1)

bias <- c(0, 0.5, 1) 

df <- as.data.frame(cbind(abundance_intercept, abundance_effect, bias)) %>%
  mutate(expected_abundance_control = exp(abundance_intercept),
         expected_abundance_restored = exp(abundance_intercept + abundance_effect + bias)) %>%
  pivot_longer(cols = c(expected_abundance_control, expected_abundance_restored), 
               names_to = "treatment", values_to = "expected_abundance")

(p <- ggplot(df) +
    theme_bw() +
    geom_point(aes(x=as.factor(bias), y=expected_abundance, colour=treatment),
               size = 12, position=position_jitterdodge(jitter.width = 0,
                                                        jitter.height = 0,
                                                        dodge.width = 1)) +
    ggtitle("") +
    scale_y_continuous(name = "Expected Abundance", limits =  c(0, 75)) +
    scale_x_discrete(name = "Bias of Restoration Effect Estimate") +
    scale_color_manual( limits = c("expected_abundance_control", "expected_abundance_restored"),
                        labels = c("control", "restored"),
                        values = c("firebrick1", "skyblue2")) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

##------------------------------------------------------------------------------
# now with real estimates

multimix <- as.data.frame(readRDS("./model_outputs/real_data/multimix.rds"))
binmix <- as.data.frame(readRDS("./model_outputs/real_data/binmix.rds"))
glm <- as.data.frame(readRDS("./model_outputs/real_data/GLM.rds"))

# parameter means
# mu_alpha0, mu_alpha1, mu_beta0, mu_beta1
# mean of eco params
mean <- c(quantile(multimix$mu_alpha0, 0.5), # mu_alpha0
           quantile(binmix$mu_alpha0, 0.5), # mu_alpha0
           quantile(glm$mu_alpha0, 0.5), # mu_alpha0
           quantile(multimix$mu_alpha1, 0.5), # mu_alpha1
           quantile(binmix$mu_alpha1, 0.5), # mu_alpha1
           quantile(glm$mu_alpha1, 0.5) # mu_alpha1
)

# confidence intervals
lower_90 <- c(quantile(multimix$mu_alpha0, 0.05), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.05), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.05), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.05), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.05), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.05) # mu_alpha1
)

upper_90 <- c(quantile(multimix$mu_alpha0, 0.95), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.95), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.95), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.95), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.95), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.95) # mu_alpha1
)

# confidence intervals
lower_50 <- c(quantile(multimix$mu_alpha0, 0.25), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.25), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.25), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.25), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.25), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.25) # mu_alpha1
)

upper_50 <- c(quantile(multimix$mu_alpha0, 0.75), # mu_alpha0
                  quantile(binmix$mu_alpha0, 0.75), # mu_alpha0
                  quantile(glm$mu_alpha0, 0.75), # mu_alpha0
                  quantile(multimix$mu_alpha1, 0.75), # mu_alpha1
                  quantile(binmix$mu_alpha1, 0.75), # mu_alpha1
                  quantile(glm$mu_alpha1, 0.75) # mu_alpha1
)


df_estimates <- as.data.frame(cbind(mean, 
                                    lower_90, upper_90,
                                    lower_50, upper_50))

model <- c("multimix", "binmix", "GLMM", "multimix", "binmix", "GLMM") 
parameter <- c("intercept", "intercept", "intercept", "effect", "effect", "effect") 

df_estimates <- cbind(df_estimates, model, parameter)

df <- df_estimates %>%
  pivot_longer(cols = c(mean, 
                        lower_90, upper_90,
                        lower_50, upper_50), names_to = "estimate") %>%
  pivot_wider(names_from = parameter, values_from = "value") %>%
  mutate(control = exp(intercept),
         restored = exp(intercept + effect),
         relative_increase = restored / control) %>%
  select(-intercept, -effect, -control, -restored) %>%
  #pivot_longer(cols = c(control, restored), names_to = "treatment") %>%
  #pivot_wider(id_cols = c(model, treatment), names_from = estimate, values_from = "value") %>%
  pivot_wider(id_cols = model, names_from = estimate, values_from = "relative_increase")

#re-order factor levels for region
df$model <- factor(df$model, levels=c('GLMM', 'binmix', 'multimix'))

(q <- ggplot(df) +
    theme_bw() +
    ggtitle("") +
    scale_y_continuous(name = "Proportional Increase in Expected \nAbundance (Restored vs. Control)", 
                       limits =  c(0, 10)) +
    scale_x_discrete(name = "") +
    scale_color_manual(    values = c("firebrick1", "skyblue2")) +
    theme(plot.title = element_text(size = 20, face = "bold"),
          legend.text=element_text(size=10),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 20, angle=45, vjust=-0.5),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) 
)

q <- q +
  
  geom_errorbar(aes(x=as.factor(model), ymin=lower_90, ymax=upper_90),
                color="black",width=0.1,size=1,alpha=0.5) +
  geom_errorbar(aes(x=as.factor(model), ymin=lower_50, ymax=upper_50),
                color="black",width=0,size=3,alpha=0.8) +
  geom_point(aes(x=as.factor(model), y=mean),
             size = 8)

q

cowplot::plot_grid(q,p, ncol = 2, 
                   labels = c('a)', 'b)'), 
                   label_size = 20)

