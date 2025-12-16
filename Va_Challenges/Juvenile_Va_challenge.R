#Survival analyses
#Project: MU42022

#Packages ----

********Cox Analysis*******
  install.packages("survival")

install.packages("survminer")

install.packages("dplyr")

install.packages("ggplot2")

install.packages("ggpubr")

install.packages("coxme")

***
  
library(survival)

library(survminer)

library(dplyr)

library(ggplot2)

library(ggpubr)

library(emmeans)

library(coxme)

#Set theme ----
theme.marissa <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16, face = "bold"))
}

theme_set(theme.marissa())

#Import data ----
spat_challenge_2_results <- read_excel("~/Documents/MSc/vibrio spat challenge/spat challenge 2 results.xlsx", 
                                       +     sheet = "R_highdoseUSE (2)")

Data <- spat_challenge_2_results

#MU42022 renaming
Data$treatment[Data$treatment == "PH"] <- "Probiotics + HT"
Data$treatment[Data$treatment == "Probiotic"] <- "Probiotics"
Data$treatment[Data$treatment == "Heat"] <- "High temperature (HT)"
Data$Treatment <- Data$treatment

Data$Family <- Data$genetics

Data$Family[Data$Family == "1"] <- "1"
Data$Family[Data$Family == "2"] <- "2"
Data$Family[Data$Family == "3"] <- "3"


#filter data
Data <- Data %>% 
  filter(!Treatment %in% c("High temperature (HT)"))

#control tank of family 4 died - remove
Data <- Data %>% 
  filter(!Family %in% c("4"))

Data$TE <- as.numeric(Data$TE)

#Surv curve ----
#Surv function can only use fixed effects
set.seed(524)
surv_object = Surv(time=Data$TE, event=Data$Outcome)
fit1 = survfit(surv_object~ Treatment, data=Data)
summary(fit1)

#plot
ggsurvplot(fit1, data=Data, 
           pval = TRUE, 
           legend = c(0.8, 0.2),
           conf.int = FALSE,
           lwd = 2.5,
           legend.title = "Treatment", 
           legend.labs = c("Control", "Roseobacter-enriched", "Roseobacter-enriched + HT"),
           font.legend = c(13, "bold", "black"),  
           xlab = "Time (days)",
           ylab = "Survival Probability",
           font.x = c(14, "bold", "black"),        # X-axis title font
           font.y = c(14, "bold", "black"),
           linetype = c("solid", "dashed", "dotted"),
           censor.size = 10,
           palette = c("darkgrey", "cornflowerblue", "orange"))

#facet by genetics (family)

ggsurvplot_facet(
  fit1,
  data = Data,
  facet.by = "Family",
  pval = TRUE,
  lwd = 1.5,
  font.legend =c(12,"plain","black"),
  palette = c("darkgrey", "cornflowerblue", "orange"),
  conf.int = TRUE,                             
  risk.table = FALSE,
  legend.labs=c("Control", "Probiotics", "Probiotics + HT"),
) + 
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )


=#Mixed effects Cox model ----
#Note, cannot test for PH assumption with mixed effects model
#May have some limitations that I need to look into compared to fixed effects model

fit_mixed <- coxme(Surv(TE, Outcome) ~ treatment + (1|genetics), data = Data)

#compare model when fixed effect added or not

summary(fit_mixed)
#check proportional hazard assumption - passed
cox.zph(fit_mixed) 

fit_strat <- coxph(Surv(TE, Outcome) ~ treatment + strata(genetics), data=Data)
summary(fit_strat)

#For simplicity, using the fixed effects model for the results, as same statistical conclusions as
#the mixed model but does not assume guassian distribution of random effect
#And can be plotted easily

#Hazard ratio
fit.coxph =coxph(Surv(TE,Outcome)~ treatment,data=Data)
summary(fit.coxph)

plot <- ggforest(fit.coxph, data=Data, main = "Hazard ratio",
                 fontsize = 1.1,
                 noDigits = 2
)

plot

#add frailty effect
fit.coxph_frailty <- coxph(Surv(TE, Outcome) ~ treatment + frailty(genetics), data = Data)
summary(fit.coxph_frailty)

#Plot mixed effects hazard ratio

# ggforest works with coxph, so you might need to use your stratified model
fit_strat <- coxph(Surv(TE, Outcome) ~ treatment + strata(Family), data = Data)

# Create forest plot
ggforest(fit_strat, data = Data)

#Using ggplot instead

# Extract coefficients
coefs <- summary(fit_strat)$coefficients

# Create data frame
hr_data <- data.frame(
  Treatment = c("Control (Reference)", "Probiotics", "Probiotics + HT"),
  HR = c(1, coefs[1, "exp(coef)"], coefs[2, "exp(coef)"]),
  Lower = c(1, 
            exp(coefs[1, "coef"] - 1.96 * coefs[1, "se(coef)"]),
            exp(coefs[2, "coef"] - 1.96 * coefs[2, "se(coef)"])),
  Upper = c(1,
            exp(coefs[1, "coef"] + 1.96 * coefs[1, "se(coef)"]),
            exp(coefs[2, "coef"] + 1.96 * coefs[2, "se(coef)"])),
  P = c(NA, coefs[1, "Pr(>|z|)"], coefs[2, "Pr(>|z|)"])
)

# Reverse order and add colors
hr_data$Treatment <- factor(hr_data$Treatment, 
                            levels = rev(hr_data$Treatment))
hr_data$Color <- factor(c("darkgray", "orange", "cornflowerblue"),
                        levels = c("darkgray", "orange", "cornflowerblue"))

# Create plot
ggplot(hr_data, aes(y = Treatment)) +
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray", linewidth = 1) +
  # Confidence interval lines
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, x = HR), 
                 height = 0, linewidth = 1) +
  # Boxes for point estimates
  geom_tile(aes(x = HR, width = 0.04, height = 0.4, fill = Color),
            color = "black", linewidth = 1) +
  scale_fill_manual(values = c("darkgray" = "darkgray", 
                               "orange" = "orange", 
                               "cornflowerblue" = "cornflowerblue"),
                    guide = "none") +
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.7, 1.0, 1.5)) +
  scale_y_discrete(labels = c("Control (Reference)" = "Control (reference)",
                              "Probiotics" = "Roseobacter-enriched",
                              "Probiotics + HT" = "Roseobacter-enriched + HT")) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
  ) +
  theme(
    axis.text.y = element_text(size = 14, face = "bold")
  )

ggplot(hr_data, aes(y = Treatment)) +
  # Reference line
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray", linewidth = 1) +
  # Confidence interval lines
  geom_errorbarh(aes(xmin = Lower, xmax = Upper, x = HR), 
                 height = 0, linewidth = 1) +
  # Boxes for point estimates
  geom_tile(aes(x = HR, width = 0.04, height = 0.4, fill = Color),
            color = "black", linewidth = 1) +
  scale_fill_manual(values = c("darkgray" = "darkgray", 
                               "orange" = "orange", 
                               "cornflowerblue" = "cornflowerblue"),
                    guide = "none") +
  # Reverse log-scale x-axis
  scale_x_continuous(trans = "log", breaks = c(0.5, 0.7, 1.0, 1.5), 
                     limits = c(1.5, 0.5)) +
  scale_y_discrete(labels = c("Control (Reference)" = "Control (reference)",
                              "Probiotics" = "Roseobacter-enriched",
                              "Probiotics + HT" = "Roseobacter-enriched + HT")) +
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL
  ) +
  theme(
    axis.text.y = element_text(size = 14, face = "bold")
  )
  