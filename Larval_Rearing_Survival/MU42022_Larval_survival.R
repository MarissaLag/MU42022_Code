#Load packages ----
#library
#install.packages("dplyr")
library(dplyr)
#install.packages("datarium")
library("datarium")
#install.packages("tidyverse")
library(tidyverse)
#install.packages("ggpubr")
library(ggpubr)
library(rstatix)
#install.packages("ggResidpanel")
library(ggResidpanel)
#install.packages("DHARMa")
library(DHARMa)
#install.packages("lme4")
library(lme4)
#install.packages("fitdistrplus")
library(fitdistrplus)
library(ggplot2)
#install.packages("hrbrthemes")
library(hrbrthemes)
library(tidyr)
#install.packages("viridis")
library(viridis)
library(car)
#install.packages("agricolae")
library(agricolae)
#install.packages("mgcv")
library(mgcv)
#install.packages("glmmTMB")
library(glmmTMB)
#install.packages("mgcViz")
library(mgcViz)
#install.packages("gamm4")
library(gamm4)



#set theme
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

#Experimental notes ----
#Notes about MU42022 data: Repeated measures (tank sampled repeatedly), 
#issues with colinearity between factors (Tank, Treatment/Family) due to experimental set-up and
#lack of true tank replication - i.e., tank cannot be used as random factor - family could. 

set.seed(765)

#Load Data and Fix ----
MU42022_counts_R_USE <- read.delim2("~/Documents/MSc/MU42022 larval counts/MU42022_counts_R_USE.txt")
Data <- MU42022_counts_R_USE
str(Data)

#Data filtering steps****
Data <- Data %>%
  filter(!Treatment %in% c("Antibiotics", "Antibiotics + HT"))

#Since tanks may be stocked differently, remove day 0
Data <- Data %>%
  filter(!Day %in% c("0"))

#Fix data structure****
Data$Day <- as.numeric(Data$Day)
Data$Treatment <- as.factor(Data$Treatment)
Data$Tank <- as.factor(Data$Tank)
Data$Family <- as.factor(Data$Family)
# Remove commas and convert to numeric
Data$Count <- as.numeric(gsub(",", "", Data$Count))
str(Data)

#Set references for treatment and family****
Data$Treatment <- relevel(Data$Treatment, ref = "Control")
#Use Family 3 as baseline as 1) all tanks start with similar Numbers and 2) HT survived
Data$Family <- relevel(Data$Family, ref = "3")


#GAM interaction ----
#GAM with Treatment*Familyinteraction***
Data$Treatment <- factor(Data$Treatment, 
                         levels = c("Control", "Probiotics", "Probiotics + HT", "High temperature (HT)"))
model_gam <- gam(Count ~ s(Day, by=Treatment,k=4) + s(Day, by=Family,k=4) + Family*Treatment, data = Data)
summary(model_gam) 
plot(model_gam, pages = 1, shade = TRUE)
gam.check(model_gam) 
concurvity(model_gam, full = TRUE) #issue
AIC(model_gam)

#Dharma Model check fit****
simulation_output <- simulateResiduals(model_gam) 
plot(simulation_output)

#Visualize smooth interactions*****
#2d
vis.gam(model_gam, 
        type = 'response', 
        plot.type = 'contour', 
        color = "heat",         
        contour.col = "black")  
#3d 
vis.gam(model_gam, 
        type = 'response', 
        plot.type = 'persp', 
        color = "heat",         
        contour.col = "black",
        axes = TRUE,
        ticktype = "detailed",
        theta = 0, phi = 50
)

#Random factor GAM ----
#Make family random factor instead of fixed
#Issues with model fit when family added as random factor - fixed with tweedie dist'd

#Getting error in below: model.matrix.formula(form, data = if (is.list(data)) data[all.vars(reformulate(names(data))) %in% data must be a data.frame
Data <- na.omit(Data)
mod <- gam(Count ~ s(Day, by = Treatment, k = 4) + Treatment + s(Family, bs = "re"), 
           data = Data)
summary(mod)

#Using gamm4 package as mgcv with Family, bs = "re" function getting errors
#Errors in mgcv package fixed
Mod4 <- gamm4(Count ~ s(Day, by = Treatment, k = 4) + Treatment, 
              random = ~(1 | Family),  
              data = Data,
              family = gaussian(link = "identity"))

summary(Mod4$mer) 
summary(Mod4$gam)
plot(Mod4$gam, pages = 1, shade = TRUE)
plot(Mod4$mer)
gam.check(Mod4$gam)

#DHARma package
simulation_output <- simulateResiduals(Mod4$mer)
# Plot the residuals
plot(simulation_output)

#Raw data visualization
hist(Data$Count) #Data zero-inflated and right skewed

ggplot(Data, aes(x = Count)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Histogram of Count Data",
       x = "Count",
       y = "Frequency")


#Random factor tweedie dist'd ----
#Try random factor with new dist'd 

Data <- na.omit(Data)

#Use the tweedie distributin for heavily tailed data with some zeros: https://www.statisticshowto.com/tweedie-distribution/
# Mod4 <- gamm4(Count ~ s(Day, by = Treatment, k = 4) + Treatment, 
#               random = ~(1 | Family),  
#               data = Data,
#               family = tw) #gamm4 does not have tw distribution as option

#Watch - if getting data.fram error - unload all other packages (other packages have same gam fxn)
Data$Treatment <- factor(Data$Treatment, 
                         levels = c("Control", "Probiotics", "Probiotics + HT", "High temperature (HT)"))

Mod5 <- mgcv::gam(Count ~ s(Day, by = Treatment, k = 4) + Treatment + 
                    s(Family, bs = "re"),  # Random effect for Family
                  family = tw(), 
                  data = Data)

vis.gam(Mod5, 
        type = 'response', 
        plot.type = 'contour', 
        color = "heat",         
        contour.col = "black")  

#Use your brain to assess fit
qqnorm(residuals(Mod5))
plot(fitted(Mod5)~residuals(Mod5))

#Use Dharma
simulation_output <- simulateResiduals(Mod5)
plot(simulation_output)

#Model results
summary(Mod5)
plot(Mod5, pages = 1, shade = TRUE)
summary(Mod5)
gam.check(Mod5)


#GLMM ----
mod_lmm <- lmer(Count ~ Day * Treatment + (1 | Family), 
                data = Data)

summary(mod_lmm)
plot(mod_lmm)

qqnorm(residuals(mod_lmm))
plot(fitted(mod_lmm)~residuals(mod_lmm))


library(glmmTMB)

mod_glmm_tweedie <- glmmTMB(Count ~ Treatment + (1 | Family), 
                            data = Data, 
                            family = tweedie())
summary(mod_glmm_tweedie)

qqnorm(residuals(mod_glmm_tweedie))
plot(fitted(Mod5)~residuals(mod_glmm_tweedie))

#Use Dharma
simulation_output <- simulateResiduals(mod_glmm_tweedie)
plot(simulation_output)



#Plotting ----

#Plot predicted response with confidence intervals 
# Create a new data frame for predictions
new_data <- expand.grid(Day = seq(min(Data$Day), max(Data$Day), length.out = 100),
                        Treatment = levels(Data$Treatment),
                        Family = levels(Data$Family))

# Get predictions and standard errors
predictions <- predict(Mod5, newdata = new_data, se.fit = TRUE)

predictions <- predict(model_gam, newdata = new_data, se.fit = TRUE)

# Add predicted counts and confidence intervals to new_data
new_data$Count <- predictions$fit
new_data$Upper <- predictions$fit + 1.96 * predictions$se.fit  # Upper CI
new_data$Lower <- predictions$fit - 1.96 * predictions$se.fit  # Lower CI

#If using tweedie distribution - remove log transformation
predictions <- predict(Mod5, newdata = new_data, se.fit = TRUE)

# Transform predictions back to the original scale
new_data$fit <- exp(predictions$fit)
new_data$upper <- exp(predictions$fit + 1.96 * predictions$se.fit)  # Upper CI
new_data$lower <- exp(predictions$fit - 1.96 * predictions$se.fit)  # Lower CI



# Plot with confidence intervals
ggplot(new_data, aes(x = Day, y = fit, color = Treatment)) +
  geom_line(size = 0.5) +  # Predicted values
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Treatment), alpha = 0.2, color = NA) +
  labs(title = "",
       x = "Day",
       y = "Predicted Larval Survival") +
  scale_color_manual(values = c("Control" = "darkgrey", 
                                "High temperature (HT)" = "red",
                                "Probiotics" = "skyblue",
                                "Probiotics + HT" = "orange")) +
  scale_fill_manual(values = c("Control" = "darkgrey", 
                               "High temperature (HT)" = "red",
                               "Probiotics" = "skyblue",
                               "Probiotics + HT" = "orange")) +
  facet_wrap(~Treatment) +
  theme(legend.position = "none")




# Plot with separate panels for each treatment*Family
plot1 <- ggplot(Data, aes(x = Day, y = Count, color = Treatment, linetype = Family)) +
  geom_line() + 
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white"),  
        strip.text = element_text(size = 14, face = "bold")) + 
  facet_wrap(~ Treatment) 

# Summarize the data
summary_data <- Data %>%
  group_by(Day, Treatment, Family) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    se_count = sd(Count, na.rm = TRUE) / sqrt(n())
  )

#order 
summary_data$Family <- factor(summary_data$Family,
                              levels = c("1", "2", "3", "4"))

#rename treatments
summary_data <- summary_data %>%
  mutate(Treatment = fct_recode(Treatment,
                                "Roseobacter-enriched" = "Probiotics",
                                "Roseobacter-enriched + HT" = "Probiotics + HT"))


ggplot(summary_data, aes(x = Day, y = mean_count, color = Treatment, group = Treatment)) +
  geom_line(size = 1.5) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.3) +
  scale_color_manual(
    values = c(
      "Control" = "darkgrey",
      "Roseobacter-enriched" = "cornflowerblue",
      "Roseobacter-enriched + HT" = "orange",
      "High temperature (HT)" = "red"
    )
  ) +
  scale_x_continuous(breaks = seq(min(summary_data$Day), max(summary_data$Day), by = 2)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
  ) +
  labs(
    title = "",
    x = "Days post-fertilization",
    y = "Count",
    color = "Treatment"
  ) +
  facet_wrap(~Family, labeller = as_labeller(c(
    "1" = "Family 1",
    "2" = "Family 2",
    "3" = "Family 3",
    "4" = "Family 4"
  )))


#Or plot by treatment
summary_data_2 <- Data %>%
  group_by(Day, Treatment) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    se_count = sd(Count, na.rm = TRUE) / sqrt(n())
  )

ggplot(summary_data_2, aes(x = Day, y = mean_count, color = Treatment, group = Treatment)) +
  geom_line(size = 1.5) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.3) +
  scale_color_manual(
    values = c(
      "Control" = "darkgrey",
      "Probiotics" = "cornflowerblue",
      "Probiotics + HT" = "orange"
    )
  ) +
  scale_x_continuous(breaks = seq(min(summary_data$Day), max(summary_data$Day), by = 2)) +  # Increment x-axis by 2
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 11, face = "bold")
  ) +
  labs(
    title = "",
    x = "Day",
    y = "Count",
    color = "Treatment"
  )


#Or include tank info
summary_data_2 <- Data %>%
  group_by(Day, Treatment,Tank) %>%
  summarise(
    mean_count = mean(Count, na.rm = TRUE),
    se_count = sd(Count, na.rm = TRUE) / sqrt(n())
  )


# Create the plot
ggplot(summary_data_2, aes(x = Day, y = mean_count, color = Treatment, fill = Treatment, group = Treatment)) +
  geom_line(size = 1.5) +  # Line for mean counts
  geom_point(size = 3) +  # Points for mean counts
  geom_errorbar(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), width = 0.3) +  # Error bars
  geom_ribbon(aes(ymin = mean_count - se_count, ymax = mean_count + se_count), alpha = 0.2, color = NA) +  # Shaded area for treatment
  scale_color_manual(
    values = c(
      "Control" = "darkgrey",
      "Probiotics" = "cornflowerblue",
      "Probiotics + HT" = "orange"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "darkgrey",
      "Probiotics" = "cornflowerblue",
      "Probiotics + HT" = "orange"
    )
  ) +
  scale_x_continuous(breaks = seq(min(summary_data$Day), max(summary_data$Day), by = 2)) +  # Increment x-axis by 2
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),  # Remove major gridlines
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14, face = "bold")
  ) +
  labs(
    title = "",
    x = "Day",
    y = "Count",
    color = "Treatment",
    fill = "Treatment"
  ) +
  facet_wrap(~Tank)  # Create a facet for each tank