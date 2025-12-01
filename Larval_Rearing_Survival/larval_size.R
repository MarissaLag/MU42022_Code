#Analyzing larval sizes (deep bay) MU42022

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

#Load data
Deep.Bay.Sizes <- read.delim("~/Documents/MSc/MU42022 larval counts/Deep Bay Sizes.txt")

Data <- Deep.Bay.Sizes
Data$Family <- as.factor(Data$Family)
Data$Treatment <- as.factor(Data$Treatment)
Data$Tank.. <- as.factor(Data$Tank..)
Data$Header <- as.factor(Data$Header)

# Remove NA values
Data_clean <- Data %>% filter(!is.na(Size))

# Calculate average Size per Treatment per Day
Data_summary <- Data_clean %>%
  group_by(Treatment, Day) %>%
  summarise(
    Avg_Size = mean(Size, na.rm = TRUE),
    SD = sd(Size, na.rm = TRUE),
    .groups = 'drop'
  )

#Filter out certain treatments
Data_summary <- Data_summary %>%
  filter(!Treatment %in% c("Antibiotics", "Antibiotics + HT"))

#Colors
custom_colors <- c(
  "Control" = "grey",
  "Probiotics" = "skyblue",
  "Probiotics + HT" = "orange",
  "High temperature (HT)" = "red"
)


ggplot(Data_summary, aes(x = Day, y = Avg_Size, color = Treatment, group = Treatment)) +
  geom_line(size = 1.5) +  # Custom line colors
  geom_point(size = 2.5) + # Custom point colors
  geom_errorbar(aes(ymin = Avg_Size - SD, ymax = Avg_Size + SD), 
                width = 0.2, alpha = 0.6) + # Error bars
  scale_color_manual(values = custom_colors) +  
  scale_x_continuous(breaks = unique(Data_summary$Day)) + 
  labs(title = "",
       x = "Day",
       y = "Average Size") +
  theme(
    legend.title = element_blank(),
    text = element_text(size = 14)
  )

#Box plot per day

library(ggplot2)

ggplot(Data_filtered, aes(x = Treatment, y = Size, fill = Treatment)) +
  geom_boxplot(outlier.size = 2, alpha = 0.6) +  # Box plot with custom outlier size
  scale_fill_manual(values = custom_colors) +  # Custom colors for the fill
  labs(title = "",
       x = "Treatment",
       y = "Larval shell width (microns)") +
  theme(
    axis.text.x = element_blank(),
    text = element_text(size = 14),
    axis.ticks = element_blank(),
    strip.text = element_text(size = 14, face = "bold")  # Make facet labels bigger and bold
  ) +
  facet_grid(~ Day, labeller = as_labeller(c(`3` = "Day 3", 
                                             `6` = "Day 6", 
                                             `8` = "Day 8", 
                                             `10` = "Day 10", 
                                             `13` = "Day 13"))) 


#GAM

Data_filtered <- Data_clean %>%
  filter(!Treatment %in% c("Antibiotics", "Antibiotics + HT"))

Data_filtered$Treatment <- droplevels(Data_filtered$Treatment)

new_Model <- mgcv::gam(Size ~ s(Day, by = Treatment, k = 3) + Treatment +
                    s(Family, bs = "re"),
                 family = tw(), 
                  data = Data_filtered)

summary(new_Model)
plot(new_Model)
gam.check(new_Model)
qqnorm(residuals(Model)) 
plot(fitted(Model)~residuals(Model))


#DHARma package
simulation_output <- simulateResiduals(new_Model)
# Plot the residuals
plot(simulation_output)

vis.gam(new_Model, 
        type = 'response', 
        plot.type = 'contour', 
        color = "heat",         
        contour.col = "black")  


library(mgcv)

Mod <- gamm(Size ~ s(Day, by = Treatment, k = 4) + Treatment + 
              s(Day, by = Family, k = 4) + Family,
            random = list(Header = ~1),  # Header as a random effect
            data = Data)


Mod_glmer <- glmer(Size ~ Day * Treatment * Family + (1 | Header),
             data = Data, 
             family = Gamma(link = "identity"))
summary(Mod_glmer)
plot(Mod_glmer)


model  <- lm(Size ~ Treatment*Family*Day, data = Data)

# Create a QQ plot of residuals
resid_panel(model)
#DHARma package
simulation_output <- simulateResiduals(model)
# Plot the residuals
plot(simulation_output)



Model <- lme(Size ~ Day * Treatment, 
             random = ~ 1 | Family,
             data = Data_clean)
Model <- lme(log(Size) ~ Day * Treatment, 
             random = ~ 1 | Family,
             data = Data_clean)

summary(Model)

resid_panel(Model)
#DHARma package
simulation_output <- simulateResiduals(Model)
# Plot the residuals
plot(simulation_output)



#Test each day

library(dplyr)
library(broom)
library(car)

Data_filtered$Day <- as.factor(Data_filtered$Day)

# Create a function to perform ANOVA for each day
anova_results <- Data_filtered %>%
  group_by(Day) %>%
  do(tidy(aov(Size ~ Treatment + Family, data = .)))

# Display the results
print(anova_results)

# Fit ANOVA model for Day 3 as an example
model_day3 <- aov(Size ~ Treatment + Family, data = Data_filtered %>% filter(Day == 3))
qqnorm(residuals(model_day3))
qqline(residuals(model_day3))

#check assumptions met each day
# Function to run Levene's Test and Shapiro-Wilk Test
run_tests <- function(data) {
  # Levene's Test for homogeneity of variances
  levene_result <- leveneTest(Size ~ Treatment, data = data)
  
  # Shapiro-Wilk Test for normality
  shapiro_result <- shapiro.test(data$Size)
  
  # Return a tidy summary
  return(tibble(
    levene_p_value = levene_result$`Pr(>F)`[1],  # Extract p-value from Levene's Test
    shapiro_p_value = shapiro_result$p.value  # Extract p-value from Shapiro-Wilk Test
  ))
}

# Run tests for each day and combine results
test_results <- Data_filtered %>%
  group_by(Day) %>%
  do(run_tests(.)) %>%
  ungroup()  # Remove grouping

# Display the results
print(test_results)

#Results - some days are not normally dist'd but all pass homo variance
#proceeding with kruskal wallis

# Function to run Kruskal-Wallis Test
run_kruskal_test <- function(data) {
  kruskal_result <- kruskal.test(Size ~ Treatment, data = data)
  
  # Return a tidy summary
  return(tibble(
    p_value = kruskal_result$p.value,  # Extract p-value from the test
    statistic = kruskal_result$statistic  # Extract statistic
  ))
}

# Run Kruskal-Wallis test for each day and combine results
kruskal_results <- Data_filtered %>%
  group_by(Day) %>%
  do(run_kruskal_test(.)) %>%
  ungroup()  # Remove grouping

# Display the results
print(kruskal_results)

#Only day 6 significant 
#run post hoc
library(dunn.test)
# Filter the data for Day 6
data_day6 <- Data_filtered %>% filter(Day == 6)
# Perform Dunn's Test
dunn_test_results <- dunn.test(
  x = data_day6$Size,  # The response variable
  g = data_day6$Treatment,  # The grouping factor
  kw = TRUE,  # Indicate that the Kruskal-Wallis test was performed
  label = TRUE,  # Show labels in output
  wrap = TRUE,  # Wrap long labels
  r = TRUE  # Return the rank of the data
)


# Load necessary libraries
library(lme4)

# Fit a GLMM with Gamma distribution
model <- glmer(Size ~ Treatment + (1 | Family), 
               data = Data_filtered %>% filter(Day == 6), 
               family = gaussian(link = "log"))

# View the summary of the model
summary(model)
resid_panel(model)
#DHARma package
simulation_output <- simulateResiduals(model)
# Plot the residuals
plot(simulation_output)

#Look at histogram
# Filter the data for Day 6
data_day6 <- Data_filtered[Data_filtered$Day == 6, ]

# Create a histogram of Size for Day 6 using base R
hist(data_day6$Size,
     breaks = 10,               
     col = "skyblue",       
     border = "black",        
     main = "Histogram of Size on Day 6",
     xlab = "Size",           
     ylab = "Frequency",        
     xlim = c(min(data_day6$Size), max(data_day6$Size)))


# Calculate mean size and standard deviation for each treatment on Day 6
mean_sd_day6 <- Data_filtered %>%
  filter(Day == 6) %>%
  group_by(Treatment) %>%
  summarize(
    mean_size = mean(Size, na.rm = TRUE),        # Calculate mean
    sd_size = sd(Size, na.rm = TRUE)             # Calculate standard deviation
  )
