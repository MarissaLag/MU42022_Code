#Alpha diveristy anlayses
#Sept 28th, 2023
#Project: MU42022

#source: https://microbiome.github.io/tutorials/PlotDiversity.html
#And code adapted from code from Andrew Loudon

#Packages ----

install.packages(c(
  "microbiome",
  "ggpubr",
  "knitr",
  "dplyr",
  "dunn.test",
  "gridExtra"
))

library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(phyloseq)
library(dunn.test)
library(gridExtra)

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

#Load data ----

setwd("C:/Users/maris/OneDrive/Documents/GitHub/mb2021_phyloseq")


Marissa_MU42022_rare <- readRDS("~/GitHub/mb2021_phyloseq/MU42022_filtered_Oct92024.rds")

pseq <- MU42022_filtered_Oct92024


#Filter data ----
#for treatment comparisons (below) need to remove "NA" values (algae)

pseq <- subset_samples(pseq, !Genetics %in% "4")
pseq <- subset_samples(pseq, !Organism %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")

#Select time point
#pseq<- subset_samples(pseq, Age %in% "Day 01")

#MU42022 - if want to test without supplemented bacteria
#taxa_to_remove <- c("ASV7", "ASV18")
#taxa_to_keep <- !(taxa_names(pseq) %in% taxa_to_remove)
#pseq <- prune_taxa(taxa_to_keep, pseq)


#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

#objects

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data

#Change chr to factor ----

Metadata = pseq@sam_data

Metadata$Treatment <- as.factor(Metadata$Treatment)

Metadata$Age <- as.factor(Metadata$Age)

Metadata$Genetics<- as.factor(Metadata$Genetics)

#alpha diversity ----

set.seed(678)

ps1 <- prune_taxa(taxa_sums(pseq) > 0, pseq)

ps1@sam_data$Genetics <- as.factor(ps1@sam_data$Genetics)
ps1@sam_data$Treatment <- as.factor(ps1@sam_data$Treatment)

tab <- microbiome::alpha(ps1, index = "all")
kable(head(tab))

#below code require "Treatment" to be a factor rather than a character
#check if your variable is a character or factor

str(Metadata)

#look for meta

ps1.meta <- meta(ps1)
kable(head(ps1.meta))

#combine meta and alpha diversity - add whatever diversity index

ps1.meta$Shannon <- tab$diversity_shannon 
ps1.meta$InverseSimpson <- tab$diversity_inverse_simpson
ps1.meta$chao1 <- tab$chao1

#test for normality ----
#Shapiro-wilk normality test

shapiro.test(ps1.meta$Shannon)
shapiro.test(ps1.meta$InverseSimpson)
shapiro.test(ps1.meta$chao1)

#for normal, run 3-way ANOVA

model <- aov(InverseSimpson ~Treatment, data=ps1.meta)

summary(model)

model <- aov(Shannon ~ Treatment, data=ps1.meta)

summary(model)

model <- aov(chao1 ~  Treatment, data=ps1.meta)

summary(model)

#run Tukey post-Hoc test to find where signif differences are

TukeyHSD(model, conf.level=.95)

#If giving error above, convert covariates to factors
ps1.meta$Treatment <- as.factor(ps1.meta$Treatment)
ps1.meta$Genetics <- as.factor(ps1.meta$Genetics)

#Age, Genetics, and Age*Genetics signif (p < 0.05)

#for non-normal - run Kruskal-Wallis test for 2 or more factors
#run Kruskal-Wallis test on chao1 and InvSimpson diversity

kruskal.test(Shannon ~ Genetics, data = ps1.meta)

kruskal.test(Shannon ~ Treatment, data = ps1.meta)

kruskal.test(Shannon ~ Age, data = ps1.meta)

kruskal.test(InverseSimpson ~ Treatment, data = ps1.meta)

kruskal.test(InverseSimpson ~ Age, data = ps1.meta)

kruskal.test(InverseSimpson ~ Genetics, data = ps1.meta)

kruskal.test(chao1 ~ Genetics, data = ps1.meta)

kruskal.test(chao1 ~ Age, data = ps1.meta)

kruskal.test(chao1 ~ Treatment, data = ps1.meta)

#On signif kruskal-wallis terms run Pairwise significance test as Posthoc test

res <- wilcox.test(chao1 ~ Treatment, data = ps1.meta,
                   exact = FALSE)

#More than 3 groups use Dunn's test
dunn.test(ps1.meta$chao1, ps1.meta$Age, method = "bonferroni")
dunn.test(ps1.meta$InverseSimpson, ps1.meta$Age, method = "bonferroni")
dunn.test(ps1.meta$Shannon, ps1.meta$Age, method = "bonferroni")

#Line plots of alpha diversity
# Summarize data: Mean and CI for Shannon diversity at each Age point
summary_data <- ps1.meta %>%
  group_by(Treatment) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    n = n(),
    ci_lower = mean_shannon - qt(0.975, n - 1) * (sd_shannon / sqrt(n)),
    ci_upper = mean_shannon + qt(0.975, n - 1) * (sd_shannon / sqrt(n))
  )

summary_data <- ps1.meta %>%
  group_by(Treatment) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    n = n(),
    ci_lower = mean_shannon - qt(0.975, n - 1) * (sd_shannon / sqrt(n)),
    ci_upper = mean_shannon + qt(0.975, n - 1) * (sd_shannon / sqrt(n))
  )

# Create the line plot
ggplot(summary_data, aes(x = Treatment, y = mean_shannon, group = 1)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "blue") +
  labs(
    title = "",
    x = "",
    y = "Mean Shannon Diversity"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Define custom colors for the treatments
custom_colors <- c("Control" = "darkgrey", 
                   "Probiotics" = "cornflowerblue", 
                   "Probiotics + HT" = "orange",
                   "High temperature (HT)" = "red")

custom_colors <- c("Control" = "darkgrey", 
                   "Probiotics" = "orange", 
                   "Probiotics + HT" = "forestgreen")

# Create the line graph
ggplot(summary_data, aes(x = Treatment, y = mean_shannon, color = Treatment, group = Treatment)) +
  geom_line(size = 1, linetype = "dashed") +  # Line for each treatment
  geom_point(size = 4) + # Points for each mean
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Treatment), 
              alpha = 0.2, color = NA) +  # Confidence intervals as ribbons
  labs(
    title = "",
    x = "",
    y = "Mean Shannon Diversity"
  ) +
  scale_color_manual(values = custom_colors) + 
  scale_fill_manual(values = custom_colors) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15, face = "bold"),
        axis.text.y = element_text(size=15, face = "bold"),
        axis.title.y = element_text(size=17))  

####Violin plots ----

p1 <- ggviolin(ps1.meta, 
               x = "Microbial.Source", 
               y = "Shannon", 
               add = "boxplot", 
               fill = "Microbial.Source", 
               ggtheme = theme_pubr()) + 
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange"))

#p1 <- p1 + scale_x_discrete(labels = c("Control", "Probiotics", "Probiotics + HT", "High temperature (HT)", "Algae"))

p1 <- p1 + theme(legend.position = "none") + 
  facet_wrap(~Age) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


print(p1)

# Create box plot
ps1.meta$Treatment <- factor(ps1.meta$Treatment, levels = c("Control", "High salinity", "Low salinity"))

color_palette <- c("grey", "#FF7F00", "#4DAF4A")

color_palette <- c("grey", "#FFB74D", "#A5D6A7")


color_palette <- c("#E41A1C", "#4DAF4A", "#377EB8")

color_palette <- c("#D85A5A", "#A8E0A8", "#A0C8E9")

color_palette <- c("#D85A5A", "#A8E0A8", "#5B9BD6")


p3 <- ggboxplot(ps1.meta, x = "Treatment", y = "InverseSimpson", 
                fill = "Treatment", 
                palette = color_palette,
                font.label = 25,
                outlier.shape = NA,
                size = 0.7) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank(),
        strip.text = element_text(size = 14, face = "bold")) +
  theme(axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 13, face = "bold"))


# Print the plot
print(p2)

p1 <- p1 + stat_compare_means(comparisons = Treatment.pairs) 
print(p1)

#align plots 

grid.arrange(p1, p2, p3, ncol = 3, nrow = 1)

#Kruskal test ----

Age <- levels(ps1.meta$Age)

Age <- levels(Metadata$Age.fact)
print(Age) 

Age.pairs <- combn(seq_along(Age), 2, simplify = FALSE, FUN = function(i)Age[i])
print(Age.pairs)

Treatment <- levels(ps1.meta$Treatment)

Treatment.pairs <- combn(seq_along(Treatment), 2, simplify = FALSE, FUN = function(i)Treatment[i])
print(Treatment.pairs)

####Violin plots ----

p1 <- ggviolin(ps1.meta, x = "Treatment", y = "InverseSimpson", add = "boxplot", 
               fill = "Age", 
               palette = c("#a6cee3", "#b2df8a", "#fdbf6f", "pink", "plum4"),
               ggtheme = theme_pubr(),
               font.label = 14)

p1 <- ggviolin(ps1.meta, x = "Treatment", y = "chao1", add = "boxplot", 
               fill = "Treatment", 
               palette = c("grey", "#fdbf6f","#a6cee3","#b2df8a", "pink", "plum4"),
               ggtheme = theme_pubr(),
               font.label = 14)

p1 <- p1 + scale_x_discrete(limits = c("Day 01", "Day 03", "Day 06", "Day 15", "Spat"))

p1 <- p1 + scale_x_discrete(limits = c("1 dpf", "18 dpf", "Spat"))

p1 <- p1 + facet_grid(~Age)

print(p1)

p1 <- p1 + stat_compare_means(comparisons = Treatment.pairs)

print(p1)


#other alpha diversity indices ----

tab <-microbiome::alpha(pseq, index = "all")
kable(head(tab))

tab <- richness(pseq)
kable(head(tab))

tab <- dominance(pseq, index = "all")
kable(head(tab))

dominant(pseq)

tab <- rarity(pseq, index = "all")
kable(head(tab))

tab <- coverage(pseq, threshold = 0.5)
kable(head(tab))

tab <- core_abundance(pseq, detection = .1/100, prevalence = 90/100)

tab <- inequality(pseq)


tab <- evenness(pseq, "all")
kable(head(tab))


#Visualization ----

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Treatment"
)

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Age"
)

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

p.shannon <- boxplot_alpha(pseq, 
                           index = "shannon",
                           x_var = "Sample.type"
)

p.shannon <- p.shannon + theme_classic() + 
  labs(x="", y="Shannon diversity") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
        legend.text = element_text(size=12),
        legend.title = element_text(size=16))

fill.colors = c(Control="cyan4",High.temperature = "deeppink4", Probiotics="forestgreen", Probiotics+HT= "slateblue2", NA="goldenrod1"))


palette = c("#a6cee3", "#b2df8a", "#fdbf6f", "pink", "plum4"),

p.shannon 


#Test significance (alpha diversity) ----

d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon
# Split the values by group
spl <- split(d$diversity, d$Treatment)
# Kolmogorov-Smironv test
pv <- ks.test(spl$Control, spl$"Probiotics + HT", spl$"High temperature", spl$"Probiotics", spl$"NA")$p.value

pv <- ks.test(spl$"Control", spl$"Low salinity", spl$"High salinity")$p.value


# Adjust the p-value
padj <- p.adjust(pv)

#p value = 0.419

d <- meta(pseq)
d$diversity <- microbiome::diversity(pseq, "shannon")$shannon
# Split the values by group
spl <- split(d$diversity, d$Treatment)
# Kolmogorov-Smironv test
pv <- ks.test(spl$"Probiotics", spl$"High temperature")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#p value = 0.064 * most signif difference between treatments (alpha diversity)


# Split the values by group
spl <- split(d$diversity, d$Sample.type)
# Kolmogorov-Smironv test
pv <- ks.test(spl$"Algae", spl$"Larvae", spl$"Spat")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#pvalue = 0.01268 (algae, larvae, spat)

pv <- ks.test(spl$"Larvae", spl$"Spat")$p.value
# Adjust the p-value
padj <- p.adjust(pv)

#p value = 7.05e-14 (larvae and spat)


#Day 1 only ----
#Since high temp treatment is missing spat samples = diversity skewed lower - look at day 1 diversity between treatments

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 01")

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "1 dpf")

p.shannon1 <- boxplot_alpha(pseq_Day1, 
                            index = "shannon",
                            x_var = "Treatment"
)

p.shannon1 <- p.shannon1 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 1") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon1 <- p.shannon1 + theme(legend.position = "none")

p.shannon1

#Day 3 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 03")

p.shannon3 <- boxplot_alpha(pseq_Day1, 
                            index = "shannon",
                            x_var = "Treatment"
)

p.shannon3 <- p.shannon3 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 3") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon3 <- p.shannon3 + theme(legend.position = "none")

p.shannon3

#Day 6 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Day 06")

p.shannon6 <- boxplot_alpha(pseq_Day1, 
                            index = "shannon",
                            x_var = "Treatment"
)

p.shannon6 <- p.shannon6 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 6") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon6 <- p.shannon6 + theme(legend.position = "none")

p.shannon6

#Day 18 only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day18 <- subset_samples(pseq_filtered, Age %in% "18 dpf")

p.shannon18 <- boxplot_alpha(pseq_Day18, 
                             index = "shannon",
                             x_var = "Treatment"
)

p.shannon18 <- p.shannon18 + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Day 15") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))

p.shannon18 <- p.shannon18 + theme(legend.position = "none")

p.shannon18

#Spat only ----

pseq_filtered <- subset_samples(pseq, !Sample.type %in% c("Algae"))

pseq_Day1 <- subset_samples(pseq_filtered, Age %in% "Spat")

p.shannons <- boxplot_alpha(pseq_filtered, 
                            index = "shannon",
                            x_var = "Treatment"
)

p.shannons <- p.shannons + theme_classic() + 
  labs(x="", y="Shannon diversity", title = "Spat") +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size = 12, hjust = 0.5))


p.shannons <- p.shannons + theme(legend.position = "none") +
  facet_wrap(~Age)

p.shannons

#Align plots into grid ----

install.packages("gridExtra")
library(gridExtra)

grid.arrange(p.shannon1, p.shannon3, p.shannon6, p.shannon15, p.shannons, ncol = 2)
