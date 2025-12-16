#PCOa plots 
#Project: MU42022

##Packages ----
install.packages("devtools")
library(devtools)
install.packages("microbiome")
library(microbiome)
install.packages("ggalt")
library(ggalt)
install.packages("gridExtra")
library(gridExtra)
install.packages("ggpubr")
library(ggpubr)
install.packages("phyloseq")
library(phyloseq)
install.packages("ggplot2")
library(ggplot2)

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
MU42022_filtered_Oct92024 <- readRDS("~/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")

#Renaming Age
MU42022_filtered_Oct92024@sam_data$Age <- dplyr::recode(
  MU42022_filtered_Oct92024@sam_data$Age %>% as.character(),
  "Day 01" = "1 dpf",
  "Day 03" = "3 dpf",
  "Day 06" = "6 dpf",
  "Day 15" = "15 dpf",
  "Spat"   = "90 dpf"
)

MU42022_filtered_Oct92024@sam_data$Treatment <- dplyr::recode(
  MU42022_filtered_Oct92024@sam_data$Treatment %>% as.character(),
  "Control" = "Control",
  "Probiotics" = "Roseobacter-enriched",
  "Probiotics + HT" = "Roseobacter-enriched + HT",
)

#plot MDS/PcoA ----
#Create PCO for each time-point then combine plots **must be done for correct ordination
#note: different projects have different time-points
#note: if a treatment is lost at a certain time-point, adjust colours manually 

set.seed(4235421)

#MU42022 ----
#All time-points
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

#order by age
pseq.rel@sam_data$Age <- factor(
  pseq.rel@sam_data$Age,
  levels = c("1 dpf", "3 dpf", "6 dpf", "15 dpf", "90 dpf")
)

p_legend <- plot_ordination(pseq.rel, ord, color = "Treatment", shape = "Age") +
  geom_point(aes(fill = Treatment), size = 6) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")

p_legend

#1 dpf
#filtering
#important to note that the PB sequence was not removed
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature")) #Add/remove "red" if HT included
pseq <- subset_samples(pseq, Age %in% c("1 dpf"))
pseq.rel <- microbiome::transform(pseq, "compositional")

#ordination
ord <- ordinate(pseq.rel, "MDS", "bray")

#plot
p1 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  ggtitle("1 dpf") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
  theme(axis.line = element_line(size = 1, colour = "black")
  )
#geom_encircle(aes(fill = Treatment), expand = 0.2, alpha = 0.2)
p1

#3 dpf
#filtering
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("3 dpf"))
pseq.rel <- microbiome::transform(pseq, "compositional")
#ordination
ord <- ordinate(pseq.rel, "MDS", "bray")
#plot
p3 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, color = "black", stroke = 1) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  ggtitle("3 dpf") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
p3


#6 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("6 dpf"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p6 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, color = "black", stroke = 1) +
  scale_colour_manual(values = c("darkgrey",  "cornflowerblue", "orange")) +
  scale_fill_manual(values = c("darkgrey",  "cornflowerblue", "orange")) + # Matching colors for ellipses and points
  ggtitle("6 dpf") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

p6

#15 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% c("High temperature"))
pseq <- subset_samples(pseq, Age %in% c("15 dpf"))
pseq.rel <- microbiome::transform(pseq, "compositional")

ord <- ordinate(pseq.rel, "MDS", "bray")

p15 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  ggtitle("15 dpf") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

p15

#90 dpf
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% c("90 dpf"))
pseq.rel <- microbiome::transform(pseq, "compositional")
ord <- ordinate(pseq.rel, "MDS", "bray")

p90 <- plot_ordination(pseq.rel, ord, color = "Treatment") +
  geom_point(aes(fill = Treatment), size = 6, shape = 21, color = "black", stroke = 1) +
  scale_fill_manual(values = c("darkgrey", "cornflowerblue", "orange")) +
  ggtitle("90 dpf") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

p90

#arrange PCO plots -----

grid.arrange(p1, p3, p6, p15, p90, ncol = 3)
ggarrange(p1, p3, p6, p15, p90, common.legend = TRUE, legend="bottom") 
