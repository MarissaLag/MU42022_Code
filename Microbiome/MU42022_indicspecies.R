##Code for running indicspecies analysis and Roseoabcter phylogeny
## Source:  https://jkzorz.github.io/2019/07/02/Indicator-species-analysis.html 
#Project: MU42022

#Packages ----
#install.packages("indicspecies")
library(indicspecies)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("microbiome")
library(microbiome)
#install.packages("phyloseq")
library(phyloseq)
#install.packages("RColorBrewer")
library(RColorBrewer)
#install.packages("tidyr")
library(tidyr)

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
MU42022_filtered_Oct92024 <- readRDS("~/GitHub/Old_RDS/MU42022_filtered_Oct92024.rds")
pseq <- MU42022_filtered_Oct92024

#Load objects ----

OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data
Tree = pseq@phy_tree

#Extract abundance matrix ----
#from the phyloseq object using phyloseq

OTU1 = as(OTU, "matrix")
write.csv(OTU1, file="Data_fram_1.cvs",row.names=TRUE)
write.table(OTU1,file="data_table_MU42022.csv",sep=",",dec = ".")

####Format to example data and reload below for actual test 

#reload edited table
pc_FUN <- data_table_MU42022

#Optional: filtering step if needed
pc_FUN <- pc_FUN[pc_FUN$Sample_type == "Spat", ] #filter for juvenile (Spat) samples only

#remove NAs if above not working
pc_FUN_clean <- pc_FUN %>%
  filter(!if_all(everything(), is.na))
pc_FUN <-pc_FUN_clean

#If present, filter NAs in column 1
pc_FUN <- pc_FUN[-1, ]

#Remove ASVs with zero abundance
metadata_cols <- 7  # Example: 7 metadata columns before ASVs

# Split metadata and ASV columns
meta <- pc_FUN[, 1:metadata_cols]
asvs <- pc_FUN[, (metadata_cols + 1):ncol(pc_FUN)]

# Remove ASVs with all zero values
asvs_filtered <- asvs[, colSums(asvs) > 0]

# Combine metadata and filtered ASVs
pc_FUN_filtered <- bind_cols(meta, asvs_filtered)
pc_FUN <- pc_FUN_filtered

#Test ASVs ----
#Inverse data
funi_df<- t(pc_FUN) 
dim(pc_FUN)

#Make matrix
matrix_F = pc_FUN[ ,8:586]

### Make the equation. Saying we want to examine specific column of metadata
time_a_F = pc_FUN$Treatment

### Run test 
set.seed(123)
inv_F = multipatt(matrix_F, time_a_F, func = "r.g", 
                  control = how(nperm=9999))
results <- summary(inv_F)

#Save results
write.csv(inv_F, "INVALsummary_results.csv", row.names = TRUE)

#Make horizontal box plots of signif ASVs ----
#Data
pseq <- MU42022_filtered_Oct92024

#Renaming Age
MU42022_filtered_Oct92024@sam_data$Age <- dplyr::recode(
  MU42022_filtered_Oct92024@sam_data$Age,
  "Day 01" = "1 dpf",
  "Day 03" = "3 dpf",
  "Day 06" = "6 dpf",
  "Day 15" = "15 dpf",
  "Spat"   = "90 dpf"
)

#Renaming Treatment
MU42022_filtered_Oct92024@sam_data$Treatment <- dplyr::recode(
  MU42022_filtered_Oct92024@sam_data$Treatment %>% as.character(),
  "Control" = "Control",
  "Probiotics" = "Roseobacter-enriched",
  "Probiotics + HT" = "Roseobacter-enriched + HT",
)

#Filter
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, Age %in% "90 dpf")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- microbiome::transform(pseq, "compositional")
ps <- psmelt(pseq)

# List of ASVs you're interested in
#MU42022 spat inval results - PB and PBH treatments shared
selected_asvs <- c("ASV88", "ASV178")
#Only PB treatment signif
selected_asvs <- c("ASV11", "ASV198", "ASV201", "ASV471", "ASV613")

# Subset the phyloseq object to keep only the selected ASVs
filtered_phyloseq <- prune_taxa(selected_asvs, pseq)

# Melt the phyloseq object to a long data frame
melted_phyloseq <- psmelt(filtered_phyloseq)

#Calculate log-fold change relative to the Control treatment
average_abundance <- melted_phyloseq %>%
  group_by(Treatment, OTU) %>%
  summarise(Average_Abundance = mean(Abundance)) %>%
  pivot_wider(names_from = Treatment, values_from = Average_Abundance, names_prefix = "Treatment_") %>%
  mutate(across(starts_with("Treatment_"), log, .names = "Log_{.col}")) %>%
  mutate(across(starts_with("Log_Treatment_"), 
                ~ . - Log_Treatment_Control, 
                .names = "Log_Fold_Change_{.col}")) %>%
  select(OTU, starts_with("Log_Fold_Change_"))


average_abundance_long <- average_abundance %>%
  pivot_longer(
    cols = starts_with("Log_Fold_Change_"),
    names_to = "Treatment",
    values_to = "Log_Fold_Change"
  ) %>%
  mutate(Treatment = gsub("Log_Fold_Change_Log_Treatment_", "", Treatment))

average_abundance_long <- average_abundance_long %>%
  filter(Treatment != "Control")



#Update names to include genus
# Update OTU names
average_abundance_long$OTU <- ifelse(
  average_abundance_long$OTU == "ASV88", "ASV88 - Loktanella",
  ifelse(average_abundance_long$OTU == "ASV178", "ASV178 - Loktanella", average_abundance_long$OTU)
)

average_abundance_long$OTU <- dplyr::recode(
  average_abundance_long$OTU,
  "ASV11"  = "ASV11 - unknown Roseobacter",
  "ASV198" = "ASV198 - unknown Roseobacter",
  "ASV201" = "ASV201 - Sulfitobacter",
  "ASV471" = "ASV471 - Sulfitobacter",
  "ASV613" = "ASV613 - unknown Roseobacter"
)


# Plotting
p1 <- ggplot(average_abundance_long, aes(x = Log_Fold_Change, y = OTU, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") + # Adds black borders for clarity
  labs(
    title = "",
    x = "",
    y = ""
  ) +
  scale_fill_manual(values = c("Roseobacter-enriched" = "cornflowerblue", "Roseobacter-enriched + HT" = "orange")) +
  theme(
    legend.position = "top",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

p2 <- ggplot(average_abundance_long, aes(x = Log_Fold_Change, y = OTU, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.9), color = "black") + # Adds black borders for clarity
  labs(
    title = "",
    x = "Log-fold Change in Relative Abundance",
    y = ""
  ) +
  scale_fill_manual(values = c("Roseobacter-enriched" = "cornflowerblue", "Roseobacter-enriched + HT" = "orange")) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14)
  )

p1 / p2

#Plot supplemented Roseobacter (ASV7/18) overtime ----
pseq <- MU42022_filtered_Oct92024
pseq <- subset_samples(pseq, !Genetics %in% c("4"))
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")
pseq <- subset_samples(pseq, !Treatment %in% "High temperature")
pseq <- microbiome::transform(pseq, "compositional")

asv_ids <- c("ASV7", "ASV18")
pseq_filtered <- prune_taxa(taxa_names(pseq) %in% asv_ids, pseq)
pseq_filtered <- psmelt(pseq_filtered)

#Average abundace for each treatmentxage group
pseq_avg <- pseq_filtered %>%
  group_by(Genus, Treatment, Age) %>%
  summarise(Avg_Abundance = mean(Abundance), .groups = 'drop')

ggplot(pseq_avg, aes(x = Age, y = Avg_Abundance, fill = Treatment)) +
  # Black outline line
  geom_line(aes(group = Treatment), color = "black", linewidth = 2.5) +
  # Colored line on top
  geom_line(aes(group = Treatment, color = Treatment), linewidth = 2) +
  # Points with black outline
  geom_point(size = 5, shape = 21, color = "black", stroke = 0.9, aes(fill = Treatment)) +
  scale_fill_manual(values = c(
    "Control" = "darkgrey",
    "Roseobacter-enriched" = "cornflowerblue",
    "Roseobacter-enriched + HT" = "orange"
  )) +
  scale_color_manual(values = c(
    "Control" = "darkgrey",
    "Roseobacter-enriched" = "cornflowerblue",
    "Roseobacter-enriched + HT" = "orange"
  )) +
  labs(x = "", y = "Celeribacter baekdonensis Relative Abundance") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.line = element_line(size = 0.9, color = "black"),
    axis.ticks = element_line(size = 0.9, color = "black")
  )


