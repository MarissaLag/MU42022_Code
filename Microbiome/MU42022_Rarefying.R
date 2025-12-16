#Rarefication for 16S rRNA sequences

#Load packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("microbiome")

install.packages("phyloseq")
install.packages("devtools")

library("devtools")
library(phyloseq)
library(microbiome)

#Load un-rarefied data ----
pseq<- readRDS("Marissa_.rds")

#Sample filtering
#Note: Analyzing algal samples seperately because very low prev data (few microbes associated with algae)
pseq <- subset_samples(pseq, !Sample.type %in% "Algae")

#create objects
OTU = pseq@otu_table
Tax = pseq@tax_table
Metadata = pseq@sam_data

#Sanity check
#check if any OTUs are not present in any samples (want false)
any(taxa_sums(pseq) == 0)

#if true

pseq_filtered <- prune_taxa(taxa_sums(pseq) > 0, pseq)
any(taxa_sums(pseq_filtered) == 0)

pseq <- pseq_filtered

#Remove chloroplast, mito, archaea, chimera ----

pseq <- subset_taxa(pseq,Class!="c__Chloroplast")

pseq <- subset_taxa(pseq,Order!="o__Mitochondria")
pseq <- subset_taxa(pseq,Family!="Mitochondria")

pseq <- subset_taxa(pseq,Kingdom!="k__Archaea")
pseq <- subset_taxa(pseq,Kingdom!="Archaea")

pseq <- subset_taxa(pseq,Order!="Chloroplast")

#Check if any chloro, ,mito, or achaeae
tax_levels <- colnames(pseq@tax_table)
chloroplast_archaea_mitochondria <- c("chloroplast", "archaea", "mitochondria")

any_contain <- sapply(chloroplast_archaea_mitochondria, function(term) {
  any(grepl(term, tolower(pseq@tax_table[, tax_levels]), ignore.case = TRUE))
})

any(any_contain)

#If true still contains bad stuff
#if false, move on
rank_names(pseq)

# Create table, number of features for each phyla
#Checking for NA phyla - these are often sequence artifacts and should be removed
table(tax_table(pseq)[, "Phylum"], exclude = NULL)

#If need to remove: 
#pseq <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#Remove low prev ----
plot(sort(taxa_sums(x2), TRUE), type="h", ylim=c(0, 10000))
# 

x1 = prune_taxa(taxa_sums(pseq) > 100, pseq) 
x2 = prune_taxa(taxa_sums(pseq) > 200, pseq)
x3 = prune_taxa(taxa_sums(pseq) > 800, pseq)

#used x2
summarize_phyloseq(x2)
summarize_phyloseq(x3)

pseq <- x2

#Rarefy read depth ----
##rarify data to make sequences an even read depth
Rare <-rarefy_even_depth(pseq, sample.size= 5128, rngseed = TRUE)

sample_depths <- sample_sums(Rare)

print(sample_depths)

pseq <- Rare

#Save rds
saveRDS(pseq, "/Users/maris/Documents/GitHub/Phyloseq and microbiome analysis/Old RDS files/MU42022_filtered_Oct92024.rds")


