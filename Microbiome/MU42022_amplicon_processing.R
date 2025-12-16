# https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html

#Packages ----
setwd("/Users/greent/Desktop/Amplicon_Seq_data_Aug_2023/Marissa_Seqs")

# Clear workspace 
# rm(list=ls())

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
BiocManager::install("phangorn")
BiocManager::install("dada2")
BiocManager::install("BiocStyle")
install.packages("devtools")
install.packages("seqinr")
install.packages("Biostrings")

#Load 
library("dada2")
library("seqinr")
library("knitr")
library("BiocStyle")
library("devtools")
library("seqinr")
library("phangorn")
library("Biostrings")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

#Set seed ----
set.seed(100)

#Load seqs ----
####Tell R where the data is...
miseq_path <- "R1R2_SMK_2024/"
list.files(miseq_path)

# Sort ensures forward/reverse reads are in same order. notice the pattern (two different reads, Forward and Reverse)
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
#fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
#fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)

#Did it work? 
fnFs[1:3]
fnRs[1:3]

#Seq quality ----
#Quality of reads: Most Illumina sequencing data shows a trend of decreasing average quality towards the end of sequencing reads.
plotQualityProfile(fnFs[1:10])
plotQualityProfile(fnRs[50:59])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

length(fnFs)
length(fnRs)
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250, 180),
                     maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft = 10,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Plot quality again
plotQualityProfile(filtFs[148:149])
plotQualityProfile(filtRs[70:80])

#Make ASVs ----
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

#If getting error "files do not exist" from above, run this code:
file_check <- file.exists(filtFs) 
false_indices <- which(!file_check)
missing_files <- filtFs[false_indices]
print(missing_files)

#learn error rates - long step
errF <- learnErrors(filtFs, multithread=TRUE)

#Check ASVs - if many unknowns, you may want to change your trimming parameters (make them larger) to provide more info
#If using just forward reads
mergers <- mergePairs(dadaFs, derepFs) 

seqtabAll <- makeSequenceTable(dadaFs[!grepl("Mock", names(dadaFs))])
table(nchar(getSequences(seqtabAll)))

seqtabNoC <- removeBimeraDenovo(seqtabAll)

#ASVs to Taxa info ----
# now replace the long ASV names (the actual sequences) with human-readable names
#save the new names and sequences as a .fasta file in your project working directory, and save a table that shows the mapping of sequences to new ASV names
my_otu_table <- t(as.data.frame(seqtabNoC)) #transposed (OTUs are rows) data frame. unclassing the otu_table() output avoids type/class errors later on
ASV.seq <- as.character(unclass(row.names(my_otu_table))) #store sequences in character vector
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') #create new names

write.table(cbind(ASV.num, ASV.seq), "sequence_ASVname_mapping_SMK.txt", sep="\t", quote=F, row.names=F, col.names=F)

write.fasta(sequences=as.list(ASV.seq), names=ASV.num, "16s_ASV_sequences_all_SMK.fasta") #save sequences with new names in fasta format

#Phylogenetic tree ----

###Try 
Object1<- cbind(ASV.num, taxTab)

#IMPORTANT: sanity checks
colnames(seqtabNoC) == ASV.seq #only proceed if this tests as true for all elements -true
row.names(taxTab) == ASV.seq #only proceed if this tests as true for all elements -true

#rename your ASVs in the taxonomy table and sequence table objects
colnames(seqtabNoC) <- ASV.num
row.names(taxTab) <- ASV.num

#re-save sequence and taxonomy tables with updated names
write.table(data.frame("row_names"=rownames(seqtabNoC),seqtabNoC),"sequence_table.16s.all_merged_SMK.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(taxTab),taxTab),"taxonomy_table.16s_all_merged_SMK.txt", row.names=FALSE, quote=F, sep="\t")


#### Phylogenetic tree 
Fast1<- readDNAStringSet(file= "./16s_ASV_sequences_SMK.fasta",format = "fasta")

seqs <- getSequences(Fast1)
names(Fast1) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

plot(fitGTR)

#RDS phyloseq object ----
meta<-import_qiime_sample_data("Meta_all.txt")
ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), tax_table(taxTab))
ps1 <-merge_phyloseq(ps,meta)
ps1

