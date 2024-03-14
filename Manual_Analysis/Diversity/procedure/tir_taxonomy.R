#code to run ANCOM
#load libraries
library(readr)
library(tidyverse)
library(ape)
library(nlme)
library(compositions)
library(phyloseq)
library(vegan)
library(BiocManager)
library(ANCOMBC)
library(ComplexHeatmap)
library(microbiome)
library(ggplot2)
#source("programs/ancom.R")
#library(DeSeq2)

#read in data
asv_data <- read.table(file="input/tir_feature_tables/TIR_skeleton_feature_table_filtered.tsv", sep ='\t', header = T, check.names=FALSE, row.names = 1, skip = 1, comment.char = "")
#want to order data by column interested in for heatmap
metadata <- read_tsv(file="input/GCMP_TIR_genomes_mapping.txt")
metadata.df<-data.frame(metadata)
rownames(metadata.df)<-metadata.df$SampleID

taxonomy <- read.table(file="input/taxonomy.tsv", sep = "\t", header = T, row.names=1)
phylo_tree <- "input/physeq.noncton-rooted-tree.nwk"
#clean the taxonomy file
tax <- taxonomy %>%
  select(Taxon) %>%
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";")
#a warning comes up about missing peaces but this seems to be unclassified bacteria
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2],"D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        Species = str_replace(tax[,7], "D_6__",""),
                        stringAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep=" ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified",tax.clean$Genus[i], sep = " ")
  }
}


#prepare data for the phyloseq object
asv = otu_table(as.matrix(asv_data), taxa_are_rows = TRUE)
tax = tax_table(as.matrix(tax.clean))
sample <- sample_data(metadata.df)
tree <- read_tree(phylo_tree)

#merge files into phyloseq object
#this takes a few minutes
#merge asv and tax file first followed by adding them to the metadata and tree
asv.tax <-phyloseq(asv,tax)
ps_new <- merge_phyloseq(asv.tax,sample,tree)

#choose a level to analyze the data
level="Class"

# Aggregate to phylum level
phylum_data = aggregate_taxa(ps_new, level)

#Differencial abundance using ANCOMBC if you want to run as categorical
#sample_data(phylum_data)$IL1R <- as.factor(sample_data(phylum_data)$IL1R)
#ps.taxa <-tax_glom(phylum_data, taxrank = level) #took out NArm=FALSE


#run ANCOMBC
#use the direct aggregate taxa script to run continuous data
#lets see if this runs without further subsetting
#for p_adj_method there a a few options: "holm" (default), "hochberg", "hommel",
#"bonferroni", "BH", "BY", "fdr", "none"
out <- ancombc(phyloseq = phylum_data, formula = "TLR", p_adj_method = "fdr", zero_cut =0.90, 
               lib_cut = 0, group = "TLR", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5,
               max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
res <- out$res
res_global <- out$res_global
res.global.df <- data.frame(res_global)
write.csv(res.global.df, file="coral_output/TIR_Taxonomy/skeleton/res_global_class_TLR_skeleton.csv", row.names = TRUE)

tab_coef = res$beta
tab_coef.df<-data.frame(tab_coef)
write.csv(tab_coef.df, file="coral_output/TIR_Taxonomy/skeleton/coef_class_TLR_skeleton.csv", row.names = TRUE)

tab_se = res$se

tab_w =res$W
tab_w.df<-data.frame(tab_w)
write.csv(tab_w.df, file="coral_output/TIR_Taxonomy/skeleton/w_class_TLR_skeleton.csv", row.names = TRUE)

tab_p = res$p_val
tab_p.df<-data.frame(tab_p)
write.csv(tab_p.df, file="coral_output/TIR_Taxonomy/skeleton/p_val_class_TLR_skeleton.csv", row.names = TRUE)

tab_q = res$q
tab_q.df<-data.frame(tab_q)
write.csv(tab_q.df, file="coral_output/TIR_Taxonomy/skeleton/q_val_class_TLR_skeleton.csv", row.names = TRUE)

#create a heatmap of the taxon data to visualize the differences
theme_set(theme_bw())
ps.taxa <-tax_glom(phylum_data, taxrank = level) #took out NArm=FALSE
prune <- prune_taxa(names(sort(taxa_sums(ps.taxa),TRUE)[1:300]), ps.taxa)
#sample_order<- sample_names(ps.taxa)
#sample.order <- c("3","3","3","3","3","3","3","5","5","5","6","6","6","6","6","6","6","6","6","6","6","6","6","6","6","7","7","7","7","7","7","7","8","8","8","8","8","8","8","10","10","11","11","11","11","11","11","11","11","11","11","11","11","11","11","11","11","14")
#Sample.order <- c(3:14)
#method may be the same plot ordination function: perhaps try "PCoA" but
#needs to called RDA

heatmap <- plot_heatmap(prune, method="RDA", sample.label="TLR", sample.order=metadata.df$SampleID, taxa.label="Class")
ggsave(filename ="coral_output/TIR_Taxonomy/skeleton/heatmap_class_TLR_skeleton.pdf", plot = heatmap, width=11, height=8) #9x8 good for class IL1R