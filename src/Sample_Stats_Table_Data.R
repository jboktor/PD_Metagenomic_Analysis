#################################### Fig 1 Table stats #################################### 

##### Alpha Diversity Boxplots Script

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(ggdendro);library(plyr);library(fastcluster);library(dendextend);
library(grid);library(devtools);library(RColorBrewer);library(ggfortify);library(vegan);library(vegan3d);library(MASS);
library(compositions);library(zCompositions);library(phyloseq);library(gplots);library(ape);library(lme4);library(phangorn);
library(plotly);library(VennDiagram);library(ggvegan);library(Biobase);library(BiocInstaller);library(viridis);
library("foreach"); packageVersion("foreach");library("doParallel"); packageVersion("doParallel");library(jakR);
library(ggbeeswarm);library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future)

load("Metaphlan2_PhyloseqObj.RData")

pd_dat = subset_samples(dat, donor_group == "PD")
hc_dat = subset_samples(dat, donor_group == "HC")
pc_dat = subset_samples(dat, donor_group == "PC")

# pull metadata
pd_dat_meta <- meta(pd_dat)
hc_dat_meta <- meta(hc_dat)
pc_dat_meta <- meta(pc_dat)
pd_dat_meta[pd_dat_meta=="not provided"] <-  NA
hc_dat_meta[hc_dat_meta=="not provided"] <-  NA
pc_dat_meta[pc_dat_meta=="not provided"] <-  NA
#Setting up AGE Cols
pd_dat_meta$host_age <- as.numeric(pd_dat_meta$host_age)
hc_dat_meta$host_age <- as.numeric(hc_dat_meta$host_age)
pc_dat_meta$host_age <- as.numeric(pc_dat_meta$host_age)
#Setting up BMI Cols
pd_dat_meta$host_body_mass_index <- as.numeric(pd_dat_meta$host_body_mass_index)
hc_dat_meta$host_body_mass_index <- as.numeric(hc_dat_meta$host_body_mass_index)
pc_dat_meta$host_body_mass_index <- as.numeric(pc_dat_meta$host_body_mass_index)
#Setting up BSS Cols
pd_dat_meta$bristol_stool_scale <- as.numeric(pd_dat_meta$bristol_stool_scale)
hc_dat_meta$bristol_stool_scale <- as.numeric(hc_dat_meta$bristol_stool_scale)
pc_dat_meta$bristol_stool_scale <- as.numeric(pc_dat_meta$bristol_stool_scale)


# Age Mean
mean(na.omit(pd_dat_meta$host_age))
mean(na.omit(hc_dat_meta$host_age))
mean(na.omit(pc_dat_meta$host_age))
# Age Std Dev
sd(na.omit(pd_dat_meta$host_age))
sd(na.omit(hc_dat_meta$host_age))
sd(na.omit(pc_dat_meta$host_age))

# BMI Mean
mean(na.omit(pd_dat_meta$host_body_mass_index))
mean(na.omit(hc_dat_meta$host_body_mass_index))
mean(na.omit(pc_dat_meta$host_body_mass_index))
# BMI Std Dev
sd(na.omit(pd_dat_meta$host_body_mass_index))
sd(na.omit(hc_dat_meta$host_body_mass_index))
sd(na.omit(pc_dat_meta$host_body_mass_index))

# Sex Count
count(na.omit(pd_dat_meta$sex))
count(na.omit(hc_dat_meta$sex))
count(na.omit(pc_dat_meta$sex))

# BSS Mean
mean(na.omit(pd_dat_meta$bristol_stool_scale))
mean(na.omit(hc_dat_meta$bristol_stool_scale))
mean(na.omit(pc_dat_meta$bristol_stool_scale))
# BSS Std Dev
sd(na.omit(pd_dat_meta$bristol_stool_scale))
sd(na.omit(hc_dat_meta$bristol_stool_scale))
sd(na.omit(pc_dat_meta$bristol_stool_scale))


## READS 
QC_data <- read.table(file = "qc_counts_pairs_table.tsv")
redz <- mutate(QC_data, group = if_else(grepl("HC", Helper), "HC", 
                                        if_else(grepl("PC", Helper), "PC",
                                                "PD")))
redz_pd <- filter(redz, group=="PD")
redz_pc<- filter(redz, group=="PC")
redz_hc<- filter(redz, group=="HC")

# Mean Reads per group
mean(redz_pd$Raw)/1000000
mean(redz_pc$Raw)/1000000
mean(redz_hc$Raw)/1000000

# Mean Reads per group
sd(redz_pd$Raw)/1000000
sd(redz_pc$Raw)/1000000
sd(redz_hc$Raw)/1000000

