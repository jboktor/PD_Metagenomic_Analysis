## PERMANOVA script

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot)


load("files/Species_PhyloseqObj.RData")
load("files/KOs.all_PhyloseqObj.RData")

dat <- microbiome::transform(dat, "compositional")

dat_clr <- microbiome::transform(dat, "clr")
sp <- t(abundances(dat_clr))

# Run Metadata pre-processing function
source("src/Metadata_prep_funcs.R")
process_meta(dat)
trim_meta(env)

env.pa <- dplyr::select(env, -matches(other_metadata))

## Adding Alpha Diversity to PERMANOVA
env.pa$Shannon <- alpha(abundances(dat), 'shannon')$diversity_shannon


####################  Prep Permanova  #################### 
distlist <- c("euclidean")
sigmetadata <- vector("list", length(distlist))
names(sigmetadata) = distlist

for (j in distlist) {
  sigmetadata[[j]]$siglst <- c()
  sigmetadata[[j]]$almostSiglst <- c()
}
siglst <- c()
env.pa_perm <- tibble()

#################### PERMANOVA LOOP  #################### 
# Note with permuations = 9999, will take 5-10 mins to run

# Runs through metadata columns, if NA detected, 
# sample is omitted from metdata & same sample from abundnace table
for (i in 1:length(env.pa)) { 
  a <- env.pa[,i]
  a.narm <- na.omit(a)
  if (any(is.na(a))) {
    sp.narm <- sp[-attr(a.narm, "na.action"), ]
  } else {
    sp.narm <- sp
  }
  for (j in distlist) {
    meta_ano = adonis(dist(sp.narm, method = "euclidean") ~ a.narm, permutations = 9999)
    row2add <- cbind(colnames(env.pa[i]),meta_ano$aov.tab[1,],length(a.narm))
    env.pa_perm <- rbind(env.pa_perm, row2add)
    sig <- meta_ano$aov.tab$`Pr(>F)`[1]
    if (sig < .1 & sig > .05) {
      mssage <- paste(colnames(env.pa[i]), ": near significant metadata column (p < 0.1) with", j, "Distance")
      print(mssage)
      sigmetadata[[j]]$almostSiglst <- c(sigmetadata[[j]]$almostSiglst, colnames(env.pa[i]))
    }
    if (sig < .05) {
      mssage <- paste(colnames(env.pa[i]), "is a significant metadata column (p < 0.05) with", j, "Distance")
      print(mssage)
      sigmetadata[[j]]$siglst <- c(sigmetadata[[j]]$siglst, colnames(env.pa[i]))
    }
  }
}


####################  VARS  - For Plotting & Confounder Analysis #################### 
# 
permdf <- env.pa_perm
permdf$FDR <- p.adjust(permdf$`Pr(>F)`, method = 'BH', n = nrow(permdf)) # Multiple comparison correction using Benjamini Hochberg Method
rownames(permdf) <- NULL
colnames(permdf)[1] <- "vars"
colnames(permdf)[8] <- "n_meta"
permdf$R2 <- permdf$R2*100

## Add Metadata catagory column
permdf <- mutate(permdf, metacat = if_else(permdf$vars %in% grouping, "Disease-Status/Grouping",
                                           if_else(permdf$vars %in% Anthro, "Anthropomorphic",
                                                   if_else(permdf$vars %in% GI, "GI distress",
                                                           if_else(permdf$vars %in% nonstarters, "Microbiota disruption",
                                                                   if_else(permdf$vars %in% allergy, "Allergies",
                                                                           if_else(permdf$vars %in% smoke, "Smoking stats",
                                                                                   if_else(permdf$vars %in% general_drugs, "Other medication",
                                                                                           if_else(permdf$vars %in% pd_drugs, "PD medication",
                                                                                                   if_else(permdf$vars %in% supplements, "Supplements",
                                                                                                           if_else(permdf$vars %in% diet, "Diet",
                                                                                                                   if_else(permdf$vars == "Shannon", "Shannon Diversity", 
                                                                                                                           "AA_notmatched"))))))))))))



write.csv(permdf, file = 'files/permanova_data.csv')

####################  FILTERING FOR FDR SIGNIFICANT 

permdfsig <- filter(permdf, FDR <= 0.05)


#################### FILTERING FOR CORE METADATA - n > 90 (Metadata present with most samples)

permdfcore <- filter(permdf, n_meta >= 90)
permdfcoresig <- filter(permdfcore, FDR <= 0.5)


#################### FILTERING FOR Essential METADATA

essential <- c("Disease-Status/Grouping", "Anthropomorphic", "GI distress", "Microbiota disruption", "Shannon Diversity")
permdfess <- filter(permdf, metacat %in% essential)


#################### FILTERING FOR DIET METADATA

permdf_diet <- filter(permdf, metacat == "Diet")


#################### FILTERING FOR Supplement METADATA

permdf_supp <- filter(permdf, metacat == "Supplements")


#################### FILTERING FOR General Drug use METADATA

permdf_gendrug <- filter(permdf, metacat == "Other medication")


#################### FILTERING FOR PD Drug use METADATA

permdf_pddrug <- filter(permdf, metacat == "PD medication")




