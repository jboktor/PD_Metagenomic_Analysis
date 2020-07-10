# Alpha Diversity Exploration 


library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot);
library(EnvStats);library(sjlabelled);library(sjmisc);library(sjPlot);library(nlme)

rm(list = ls())

######## Load Data & functions
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")


################################################################################# 
########################### SWAP INPUT LEVEL HERE: ########################### 
################################################################################# 

# Phylo_Objects$Species
obj <- "Species"
dat_obj <- Phylo_Objects[[obj]]


alpha(abundances(dat_obj), 'observed')
alpha(abundances(dat_obj), 'shannon')
evenness(abundances(dat_obj), 'simpson')

################################################################################# 

### WARNING : OTU TABLE CONVERTED TO BINARY - DO NOT USE DOWNSTREAM FOR other analyses
dat_alpha <- dat_obj 
otu_table(dat_alpha) <- ((otu_table(dat_alpha) >  0) + 0)

## Calculate Alpha Diversity Metrics
tab <- microbiome::alpha(dat_alpha, index = "all")
# tab <- microbiome::alpha(dat_alpha, index= c("observed" ,'shannon', "rarity_log_modulo_skewness", "diversity_gini_simpson", 
#                                              "diversity_inverse_simpson", "diversity_coverage", "dominance_simpson", "dominance_core_abundance",
#                                              "dominance_gini"))


# Run Metadata pre-processing function
process_meta(dat_alpha)
env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))


# Add identifier cols to tab df
tab$donor_group = env$donor_group


tabm <- melt(tab)

ggplot(tabm, aes(x = donor_group, y = value)) + theme_minimal() + 
  geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
  geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
  geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
  facet_wrap(~variable, scales = "free") + 
  theme(axis.title.x=element_blank(),
        legend.position = "none") 



# Distribution sanity check
histogram(melt(abundances(dat_obj))$value)
abund.melt <- melt(abundances(dat_obj))
abund.melt <- mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC", if_else(grepl("PC", Var2), "PC", "PD")))



