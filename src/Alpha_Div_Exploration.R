# Alpha Diversity EDA 

rm(list = ls())

######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")


################################################################################# 
########################### SWAP INPUT LEVEL HERE: ########################### 
################################################################################# 

# Phylo_Objects$Species
obj <- "Species"
dat_obj <- Phylo_Objects[[obj]]


# alpha(abundances(dat_obj), 'observed')
# alpha(abundances(dat_obj), 'shannon')
# evenness(abundances(dat_obj), 'simpson')

################################################################################# 

#-------------------------------------------
# Get PseudoCounts for Alpha Diversity

reads <- load_reads()
reads$id <- gsub("_", ".", reads$id)
pseudo_counts <- PseudoCounts(dat_obj, reads) %>% 
  round(digits = 0)
## Calculate Alpha Diversity Metrics
tab <- microbiome::alpha(pseudo_counts, index = "all")

#-------------------------------------------
# Run Metadata pre-processing function
process_meta(dat)
env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))

#-------------------------------------------
# Add identifier cols to tab df
tab$donor_group = env$donor_group
tabm <- melt(tab)

#-------------------------------------------
ggplot(tabm, aes(x = donor_group, y = value)) + theme_minimal() + 
  geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
  geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
  geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
  facet_wrap(~variable, scales = "free") + 
  theme(axis.title.x=element_blank(),
        legend.position = "none") 


# # Distribution sanity check
# distribution_sanity()
# # histogram(melt(abundances(dat_obj))$value)
# abund.melt <- melt(abundances(dat_obj))
# abund.melt <- mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC", if_else(grepl("PC", Var2), "PC", "PD")))




