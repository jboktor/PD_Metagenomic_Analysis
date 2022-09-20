# Caltech - Mazmanian Lab
# Joe Boktor
# October 2021
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")

phyloseq_objs_rare <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean_rarefied.rds")
phyloseq_objs_rare_TBC <- phyloseq_objs_rare %>% 
  purrr::map(~.x %>% subset_samples(cohort == "TBC" ))
phyloseq_objs_rare_Rush <- phyloseq_objs_rare %>% 
  purrr::map(~.x %>% subset_samples(cohort == "Rush" ))



# TBC Species ADONIS stats
obj_processed <- phyloseq_objs_rare_TBC[["Species"]] %>% 
  microbiome::transform("compositional") %>% 
  microbiome::transform("clr")

sp <- t(abundances(obj_processed))
env <- process_meta(obj_processed, cohort = "Merged")
a <- env[ ,"description"]
permanova_tbc <- adonis(vegdist(sp, method = "euclidean") ~ a, permutations = 99999)


# RUMC Species ADONIS stats
obj_processed <- phyloseq_objs_rare_Rush[["Species"]] %>% 
  microbiome::transform("compositional") %>% 
  microbiome::transform("clr")

sp <- t(abundances(obj_processed))
env <- process_meta(obj_processed, cohort = "Merged")
a <- env[ ,"description"]
permanova_rush <- adonis(vegdist(sp, method = "euclidean") ~ a, permutations = 99999)

