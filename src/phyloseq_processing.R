# Caltech - Mazmanian Lab
# Joe Boktor
# October 2021

source("src/load_packages.R")
load("EDA_App/low_quality_samples.RData")
phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj.rds")

filter_samples <- function(dat) {
  dat %>%
    subset_samples(antibiotics != "Yes") %>%
    subset_samples(donor_id %nin% low_qc[[1]])
}

# Remove antibiotic using and low read samples
phyloseq_objs_clean <- phyloseq_objs %>% map(filter_samples)
saveRDS(phyloseq_objs_clean, "files/Phyloseq_Merged/PhyloseqObj_clean.rds")
saveRDS(phyloseq_objs_clean, file = "EDA_App/PhyloseqObj_clean.rds")


phyloseq_objs$Species %>% 
  subset_samples(antibiotics != "Yes")
