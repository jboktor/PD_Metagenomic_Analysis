# Machine Learning - Data Visualization

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/ml_models.R")
source("src/ml_plots.R")
load("files/low_quality_samples.RData")

# remove_dats()
# load_data("Merged_ML")

LOSO_KOs <-
  read.csv(file = "files/Machine_Learning_Models/KOs/LOSO_performance.csv", header = TRUE)
