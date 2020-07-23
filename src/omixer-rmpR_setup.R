# Gut-Brain Module (GBM) and Gut-Metabolic Module (GMM) Analysis

############  Prep Metadata  ############
library(microbiome); library(phyloseq); library(Maaslin2); library(dplyr); library(tidyr); 
library(reshape2); library(lattice);library(EnvStats);
library(omixerRpm)

rm(list = ls())

######## Load Data & functions
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")

input_kos <- microbiome::abundances(dat.KOs.slim) %>% 
  data.frame() %>% 
  rownames_to_column(var = "entry")
# listDB()


# Run the module mapping on the loaded table.
gmm <- rpm(input_kos, minimum.coverage=0.3, annotation = 1, 
           module.db	= loadDB("GMMs.v1.07"))
gmm.abundnace <- gmm@abundance
gmm.annotation <- gmm@annotation
gmm.coverage <- gmm@coverage

gbm <- rpm(input_kos, minimum.coverage=0.3, annotation = 1, 
            module.db	= loadDB("GBMs.v1.0"))
gbm.abundnace <- gbm@abundance
gbm.annotation <- gbm@annotation
gbm.coverage <- gbm@coverage

# Load the default mapping database
db_gmm <- loadDB("GMMs.v1.07")
db_gbm <- loadDB("GBMs.v1.0")


# Create dataframes with abundnces and modules along with translated name
####### GBM ####### 
translated_annotations_gbm <- c()
modname_list_gbm <- c()
for (mod in gbm@annotation){
  for (modname in mod){
    modname_list_gbm <- c(modname_list_gbm, modname)
    print(modname)
    print(as.character(omixerRpm::getNames(db_gbm, modname) ))
    translated_annotations_gbm <- c(translated_annotations_gbm, omixerRpm::getNames(db_gbm, modname) )
  }
}
df.gbm <- data.frame(row.names = translated_annotations_gbm,
                     "module"=modname_list_gbm,
                     gbm.abundnace)

# Distribtution Sanity Check
source("src/miscellaneous_funcs.R")
distribution_sanity(df.gbm)

####### GMM ####### 
translated_annotations_gmm <- c()
modname_list_gmm <- c()
for (mod in gmm@annotation){
  for (modname in mod){
    modname_list_gmm <- c(modname_list_gmm, modname)
    print(modname)
    print(as.character(omixerRpm::getNames(db_gmm, modname) ))
    translated_annotations_gmm <- c(translated_annotations_gmm, omixerRpm::getNames(db_gmm, modname) )
  }
}
df.gmm <- data.frame(row.names=translated_annotations_gmm,
                     "module"=modname_list_gmm,
                     gmm.abundnace)

# Distribtution Sanity Check
distribution_sanity(df.gmm)


# Split R objects for metadata
dat_pdhc = subset_samples(dat, Paired !="No")
dat_pdpc = subset_samples(dat, donor_group !="HC")

# GMM abundance tables
df.gmm_pdpc <- df.gmm[, which(!grepl("HC", colnames(df.gmm)))]
df.gmm_pdhc <- df.gmm[, which(!grepl("PC", colnames(df.gmm)))]

# GBM abundance tables
df.gbm_pdpc <- df.gbm[, which(!grepl("HC", colnames(df.gbm)))]
df.gbm_pdhc <- df.gbm[, which(!grepl("PC", colnames(df.gbm)))]


######## Format Metadata 
# Run Metadata pre-processing function
process_meta(dat)
df_input_metadata <- env %>% column_to_rownames(var = "donor_id")

process_meta(dat_pdpc)
df_input_metadata_pdpc <- env %>% column_to_rownames(var = "donor_id")
df_input_metadata_pdpc$description <- factor(df_input_metadata_pdpc$description, 
                                             levels = c("PD Patient", "Population Control"))

process_meta(dat_pdhc)
df_input_metadata_pdhc <- env %>% column_to_rownames(var = "donor_id")
df_input_metadata_pdhc$description <- factor(df_input_metadata_pdhc$description, 
                                             levels = c("PD Patient", "Household Control"))


################################################################################### 
###################                 GMM MODELS                  ###################
################################################################################### 
#set file path 
wkd <- getwd()

############  PD v PC - GMMs ############

fit_data = Maaslin2(
  input_data = df.gmm_pdpc[,-1],
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/GMM_PDvPC_maaslin2_output"), 
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
    transform = "AST",
  cores = 1
)

############  PD v HC Paired - GMMs ############

fit_data = Maaslin2(
  input_data = df.gmm_pdhc[,-1], 
  input_metadata = df_input_metadata_pdhc, 
  output = paste0(wkd, "/data/MaAsLin2_Analysis/GMM_PDvHC_maaslin2_output"), 
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

################################################################################### 
###################                 GBM MODELS                  ###################
################################################################################### 

############  PD v PC - GBMs ############

fit_data = Maaslin2(
  input_data = df.gbm_pdpc[,-1],
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/GBM_PDvPC_maaslin2_output"), 
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)
############  PD v HC Paired - GBMs ############

fit_data = Maaslin2(
  input_data = df.gbm_pdhc[,-1],
  input_metadata = df_input_metadata_pdhc, 
  output = paste0(wkd, "/data/MaAsLin2_Analysis/GBM_PDvHC_maaslin2_output"), 
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)



#################
