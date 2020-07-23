# Differential Abundance of Features - MaAsLiN2 Script

######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")


#set file path 
wkd <- getwd()


## Add prevalence shift on all samples prior to partioning 
dat <- core(dat, detection = 0, prevalence = 0.1)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.3) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.3) %>%  transform("compositional")


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
###################                    MODELS                  ###################
################################################################################### 


################################################################################### 
###############################         TAXA        ############################## 
 

############  PD v PC - SPECIES ############

fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Species_PDvPC_maaslin2_output"), 
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

############  PD v HC Paired  - SPECIES ############

fit_data = Maaslin2(
  input_data = df_input_data_pdhc, 
  input_metadata = df_input_metadata_pdhc, 
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Species_PDvHC_maaslin2_output"), 
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

######################## Simplified Model for Confounder Analysis  ########################   

############  Taxa Data  ############ All Groups

# fit_data = Maaslin2(
#   input_data = df_input_data, 
#   input_metadata = df_input_metadata, 
#   output =  paste0(wkd, "/data/MaAsLin2_Analysis/taxa_allsamples_PDeffect_maaslin2_output"),
#   min_prevalence = 0,
#   random_effects = c("Paired"),
#   fixed_effects = c("PD"),
#   analysis_method = "LM",
#   normalization = "NONE",
#   transform = "AST",
#   cores = 1
# )


###################   TAXA ANALYSIS  - GENUS  ################### 
## Add prevalence shift on all samples prior to partioning 
dat.genus <- core(dat.genus, detection = 0, prevalence = 0.1)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.genus, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.1) %>% transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.genus, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.1) %>% transform("compositional")


#############  PD v PC - GENUS ############

fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Genus_PDvPC_maaslin2_output"), 
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

############  PD v HC Paired  - GENUS ############

fit_data = Maaslin2(
  input_data = df_input_data_pdhc, 
  input_metadata = df_input_metadata_pdhc, 
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Genus_PDvHC_maaslin2_output"), 
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)




###################   TAXA ANALYSIS  - PHYLUM  ################### 

## Add prevalence shift on all samples prior to partioning 
dat.phylum <- core(dat.phylum, detection = 0, prevalence = 0.1)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.phylum, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>% abundances() %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.phylum, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>% abundances() %>%  transform("compositional")


#############  PD v PC - PHYLUM ############

fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Phylum_PDvPC_maaslin2_output"), 
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

############  PD v HC Paired  - PHYLUM ############

fit_data = Maaslin2(
  input_data = df_input_data_pdhc, 
  input_metadata = df_input_metadata_pdhc, 
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Phylum_PDvHC_maaslin2_output"), 
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)





######################################################################################
###############################        PATHWAYS        ###############################

############  Prep Data  ############
## Add prevalence shift on all samples prior to partioning 
dat.path <- core(dat.path, detection = 0, prevalence = 0.1)
taxa_names(dat.path) <- prep_pathway_names(dat.path)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.path, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <-  LowVarianceFilter(dat_pdpc, filter.percent = 0.3) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.path, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.3) %>%  transform("compositional")


############  Pathway Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Pathways_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  Pathway Data  ############ PD v HC Paired
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Pathways_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)


############  Prep Data PATHWAYS SLIM  ############
## Add prevalence shift on all samples prior to partioning 
dat.path.slim <- core(dat.path.slim, detection = 0, prevalence = 0.1)
taxa_names(dat.path.slim) <- prep_pathway_names(dat.path.slim)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.path.slim, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.1) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.path.slim, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc) #Trim 10%
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.1) %>%  transform("compositional")


############  Pathway Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Pathways.slim_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  Pathway Data  ############ PD v HC Paired
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Pathways.slim_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)



######################################################################################
###############################        ENZYMES        ###############################

############  Prep Data  ############
## Add prevalence shift on all samples prior to partioning 
dat.ec <- core(dat.ec, detection = 0, prevalence = 0.1)
taxa_names(dat.ec) <- paste0("ENZYME_", taxa_names(dat.ec))

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.ec, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.7) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.ec, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.7) %>%  transform("compositional")


############  ENZYME Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Enzymes_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  ENZYME Data  ############ PD v HC Paired - Takes a really long time
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Enzymes_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)


############  Prep Data ENZYMES - SLIM ############

## Add prevalence shift on all samples prior to partioning 
dat.ec.slim <- core(dat.ec.slim, detection = 0, prevalence = 0.1)
taxa_names(dat.ec.slim) <- paste0("ENZYME_", taxa_names(dat.ec.slim))

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.ec.slim, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.4) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.ec.slim, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.4) %>%  transform("compositional")


############  ENZYME Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Enzymes.slim_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  ENZYME Data  ############ PD v HC Paired - Takes a really long time
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/Enzymes.slim_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)




######################################################################################
###############################        KOs        ###############################

############  Prep Data  ############
## Add prevalence shift on all samples prior to partioning 
dat.KOs <- core(dat.KOs, detection = 0, prevalence = 0.1)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.KOs, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.7) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.KOs, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.7) %>%  transform("compositional")


############  Gene Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/KOs_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  Gene Data  ############ PD v HC Paired
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/KOs_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)



## KOs only with no stratification
############  Prep Data  ############
## Add prevalence shift on all samples prior to partioning 
dat.KOs.slim <- core(dat.KOs.slim, detection = 0, prevalence = 0.1)

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc = subset_samples(dat.KOs.slim, donor_group !="HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- LowVarianceFilter(dat_pdpc, filter.percent = 0.5) %>%  transform("compositional")

# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc = subset_samples(dat.KOs.slim, Paired !="No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- LowVarianceFilter(dat_pdhc, filter.percent = 0.5) %>%  transform("compositional")


############  Gene Data  ############ PD v PC
fit_data = Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/KOs.slim_PDvPC_maaslin2_output"),
  min_prevalence = 0,
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)

############  Gene Data  ############ PD v HC Paired
fit_data = Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/KOs.slim_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  plot_scatter = F,
  cores = 1)


# sessionInfo()
