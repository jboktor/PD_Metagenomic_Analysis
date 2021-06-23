# Feature Correlation Analysis

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
load("files/low_quality_samples.RData")
load_data("Merged")

#_______________________________________________________________________________
#####                      Correlation Analysis                           ##### 
#_______________________________________________________________________________

# Metadata Selection
corr.meta.clinical <-
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  process_meta(cohort = "Merged") %>%
  dplyr::select(contains(
    c(motor_severity_scores_summary,
      clinical_variables
    )
  ))
corr.meta.diet <-
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  process_meta(cohort = "Merged") %>%
  dplyr::select(contains(
    c(diet
    )
  ))

#_______________________________________________________________________________
#####                              Species                               ##### 
#_______________________________________________________________________________

corr.abund.species <- 
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  core(detection = 0, prevalence = 1/234) %>% 
  abundances() %>% 
  t() %>%
  as.data.frame()

# Correlation calculations
corrs.clinical.species <-
  corr_loop_parallel(metadata = corr.meta.clinical, abundance = corr.abund.species, obj.name = "Species")
corrs.diet.species <-
  corr_loop_parallel(metadata = corr.meta.diet, abundance = corr.abund.species, obj.name = "Species")

# Scatter plots
# top_n_scatterplots(dat = dat.species, obj.name = "Species", 
#                    df.cors = corrs.clinical.species, metaclass = "Clinical")
# top_n_scatterplots(dat = dat.species, obj.name = "Species", 
#                    df.cors = corrs.diet.species, metaclass = "Diet")

# Heatmaps
h1 <- corr_heatmap(corrs.clinical.species)
h2 <- corr_heatmap(corrs.diet.species)
esedhh

#_______________________________________________________________________________
#####                              Pathways                               ##### 
#_______________________________________________________________________________

corr.abund.pathways <-
  corr_abund_prep(obj.name = "Pathways.slim",
                  obj = dat.path.slim,
                  cohort = "Merged")

# Correlation calculations
corrs.clinical.pathways.slim <-
  corr_loop_parallel(metadata = corr.meta.clinical, abundance = corr.abund.pathways, obj.name = "Pathways")
corrs.diet.pathways.slim <-
  corr_loop_parallel(metadata = corr.meta.diet, abundance = corr.abund.pathways, obj.name = "Pathways")

# Scatter plots
# top_n_scatterplots(dat = dat.path.slim, obj.name = "Pathways.slim", 
#                    df.cors = corrs.clinical.pathways.slim, metaclass = "Clinical")
# top_n_scatterplots(dat = dat.path.slim, obj.name = "Pathways.slim", 
#                    df.cors = corrs.diet.pathways.slim, metaclass = "Diet")

# Heatmaps
h3 <- corr_heatmap(corrs.clinical.pathways.slim)
h4 <- corr_heatmap(corrs.diet.pathways.slim)
ggsave(h3, filename = "data/Correlations/Pathways.slim/Clinical/Heatmap_Pathways.slim_clinical.svg",
       height = 10, width = 13)
ggsave(h4, filename = "data/Correlations/Pathways.slim/Diet/Heatmap_Pathways.slim_diet.svg",
       height = 9, width = 11)

#_______________________________________________________________________________
#####                              KOs                               ##### 
#_______________________________________________________________________________

corr.abund.KOs <-
  corr_abund_prep(obj.name = "KOs.slim",
                  obj = dat.KOs.slim,
                  cohort = "Merged")

# Correlation calculations
corrs.clinical.KOs.slim <-
  corr_loop_parallel(metadata = corr.meta.clinical, abundance = corr.abund.KOs, obj.name = "KOs")
corrs.diet.KOs.slim <-
  corr_loop_parallel(metadata = corr.meta.diet, abundance = corr.abund.KOs, obj.name = "KOs")

# Scatter plots
# top_n_scatterplots(dat = dat.KOs.slim, obj.name = "KOs.slim", 
#                    df.cors = corrs.clinical.KOs.slim, metaclass = "Clinical")
# top_n_scatterplots(dat = dat.KOs.slim, obj.name = "KOs.slim", 
#                    df.cors = corrs.diet.KOs.slim, metaclass = "Diet")

# Heatmaps
h5 <- corr_heatmap(corrs.clinical.KOs.slim)
h6 <- corr_heatmap(corrs.diet.KOs.slim)
ggsave(h5, filename = "data/Correlations/KOs.slim//Clinical/Heatmap_KOs.slim_clinical.svg", 
       height = 10, width = 12)
ggsave(h6, filename = "data/Correlations/KOs.slim/Diet/Heatmap_KOs.slim_diet.svg", 
       height = 10, width = 13)



#_______________________________________________________________________________
#####                              GOs                               ##### 
#_______________________________________________________________________________

corr.abund.GOs <-
  corr_abund_prep(obj.name = "GOs.slim",
                  obj = dat.GOs.slim,
                  cohort = "Merged")

# Correlation calculations
corrs.clinical.GOs.slim <-
  corr_loop_parallel(metadata = corr.meta.clinical, abundance = corr.abund.GOs, obj.name = "GOs")
corrs.diet.GOs.slim <-
  corr_loop_parallel(metadata = corr.meta.diet, abundance = corr.abund.GOs, obj.name = "GOs")

# Scatter plots
# top_n_scatterplots(dat = dat.GOs.slim, obj.name = "GOs.slim",
#                    df.cors = corrs.clinical.GOs.slim, metaclass = "Clinical")
# top_n_scatterplots(dat = dat.GOs.slim, obj.name = "GOs.slim",
#                    df.cors = corrs.diet.GOs.slim, metaclass = "Diet")

# Heatmaps
h7 <- corr_heatmap(corrs.clinical.GOs.slim)
h8 <- corr_heatmap(corrs.diet.GOs.slim)
ggsave(h7, filename = "data/Correlations/GOs.slim/Clinical/Heatmap_GOs.slim_clinical.svg", 
       height = 10, width = 13)
ggsave(h8, filename = "data/Correlations/GOs.slim/Diet/Heatmap_GOs.slim_diet.svg", 
       height = 10, width = 13)



#_______________________________________________________________________________
#####                              PFAMs                               ##### 
#_______________________________________________________________________________

corr.abund.PFAMs <-
  corr_abund_prep(obj.name = "PFAMs.slim",
                  obj = dat.PFAMs.slim,
                  cohort = "Merged")

# Correlation calculations
corrs.clinical.PFAMs.slim <-
  corr_loop_parallel(metadata = corr.meta.clinical, abundance = corr.abund.PFAMs, obj.name = "PFAMs")
corrs.diet.PFAMs.slim <-
  corr_loop_parallel(metadata = corr.meta.diet, abundance = corr.abund.PFAMs, obj.name = "PFAMs")

# Scatter plots
# top_n_scatterplots(dat = dat.PFAMs.slim, obj.name = "PFAMs.slim", 
#                    df.cors = corrs.clinical.PFAMs.slim, metaclass = "Clinical")
# top_n_scatterplots(dat = dat.PFAMs.slim, obj.name = "PFAMs.slim", 
#                    df.cors = corrs.diet.PFAMs.slim, metaclass = "Diet")

# Heatmaps
h9 <- corr_heatmap(corrs.clinical.PFAMs.slim)
h10 <- corr_heatmap(corrs.diet.PFAMs.slim)
ggsave(h9, filename = "data/Correlations/PFAMs.slim/Clinical/Heatmap_PFAMs.slim_clinical.svg", 
       height = 10, width = 12)
ggsave(h10, filename = "data/Correlations/PFAMs.slim/Diet/Heatmap_PFAMs.slim_diet.svg", 
       height = 10, width = 13)

# corr_xy(obj = dat.PFAMs.slim, corr_obj = corrs.clinical.PFAMs.slim,
#         feature_var = "PF12675: Protein of unknown function (DUF3795)",
#         metadata_var = "family_history_pd_degree_relative")


#_______________________________________________________________________________
#####                      Saving as Excel files                           ##### 
#_______________________________________________________________________________
clinical_corrs <-
  list(Species = corrs.clinical.species,
       Pathways = corrs.clinical.pathways.slim,
       Kegg_Orthology = corrs.clinical.KOs.slim,
       Gene_Ontology = corrs.clinical.GOs.slim,
       PFAMs = corrs.clinical.PFAMs.slim)

dietary_corrs <-
  list(Species = corrs.diet.species,
       Pathways = corrs.diet.pathways.slim,
       Kegg_Orthology = corrs.diet.KOs.slim,
       Gene_Ontology = corrs.diet.GOs.slim,
       PFAMs = corrs.diet.PFAMs.slim)

# openxlsx::write.xlsx(clinical_corrs, file = 'files/Correlations/Clinical_Correlations.xlsx')
# openxlsx::write.xlsx(dietary_corrs, file = 'files/Correlations/Dietary_Correlations.xlsx')
# save(clinical_corrs, file = 'files/Correlations/Clinical_Correlations.RData')
# save(dietary_corrs, file = 'files/Correlations/Dietary_Correlations.RData')


#' #_______________________________________________________________________________
#' #####                      Merge as one large analysis                           ##### 
#' #_______________________________________________________________________________\
#' 
#' all_corrs <- bind_rows(corrs.clinical.species, 
#'                           corrs.diet.species,
#'                           corrs.clinical.pathways.slim,
#'                           corrs.diet.pathways.slim,
#'                           corrs.clinical.KOs.slim,
#'                           corrs.diet.KOs.slim,
#'                           corrs.clinical.PFAMs.slim,
#'                           corrs.diet.PFAMs.slim,
#'                           corrs.clinical.GOs.slim,
#'                           corrs.diet.GOs.slim)
#' 
#' write.csv(all_corrs, file = 'files/correlations_analysis.csv')

