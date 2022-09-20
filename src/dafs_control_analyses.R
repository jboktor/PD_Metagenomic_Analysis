# Caltech - Mazmanian Lab
# Joe Boktor
# October 2021

# Load Data & functions
rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
wkd <- getwd()
phylo_objects <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")

# Testing levels and variance thresholds determined previously
test_list <- c(
  "Species" = 0.1,
  "Genus" = 0,
  "Phylum" = 0,
  "Pathways.slim" = 0.2,
  "Enzymes.slim" = 0.3,
  "KOs.slim" = 0.2,
  "GOs.slim" = 0.3,
  "Pfams.slim" = 0.2,
  "eggNOGs.slim" = 0.7
)

df_input_metadata <- process_meta(phylo_objects$Species, cohort = "Merged") %>% 
  mutate(bristol_stool_scale = as.numeric(bristol_stool_scale)) 

general_confounder_list <- c(
  "bristol_stool_scale",
  "nonsteroidal_anti_inflammatory",
  "ssri_antidepressants",
  "vitamin_C",
  "calcium",
  "proton_pump_inhibitors",
  "laxatives")

pd_meds <- c(
  "MAO_B_inhibitors",
  "levodopa",
  "carbidopa",
  "dopamine_agonists",
  "amantadine"
)

#________________________________________________________________________________
#                          Model Fitting ----
#________________________________________________________________________________

for (level in names(test_list)) {
  
  df_input_data <- phylo_objects[[level]] %>%
    microbiome::transform("compositional") %>%
    variance_filter(filter.percent = test_list[level])
  
  for (confounding_var in general_confounder_list) {
    
    cat("Now Analyzing: ", level, "for ", confounding_var, "\n")
    df_input_metadata[[confounding_var]] %>% table() %>% print
    outfile <- paste0("/data/MaAsLin2_Analysis/Confounder_Testing/general_variables/", 
                      level, "__", confounding_var, "_maaslin2_output")
    
    
    fit_data = Maaslin2(
      input_data = df_input_data,
      input_metadata = df_input_metadata,
      output = paste0(wkd, outfile), 
      random_effects = "cohort",
      fixed_effects =  c("description", "host_age_factor", "sex", "host_body_mass_index", confounding_var),
      reference = c("description,PD Patient"),
      min_prevalence = 0,
      analysis_method = "LM",
      normalization = "NONE",
      transform = "AST",
      cores = 8,
      plot_scatter = F
    )
    
  }
}


for (level in names(test_list)) {
  
  df_input_data <- phylo_objects[[level]] %>%
    subset_samples(PD == "Yes") %>% 
    microbiome::transform("compositional") %>%
    variance_filter(filter.percent = test_list[level])
  
  for (confounding_var in pd_meds) {
    
    cat("Now Analyzing: ", level, "for ", confounding_var, "\n")
    df_input_metadata[[confounding_var]] %>% table() %>% print
    outfile <- paste0("/data/MaAsLin2_Analysis/Confounder_Testing/PD_medications/", 
                      level, "__", confounding_var, "_maaslin2_output")
    
    
    fit_data = Maaslin2(
      input_data = df_input_data,
      input_metadata = df_input_metadata,
      output = paste0(wkd, outfile), 
      random_effects = "cohort",
      fixed_effects =  c("host_age_factor", "sex", "host_body_mass_index", confounding_var),
      min_prevalence = 0,
      analysis_method = "LM",
      normalization = "NONE",
      transform = "AST",
      cores = 8,
      plot_scatter = F
    )
    
  }
}

#_______________________________________________________________________________
# Summary Table for Associations




disease_stats <- tibble()
# shared_features <- vector(mode = "list")
reportdir <- paste0("data/MaAsLin2_Analysis/Merged/")
filepaths <- list.files(path = reportdir)
for (report in filepaths) {
  maaslin_data <-
    read_tsv(
      paste0(reportdir, report, "/all_results.tsv"),
      col_names = T
    ) %>% 
    mutate(feature_level = str_split(report, "_")[[1]][1])
  disease_stats <- bind_rows(disease_stats, maaslin_data)
}

disease_stats_sig <- disease_stats %>% 
  filter(metadata == "description", 
        qval <= 0.1) %>% 
  decode_rfriendly_rows(passed_column = "feature") %>%
  dplyr::select(-feature) %>%
  dplyr::rename("feature" = "fullnames")


general_variables_stats <- tibble()
# shared_features <- vector(mode = "list")
reportdir <- paste0("data/MaAsLin2_Analysis/Confounder_Testing/general_variables/")
filepaths <- list.files(path = reportdir)
for (report in filepaths) {
  maaslin_data <-
    read_tsv(
      paste0(reportdir, report, "/all_results.tsv"),
      col_names = T
    ) %>% 
    mutate(feature_level = str_split(report, "__")[[1]][1])
  general_variables_stats <- bind_rows(general_variables_stats, maaslin_data)
}
medications_stats <- tibble()
# shared_features <- vector(mode = "list")
reportdir <- paste0("data/MaAsLin2_Analysis/Confounder_Testing/PD_medications/")
filepaths <- list.files(path = reportdir)
for (report in filepaths) {
  maaslin_data <-
    read_tsv(
      paste0(reportdir, report, "/all_results.tsv"),
      col_names = T
    ) %>% 
    mutate(feature_level = str_split(report, "__")[[1]][1])
  medications_stats <- bind_rows(medications_stats, maaslin_data)
}

confounder_stats <- bind_rows(medications_stats, general_variables_stats)

confounder_stats$metadata %>% table()

confounder_stats_sig <- confounder_stats %>% 
  filter(metadata %nin% c("description", "bristol_stool_scale")) %>% 
  filter(qval <= 0.1) %>% 
  decode_rfriendly_rows(passed_column = "feature") %>%
  dplyr::select(-feature) %>%
  dplyr::rename("feature" = "fullnames")

  confounder_list <- vector(mode = "list")
  confounder_list$disease_confound_overlap <- disease_stats_sig %>% filter(feature %in% confounder_stats_sig$feature)
  for (confnd in unique(confounder_stats_sig$metadata)){
    confounder_list[[confnd]] <- confounder_stats_sig %>% 
      filter(metadata == confnd) %>% 
      arrange(feature_level, pval)
  }
  openxlsx::write.xlsx(confounder_list, file = 'files/Supplementary Tables/Table_S7 Confounder Estimation.xlsx', overwrite = T)
  


































