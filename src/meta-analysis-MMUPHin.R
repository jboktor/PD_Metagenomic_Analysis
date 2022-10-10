# Meta-Analysis with MMUPHin

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
library(MMUPHin)

#_______________________________________________________________________________
# Utility Functions ----

format_metadata <- function(ps) {
  df <- meta(ps) %>%
    as.data.frame() %>%
    select(PD,
           paired,
           sex,
           host_age,
           host_body_mass_index,
           country,
           cohort) %>%
    mutate_if(is.character, na_if, "not provided") %>%
    dplyr::mutate(
      PD = factor(PD, levels = c("No", "Yes")),
      sex = factor(sex, levels = c("male", "female")),
      country = if_else(cohort %in% c("TBC", "Rush"), "USA", country),
      country = factor(country, levels = c("USA", "Germany", "China")),
      paired = factor(paired),
      host_age = as.numeric(host_age),
      host_body_mass_index = as.numeric(host_body_mass_index)
    )
  
  if (!all(is.na(df$host_age))) {
    df %<>% mutate(host_age_factor =
                     as.numeric(
                       cut(
                         host_age,
                         breaks = quantile(host_age, na.rm = T),
                         labels = c(1, 2, 3, 4),
                         include.lowest = T
                       )
                     ))
  }
  if (!all(is.na(df$host_body_mass_index))) {
    df %<>% mutate(
      host_body_mass_index_factor = case_when(
        host_body_mass_index < 18.5 ~ 0,
        host_body_mass_index >= 35 ~ 4,
        host_body_mass_index < 35 &
          host_body_mass_index >= 30 ~ 3,
        host_body_mass_index < 30 &
          host_body_mass_index >= 25 ~ 2,
        host_body_mass_index < 25 &
          host_body_mass_index >= 18.5 ~ 1,
        TRUE ~ host_body_mass_index
      )
    )
  }
  return(df)
}
#_______________________________________________________________________________

ps <- readRDS(file = "files/Phyloseq_all-cohorts_clean_slim.rds")
cohorts <- c("TBC", "Rush", "Shanghai", "Bonn")
df_meta <- format_metadata(ps$Species) %>% 
  mutate(cohort = factor(cohort, levels = c("TBC", "Rush", "Bonn", "Shanghai")))
MMUPHin_models <- list()
MMUPHin_dir <- "data/MMUPHin_Meta-Analysis"

for (level in c("Species", "Pathways.slim", "KOs.slim", "Pfams.slim", "GOs.slim")) {
  output_dir <- glue("{MMUPHin_dir}/{level}")
  dir.create(output_dir, recursive = T)
  
  df_abd <- ps[[level]] %>% 
    microbiome::transform('compositional') %>% 
    core(detection = 0, prevalence = 0.1) %>% 
    decoded_abundance()
  
  MMUPHin_models[[glue("{level}_adj")]] <- lm_meta(
    feature_abd = df_abd,
    batch = "cohort",
    exposure = "PD",
    covariates = c("host_age", "sex", "host_body_mass_index"),
    covariates_random = c("cohort", "paired"),
    data = df_meta,
    control = c(
      normalization = "NONE",
      transform = "AST",
      analysis_method = "LM",
      rma_method = "REML",
      output = glue("{output_dir}_adj"),
      verbose = TRUE
    )
  )
}

saveRDS(MMUPHin_models, 
        glue("{MMUPHin_dir}/results_{Sys.Date()}.rds") )
