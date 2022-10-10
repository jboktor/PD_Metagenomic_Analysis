## PERMANOVA script

rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_slim_clean_rarefied.rds")

#-------------------------------------------------------------------------------
#####                       PERMANOVA Analysis                             #####
#-------------------------------------------------------------------------------

obj_names <-
  c(
    "Species",
    "Pathways.slim",
    "Enzymes.slim",
    "KOs.slim",
    "Pfams.slim",
    "eggNOGs.slim"
  )

start_time <- Sys.time()
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

permanova_data <- tibble()

study_set <- c("Merged", "Rush", "TBC") 

for (study in study_set) {
  for (level in obj_names) {
    
    if (study == "Merged") { 
      ps <- phyloseq_objs[[level]]
    } else if (study == "Rush") {
      ps <- phyloseq_objs[[level]] %>% 
        subset_samples(cohort == "Rush")
    } else if (study == "TBC") {
      ps <- phyloseq_objs[[level]] %>% 
        subset_samples(cohort == "TBC")
    }
    
    obj_processed <- ps %>%
      microbiome::transform("compositional") %>%
      microbiome::transform("clr")
    sp <- t(abundances(obj_processed))
    # Run Metadata pre-processing function
    env <- process_meta(obj_processed, cohort = study)
    env <- trim_meta(env, min_ratio = 2/nsamples(ps) )
    # remove ID variable and duplicate grouping variable
    env.pa <- dplyr::select(env, -matches(other_metadata))
    ## Adding Alpha Diversity to PERMANOVA
    alpha_diversity <- obj_processed %>%
      abundances() %>%
      microbiome::alpha("shannon")
    env.pa$diversity_shannon <- alpha_diversity$diversity_shannon
    
    permanova_parallel <- foreach(
      var = 1:length(env.pa),
      .packages = c("phyloseq", "magrittr", "microbiome", "Hmisc", "dplyr", "vegan"),
      .combine = "rbind"
    ) %dopar% {
      a <- env.pa[, var]
      a.narm <- na.omit(a)
      
      if (any(is.na(a))) {
        sp.narm <- sp[-attr(a.narm, "na.action"), ]
      } else {
        sp.narm <- sp
      }
      
      meta_ano <- adonis(
        vegdist(sp.narm, method = "euclidean") ~ a.narm,
        permutations = 9999
      )
      
      cbind(
        colnames(env.pa[var]),
        meta_ano$aov.tab[1, ],
        length(a.narm),
        level,
        "Aitchisons",
        study
      )
    }
    permanova_data <- bind_rows(permanova_data, permanova_parallel)
  }
}

stopCluster(cl)
end_time <- Sys.time()
cat(
  "PERMANOVA calculated in : ",
  end_time - start_time,
  attr(end_time - start_time, "units"), "\n"
)

permanova_aitch <- permanova_data %>%
  remove_rownames() %>%
  dplyr::rename(
    "p_value" = "Pr(>F)",
    "vars" = "colnames(env.pa[var])",
    "n_meta" = "length(a.narm)",
    "data_type" = "level",
    "distance" = '"Aitchisons"'
  )

permdf <-
  permanova_aitch %>%
  dplyr::mutate(vars = as.character(vars)) %>%
  dplyr::mutate(R2 = R2 * 100) %>%
  group_by(data_type, study) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  ungroup()


## Add Metadata category column
permdf <-
  mutate(permdf,
    metacat = if_else(
      permdf$vars %in% grouping,
      "Disease-Status/Grouping",
      if_else(
        permdf$vars %in% Anthro,
        "Anthropomorphic",
        if_else(
          permdf$vars %in% GI,
          "GI distress",
          if_else(
            permdf$vars %in% nonstarters,
            "Microbiota disruption",
            if_else(
              permdf$vars %in% allergy,
              "Allergies",
              if_else(
                permdf$vars %in% smoke,
                "Smoking",
                if_else(
                  permdf$vars %in% general_drugs,
                  "Other medication",
                  if_else(
                    permdf$vars %in% pd_drugs,
                    "PD medication",
                    if_else(
                      permdf$vars %in% supplements,
                      "Supplements",
                      if_else(
                        permdf$vars %in% diet,
                        "Diet",
                        if_else(
                          permdf$vars %in% motor_severity_scores_summary,
                          "Motor Symptom Severity",
                          if_else(
                            permdf$vars %in% clinical_variables,
                            "Clinical Variables",
                            if_else(
                              permdf$vars %in% environmental,
                              "Environment",
                              if_else(permdf$vars == "Shannon", "Shannon Diversity",
                                "AA_notmatched"
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )

write.csv(permdf, file = paste0("files/permanova_analysis_", Sys.Date(), ".csv"))
openxlsx::write.xlsx(permdf, file = "files/Supplementary Tables/Table_S2_permanova_analysis.xlsx", overwrite = T)
saveRDS(permdf, paste0("files/permanova_analysis_", Sys.Date(), ".rds"))
