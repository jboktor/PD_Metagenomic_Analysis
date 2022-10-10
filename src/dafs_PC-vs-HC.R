# Comparing PC vs HC as control groups

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
load("files/low_quality_samples.RData")
wkd <- getwd()

phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_slim_clean.rds")

levels <- c("Species", "KOs.slim")
cohort <- "Merged"

variance_thresholds <- list(
  "Species" = 0.1,
  "KOs.slim" = 0.2
)

for (level in levels) {
  ## Filter low QC samples and trim low prevalence feature
  dat.object <- maaslin_prep(phyloseq_objs[[level]])
  
  # PC v HC abundance data
  dat_pchc <- subset_samples(dat.object, donor_group != "PD") %>% 
    microbiome::transform("compositional") %>%
    core(detection = 0, prevalence = 0.1)
  
  # Plot Variance Estimate
  variance_plot(dat_pchc)
  df_input_data_pchc <- dat_pchc %>%
    microbiome::transform("compositional") %>%
    variance_filter(filter.percent = variance_thresholds[[level]])

  df_input_metadata_pchc <-
    process_meta(dat_pchc, cohort = cohort) %>%
    dplyr::mutate(description = factor(description, levels = c("Household Control", "Population Control")))
  
  # PC v HC
  fit_data <- Maaslin2(
    input_data = df_input_data_pchc,
    input_metadata = df_input_metadata_pchc,
    output = paste0(wkd, "/data/MaAsLin2_Analysis/", cohort, "/", level, "_HCvPC_maaslin2_output"),
    random_effects = c("cohort", "quadrant_of_residence"),
    fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
    min_prevalence = 0,
    analysis_method = "LM",
    normalization = "NONE",
    transform = "AST",
    cores = 12,
    plot_scatter = TRUE, correction = 
  )
}

#_______________________________________________________________________________

cohort <- "Merged"
maaslin_results <- list()
group_comparisons <-
  list("PDvPC" = "Population Control",
       "PDvHC" = "Household Control",
       "HCvPC" = "Population Control")

for (level in levels) {
  for (comp in names(group_comparisons)) {
    f <- glue("data/MaAsLin2_Analysis/{cohort}/{level}_{comp}_maaslin2_output/all_results.tsv")
    g <- glue("{comp}_{level}")
    maaslin_results[[g]] <- read_tsv(f, col_names = T) %>%
      filter(value == group_comparisons[[comp]]) %>% 
      mutate(comparison = comp,
             data_level = level)
  }
}

maaslin_results_df <- maaslin_results %>% bind_rows()

# Calculating Group-wise AUROC values -----
auroc_data <- tibble()
for (level in levels) {
  ps <- phyloseq_objs[[level]]
  df_pc <- subset_samples(ps, donor_group == "PC") %>% abundances()
  df_hc <- subset_samples(ps, donor_group == "HC") %>% abundances()
  df_pd <- subset_samples(ps, donor_group == "PD") %>% abundances()
  
  auroc_data %<>%
    bind_rows(
      calculate_auroc(case_df = df_pc, control_df = df_pd) %>%
        mutate(comparison = "PDvPC", data_level = level),
      calculate_auroc(case_df =  df_hc, control_df = df_pd) %>%
        mutate(comparison = "PDvHC", data_level = level),
      calculate_auroc(case_df = df_pc, control_df = df_hc) %>%
        mutate(comparison = "HCvPC", data_level = level)
    )
}

library(ggside)
library(plotly)

group_comp_data <- maaslin_results_df %>% 
  left_join(auroc_data) %>%
  mutate(data_level = factor(data_level, levels = c("Species", "KOs.slim"))) %>% 
  decode_rfriendly_rows(passed_column = "feature") %>% 
  relocate(fullnames) %>% 
  select(-feature) %>% 
  dplyr::rename(feature = fullnames)


#____________________________________________________________________________
# Save PD vs HC data ----

supp_loc <- "files/Supplementary Tables"
openxlsx::write.xlsx(group_comp_data,
                     file = glue("{supp_loc}/Table_S9 Household-vs-Population-Control-Comparisons_{Sys.Date()}.xlsx"))

#____________________________________________________________________________
# plot group comparisons ----
p <-
group_comp_data %>% 
  arrange(desc(qval)) %>% 
  ggplot(aes(x=auroc, y=coef )) +
  geom_vline(xintercept = c(0.4, 0.6), linetype = "dashed") +
  geom_point(aes(color = qval, group = feature), alpha = 1) +
  labs(x = "AUROC",
       y = "MaAsLin2 Model Coefficient", 
       color = "MaAsLin2\nq-value") +
  theme_bw() +
  facet_grid(data_level ~ comparison, scales = "free_y") +
  theme(panel.grid = element_blank())
p


ggsave(p, filename = glue("data/DAF_Analysis/Merged/Figure_2S_{Sys.Date()}_HC-vs-PC_covar_adj.svg"),
       width = 8, height = 3.5)
