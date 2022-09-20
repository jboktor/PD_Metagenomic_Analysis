# Dataset summary statistics 

source("src/load_packages.R")
source("src/metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")
# load("files/Phyloseq_Merged/Species_PhyloseqObj.RData")
phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")

# For table summary vis
library(gtsummary) 
library(webshot)

dat <- phyloseq_objs[["Species"]]
#_______________________________________________________________________________
# Samples prior to filtering
env <- process_meta(dat, cohort = "Merged")
reads <- meta(dat) %>% select(total_reads, donor_id)
study_metadata <- env %>% left_join(reads, by = "donor_id") %>% 
  mutate(cohort_donor_group = paste0(cohort, "-", donor_group))

# number of pairs after removing low quality samples and abx usage
study_metadata %>% 
  select(paired, cohort) %>% 
  filter(paired != "No") %>% 
  get_dupes(paired) %>% 
  group_by(dupe_count, cohort) %>% 
  select(-paired) %>% 
  table() / 2

# display analyzed sample #'s
study_metadata %>% 
  select(donor_group, cohort) %>% 
  tbl_summary(by = cohort) 

#_______________________________________________________________________________
# Summary tabels ----

cont_pd_vars <-
  c(
    "disease_duration",
    "age_of_onset",
    "bristol_stool_scale",
    "updrs_total",
    "mds_updrs_3_total",
    "mds_updrs_survey_total",
    "hy_stage",
    "smell_score_upsit",
    "olfactory_diagnosis"
  )

med_pd_vars <-
  c(
    "levodopa",
    "carbidopa",
    "dopamine_agonists",
    "MAO_B_inhibitors",
    "rasagiline",
    "selegiline"
  )

# Clinical PD variables ----
study_metadata %>% 
  filter(PD == "Yes") %>%
  select(all_of(cont_pd_vars)) %>% 
  dplyr::mutate_at(vars(all_of(cont_pd_vars)), ~ na_if(., "not provided")) %>%
  mutate_at(vars(all_of(cont_pd_vars)), as.character) %>%
   mutate_at(vars(all_of(cont_pd_vars)), as.numeric) %>%
   tbl_summary(
     missing = "no",
     type = list(cont_pd_vars ~ "continuous")
   ) %>% add_n() %>% bold_labels() %>% 
  as_gt() %>%
  gt::gtsave(filename = "files/Supplementary Tables/Table_S1_Clinical_variables.png")
 

# PD medication variables ----
study_metadata %>% 
  filter(PD == "Yes") %>%
  select(all_of(med_pd_vars)) %>% 
  dplyr::mutate_all(~ na_if(., "not provided")) %>% 
  mutate_at(vars(all_of(med_pd_vars)), as.character ) %>%
  mutate_at(vars(all_of(med_pd_vars)), ~factor(., labels = c("No", "Yes") ) ) %>%
  tbl_summary(missing = "no") %>% 
  add_n() %>% bold_labels() %>%
  as_gt() %>% 
  gt::gtsave(filename = "files/Supplementary Tables/Table_S1_PD_medication_usage.png")



# all donor demographic variables ----
general_vars <-
  c("host_age",
    "sex",
    "host_body_mass_index",
    "total_reads",
    "bristol_stool_scale")
general_vars_cont <-
  c("host_age",
    "host_body_mass_index",
    "total_reads",
    "bristol_stool_scale")
 
meta(dat) %>%
  mutate(cohort_donor_group = paste0(cohort, "-", donor_group)) %>%
  select(all_of(general_vars), cohort_donor_group) %>%
  mutate_at(vars(all_of(general_vars_cont)), as.character) %>%
  mutate_at(vars(all_of(general_vars_cont)), as.numeric) %>%
  tbl_summary(
    by = cohort_donor_group,
    missing = "no",
    type = list(general_vars_cont ~ "continuous")
  ) %>%
  add_p() %>% add_n() %>% bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "files/Supplementary Tables/Table_S1_Demographic_variables.png")

 
 #_______________________________________________________________________________
 
## Age, BMI, READS, N,
summary_A <- 
  env %>%
  group_by(donor_group, cohort) %>%
  summarise_at(vars(cont_vars),
               funs(
                 mean(., na.rm = T),
                 # n = sum(!is.na(.)),
                 se = sd(., na.rm = T) / sqrt(sum(!is.na(.)))
               ))



#_______________________________________________________________________________
# Function to pull counts and ratio of Yes/No variables
#_______________________________________________________________________________

prop_summary <- function(metvar) {
  metvar_col <- rlang::sym(metvar)
  var <-
    group_by(env, donor_group, cohort,!!metvar_col) %>%
    dplyr::summarise(n = n()) %>%
    # na.omit() %>%
    dplyr::mutate(percentage = n / sum(n) * 100) %>%
    dplyr::rename(response = !!metvar) %>%
    dplyr::mutate(metavar = metvar)
  return(var)
}

#_______________________________________________________________________________
# Non-PD Medications
#_______________________________________________________________________________
# prop_vars <-
#   c(
#     "antibiotics",
#     "laxatives",
#     "statins",
#     "proton_pump_inhibitors",
#     "ssri_antidepressants"
#   )

antibiotics <- prop_summary("antibiotics")
laxatives <- prop_summary("laxatives")
suppositories <- prop_summary("suppositories")
statins <- prop_summary("statins")
proton_pump_inhibitors <- prop_summary("proton_pump_inhibitors")
ssri_antidepressants <- prop_summary("ssri_antidepressants")
amantadine <- prop_summary("amantadine")

prop_df <-
  antibiotics %>%
  full_join(laxatives) %>%
  full_join(suppositories) %>%
  full_join(statins) %>%
  full_join(proton_pump_inhibitors) %>%
  full_join(ssri_antidepressants)
core_stats <- 
  prop_df %>% 
  dplyr::mutate(response = replace_na(response, "not provided")) %>% 
  ggplot(aes(x = percentage, y = metavar, fill = response)) + 
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(cohort~donor_group) +
  theme_bw() +
  labs(x = "Percentage %") +
  geom_text(aes(label=n), size=3, color = "white", #nudge_y = -1,
            position = position_stack(vjust = 0.6)) +
  scale_fill_manual(values =
                      c("Yes" = "#800000",
                        "No" = "#03353E",
                        "not provided" = "#959595")) +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank())
core_stats

# ggsave(core_stats, filename = "data/Community_Composition/Quality_Control/core_stats.svg",
#        width = 10, height = 6)

#_______________________________________________________________________________
# PD Medications
#_______________________________________________________________________________
levodopa <- prop_summary("levodopa")
carbidopa <- prop_summary("carbidopa")
dopamine_agonists <- prop_summary("dopamine_agonists")
MAO_B_inhibitors <- prop_summary("MAO_B_inhibitors")
rasagiline <- prop_summary("rasagiline")
selegiline <- prop_summary("selegiline")
mirapex <- prop_summary("mirapex")


pd_drugs <-
  levodopa %>%
  full_join(carbidopa) %>%
  full_join(dopamine_agonists) %>%
  full_join(MAO_B_inhibitors) %>%
  full_join(amantadine)

pd_drug_stats <- 
  pd_drugs %>% 
  dplyr::mutate(response = replace_na(response, value = "not provided")) %>% 
  ggplot(aes(x = percentage, y = metavar, fill = response)) + 
  geom_bar(stat = "identity", width = 0.6) +
  facet_grid(cohort~donor_group) +
  theme_bw() +
  labs(x = "Percentage %") + 
  geom_text(aes(label=n), size=3, color = "white", #nudge_y = -1,
            position = position_stack(vjust = 0.6)) +
  scale_fill_manual(values =
                      c("Yes" = "#800000",
                        "No" = "#03353E",
                        "not provided" = "#959595")) +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank())
pd_drug_stats

# ggsave(pd_drug_stats, filename = "data/Community_Composition/Quality_Control/pd_drug_stats.svg",
#        width = 10, height = 6)




#_______________________________________________________________________________
# Pre-processed BMI values
#_______________________________________________________________________________
env <- 
  dat %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  meta() %>%
  left_join(reads, by = "donor_id")

## BMI
env$host_body_mass_index <- as.character(env$host_body_mass_index)
env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
env %>%
  group_by(donor_group, cohort) %>%
  dplyr::summarise(
    count = n(),
    mean = mean(host_body_mass_index, na.rm = TRUE),
    sd = sd(host_body_mass_index, na.rm = TRUE),
    se = sd(host_body_mass_index, na.rm = T) / sqrt(sum(!is.na(
      host_body_mass_index
    )))
  )

# Sample sex count
env %>% 
  group_by(donor_group, cohort, sex) %>% 
  dplyr::summarise(count = n())



