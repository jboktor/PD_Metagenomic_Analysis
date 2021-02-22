# Dataset summary statistics 

source("src/load_packages.R")
source("src/metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")
load("files/Phyloseq_Merged/Species_PhyloseqObj.RData")

#### Summary Stats ####
## Prep Metadata
dat <- dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]])
reads <- load_reads("Merged") %>% 
  dplyr::select(-cohort)
env <- 
  process_meta(dat, cohort = "Merged") %>%
  left_join(reads, by = "donor_id")
cont_vars <-
  c("host_age", "host_body_mass_index", "clean_total_reads")
prop_vars <-
  c(
    "antibiotics",
    "laxatives",
    "statins",
    "proton_pump_inhibitors",
    "ssri_antidepressants"
  )


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


#--------------------------------------------------------------------------------
# Function to pull counts and ratio of Yes/No variables
#--------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------
# Non-PD Medications
#--------------------------------------------------------------------------------
antibiotics <- prop_summary("antibiotics")
laxatives <- prop_summary("laxatives")
statins <- prop_summary("statins")
proton_pump_inhibitors <- prop_summary("proton_pump_inhibitors")
ssri_antidepressants <- prop_summary("ssri_antidepressants")
amantadine <- prop_summary("amantadine")

prop_df <-
  antibiotics %>%
  full_join(laxatives) %>%
  full_join(statins) %>%
  full_join(proton_pump_inhibitors) %>%
  full_join(ssri_antidepressants)
core_stats <- 
  prop_df %>% 
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
core_stats

# ggsave(core_stats, filename = "data/Community_Composition/Quality_Control/core_stats.svg",
#        width = 10, height = 6)


#--------------------------------------------------------------------------------
# PD Medications
#--------------------------------------------------------------------------------
levodopa <- prop_summary("levodopa")
carbidopa <- prop_summary("carbidopa")
dopamine_agonists <- prop_summary("dopamine_agonists")
MAO_B_inhibitors <- prop_summary("MAO_B_inhibitors")
rasagiline <- prop_summary("rasagiline")
selegiline <- prop_summary("selegiline")
selegiline <- prop_summary("mirapex")


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




# -------------------------------------------------------------------------------
# Pre-processed BMI values
# -------------------------------------------------------------------------------
env <- 
  meta(dat) %>%
  # process_meta(dat, cohort = "Merged") %>%
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



