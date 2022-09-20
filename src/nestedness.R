# NESTEDNESS ANALYSIS

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
load("files/low_quality_samples.RData")

remove_dats()
load_all_cohorts()

datObj <- dat.genus %>% 
  subset_samples(donor_id %ni% low_qc[[1]])
rich <- 
  microbiome::alpha(abundances(datObj), 'observed') %>% 
  rownames_to_column(var = "donor_id") 
meta_all <- process_meta(datObj, cohort = "Merged")
prevalence <-
  tibble::enframe(prevalence(datObj)) %>% 
  dplyr::rename(features = name, prevalence = value)
taxon_fill <- tax_table(datObj)[, "Order"] %>% 
  as.data.frame() %>% 
  rownames_to_column("features") %>% 
  dplyr::rename(higher_lev = Order)

richness_df <- 
  datObj %>% 
  abundances %>% 
  binarize() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "donor_id") %>% 
  left_join(rich, by = "donor_id") %>% 
  left_join(meta_all, by = "donor_id") %>% 
  pivot_longer(!c(donor_id, observed, colnames(meta_all)), 
               names_to = "features",
               values_to = "detection") %>% 
  left_join(prevalence, by = "features") %>% 
  left_join(taxon_fill, by = "features") %>% 
  mutate(detection = as.factor(detection)) %>% 
  dplyr::relocate(features:higher_lev, .after = donor_group)


richA <- 
  richness_df %>% 
  mutate(donor_id = fct_reorder(donor_id, observed)) %>%
  mutate(features = fct_reorder(features, prevalence)) %>%
  ggplot(aes(x=donor_id, y=features, fill=detection)) +
  geom_tile() +
  facet_wrap(cohort~donor_group, scales = "free") +
  theme_classic() +
  labs(x = "Samples (ranked by richness)", y = "Genera (ranked by prevalence)") +
  scale_fill_manual(values = c("0" = "#f1f1f1", "1" = "#434343")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none")

ggsave(richA,
       filename = "data/Community_Composition/Nestedness/Nestedness_cohort_group.png",
       width = 8,
       height = 8)


# Plot all samples
richB <- 
  richness_df %>% 
  mutate(donor_id = fct_reorder(donor_id, observed)) %>%
  mutate(features = fct_reorder(features, prevalence)) %>%
  ggplot(aes(x=donor_id, y=features, fill=detection)) +
  geom_tile() +
  theme_classic() +
  labs(x = "Samples (ranked by richness)", y = "Species (ranked by prevalence)") +
  scale_fill_manual(values = c("0" = "#f1f1f1", "1" = "#434343")) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none") 

ggsave(richB,
       filename = "data/Community_Composition/Nestedness/Nestedness_allsamples.png",
       width = 7,
       height = 4)









