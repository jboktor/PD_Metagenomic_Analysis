rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")

#----------------------------------------------
#  Plot Seq Depth Distributions by group ----
#----------------------------------------------

cols.pdpchc <- c(
  "PD" = "#bfbfbf",
  "PC" = "#ed7d31",
  "HC" = "#5b9bd5"
)
cols.pdpchc_dark <- c(
  "PD" = "#494949",
  "PC" = "#ed7d31",
  "HC" = "#5b9bd5"
)
cols.pdpchc.rim <- c(
  "PD" = "#494949",
  "PC" = "#c15811",
  "HC" = "#2e75b5"
)

negative_controls <- c(
  "S00A4-ATCC_MSA_1003_S96",
  "S00A4-neg2_S119",
  "S00A4-neg3_S125",
  "S00A4-neg_S118",
  "S00A4NegExt_P00A4_S94",
  "S00A4NegH2O_P00A4_S95",
  "S00A4_stagPos_S117",
  "BLANK"
)

TBC_keys <- read.csv(file = "files/metadata_keys.csv", header = TRUE) %>%
  dplyr::select(c(MBI_Sample_ID, id)) %>%
  mutate(id = gsub("_", ".", id)) %>%
  mutate(MBI_Sample_ID = as.character(MBI_Sample_ID)) %>%
  mutate(id = as.character(id)) %>%
  dplyr::rename(`# samples` = MBI_Sample_ID)

RUSH_keys <- read.csv(file = "files/metadata_phyloseq_RUSH.csv", header = TRUE) %>%
  dplyr::filter(study_group == "PD") %>%
  dplyr::select(donor_id, host_subject_id) %>%
  dplyr::mutate(donor_id = as.character(donor_id)) %>%
  dplyr::mutate(host_subject_id = as.character(host_subject_id)) %>%
  dplyr::rename(`# samples` = host_subject_id)

knead_reads_TBC <-
  read_tsv(
    "files/TBC_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
    col_names = T
  ) %>%
  dplyr::rename("# samples" = "Sample")
func_reads_TBC <-
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
reads_TBC <-
  func_reads_TBC %>%
  left_join(knead_reads_TBC, by = "# samples") %>%
  dplyr::filter(`# samples` %ni% negative_controls) %>%
  dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>%
  left_join(TBC_keys, by = "# samples") %>%
  dplyr::rename("donor_id" = "id", "clean_total_reads" = "total reads") %>%
  dplyr::mutate(cohort = "TBC")

knead_reads_RUSH <-
  read_tsv(
    "files/RUSH_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
    col_names = T
  ) %>%
  dplyr::rename("# samples" = "Sample")
func_reads_RUSH <-
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
reads_RUSH <-
  func_reads_RUSH %>%
  left_join(knead_reads_RUSH, by = "# samples") %>%
  dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>%
  dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>%
  dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>%
  left_join(RUSH_keys, by = "# samples") %>%
  dplyr::rename("clean_total_reads" = "total reads") %>%
  dplyr::mutate(cohort = "RUSH")

reads <- rbind(reads_TBC, reads_RUSH)
df.reads <- group_col_from_ids(reads, reads$donor_id)
df.reads$group <- factor(df.reads$group, levels = c("PC", "PD", "HC"))

#--------------------------------------------------------------------
#           Determine Samples with low quality reads ----
#--------------------------------------------------------------------

low_qc <-
  df.reads %>%
  filter(clean_total_reads < 1000000 | is.na(clean_total_reads)) %>%
  select(donor_id)

save(low_qc, file = "files/low_quality_samples.RData")

#--------------------------------------------------------------------
#        Remove low-read & ABX usage samples ----
#--------------------------------------------------------------------

ps <- readRDS(file = "files/Phyloseq_Merged/PhyloseqObj_slim.rds")
ps_clean <- ps %>% purrr::map(~.x %>% 
                                subset_samples(donor_id %ni% low_qc[[1]]) %>% 
                                subset_samples(antibiotics != "Yes"))
saveRDS(ps_clean, file = "files/Phyloseq_Merged/PhyloseqObj_slim_clean.rds")

#--------------------------------------------------------------------
#                          Plotting
#--------------------------------------------------------------------

# Total Clean Reads Per Sample - Distributions
histo_plot <-
  ggplot(df.reads, aes(x = clean_total_reads, color = group), alpha = 0.3) +
  theme_bw() +
  geom_density() +
  facet_wrap(~cohort) +
  scale_color_manual(values = cols.pdpchc) +
  theme(
    axis.title.x = element_blank(),
    legend.position = c(0.9, 0.5),
    panel.grid = element_blank()
  )

ecdf_plot <-
  ggplot(df.reads, aes(x = clean_total_reads, colour = group)) +
  stat_ecdf(geom = "step", pad = FALSE, alpha = 0.5) +
  stat_ecdf(geom = "point", pad = FALSE, alpha = 0.9, size = 0.75) +
  theme_bw() +
  facet_wrap(~cohort) +
  labs(x = "Clean Reads", y = "ECDF") +
  scale_color_manual(values = cols.pdpchc, guide = FALSE) +
  theme(
    legend.position = c(0.9, 0.5),
    panel.grid = element_blank()
  )
c1 <- cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align = "v")
c1


# Total Clean Reads Per Sample - Boxplot
ggplot(df.reads, aes(x = group, y = clean_total_reads)) +
  geom_point(aes(fill = group, color = group), position = position_jitterdodge(dodge.width = 1), shape = 21, size = 1.25, alpha = 1) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.3, width = 0.45) +
  theme_minimal() +
  labs(y = "Total Clean Reads per Sample") +
  scale_fill_manual(values = cols.pdpchc) +
  scale_color_manual(values = cols.pdpchc.rim) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )

# Total Clean Reads Per Sample - stratified by cohort
ggplot(df.reads, aes(x = group, y = clean_total_reads)) +
  geom_point(aes(fill = group, color = group), position = position_jitterdodge(dodge.width = 1), shape = 21, size = 1.25, alpha = 1) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.3, width = 0.45) +
  theme_minimal() +
  labs(y = "Total Clean Reads per Sample") +
  scale_fill_manual(values = cols.pdpchc) +
  scale_color_manual(values = cols.pdpchc.rim) +
  facet_wrap(~cohort) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank()
  )

#_______________________________________________________________________________
# Exploring sample richness across data-levels


get_list_colsums <- function(x, y, ...) {
  x  %>%
    abundances() %>% 
    colSums() %>% 
    as.data.frame() %>%
    rownames_to_column() %>%
    mutate(level = y)
}

read_summary_stats <- 
  map2_dfr(ps_clean, names(ps), get_list_colsums) %>% 
  dplyr::rename("sample_sum" = ".",
                "donor_id" = "rowname")

# Number of samples at each level
read_summary_stats %>% count(level)

process_meta(ps_clean$Species, cohort = "Merged")

# Number of samples at each level
read_summary_stats %>% 
  mutate(level = factor(level, levels = names(ps_clean))) %>% 
  left_join(env) %>% 
  ggplot(aes(x=donor_group, y=sample_sum)) +
  geom_point(aes(color = cohort),
    position = position_jitter(height = 0, width = 0.1),
    alpha = 0.7) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.1) +
  facet_wrap(~level, scales = "free") +
  scale_y_log10() 

