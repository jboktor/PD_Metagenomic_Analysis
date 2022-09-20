# Contamination analysisi

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

knead_RUSH <-
  read.csv(file = "files/RUSH_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv", header = TRUE)
knead_TBC <-
  read.csv(file = "files/TBC_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv", header = TRUE)

# trim_cols <- function(x, cohort) {
#   if (cohort == "TBC"){
#     colnames(x) <- substr(colnames(x), 1, 10); x
#   } else if (cohort == "RUSH") {
#     colnames(x) <- substr(colnames(x), 1, 6); x
#   }
# }

knead_RUSH <- read_tsv(
  "files/RUSH_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
  col_names = T
) %>%
  janitor::clean_names()
knead_TBC <- read_tsv(
  "files/TBC_biobakery_output_slim/kneaddata/merged/kneaddata_read_count_table.tsv",
  col_names = T
) %>%
  janitor::clean_names()

knead_RUSH.human <- knead_RUSH %>%
  mutate(
    human_reads =
      (trimmed_pair1 - decontaminated_hg37dec_v0_1_pair1) +
        (trimmed_pair2 - decontaminated_hg37dec_v0_1_pair2) +
        (trimmed_orphan1 - decontaminated_hg37dec_v0_1_orphan1) +
        (trimmed_orphan2 - decontaminated_hg37dec_v0_1_orphan2),
    total_reads = final_pair1 + final_pair2 + final_orphan1 + final_orphan2,
    human_vs_microbial = human_reads * 100 / total_reads
  ) %>%
  mutate(sample = substr(sample, 1, 6)) %>%
  filter(!grepl("BLANK", sample)) %>%
  mutate(cohort = "RUSH")


mean(knead_RUSH.human$human_reads)
mean(knead_RUSH.human$total_reads)
mean(knead_RUSH.human$human_vs_microbial)

knead_TBC.human <- knead_TBC %>%
  mutate(
    human_reads =
      (trimmed_pair1 - decontaminated_hg37dec_v0_1_pair1) +
        (trimmed_pair2 - decontaminated_hg37dec_v0_1_pair2) +
        (trimmed_orphan1 - decontaminated_hg37dec_v0_1_orphan1) +
        (trimmed_orphan2 - decontaminated_hg37dec_v0_1_orphan2),
    total_reads = final_pair1 + final_pair2 + final_orphan1 + final_orphan2,
    human_vs_microbial = human_reads * 100 / total_reads
  ) %>%
  filter(sample %ni% negative_controls) %>%
  mutate(cohort = "TBC")

mean(knead_TBC.human$human_reads)
mean(knead_TBC.human$total_reads)
mean(knead_TBC.human$human_vs_microbial)




# bind_rows(knead_TBC.human, knead_RUSH.human) %>%
knead_TBC.human %>%
  filter(total_reads > 1000000) %>%
  ggplot(aes(x = total_reads, y = log10(human_reads))) +
  geom_point(aes(color = cohort)) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  scale_color_d3()

# bind_rows(knead_TBC.human, knead_RUSH.human) %>%
knead_TBC.human %>%
  filter(total_reads > 1000000) %>%
  ggplot(aes(x = total_reads, y = log10(human_vs_microbial))) +
  geom_point(aes(color = cohort)) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  scale_color_d3()
