# Functional Enrichment Analysis

# rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")
load("files/low_quality_samples.RData")

cohort <- "Merged"
obj.name <- "KOs.slim"

### Read-in MaAsLin2 output
Maas.pd.pc <- read_tsv(paste0(
  "data/MaAsLin2_Analysis/", cohort, "/",
  obj.name, "_PDvPC_maaslin2_output/all_results.tsv"
), col_names = T) %>%
  filter(value == "Population Control") %>%
  decode_rfriendly_rows(passed_column = "feature") %>%
  dplyr::select(-feature) %>%
  mutate(fullnames = substr(fullnames, 1, 6)) %>%
  dplyr::rename("KOs" = "fullnames")

Maas.pd.hc <- read_tsv(paste0(
  "data/MaAsLin2_Analysis/", cohort, "/",
  obj.name, "_PDvHC_maaslin2_output/all_results.tsv"
), col_names = T) %>%
  filter(value == "Household Control") %>%
  decode_rfriendly_rows(passed_column = "feature") %>%
  dplyr::select(-feature) %>%
  mutate(fullnames = substr(fullnames, 1, 6)) %>%
  dplyr::rename("KOs" = "fullnames")

KO_BRITE <- fromJSON("files/enrichment_analysis/ko00001_BRITE.json") %>%
  as.data.frame() %>%
  unnest(children.children, names_repair = "unique") %>%
  unnest(children, names_repair = "unique") %>%
  unnest(children, names_repair = "unique") %>%
  dplyr::rename(
    hierarchy = children.name,
    level_1 = name...3,
    level_2 = name...4,
    level_3 = name...5
  ) %>%
  mutate(KOs = substr(level_3, 1, 6))

annot_KOs_pdpc <- KO_BRITE %>%
  left_join(Maas.pd.pc, by = "KOs")
annot_KOs_pdhc <- KO_BRITE %>%
  left_join(Maas.pd.hc, by = "KOs")

KO_BRITE_keys <- KO_BRITE %>%
  select(hierarchy, level_1, level_2) %>%
  distinct() %>%
  dplyr::rename("path" = "level_2")


# _______________________________________________________________________________
#                    Pathway Enrichment Analysis            ----
# _______________________________________________________________________________
#                                PD vs PC (PD Enriched)
# _______________________________________________________________________________

enrichment_pdpc.PD.UP <- tibble()
for (path in unique(annot_KOs_pdpc$level_2)) {
  total_lib <- length(unique(annot_KOs_pdpc$KOs))
  total_sig_lib <- annot_KOs_pdpc %>%
    filter(qval <= 0.25, coef < 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_lib <- annot_KOs_pdpc %>%
    filter(level_2 == path) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_sig_lib <- annot_KOs_pdpc %>%
    filter(level_2 == path, qval <= 0.25, coef < 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()

  enrichment_value <-
    enrichment_formula(
      N = total_lib,
      n = total_sig_lib,
      m = path_lib,
      k = path_sig_lib
    )

  row2add <-
    cbind(
      "comparison" = "PD/PC",
      "direction" = "PD Enriched",
      "path" = path,
      "N" = total_lib,
      "n" = total_sig_lib,
      "m" = path_lib,
      "k" = path_sig_lib,
      "enrichment" = print(enrichment_value)
    )

  enrichment_pdpc.PD.UP <- rbind(enrichment_pdpc.PD.UP, row2add)
  cat(path, "enrichment_value: ", enrichment_value, "\n\n")
}

statvars <- c("N", "n", "m", "k", "enrichment")
enrichment_pdpc.PD.UP <-
  enrichment_pdpc.PD.UP %>%
  mutate(
    across(all_of(statvars), as.character),
    across(all_of(statvars), as.numeric)
  )

# _______________________________________________________________________________
#                                PD vs PC (PD DEPLTED)
# _______________________________________________________________________________

enrichment_pdpc.PD.DOWN <- tibble()
for (path in unique(annot_KOs_pdpc$level_2)) {
  total_lib <- length(unique(annot_KOs_pdpc$KOs))
  total_sig_lib <-
    annot_KOs_pdpc %>%
    filter(qval <= 0.25, coef > 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_lib <-
    annot_KOs_pdpc %>%
    filter(level_2 == path) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_sig_lib <-
    annot_KOs_pdpc %>%
    filter(level_2 == path, qval <= 0.25, coef > 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()

  enrichment_value <-
    enrichment_formula(
      N = total_lib,
      n = total_sig_lib,
      m = path_lib,
      k = path_sig_lib
    )

  row2add <-
    cbind(
      "comparison" = "PD/PC",
      "direction" = "PD Depleted",
      "path" = path,
      "N" = total_lib,
      "n" = total_sig_lib,
      "m" = path_lib,
      "k" = path_sig_lib,
      "enrichment" = print(enrichment_value)
    )

  enrichment_pdpc.PD.DOWN <- rbind(enrichment_pdpc.PD.DOWN, row2add)
  cat(path, "enrichment_value: ", enrichment_value, "\n\n")
}

statvars <- c("N", "n", "m", "k", "enrichment")
enrichment_pdpc.PD.DOWN <-
  enrichment_pdpc.PD.DOWN %>%
  mutate(
    across(all_of(statvars), as.character),
    across(all_of(statvars), as.numeric)
  )
# _______________________________________________________________________________
# _______________________________________________________________________________
#                                PD vs HC (PD Enriched)
# _______________________________________________________________________________


enrichment_pdhc.PD.UP <- tibble()
for (path in unique(annot_KOs_pdhc$level_2)) {
  total_lib <- length(unique(annot_KOs_pdhc$KOs))
  total_sig_lib <- annot_KOs_pdhc %>%
    filter(qval <= 0.25, coef < 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_lib <- annot_KOs_pdhc %>%
    filter(level_2 == path) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_sig_lib <- annot_KOs_pdhc %>%
    filter(level_2 == path, qval <= 0.25, coef < 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()

  enrichment_value <-
    enrichment_formula(
      N = total_lib,
      n = total_sig_lib,
      m = path_lib,
      k = path_sig_lib
    )

  row2add <-
    cbind(
      "comparison" = "PD/HC",
      "direction" = "PD Enriched",
      "path" = path,
      "N" = total_lib,
      "n" = total_sig_lib,
      "m" = path_lib,
      "k" = path_sig_lib,
      "enrichment" = print(enrichment_value)
    )

  enrichment_pdhc.PD.UP <- rbind(enrichment_pdhc.PD.UP, row2add)
  cat(path, "enrichment_value: ", enrichment_value, "\n\n")
}

enrichment_pdhc.PD.UP <-
  enrichment_pdhc.PD.UP %>%
  mutate(
    across(all_of(statvars), as.character),
    across(all_of(statvars), as.numeric)
  )



# _______________________________________________________________________________
#                                PD vs HC (PD Depleted)
# _______________________________________________________________________________

enrichment_pdhc.PD.DOWN <- tibble()
for (path in unique(annot_KOs_pdhc$level_2)) {
  total_lib <- length(unique(annot_KOs_pdhc$KOs))
  total_sig_lib <- annot_KOs_pdhc %>%
    filter(qval <= 0.25, coef > 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_lib <- annot_KOs_pdhc %>%
    filter(level_2 == path) %>%
    pull(KOs) %>%
    unique() %>%
    length()
  path_sig_lib <- annot_KOs_pdhc %>%
    filter(level_2 == path, qval <= 0.25, coef > 0) %>%
    pull(KOs) %>%
    unique() %>%
    length()

  enrichment_value <-
    enrichment_formula(
      N = total_lib,
      n = total_sig_lib,
      m = path_lib,
      k = path_sig_lib
    )

  row2add <-
    cbind(
      "comparison" = "PD/HC",
      "direction" = "PD Depleted",
      "path" = path,
      "N" = total_lib,
      "n" = total_sig_lib,
      "m" = path_lib,
      "k" = path_sig_lib,
      "enrichment" = print(enrichment_value)
    )

  enrichment_pdhc.PD.DOWN <- rbind(enrichment_pdhc.PD.DOWN, row2add)
  cat(path, "enrichment_value: ", enrichment_value, "\n\n")
}

enrichment_pdhc.PD.DOWN <-
  enrichment_pdhc.PD.DOWN %>%
  mutate(
    across(all_of(statvars), as.character),
    across(all_of(statvars), as.numeric)
  )

# ______________________________________________________________________________

enrichment_pdpc.PD.DOWN %>% glimpse()
enrichment_pdpc.PD.UP %>% glimpse()
enrichment_pdhc.PD.DOWN %>% glimpse()
enrichment_pdhc.PD.UP %>% glimpse()

#### Plotting ----

supplement_table <-
  bind_rows(
    enrichment_pdpc.PD.DOWN,
    enrichment_pdpc.PD.UP,
    enrichment_pdhc.PD.DOWN,
    enrichment_pdhc.PD.UP
  ) %>%
  right_join(KO_BRITE_keys, by = "path")

saveRDS(supplement_table, glue("data/Enrichment_Analysis/enrichment_stats_{Sys.Date()}.rds"))

openxlsx::write.xlsx(supplement_table, 
                     file = glue('files/Supplementary Tables/Table_S4 Pathway Enrichment Analysis Directional_{Sys.Date()}.xlsx'))
