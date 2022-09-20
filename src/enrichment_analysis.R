# Functional Enrichment Analysis

# rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")
load("files/low_quality_samples.RData")
# load("files/Phyloseq_Merged_ML.RData") #phyloseq_objs


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
# openxlsx::write.xlsx(supplement_table, file = 'files/Supplementary Tables/Table_S4 Pathway Enrichment Analysis Directional_new.xlsx')



enrichment_pdpc_trim <- bind_rows(enrichment_pdpc.PD.UP, enrichment_pdpc.PD.DOWN) %>%
  dplyr::rename("enrichment_pdpc" = "enrichment", "k_pdpc" = "k") %>%
  select(path, enrichment_pdpc, k_pdpc, direction)
enrichment_pdhc_trim <- bind_rows(enrichment_pdhc.PD.UP, enrichment_pdhc.PD.DOWN) %>%
  dplyr::rename("enrichment_pdhc" = "enrichment", "k_pdhc" = "k") %>%
  select(path, enrichment_pdhc, k_pdhc, direction)

dfplot <- enrichment_pdpc_trim %>%
  full_join(enrichment_pdhc_trim, by = c("path", "direction")) %>%
  dplyr::mutate(km = rowMeans(select(., c(k_pdpc, k_pdhc))))

dfplot_annot <- KO_BRITE_keys %>%
  left_join(dfplot, by = "path")

dfplot_annot <- dfplot_annot %>%
  mutate(color_col = if_else(
    enrichment_pdpc > -log10(.05) |
      enrichment_pdhc > -log10(.05),
    level_1, "NS"
  ))


# dfplot_annot %>%
#   filter(hierarchy %ni% c("09150 Organismal Systems", "09160 Human Diseases")) %>%
#   ggplot(aes(x = enrichment_pdpc, y = enrichment_pdhc)) +
#   geom_point(aes(size = km), shape = 21, alpha = 0.85) +
#   # scale_fill_manual(values = colormash) +
#   geom_text_repel(aes(label = path)) +
#   facet_wrap(~hierarchy) +
#   theme_bw() +
#   theme(panel.grid = element_blank())


metabolism_fig <- dfplot_annot %>%
  filter(hierarchy == "09100 Metabolism")
metabolism_fig.s <- metabolism_fig %>%
  filter(enrichment_pdpc > -log10(.05) | enrichment_pdhc > -log10(.05))
metabolism_fig.lab <- metabolism_fig %>%
  dplyr::mutate(avg_sig = rowMeans(select(., c(enrichment_pdhc, enrichment_pdpc)))) %>%
  group_by(direction) %>%
  slice_max(order_by = avg_sig, n = 7, with_ties = FALSE)

metabolism_fig_plot <-
  metabolism_fig %>%
  ggplot() +
  geom_point(aes(y = enrichment_pdhc, x = enrichment_pdpc, size = km),
    fill = "#525252", alpha = 0.2
  ) +
  geom_point(
    data = metabolism_fig.s,
    aes(
      y = enrichment_pdhc, x = enrichment_pdpc, size = km,
      fill = level_1
    ),
    alpha = 0.9, shape = 21, stroke = 0.25
  ) +
  scale_fill_manual(values = colormash) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  # scale_y_continuous(limits = c(0, 8)) +
  geom_text_repel(
    data = metabolism_fig.lab,
    aes(y = enrichment_pdhc, x = enrichment_pdpc, label = path),
    segment.color = "#e0e0e0", seed = 42,
    segment.curvature = 0.5,
    nudge_y = 1,
    direction = "y",
    force = 6,
    size = 2.75
  ) +
  facet_wrap(~direction, scales = "free", ncol = 1) +
  labs(
    x = "Enrichment Scores (PD / PC)", y = "Enrichment Scores (PD / HC)",
    fill = "KEGG Hierarchy", size = "Number of KOs with \nq<0.25 in group"
  ) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid = element_blank())
metabolism_fig_plot

# save_me_cleanly(ggobj = metabolism_fig_plot,
#                 filename = "data/Enrichment_Analysis/Kegg_Metabolism_bubbleplot_new2",
#                 plot_w = 8, plot_h = 7, leg_w = 9, leg_h = 9)
ggsave(metabolism_fig_plot,
  filename = "data/Enrichment_Analysis/Kegg_Metabolism_bubbleplot_new2.svg",
  width = 9, height = 6
)


kegg_fig <- dfplot_annot %>%
  filter(hierarchy != "09100 Metabolism") # %>%
# filter(grepl("BR:ko", path))
kegg_fig.s <- kegg_fig %>%
  filter(enrichment_pdpc > -log10(.05) | enrichment_pdhc > -log10(.05))
kegg_fig.lab <- kegg_fig %>%
  dplyr::mutate(avg_sig = rowMeans(select(., c(enrichment_pdhc, enrichment_pdpc)))) %>%
  group_by(direction) %>%
  slice_max(order_by = avg_sig, n = 7, with_ties = FALSE)

kegg_fig_plot <-
  kegg_fig %>%
  ggplot() +
  geom_point(aes(x = enrichment_pdhc, y = enrichment_pdpc, size = km),
    fill = "#525252", alpha = 0.2
  ) +
  geom_point(
    data = kegg_fig.s,
    aes(
      x = enrichment_pdhc, y = enrichment_pdpc, size = km,
      fill = level_1
    ),
    alpha = 0.9, shape = 21, stroke = 0.25
  ) +
  scale_fill_manual(values = c(colormash, colormash[1:3])) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_legend(order = 2)
  ) +
  # scale_y_continuous(limits = c(0, 10)) +
  geom_text_repel(
    data = kegg_fig.lab,
    aes(x = enrichment_pdhc, y = enrichment_pdpc, label = path),
    segment.color = "#e0e0e0",
    segment.curvature = -0.5,
    # angle = 90,
    nudge_y = 1.3,
    direction = "x",
    # vjust = 0.1,
    force = 3,
    size = 2.75
  ) +
  facet_wrap(~direction, scales = "free", ncol = 1) +
  labs(
    y = "Enrichment Scores (PD / PC)", x = "Enrichment Scores (PD / HC)",
    fill = "KEGG Hierarchy", size = "Number of KOs with \nq<0.25 in group"
  ) +
  theme_bw() +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid = element_blank())
kegg_fig_plot

save_me_cleanly(
  ggobj = kegg_fig_plot,
  filename = "data/Enrichment_Analysis/Kegg_all_others_bubbleplot_new2",
  plot_w = 6, plot_h = 7, leg_w = 9, leg_h = 9
)

# ggsave(kegg_fig_plot, filename = "data/Enrichment_Analysis/Kegg_all_others_bubbleplot_new2.svg",
# width = 11, height = 6)







#
# #_______________________________________________________________________________
# #                    Pathway Enrichment Analysis            ----
# #_______________________________________________________________________________
# #                                PD vs PC
# #_______________________________________________________________________________
#
# enrichment_pdpc <- tibble()
# for (path in unique(annot_KOs_pdpc$level_2)) {
#
#   total_lib <- length(unique(annot_KOs_pdpc$KOs))
#   total_sig_lib <- annot_KOs_pdpc %>% filter(qval < 0.25) %>%  pull(KOs) %>% unique() %>% length()
#   path_lib <- annot_KOs_pdpc %>% filter(level_2 ==  path) %>%  pull(KOs) %>% unique() %>% length()
#   path_sig_lib <- annot_KOs_pdpc %>% filter(level_2 ==  path, qval < 0.25) %>%  pull(KOs) %>% unique() %>% length()
#
#   enrichment_value <-
#     enrichment_formula(N = total_lib,
#                        n = total_sig_lib,
#                        m = path_lib,
#                        k = path_sig_lib)
#
#   row2add <-
#     cbind(
#       "comparison" = "PD/PC",
#       "path" = path,
#       "N" = total_lib,
#       "n" = total_sig_lib,
#       "m" = path_lib,
#       "k" = path_sig_lib,
#       "enrichment" = print(enrichment_value)
#     )
#
#   enrichment_pdpc <- rbind(enrichment_pdpc, row2add)
#   cat(path, "enrichment_value: ", enrichment_value, "\n\n")
# }
#
# statvars <- c("N", "n", "m", "k", "enrichment")
# enrichment_pdpc <-
#   enrichment_pdpc %>%
#   mutate(across(all_of(statvars), as.character),
#          across(all_of(statvars), as.numeric))
#
# #_______________________________________________________________________________
# #                                PD vs HC
# #_______________________________________________________________________________
#
# enrichment_pdhc <- tibble()
# for (path in unique(annot_KOs_pdhc$level_2)) {
#
#   total_lib <- length(unique(annot_KOs_pdhc$KOs))
#   total_sig_lib <- annot_KOs_pdhc %>% filter(qval < 0.25) %>%  pull(KOs) %>% unique() %>% length()
#   path_lib <- annot_KOs_pdhc %>% filter(level_2 ==  path) %>%  pull(KOs) %>% unique() %>% length()
#   path_sig_lib <- annot_KOs_pdhc %>% filter(level_2 ==  path, qval < 0.25) %>%  pull(KOs) %>% unique() %>% length()
#
#   enrichment_value <-
#     enrichment_formula(N = total_lib,
#                        n = total_sig_lib,
#                        m = path_lib,
#                        k = path_sig_lib)
#
#   row2add <-
#     cbind(
#       "comparison" = "PD/HC",
#       "path" = path,
#       "N" = total_lib,
#       "n" = total_sig_lib,
#       "m" = path_lib,
#       "k" = path_sig_lib,
#       "enrichment" = print(enrichment_value)
#     )
#
#   enrichment_pdhc <- rbind(enrichment_pdhc, row2add)
#   cat(path, "enrichment_value: ", enrichment_value, "\n\n")
# }
#
# enrichment_pdhc <-
#   enrichment_pdhc %>%
#   mutate(across(all_of(statvars), as.character),
#          across(all_of(statvars), as.numeric))
