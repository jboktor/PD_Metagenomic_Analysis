
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")

# Reading in data ----
enrichment_stats <- readRDS("data/Enrichment_Analysis/enrichment_stats_2022-10-05.rds")

enrichment_pdpc.PD.DOWN <-
  enrichment_stats %>% filter(comparison == "PD/PC", direction == "PD Depleted")
enrichment_pdpc.PD.UP <-
  enrichment_stats %>% filter(comparison == "PD/PC", direction == "PD Enriched")
enrichment_pdhc.PD.DOWN <-
  enrichment_stats %>% filter(comparison == "PD/HC", direction == "PD Depleted")
enrichment_pdhc.PD.UP <-
  enrichment_stats %>% filter(comparison == "PD/HC", direction == "PD Enriched")

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
#                 filename = glue("data/Enrichment_Analysis/Kegg_Metabolism_bubbleplot__{Sys.Date()}"),
#                 plot_w = 8, plot_h = 7, leg_w = 9, leg_h = 9)
ggsave(metabolism_fig_plot,
       filename = glue("data/Enrichment_Analysis/Kegg_Metabolism_bubbleplot_{Sys.Date()}.svg"),
       width = 9.5, height = 6.5
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
  filename = glue("data/Enrichment_Analysis/Kegg_all_others_bubbleplot_{Sys.Date()}.svg"),
  plot_w = 6, plot_h = 7, leg_w = 9, leg_h = 9
)
ggsave(kegg_fig_plot,
       filename = glue("data/Enrichment_Analysis/Kegg_all_others_bubbleplot_{Sys.Date()}.svg"),
       width = 11,
       height = 6)
