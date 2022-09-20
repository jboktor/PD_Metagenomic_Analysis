# plotting correlation analyses

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

corrs <-
  read.csv(file = "files/correlations_analysis.csv") %>%
  column_to_rownames(var = "X") %>%
  group_by(object_name, metadata) %>%
  mutate(FDR = p.adjust(p, method = "BH")) %>%
  ungroup()

corrs.diet <- corrs %>%
  filter(metadata %in% diet)
corrs.clinical <- corrs %>%
  filter(metadata %in% c(motor_severity_scores_summary, clinical_variables))

#
# clinical_summary <-
#   corrs.clinical %>%
#   filter(FDR < 0.1) %>%
#   arrange(FDR) %>%
#   ggplot(aes(x=metadata, y=fct_reorder(feature, desc(FDR)), fill = rho)) +
#   theme_bw() +
#   geom_point(aes(size = -FDR), shape=21, stroke = 0.2
#   ) +
#   # labs(x = NULL, y = NULL, fill = "FDR",
#   #      size = expression(R^"2"~"(%)") ) +
#   scale_y_discrete(position = "right") +
#   scale_fill_viridis_c(option = "magma", begin = 0, end = 1, direction = -1,
#                        na.value = "transparent") +
#   # scale_x_discrete(labels= c("Eggnogs.slim" = "eggNOGs",
#   #                            "KOs.slim" = "KOs")) +
#   # scale_size(breaks = c(1, 5, 10, 20), labels =  c(1, 5, 10, 20)) +
#   # guides(fill = guide_colourbar(barwidth = 1, barheight = 10)) +
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         legend.position =  "left",
#         axis.text.x = element_text(angle = 45, hjust = 1))
#
# ggsave(clinical_summary, filename = "data/Correlations/clinical_data_summary.png",
#        width = 15, height = 40)





# clinical_summary <-

corrs.clinical %>%
  mutate(featurecomp = paste(metadata, feature, sep = "_")) %>%
  mutate(featurecomp = fct_reorder(featurecomp, rho)) %>%
  ggplot(aes(x = featurecomp, y = rho, color = FDR)) +
  geom_bar(stat = "identity", alpha = 0.2) +
  geom_point(size = 0.1) +
  labs(x = "Ranked Correlations") +
  facet_wrap(~metadata, scales = "free_x") +
  scale_color_viridis_c(
    option = "magma", begin = 0, end = 1, direction = -1,
    na.value = "transparent"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

corrs.clinical %>%
  ggplot(aes(y = fct_reorder(metadata, rho), x = rho)) +
  geom_violin() +
  # geom_beeswarm(aes(fill = FDR, ), size= 0.2,
  #               alpha = 0.5, cex = 0.2, groupOnX=FALSE) +
  theme_bw() +
  labs(x = "Ranked Correlations") +
  # facet_wrap(~object_name, scales = "free_y", ncol = 2) +
  scale_fill_viridis_c(
    option = "magma", begin = 0, end = 1,
    direction = -1,
    na.value = "transparent"
  ) +
  theme(axis.ticks.x = element_blank())




corrs.diet %>%
  mutate(featurecomp = paste(metadata, feature, sep = "_")) %>%
  mutate(featurecomp = fct_reorder(featurecomp, rho)) %>%
  ggplot(aes(x = featurecomp, y = rho, color = FDR)) +
  geom_bar(stat = "identity", alpha = 0.2) +
  geom_point(size = 0.1) +
  labs(x = "Ranked Correlations") +
  facet_wrap(~metadata, scales = "free_x") +
  scale_color_viridis_c(
    option = "magma", begin = 0, end = 1, direction = -1,
    na.value = "transparent"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

corrs.diet %>%
  ggplot(aes(y = fct_reorder(metadata, rho), x = rho)) +
  geom_violin() +
  # geom_beeswarm(aes(fill = FDR, ), size= 0.2,
  #               alpha = 0.5, cex = 0.2, groupOnX=FALSE) +
  theme_bw() +
  labs(x = "Ranked Correlations") +
  # facet_wrap(~object_name, scales = "free_y", ncol = 2) +
  scale_fill_viridis_c(
    option = "magma", begin = 0, end = 1,
    direction = -1,
    na.value = "transparent"
  ) +
  theme(axis.ticks.x = element_blank())


# Heatmaps
h1 <- corr_heatmap(corrs.clinical.species)
h2 <- corr_heatmap(corrs.diet.species)
ggsave(h1,
  filename = "data/Correlations/Species/Clinical/Heatmap_Species_clinical.svg",
  height = 7, width = 7
)
ggsave(h2,
  filename = "data/Correlations/Species/Diet/Heatmap_Species_diet.svg",
  height = 9, width = 9
)



# brewer.pal(name = "Paired", 5)
fillpal <- c(
  "GOs" = "#a6cee3", "KOs" = "#1f78b4", "PFAMs" = "#b2df8a",
  "Pathways" = "#33a02c", "Species" = "#fb9a99"
)

# corrs.diet.sig <-
#   corrs.diet %>%
#   dplyr::filter(FDR < 0.2) %>%
#   group_by(object_name, metadata) %>%
#   dplyr::summarise(mean = mean(rho), n = n()) %>%
#   ungroup() %>%
#   dplyr::mutate(object_name = as.factor(object_name))
#
# corrs.diet.sig.barplot <-
#   corrs.diet.sig %>%
#   ggplot(aes(x = log10(n),
#              y = fct_reorder(metadata, n),
#              fill = object_name)) +
#   geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
#   theme_bw() +
#   labs(y = NULL, x = expression(paste("log "[10]*" (FDR < 0.2)")),
#        fill = "Feature level", title = "Dietary Correlations") +
#   scale_fill_manual(values = fillpal) +
#   theme(
#     panel.grid = element_blank(),
#     strip.background = element_rect(fill = "white"),
#     axis.title.y = element_blank(),
#     plot.title = element_text(hjust = 0.5))
#
# ggsave(corrs.diet.sig.barplot, filename = "data/Correlations/Dietary_summary.svg",
#        width = 6.5, height = 5)
#
# corrs.clinical.sig <-
#   corrs.clinical %>%
#   dplyr::filter(FDR < 0.2) %>%
#   group_by(object_name, metadata) %>%
#   dplyr::summarise(mean = mean(rho), n = n()) %>%
#   ungroup() %>%
#   dplyr::mutate(object_name = as.factor(object_name))
#
# corrs.clinical.sig.barplot <-
#   corrs.clinical.sig %>%
#   ggplot(aes(x = log10(n + 1),
#              y = fct_reorder(metadata, n),
#              fill = object_name)) +
#   geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
#   theme_bw() +
#   labs(y = NULL, x = expression(paste("log"[10]*" (FDR < 0.2)")),
#        fill = "Feature level", title = "Clinical Correlations") +
#   scale_fill_manual(values = fillpal) +
#   theme(
#     panel.grid = element_blank(),
#     strip.background = element_rect(fill = "white"),
#     axis.title.y = element_blank(),
#     plot.title = element_text(hjust = 0.5))




corrs.sig <-
  corrs %>%
  dplyr::filter(FDR < 0.2) %>%
  group_by(object_name, metadata) %>%
  dplyr::summarise(mean = mean(rho), n = n()) %>%
  ungroup() %>%
  dplyr::mutate(category = if_else(
    metadata %in% diet,
    "Dietary Correlations",
    if_else(
      metadata %in% c(motor_severity_scores_summary, clinical_variables),
      "Clinical Correlations",
      "Other"
    )
  )) %>%
  dplyr::mutate(object_name = as.factor(object_name))

corrs.sig.barplot <-
  corrs.sig %>%
  ggplot(aes(
    x = log10(n + 1),
    y = fct_reorder(metadata, n),
    fill = object_name
  )) +
  geom_bar(size = 10^-1, stat = "identity", width = 0.6, color = "black") +
  theme_bw() +
  labs(
    y = NULL, x = expression(paste("log"[10] * " (FDR < 0.2)")),
    fill = "Feature level"
  ) +
  scale_fill_manual(values = fillpal) +
  facet_grid(
    row = vars(category), scales = "free_y", space = "free_y",
    # strip.position = "right", space = "free_y"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white"),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(corrs.sig.barplot,
  filename = "data/Correlations/Summary.svg",
  width = 6.5, height = 5
)
