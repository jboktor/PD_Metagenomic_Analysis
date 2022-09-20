
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
# source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_seqdepthfilter5mil.RData")
# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.RData")
# aucs <- readRDS("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_Rarefied.RData")

# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")
# openxlsx::write.xlsx(aucs, file = 'data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.xlsx')
aucs <- readRDS("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.rds")


aurocs_map <- aucs %>%
  select(feature, data_level) %>%
  distinct()

aucs_shared_loose <- aucs %>%
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    mean_auroc = mean(auroc),
    mean_auroc_cent = mean(auroc_center),
    median = median(auroc),
    n = n()
  ) %>%
  left_join(aurocs_map, by = "feature") %>%
  dplyr::group_by(data_level) %>%
  # filter(mean > 0.60 | mean < 0.40, n > 2) %>%
  slice_max(order_by = abs(mean_auroc_cent), n = 200)

aucplot_errorbars_loose <- aucs %>%
  filter(feature %in% aucs_shared_loose$feature) %>%
  # filter(data_level == "Species") %>%
  ggplot(aes(
    x = fct_reorder(feature, desc(auroc)), y = auroc, fill = study,
    ymin = ci_lower, ymax = ci_upper, group = feature
  )) +
  geom_pointrange(aes(group = feature),
    position = position_jitterdodge(jitter.height = 0),
    shape = 24, stroke = 0.2, colour = "grey"
  ) +
  facet_wrap(~data_level, scales = "free_y") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
  my_clean_theme() +
  scale_fill_d3() +
  theme(axis.text.y = element_blank()) +
  coord_flip()
aucplot_errorbars_loose
# ggsave(aucplot_errorbars_loose,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_top200_facets.png",
#        width = 8, height = 8, dpi = 600)





aucs_shared <- aucs %>%
  dplyr::mutate(auroc_center = auroc - 0.5) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    mean_auroc = mean(auroc),
    mean_auroc_cent = mean(auroc_center),
    median = median(auroc),
    n = n()
  ) %>%
  left_join(aurocs_map, by = "feature") %>%
  dplyr::group_by(data_level) %>%
  # filter(mean > 0.60 | mean < 0.40, n > 2) %>%
  slice_max(order_by = abs(mean_auroc_cent), n = 25)


aucplot_errorbars <-
  aucs %>%
  filter(feature %in% aucs_shared$feature) %>%
  mutate(feature = gsub("_", " ", feature)) %>%
  # filter(data_level == "Species") %>%
  ggplot(aes(
    x = fct_reorder(feature, auroc), y = auroc, fill = study,
    ymin = ci_lower, ymax = ci_upper, group = feature
  )) +
  geom_pointrange(aes(group = feature),
    position = position_jitterdodge(jitter.height = 0),
    shape = 24, stroke = 0.2, colour = "grey"
  ) +
  facet_grid(rows = vars(data_level), scales = "free_y", space = "free") +
  # facet_wrap(~data_level, scales = "free_y", ncol = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, colour = "grey", alpha = 0.8) +
  labs(y = "AUROC", x = NULL) +
  my_clean_theme() +
  scale_fill_d3() +
  coord_flip()
aucplot_errorbars
ggsave(aucplot_errorbars,
  filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_facets_top25_new.svg",
  width = 10, height = 40, limitsize = FALSE
)


aucplot_ridges <-
  aucs %>%
  filter(feature %in% aucs_shared$feature) %>%
  # filter(data_level == "PFAMs.slim") %>%
  ggplot(aes(x = auroc, y = fct_reorder(feature, desc(auroc)))) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey", alpha = 0.8) +
  geom_density_ridges_gradient(aes(fill = stat(x))) +
  scale_fill_viridis_c(name = "AUROC", option = "C") +
  facet_wrap(~data_level, scales = "free", ncol = 1) +
  labs(x = "AUROC", y = NULL) +
  my_clean_theme()
aucplot_ridges

# ggsave(aucplot_ridges,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Ridgeline_facets_top25.png",
#        width = 10, height = 45, dpi = 600)


aucplot_errorbars <-
  aucs %>%
  filter(feature %in% aucs_shared$feature) %>%
  # filter(data_level %in%  c("Species", "KOs.slim")) %>%
  # mutate(data_level = factor(data_level, levels = c("Species", "KOs.slim"))) %>%
  ggplot(aes(
    x = fct_reorder(feature, auroc), y = auroc, fill = study,
    ymin = ci_lower, ymax = ci_upper, group = feature
  )) +
  geom_pointrange(aes(group = feature),
    position = position_jitterdodge(jitter.height = 0),
    shape = 24, stroke = 0.1, colour = "grey"
  ) +
  facet_wrap(~data_level, scales = "free", ncol = 1) +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
  labs(y = "AUROC", x = NULL) +
  my_clean_theme() +
  scale_fill_d3() +
  coord_flip()
aucplot_errorbars
ggsave(aucplot_errorbars,
  filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_facets.svg",
  width = 10, height = 45
)

#
#
# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("Species", "KOs.slim")) %>%
#   mutate(data_level = factor(data_level, levels = c("Species", "KOs.slim"))) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.1, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_Species&KOs.svg",
#        width = 14, height = 7, dpi = 600)
#
# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("Species", "KOs.slim")) %>%
#   mutate(data_level = factor(data_level, levels = c("Species", "KOs.slim"))) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.1, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_Pathwayss&KOs.svg",
#        width = 14, height = 7, dpi = 600)
#
#
#
# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("Species")) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.2, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, colour="orange", alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_Species.svg",
#        width = 10, height = 7, dpi = 600)

# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("KOs.slim")) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.2, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_KOs.svg",
#        width = 10, height = 7)
#
# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("Pathways.slim")) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.2, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, colour="orange", alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_Pathways.svg",
#        width = 10, height = 7, dpi = 600)
#
# aucplot_errorbars <-
#   aucs %>%
#   filter(feature %in% aucs_shared$feature) %>%
#   filter(data_level %in%  c("Pathways.slim")) %>%
#   ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.2, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", nrow = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, colour="orange", alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars
# ggsave(aucplot_errorbars,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_AUROCs_top25_Pathways.svg",
#        width = 10, height = 7, dpi = 600)

# Best features of each cohort ----

# top_aucs_per_cohort_datatype <- aucs %>%
#   dplyr::mutate(auroc_center = auroc - 0.5) %>%
#   dplyr::group_by(study, data_level) %>%
#   slice_max(order_by = abs(auroc_center), n = 10)
#
#
# aucplot_errorbars_top_per_cohort <-
#   aucs %>%
#   filter(feature %in% top_aucs_per_cohort_datatype$feature) %>%
#   filter(data_level == "Pathways.slim") %>%
#   ggplot(aes(x=fct_reorder(feature, desc(auroc)), y = auroc, fill = study,
#              ymin=ci_lower, ymax=ci_upper, group = feature)) +
#   geom_pointrange(aes(group = feature),
#                   position = position_jitterdodge(jitter.height = 0),
#                   shape = 24, stroke = 0.2, colour="grey") +
#   facet_wrap(~data_level, scales = "free_y", ncol = 1) +
#   geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
#   labs(y = "AUROC", x = NULL) +
#   my_clean_theme() +
#   scale_fill_d3() +
#   coord_flip()
# aucplot_errorbars_top_per_cohort
# ggsave(aucplot_errorbars_top_per_cohort,
#        filename = "data/Machine_Learning_Analysis/feature_AUROCs/Errorbar_facets_Pathways.png",
#        width = 12, height = 7, dpi = 600)
