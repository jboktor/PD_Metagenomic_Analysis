# Machine Learning Plotting functions


cohort_comparison_bars <- function(cohort_summary){
  bars <- 
    cohort_summary %>% 
    mutate(model_perf = coalesce(mean, .estimate)) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::mutate(model_type = if_else((train == test), "10-fold CV", "Prediction")) %>%
    pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
    ggplot(aes(train, roc_auc, color = model_type, fill = test)) +
    theme_classic() +
    labs(x = "Training set", y = "AUROC", fill = "Test") +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    geom_col(position=position_dodge(), alpha = 0.7) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_tableau(palette = "Tableau 10") +
    scale_color_manual(values = c("10-fold CV" = "#000000", "Prediction" = "#b3b3b3")) +
    theme(panel.grid = element_blank())
  bars
}



heatmap_cohort_groups <- function(cohort_group_summary){
  heatmap <- 
    cohort_group_summary %>% 
    mutate(model_perf = coalesce(mean, .estimate)) %>% 
    filter(.metric == "roc_auc") %>% 
    pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
    mutate(roc_auc = if_else((train_cohort == test_cohort) & 
                               (train_group != test_group), 
                             as.double(NA), roc_auc)) %>% 
    dplyr::group_by(train) %>% 
    ggplot(aes(train, test, fill = roc_auc)) +
    theme_bw() +
    geom_tile() +
    labs(x = "Training sets", y = "Testing sets", fill = "AUROC") +
    geom_text(aes(label=round(roc_auc, digits = 2)), size=3, 
              vjust = 0.77, color = "white") +
    scale_fill_viridis_c(option = "magma", begin = .9, end = 0, 
                         na.value = "transparent") +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 15)) +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(heatmap)
  return(heatmap)
}
