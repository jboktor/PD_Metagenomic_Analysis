# Machine Learning Plotting functions


cohort_comparison_bars <- function(cohort_summary){
  bars <- 
    cohort_summary %>% 
    mutate(model_perf = coalesce(mean, .estimate)) %>% 
    filter(.metric == "roc_auc") %>% 
    dplyr::mutate(model_type = if_else((train == test), "Bootstrapped \nTraining Set", "Excluded Study\nPrediction")) %>%
    pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
    ggplot(aes(train, roc_auc, fill = model_type)) +
    theme_classic() +
    expand_limits(y=c(0, 1)) +
    labs(x = NULL, y = "AUROC", fill = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    geom_col(position=position_dodge(), color = "black", alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    scale_fill_manual(values = c("Bootstrapped \nTraining Set" = "#fafafa", "Excluded Study\nPrediction" = "#404040")) +
    theme(panel.grid = element_blank(),
          legend.position = c(0.5, 0.9),
          legend.direction = "horizontal",
          axis.text.x = element_text(angle = 45, hjust = 1))
  return(bars)
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


ml_heatmap_summary <- function(ml_summary){
  ml_summary %>% 
    mutate(model_perf = coalesce(mean, .estimate)) %>% 
    filter(.metric == "roc_auc") %>% distinct() %>% 
    pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
    ggplot(aes(x=train, y=test, fill = roc_auc)) +
    theme_bw() +
    geom_tile() +
    labs(x = "Training set", y = "Test set", fill = "AUROC") +
    geom_text(aes(label=round(roc_auc, digits = 2)), size=3, 
              vjust = 0.77, color = "white") +
    scale_fill_viridis_c(option = "magma", begin = .9, end = 0, 
                         na.value = "transparent") +
    guides(fill = guide_colourbar(barwidth = 1, barheight = 7)) +
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}
  