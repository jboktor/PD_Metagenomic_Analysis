# Machine Learning - Model Comparison Analysis

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/ml_models.R")
load("files/low_quality_samples.RData")

remove_dats()
load_all_cohorts()

obj <- subset_samples(dat.species, donor_id %ni% low_qc[[1]]) %>% 
  core(detection = 0, prevalence = 0.025) 

obj_tbc_all <- subset_samples(obj, cohort == "TBC")
obj_rush_all <- subset_samples(obj, cohort == "Rush")

cohort_summary <- lasso_cohort_summary(obj_tbc_all, obj_rush_all)

heatmap <- 
  cohort_summary %>% 
  mutate(model_perf = coalesce(mean, .estimate)) %>% 
  filter(.metric == "roc_auc") %>% 
  pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
  ggplot(aes(train, test, fill = roc_auc)) +
  geom_tile() +
  labs(x = "Training sets", y = "Testing sets", fill = "AUROC") +
  geom_text(aes(label=round(roc_auc, digits = 2)), size=3, 
            vjust = 0.77, color = "white") +
  scale_fill_viridis_c(option = "magma", begin = .9, end = 0, 
                       na.value = "transparent") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 9)) +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
heatmap
ggsave(heatmap, filename = "data/Machine_Learning_Analysis/Lasso_Species_Heatmap_Cohort.svg",
       width = 4, height = 3)

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
ggsave(bars, filename = "data/Machine_Learning_Analysis/Lasso_Species_Barplot_Cohort.svg",
       width = 3, height = 3)



cohort_group_summary <- lasso_cohort_x_group_summary(obj_tbc_all, obj_rush_all)

heatmap2 <- 
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
heatmap2
ggsave(heatmap2, filename = "data/Machine_Learning_Analysis/Lasso_Species_Heatmap_CohortGroup.svg",
       width = 5, height = 4)

bars2 <- 
  cohort_group_summary %>% 
  dplyr::mutate(model_perf = coalesce(mean, .estimate)) %>% 
  dplyr::filter(.metric == "roc_auc") %>%
  dplyr::mutate(traintest = paste(train, test, sep = "_")) %>% 
  dplyr::mutate(model_type = if_else((train_cohort == test_cohort) & (train_group == test_group), "10-fold CV" , "Prediction")) %>%
  pivot_wider(names_from = .metric, values_from = "model_perf") %>%
  dplyr::filter(if_else((train_cohort == test_cohort) & (train_group != test_group), F , T)) %>% 
  dplyr::group_by(train) %>%
  dplyr::mutate(ord = row_number()) %>%
  ggplot(aes(train, roc_auc, group = ord, color = model_type, fill = test)) +
  theme_classic() +
  labs(x = "Training set", y = "AUROC", color = "Model", fill = "Test") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  geom_col(position=position_dodge(), alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_tableau(palette = "Tableau 20") +
  scale_color_manual(values = c("10-fold CV" = "#000000", "Prediction" = "#b3b3b3")) +
  theme(panel.grid = element_blank())
bars2
ggsave(bars2, filename = "data/Machine_Learning_Analysis/Lasso_Species_Barplot_CohortGroup.svg",
       width = 5, height = 3)





#-----------------------------------------------------------------------------
#                               Pathways
#-----------------------------------------------------------------------------

obj <- subset_samples(dat.path.slim, donor_id %ni% low_qc[[1]]) %>% 
  core(detection = 0, prevalence = 0.025) 

obj_tbc_all <- subset_samples(obj, cohort == "TBC")
obj_rush_all <- subset_samples(obj, cohort == "Rush")

# feature selection via mRMR
mRMR_pathways <- prep_mRMR_input(obj)
hits_pathways <- feature_selection(mRMR_pathways, n_features = 100, n_algorithm = 1)
length(unique(hits_pathways))

ml_input_tbc <- prep_ml_input(obj_tbc_all) %>% 
  select(donor_id, PD, cohort, paired, all_of(hits_pathways))
ml_input_rush <- prep_ml_input(obj_rush_all) %>%
  select(donor_id, PD, cohort, paired, matches(hits_pathways))

# Note: some features can't be selected after identification via mRMR
cohort_summary_pathways <- lasso_cohort_summary(obj_tbc_all = obj_tbc_all, 
                                                obj_rush_all = obj_rush_all,
                                                featSelection = hits_pathways)

cohort_group_summary_pathways <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
                                                     obj_rush_all = obj_rush_all,
                                                     featSelection = hits_pathways)


barplot_pathways <- cohort_comparison_bars(cohort_summary_pathways)
heatmap_pathways <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_pathways)

ggsave(heatmap_pathways, filename = "data/Machine_Learning_Analysis/Lasso_Pathways_Heatmap_CohortGroup.svg",
       width = 5, height = 4)

# # TROUBLE SHOOTING ERROR
# ml_input_tbc <-
#   prep_ml_input(obj = obj_tbc_all) %>%
#   pivot_longer(!c(donor_id, PD, cohort, paired),
#                names_to = "features", values_to = "nums") %>%
#   dplyr::select(-c(donor_id, PD, cohort, paired, nums)) %>%
#   dplyr::distinct() %>%
#   decode_rfriendly_rows("features")
# ml_input_rush <- prep_ml_input(obj_rush_all) %>%
#   pivot_longer(!c(donor_id, PD, cohort, paired),
#                names_to = "features", values_to = "nums")
# tst <-
#   hits_pathways %>%
#   as.data.frame() %>%
#   rename("features" = ".") %>%
#   decode_rfriendly_rows("features") %>%
#   anti_join(ml_input_tbc, by = "features")
# 
# ml_input_tbc <-
#   prep_ml_input(obj_tbc_all)
# ml_input_tbc2 <- ml_input_tbc %>% 
#   dplyr::select(donor_id, PD, cohort, paired, matches(hits_pathways,ignore.case = T))
# 
# ml_input_tbc3 <-
#   prep_ml_input(obj_tbc_all)


#-----------------------------------------------------------------------------
#                                  KOs
#-----------------------------------------------------------------------------

obj <- subset_samples(dat.KOs.slim, donor_id %ni% low_qc[[1]]) %>% 
  core(detection = 0, prevalence = 0.15) 

obj_tbc_all <- subset_samples(obj, cohort == "TBC")
obj_rush_all <- subset_samples(obj, cohort == "Rush")

# feature selection via mRMR
mRMR_KOs <- prep_mRMR_input(obj)
hits_KOs <- feature_selection(df = mRMR_KOs, n_features = 100, n_algorithm = 1)
length(unique(hits_KOs))

# Note: some features can't be selected after identification via mRMR
cohort_summary_KOs <- lasso_cohort_summary(obj_tbc_all = obj_tbc_all, 
                                                obj_rush_all = obj_rush_all,
                                                featSelection = hits_KOs)

cohort_group_summary_KOs <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
                                                              obj_rush_all = obj_rush_all,
                                                              featSelection = hits_KOs)


barplot_KOs <- cohort_comparison_bars(cohort_summary_KOs)
heatmap_KOs <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_KOs)

ggsave(heatmap_KOs, filename = "data/Machine_Learning_Analysis/Lasso_KOs_Heatmap_CohortGroup.svg",
       width = 5, height = 4)



#-----------------------------------------------------------------------------
#                                  eggNOGs
#-----------------------------------------------------------------------------

# add feature selection process
obj <- subset_samples(dat.EGGNOGs.slim, donor_id %ni% low_qc[[1]]) %>% 
  core(detection = 0, prevalence = 0.15) 

obj_tbc_all <- subset_samples(obj, cohort == "TBC")
obj_rush_all <- subset_samples(obj, cohort == "Rush")

# feature selection via mRMR
mRMR_eggNOGs <- prep_mRMR_input(obj)
hits_eggNOGs <- feature_selection(df = mRMR_eggNOGs, n_features = 100, n_algorithm = 1)
length(unique(hits_eggNOGs))

ml_input_tbc <- prep_ml_input(obj_tbc_all) %>% 
  select(donor_id, PD, cohort, paired, all_of(hits_eggNOGs))
ml_input_rush <- prep_ml_input(obj_rush_all) %>%
  select(donor_id, PD, cohort, paired, matches(hits_eggNOGs))

# Note: some features can't be selected after identification via mRMR
cohort_summary_eggNOGs <- lasso_cohort_summary(obj_tbc_all = obj_tbc_all, 
                                           obj_rush_all = obj_rush_all,
                                           featSelection = hits_eggNOGs)

cohort_group_summary_eggNOGs <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
                                                         obj_rush_all = obj_rush_all,
                                                         featSelection = hits_eggNOGs)


barplot_eggNOGs <- cohort_comparison_bars(cohort_summary_eggNOGs)
heatmap_eggNOGs <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_eggNOGs)






  # -----------------------------------------------------
#   #                TUNE LASSO MODEL  
#   #-----------------------------------------------------
#   # Visualize top contributors to model performance
#   cv_lasso %>%
#     pull_workflow_fit() %>%
#     vi(lambda = best_aucroc$penalty) %>%
#     decode_rfriendly_rows("Variable") %>% select(-Variable) %>%
#     dplyr::rename(Variable = fullnames) %>%
#     mutate(Importance = abs(Importance),
#            Variable = fct_reorder(Variable, Importance)) %>%
#     slice_head(n = 20) %>%
#     ggplot(aes(Importance, Variable, fill = Sign)) +
#     geom_col(width = 0.5) +
#     theme_bw() +
#     scale_x_continuous(expand = c(0,0)) +
#     labs(y = NULL)
# 
# # Finalizing CV model
# final_lasso <-
#   wf %>%
#   finalize_workflow(best_aucroc) %>%
#   fit(ml_input)
#   
# # Visualize top contributors to model performance
# final_lasso %>% 
#   pull_workflow_fit() %>% 
#   vi(lambda = best_aucroc$penalty, ) %>% 
#   decode_rfriendly_rows("Variable") %>% select(-Variable) %>% 
#   dplyr::rename(Variable = fullnames) %>% 
#   mutate(Importance = abs(Importance),
#          Variable = fct_reorder(Variable, Importance)) %>% 
#   slice_head(n = 20) %>% 
#   ggplot(aes(Importance, Variable, fill = Sign)) +
#   geom_col(width = 0.5) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0,0)) +
#   labs(y = NULL)
# 
# train_predict <- 
#   stats::predict(final_lasso, type = "prob", new_data = ml_input, penalty = 1)
# train_probs <- 
#   predict(final_lasso, type = "prob", new_data = ml_input) %>%
#   bind_cols(obs = ml_input$PD) %>%
#   bind_cols(predict(final_lasso, new_data = ml_input))
# 
# conf_mat(train_probs, obs, .pred_class)
# autoplot(roc_curve(train_probs, obs, .pred_Yes, event_level = "second"))
# roc_auc(train_probs, obs, .pred_Yes, event_level = "second")
# 
# 
# #--------------------------------------------------------------------------------
# #                                Predictions  
# #--------------------------------------------------------------------------------
# 
# prediction_lasso <- 
#   wf %>% 
#   finalize_workflow(best_aucroc) %>% 
#   fit(ml_input_rush)
# 
# test_predict <- 
#   stats::predict(prediction_lasso, type = "prob", new_data = ml_input_rush, penalty = 1)
# test_probs <- 
#   predict(prediction_lasso, type = "prob", new_data = ml_input_rush) %>%
#   bind_cols(obs = ml_input_rush$PD) %>%
#   bind_cols(predict(prediction_lasso, new_data = ml_input_rush))
# 
# conf_mat(test_probs, obs, .pred_class)
# autoplot(roc_curve(test_probs, obs, .pred_Yes, event_level = "second"))
# roc_auc(test_probs, obs, .pred_Yes, event_level = "second")

# rocplot <-  ggplot() + 
#   style_roc() +
#   theme_bw() +
#   geom_abline(intercept = 0, slope = 1, color = "grey") +
#   geom_roc(data = test_probs, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
#            n.cuts=0, color = cols.ml[1], linealpha = 0.8) +








# 1) train on TBC, test on Rush


# 2) train on Rush, test on TBC




# # Data  prep
# set.seed(42)
# obj_split <- initial_split(ml_input, strata = cohort)
# obj_train <- training(obj_split)
# obj_test <- testing(obj_split)
# 
# ml_recipe <- recipe(PD ~ ., data = obj_train) %>% 
#   update_role(donor_id, new_role = "donor_id") %>% 
#   update_role(cohort, new_role = "cohort") %>% 
#   update_role(paired, new_role = "paired") %>%
#   step_zv(all_numeric(), -all_outcomes()) %>% 
#   step_normalize(all_numeric(), -all_outcomes())
# 
# ml_prep <- ml_recipe %>% 
#   prep()
# 
# # Predict Disease status with features
# cv_splits <- rsample::vfold_cv(ml_input, strata = PD, v = 5)
# lasso_spec <- logistic_reg(penalty = 0.1, mixture = 1) %>% 
#   set_engine("glmnet")
# 
# wf <- workflow() %>% 
#   add_recipe(ml_recipe) 
# 
# lasso_fit <- wf %>% 
#   add_model(lasso_spec) %>% 
#   fit(data = obj_train)
# 
# lasso_fit %>% 
#   pull_workflow_fit() %>% 
#   tidy()
# 
# # Tune LASSO parameters
# set.seed(42)
# obj_boot <- bootstraps(obj_train, strata = cohort)
# 
# tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
#   set_engine("glmnet")
# 
# lambda_grid <- grid_regular(penalty(),
#              levels = 50)
# 
# doParallel::registerDoParallel()
# 
# set.seed(42)
# lasso_grid <- tune_grid(
#   wf %>% add_model(tune_spec),
#   resamples = obj_boot,
#   grid = lambda_grid)
# 
# 
# # visualize model metrics of grid 
# lasso_grid %>% 
#   collect_metrics() %>% 
#   ggplot(aes(penalty, mean, color = .metric)) +
#   geom_errorbar(aes(ymin = mean - std_err,
#                     ymax = mean + std_err),
#                 alpha = 0.5) +
#   geom_line(size = 1.5, show.legend = F) +
#   facet_wrap(~.metric, scales = "free", nrow = 2) +
#   theme_bw() + 
#   scale_x_log10() +
#   theme(legend.position = "none")
# 
# # select best model parameters
# best_aucroc <- lasso_grid %>% 
#   select_best("roc_auc")
# 
# 
# # Run same workflow with tuned parameters
# final_lasso <- 
#   finalize_workflow(wf %>% add_model(tune_spec),
#                     best_aucroc)
# 
# final_lasso %>% 
#   fit(obj_train) %>% 
#   pull_workflow_fit() %>% 
#   vi(lambda = best_aucroc$penalty) %>% 
#   decode_rfriendly_rows("Variable") %>% select(-Variable) %>% 
#   dplyr::rename(Variable = fullnames) %>% 
#   mutate(Importance = abs(Importance),
#          Variable = fct_reorder(Variable, Importance)) %>% 
#   # top_n(n = 20, wt = Importance) %>% 
#   slice_head(n = 20) %>% 
#   ggplot(aes(Importance, Variable, fill = Sign)) +
#   geom_col(width = 0.5) +
#   theme_bw() +
#   scale_x_continuous(expand = c(0,0)) +
#   labs(y = NULL) 
# 
# 
# # Test data 
# last_fit(final_lasso,
#          obj_split) %>% 
#   collect_metrics()



















  