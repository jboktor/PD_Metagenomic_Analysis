# Machine Learning - Model Comparison Analysis

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/ml_models.R")
source("src/ml_plots.R")
load("files/low_quality_samples.RData")
# load("files/Phyloseq_Merged_ML.RData") #phyloseq_objs
load("files/Phyloseq_Merged_ML_Rarefied.RData") #phyloseq_objs_rare
# 
# tst <- tax_table(phyloseq_objs[["Species"]]) %>% 
#   as.data.frame()

# csv_prefix <- "data/Machine_Learning_Analysis/model_stats/"
# levlist <- c("Species" #,
#              # "Genus",
#              # "Pathways.slim",
#              # "KOs.slim",
#              # "PFAMs.slim",
#              # "eggNOGs.slim"
#              )
# ml_models <- c(
#   "lasso",
#   "randomforest"
#   # "xgboost"
#   )

# # ML Model loop including feature selection
# for (level in levlist){
#   
#   print_line()
#   print(level)
#   
#   obj <- subset_samples(phyloseq_objs_rare[[level]], donor_id %ni% low_qc[[1]]) %>% 
#     microbiome::transform("clr")
#   print(obj)
#   
#   obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
#   obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
#   obj_noRush <- subset_samples(obj, cohort != "Rush")
#   obj_Rush <- subset_samples(obj, cohort == "Rush")
#   obj_noTBC <- subset_samples(obj, cohort != "TBC")
#   obj_TBC <- subset_samples(obj, cohort == "TBC")
#   obj_noBonn <- subset_samples(obj, cohort != "Bonn")
#   obj_Bonn <- subset_samples(obj, cohort == "Bonn")
#   
#   for(model_eng in ml_models){
#     for (N in c(25, 100)){
#       
#       df_loso <-
#         ml_loso(
#           Shanghai = obj_Shanghai,
#           noShanghai = obj_noShanghai,
#           TBC = obj_TBC,
#           noTBC = obj_noTBC,
#           Rush = obj_Rush,
#           noRush = obj_noRush,
#           Bonn = obj_Bonn,
#           noBonn = obj_noBonn,
#           model_type = model_eng,
#           featSelection = F,
#           data_type = level,
#           nfeats =  N
#         )
#       
#       write.csv(df_loso, 
#                 file = paste0(csv_prefix, "LOSO_", level, "_", model_eng, "_", N, 
#                               "Rarefied_CLR.csv"))
#       df_loso.plot <- cohort_comparison_bars(df_loso)
#       print(df_loso.plot)
#       print_line()
#       ggsave(df_loso.plot,
#              filename = paste("data/Machine_Learning_Analysis/LOSO", level, 
#                               model_eng, "barplot", N, "Rarefied_CLR.png", sep = "_"),
#              width = 3.5, height = 4)
#       
#       
#       df_s2s <- 
#         ml_s2s(
#           Shanghai = obj_Shanghai,
#           TBC = obj_TBC,
#           Rush = obj_Rush,
#           Bonn = obj_Bonn,
#           model_type = model_eng,
#           featSelection = F, 
#           data_type = level,
#           nfeats = N
#         )
#       write.csv(df_s2s, 
#                 file = paste0(csv_prefix, "S2S_", level, "_", model_eng, 
#                               "_", N, ".csv"))
#       df_s2s.plot <- ml_heatmap_summary(df_s2s)
#       print(df_s2s.plot)
#       print_line()
#       ggsave(df_s2s.plot, 
#              filename = paste("data/Machine_Learning_Analysis/S2S", level, 
#                               model_eng, "heatmap ", N, "Rarefied_CLR.png", 
#                               sep = "_"),
#              width = 4, height = 3)
#     }
#   }
# }





# load("files/Phyloseq_Merged_ML.RData") #phyloseq_objs
load("files/Phyloseq_Merged_ML_Rarefied.RData") #phyloseq_objs_rare  
csv_prefix <- "data/Machine_Learning_Analysis/model_stats/"
levlist <- c(
  "Species",
  "Genus",
  "Pathways.slim",
  "Enzymes.slim",
  "KOs.slim",
  "PFAMs.slim",
  "eggNOGs.slim"
)
ml_models <- c(
  "lasso" #,
  # "randomforest" #,
  # "xgboost"
  )
  
  # ML Model loop - no feature selection 
for (level in levlist) {
  print_line()
  print(level)
  
  obj <-
    subset_samples(phyloseq_objs_rare[[level]], donor_id %ni% low_qc[[1]])
  print(obj)
  
  obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
  obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
  obj_noRush <- subset_samples(obj, cohort != "Rush")
  obj_Rush <- subset_samples(obj, cohort == "Rush")
  obj_noTBC <- subset_samples(obj, cohort != "TBC")
  obj_TBC <- subset_samples(obj, cohort == "TBC")
  obj_noBonn <- subset_samples(obj, cohort != "Bonn")
  obj_Bonn <- subset_samples(obj, cohort == "Bonn")

  for (model_eng in ml_models) {
    df_loso <-
      ml_loso(
        Shanghai = obj_Shanghai,
        noShanghai = obj_noShanghai,
        TBC = obj_TBC,
        noTBC = obj_noTBC,
        Rush = obj_Rush,
        noRush = obj_noRush,
        Bonn = obj_Bonn,
        noBonn = obj_noBonn,
        model_type = model_eng,
        featSelection = NULL,
        data_type = level
      )
    
    write.csv(df_loso,
              file = paste0(csv_prefix, "LOSO_", level, "_", model_eng, 
                            "_all_features_Rarefied.csv"))
    df_loso.plot <- cohort_comparison_bars(df_loso)
    print(df_loso.plot)
    print_line()
    ggsave(
      df_loso.plot,
      filename = paste(
        "data/Machine_Learning_Analysis/LOSO",
        level,
        model_eng,
        "barplot_all_features_Rarefied_CLR.svg",
        sep = "_"
      ),
      width = 3.5,
      height = 4
    )
    
    
    df_s2s <-
      ml_s2s(
        Shanghai = obj_Shanghai,
        TBC = obj_TBC,
        Rush = obj_Rush,
        Bonn = obj_Bonn,
        model_type = model_eng,
        featSelection = NULL,
        data_type = level
      )
    write.csv(df_s2s,
              file = paste0(csv_prefix, "S2S_", level, "_", model_eng, 
                            "_all_features_Rarefied_CLR.csv"))
    df_s2s.plot <- ml_heatmap_summary(df_s2s)
    print(df_s2s.plot)
    print_line()
    ggsave(
      df_s2s.plot,
      filename = paste(
        "data/Machine_Learning_Analysis/S2S",
        level,
        model_eng,
        "heatmap_all_features_Rarefied_CLR.svg",
        sep = "_"
      ),
      width = 4,
      height = 3
    )
  }
}





















# 
# #_______________________________________________________________________________
# #                               Species                                     ----
# #_______________________________________________________________________________
# 
# obj <- subset_samples(dat.species, donor_id %ni% low_qc[[1]]) %>% 
#   microbiome::transform('compositional')
# 
# sample_data(obj)$cohort_group <- 
#   paste(meta(obj)$cohort, meta(obj)$donor_group, sep = "_")
# 
# obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
# obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
# obj_noRush <- subset_samples(obj, cohort != "Rush")
# obj_Rush <- subset_samples(obj, cohort == "Rush")
# obj_noTBC <- subset_samples(obj, cohort != "TBC")
# obj_TBC <- subset_samples(obj, cohort == "TBC")
# obj_noBonn <- subset_samples(obj, cohort != "Bonn")
# obj_Bonn <- subset_samples(obj, cohort == "Bonn")
# 
# # S2S ------ 
# s2s_species_lasso <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "lasso", 
#     data_type = "Species"
#   )
# s2s_species_lasso.plot <- ml_heatmap_summary(s2s_species_lasso)
# s2s_species_lasso.plot
# ggsave(s2s_species_lasso.plot, filename = 
#          "data/Machine_Learning_Analysis/Lasso_Species_Heatmap_all_prev_nzv_logZ.svg",
#        width = 4, height = 3)
# 
# s2s_species_rf <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "randomforest", 
#     data_type = "Species"
#   )
# s2s_species_rf.plot <- 
#   ml_heatmap_summary(s2s_species_rf)
# ggsave(s2s_species_rf.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_Species_Heatmap_all.svg",
#        width = 4, height = 3)
# 
# s2s_species_xgb <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "xgboost", 
#     data_type = "Species"
#   )
# s2s_species_xgb.plot <- 
#   ml_heatmap_summary(s2s_species_xgb)
# ggsave(s2s_species_xgb.plot, filename = 
#          "data/Machine_Learning_Analysis/xgb_Species_Heatmap_all.svg",
#        width = 4, height = 3)
# 
# write.csv(s2s_species_lasso,
#           file = 'data/Machine_Learning_Analysis/model_stats/Species/Study2Study_transfer_performance_lasso.csv')
# write.csv(s2s_species_rf,
#           file = 'data/Machine_Learning_Analysis/model_stats/Species/Study2Study_transfer_performance_rf.csv')
# write.csv(s2s_species_xgb,
#           file = 'data/Machine_Learning_Analysis/model_stats/Species/Study2Study_transfer_performance_xgb.csv')
# 
# # LOSO ------ 
# loso_species_lasso <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush, 
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso", 
#     data_type = "Species"
#     )
# loso_species_lasso.plot <- cohort_comparison_bars(loso_species_lasso)
# loso_species_lasso.plot
# 
# ggsave(loso_species_lasso.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_Species_LOSO_barplot_prev0.1_nzv_logZ.png",
#        width = 3.5, height = 4)
# 
# loso_species_rf <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest", 
#     data_type = "Species"
#   )
# loso_species_rf.plot <- cohort_comparison_bars(loso_species_rf)
# ggsave(loso_species_rf.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_Species_LOSO_barplot.png",
#        width = 3.5, height = 4)
# 
# loso_species_xgb <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush, 
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost", 
#     data_type = "Species"
#   )
# loso_species_xgb.plot <- cohort_comparison_bars(loso_species_xgb)
# loso_species_xgb.plot
# ggsave(loso_species_xgb.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_Species_LOSO_barplot_prevalence_thres_0.01.png",
#        width = 3.5, height = 4)
# 
# 
# write.csv(loso_species_lasso, file = 'data/Machine_Learning_Analysis/model_stats/Species/LOSO_performance_lasso.csv')
# write.csv(loso_species_rf, file = 'data/Machine_Learning_Analysis/model_stats/Species/LOSO_performance_rf.csv')
# write.csv(loso_species_xgb, file = 'data/Machine_Learning_Analysis/model_stats/Species/LOSO_performance_xgb.csv')
# write.csv(loso_species_lgbm, file = 'data/Machine_Learning_Analysis/model_stats/Species/LOSO_performance_xgb.csv')
# 
# 
# #_______________________________________________________________________________
# #                               KOs                                      -----
# #_______________________________________________________________________________
# 
# obj <- subset_samples(dat.KOs.slim, donor_id %ni% low_qc[[1]])
# 
# obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
# obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
# obj_noRush <- subset_samples(obj, cohort != "Rush")
# obj_Rush <- subset_samples(obj, cohort == "Rush")
# obj_noTBC <- subset_samples(obj, cohort != "TBC")
# obj_TBC <- subset_samples(obj, cohort == "TBC")
# obj_noBonn <- subset_samples(obj, cohort != "Bonn")
# obj_Bonn <- subset_samples(obj, cohort == "Bonn")
# 
# # S2S ------ 
# s2s_KOs_lasso_f25 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     featSelection = T, 
#     nfeats = 25, 
#     data_type = "KOs.slim",
#     model_type = "lasso")
# 
# s2s_KOs_lasso_f25.plot <- ml_heatmap_summary(s2s_KOs_lasso_f25)
# s2s_KOs_lasso_f25.plot
# ggsave(s2s_KOs_lasso_f25.plot, filename = 
#          "data/Machine_Learning_Analysis/Lasso_KOs_Heatmap_n25.png",
#        width = 4, height = 3)
# 
# s2s_KOs_lasso_f100 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     featSelection = T, 
#     nfeats = 100, 
#     data_type = "KOs.slim",
#     model_type = "lasso")
# 
# s2s_KOs_lasso_f100.plot <- ml_heatmap_summary(s2s_KOs_lasso_f100)
# s2s_KOs_lasso_f100.plot
# ggsave(s2s_KOs_lasso_f100.plot, filename = 
#          "data/Machine_Learning_Analysis/Lasso_KOs_Heatmap_n100.png",
#        width = 4, height = 3)
# 
# 
# s2s_KOs_rf25 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 25
#   )
# s2s_KOs_rf25.plot <- ml_heatmap_summary(s2s_KOs_rf25)
# ggsave(s2s_KOs_rf25.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_KOs_Heatmap_n25.png",
#        width = 4, height = 3)
# 
# s2s_KOs_rf100 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 100, 
#     data_type = "KOs.slim",
#     featSelection = T
#   )
# 
# s2s_KOs_rf100.plot <- ml_heatmap_summary(s2s_KOs_rf100)
# ggsave(s2s_KOs_rf100.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_KOs_Heatmap_n500.png",
#        width = 4, height = 3)
# 
# s2s_KOs_xgb_f25 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 25
#   )
# s2s_KOs_xgb_f25.plot <- ml_heatmap_summary(s2s_KOs_xgb_f25);s2s_KOs_xgb_f25.plot
# ggsave(s2s_KOs_xgb_f25.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_KOs_Heatmap_n25.png",
#        width = 4, height = 3)
# 
# s2s_KOs_xgb_f100 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 100, 
#     data_type = "KOs.slim",
#     featSelection = T
#   )
# 
# s2s_KOs_xgb_f100.plot <- ml_heatmap_summary(s2s_KOs_xgb_f100)
# ggsave(s2s_KOs_xgb_f100.plot, filename = 
#          "data/Machine_Learning_Analysis/XGB_KOs_Heatmap_n100.png",
#        width = 4, height = 3)
# 
# 
# write.csv(s2s_KOs_lasso, 
#           file = 'data/Machine_Learning_Analysis/model_stats/KOs/Study2Study_transfer_performance_lasso.csv')
# write.csv(s2s_KOs_rf, 
#           file = 'data/Machine_Learning_Analysis/model_stats/KOs/Study2Study_transfer_performance_rf.csv')
# 
# 
# # LOSO ------ 
# loso_KOs_lasso_f25 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso",
#     featSelection = TRUE, 
#     data_type = "KOs.slim",
#     nfeats = 25
#   )
# loso_KOs_lasso_f25.plot <- cohort_comparison_bars(loso_KOs_lasso_f25)
# loso_KOs_lasso_f25.plot
# ggsave(loso_KOs_lasso_f25.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_KOs_LOSO_barplot_n25.png",
#        width = 3.5, height = 4)
# 
# loso_KOs_lasso_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso",
#     featSelection = TRUE, 
#     nfeats = 100, 
#     data_type = "KOs.slim"
#   )
# loso_KOs_lasso_f100.plot <- cohort_comparison_bars(loso_KOs_lasso_f100)
# loso_KOs_lasso_f100.plot
# ggsave(loso_KOs_lasso_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_KOs_LOSO_barplot_f100.png",
#        width = 3.5, height = 4)
# 
# 
# # loso_KOs_lasso_all <-
# #   ml_loso(
# #     Shanghai = obj_Shanghai,
# #     noShanghai = obj_noShanghai,
# #     TBC = obj_TBC,
# #     noTBC = obj_noTBC,
# #     Rush = obj_Rush,
# #     noRush = obj_noRush,
# #     Bonn = obj_Bonn, 
# #     noBonn = obj_noBonn,
# #     model_type = "lasso"
# #   )
# # loso_KOs_lasso_all.plot <- cohort_comparison_bars(loso_KOs_lasso_all)
# # loso_KOs_lasso_all.plot
# # ggsave(loso_KOs_lasso_all.plot, 
# #        filename = "data/Machine_Learning_Analysis/Lasso_KOs_LOSO_barplot_all.png", device = "png",
# #        width = 3.5, height = 4)
# 
# 
# loso_KOs_rf_f10 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 10
#   )
# loso_KOs_rf_f10.plot <- cohort_comparison_bars(loso_KOs_rf_f10)
# ggsave(loso_KOs_rf_f10.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_KOs_LOSO_barplot_f10.png",
#        width = 3.5, height = 4)
# 
# 
# loso_KOs_rf_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 100
#   )
# loso_KOs_rf_f100.plot <- cohort_comparison_bars(loso_KOs_rf_f100)
# loso_KOs_rf_f100.plot
# ggsave(loso_KOs_rf_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_KOs_LOSO_barplot_f100_SNM.png",
#        width = 3.5, height = 4)
# 
# 
# loso_KOs_rf_f1000 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# loso_KOs_rf_1000.plot <- cohort_comparison_bars(loso_KOs_rf_f1000)
# ggsave(loso_KOs_rf_f1000.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_KOs_LOSO_barplot_f1000.png",
#        width = 3.5, height = 4)
# 
# 
# #______________________________
# #          XGBoost
# #______________________________
# loso_KOs_xgb_f10 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#     featSelection = TRUE, 
#     nfeats = 10
#   )
# loso_KOs_xgb_f10.plot <- cohort_comparison_bars(loso_KOs_xgb_f10)
# ggsave(loso_KOs_xgb_f10.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_KOs_LOSO_barplot_f10.png",
#        width = 3.5, height = 4)
# 
# loso_KOs_xgb_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#     featSelection = TRUE, 
#     nfeats = 100
#   )
# loso_KOs_xgb_f100.plot <- cohort_comparison_bars(loso_KOs_xgb_f100)
# loso_KOs_xgb_f100.plot
# ggsave(loso_KOs_xgb_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_KOs_LOSO_barplot_f100_SNM.png",
#        width = 3.5, height = 4)
# 
# loso_KOs_xgb_f1000 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# loso_KOs_xgb_f1000.plot <- cohort_comparison_bars(loso_KOs_xgb_f1000)
# ggsave(loso_KOs_xgb_f1000.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_KOs_LOSO_barplot_f1000_SNM.png",
#        width = 3.5, height = 4)
# 
# write.csv(loso_KOs_lasso_f10, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_lasso_f10.csv')
# write.csv(loso_KOs_lasso_f100, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_lasso_f100.csv')
# write.csv(loso_KOs_lasso_f1000, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_lasso_f1000.csv')
# write.csv(loso_KOs_lasso_all, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_lasso_all.csv')
# 
# write.csv(loso_KOs_rf_f10, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_RF_f10.csv')
# write.csv(loso_KOs_rf_f100, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_RF_f100.csv')
# write.csv(loso_KOs_rf_f1000, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_RF_f1000.csv')
# # write.csv(loso_KOs_rf_all, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_RF_all.csv')
# 
# write.csv(loso_KOs_xgb_f10, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_xgb_f10.csv')
# write.csv(loso_KOs_xgb_f100, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_xgb_f100.csv')
# write.csv(loso_KOs_xgb_f1000, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_xgb_f1000.csv')
# # write.csv(loso_KOs_xgb_all, file = 'data/Machine_Learning_Analysis/model_stats/KOs/LOSO_performance_xgb_all.csv')
# 
# 
# rm(obj)
































































# 
# 
# #_______________________________________________________________________________
# #                               eggNOGs                                   -----
# #_______________________________________________________________________________
# 
# 
# obj <- subset_samples(dat.EGGNOGs.slim, donor_id %ni% low_qc[[1]]) %>% 
#   prevalence_filter(threshold = 0.01)
# 
# obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
# obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
# obj_noRush <- subset_samples(obj, cohort != "Rush")
# obj_Rush <- subset_samples(obj, cohort == "Rush")
# obj_noTBC <- subset_samples(obj, cohort != "TBC")
# obj_TBC <- subset_samples(obj, cohort == "TBC")
# 
# # S2S ------ 
# s2s_eggNOGs_lasso100 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     model_type = "lasso",
#     featSelection = TRUE,
#     nfeats = 100)
# s2s_eggNOGs_lasso100.plot <- 
#   ml_heatmap_summary(s2s_eggNOGs_lasso100)
# 
# 
# s2s_eggNOGs_rf25 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 25)
# 
# s2s_eggNOGs_rf25.plot <- ml_heatmap_summary(s2s_eggNOGs_rf25)
# ggsave(s2s_eggNOGs_rf25.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_eggNOGs_Heatmap_n25_prevalencefilter.png",
#        width = 4, height = 3)
# 
# 
# s2s_eggNOGs_rf100 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 100)
# 
# s2s_eggNOGs_rf100.plot <- ml_heatmap_summary(s2s_eggNOGs_rf100)
# ggsave(s2s_eggNOGs_rf100.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_eggNOGs_Heatmap_n100_prevalencefilter.png",
#        width = 4, height = 3)
# 
# 
# s2s_eggNOGs_rf500 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 500
#   )
# 
# s2s_eggNOGs_rf500.plot <- ml_heatmap_summary(s2s_eggNOGs_rf500)
# ggsave(s2s_eggNOGs_rf500.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_eggNOGs_Heatmap_n500_prevalencefilter.png",
#        width = 4, height = 3)
# 
# s2s_eggNOGs_rf1000 <-
#   ml_s2s(
#     Shanghai = obj_Shanghai,
#     TBC = obj_TBC,
#     Rush = obj_Rush,
#     Bonn = obj_Bonn, 
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# 
# s2s_eggNOGs_rf1000.plot <- ml_heatmap_summary(s2s_eggNOGs_rf1000)
# ggsave(s2s_eggNOGs_rf1000.plot, filename = 
#          "data/Machine_Learning_Analysis/RF_eggNOGs_Heatmap_n1000_prevalencefilter.png",
#        width = 4, height = 3)
# 
# ## LOSO ------ 
# loso_eggNOGs_lasso_f50 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso",
#     featSelection = TRUE, 
#     nfeats = 50
#   )
# loso_eggNOGs_lasso_f50.plot <- cohort_comparison_bars(loso_eggNOGs_lasso_f50)
# ggsave(loso_eggNOGs_lasso_f50.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_eggNOGs_LOSO_barplot_f50_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# loso_eggNOGs_lasso_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso",
#     featSelection = TRUE, 
#     nfeats = 100
#   )
# loso_eggNOGs_lasso_f100.plot <- cohort_comparison_bars(loso_eggNOGs_lasso_f100)
# ggsave(loso_eggNOGs_lasso_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_eggNOGs_LOSO_barplot_f100_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# loso_eggNOGs_lasso_f1000 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "lasso",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# loso_eggNOGs_lasso_f1000.plot <- cohort_comparison_bars(loso_eggNOGs_lasso_f1000)
# ggsave(loso_eggNOGs_lasso_f1000.plot, 
#        filename = "data/Machine_Learning_Analysis/Lasso_eggNOGs_LOSO_barplot_f1000_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# 
# 
# 
# loso_eggNOGs_rf_f10 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 10
#   )
# loso_eggNOGs_rf_f10.plot <- cohort_comparison_bars(loso_eggNOGs_rf_f10)
# ggsave(loso_eggNOGs_rf_f10.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_eggNOGs_LOSO_barplot_f10_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# 
# loso_eggNOGs_rf_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 100
#   )
# loso_eggNOGs_rf_f100.plot <- cohort_comparison_bars(loso_eggNOGs_rf_f100)
# ggsave(loso_eggNOGs_rf_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_eggNOGs_LOSO_barplot_f100_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# 
# loso_eggNOGs_rf_f1000 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "randomforest",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# loso_eggNOGs_rf_f1000.plot <- cohort_comparison_bars(loso_eggNOGs_rf_f1000)
# ggsave(loso_eggNOGs_rf_f1000.plot, 
#        filename = "data/Machine_Learning_Analysis/RF_eggNOGs_LOSO_barplot_f1000_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# #______________________________
# #______________________________
# #          XGBoost
# #______________________________
# loso_eggNOGs_xgb_f100 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#     featSelection = TRUE, 
#     nfeats = 100
#   )
# loso_eggNOGs_xgb_f100.plot <-
#   cohort_comparison_bars(loso_eggNOGs_xgb_f100)
# loso_eggNOGs_xgb_f100.plot
# 
# ggsave(loso_eggNOGs_xgb_f100.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_eggNOGs_LOSO_barplot_f100_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# 
# 
# loso_eggNOGs_xgb_f1000 <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#     featSelection = TRUE, 
#     nfeats = 1000
#   )
# loso_eggNOGs_xgb_f1000.plot <- cohort_comparison_bars(loso_eggNOGs_xgb_f1000)
# ggsave(loso_eggNOGs_xgb_f1000.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_eggNOGs_LOSO_barplot_f1000_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# loso_eggNOGs_xgb_all <-
#   ml_loso(
#     Shanghai = obj_Shanghai,
#     noShanghai = obj_noShanghai,
#     TBC = obj_TBC,
#     noTBC = obj_noTBC,
#     Rush = obj_Rush,
#     noRush = obj_noRush,
#     Bonn = obj_Bonn, 
#     noBonn = obj_noBonn,
#     model_type = "xgboost",
#   )
# 
# loso_eggNOGs_xgb_all.plot <- cohort_comparison_bars(loso_eggNOGs_xgb_all)
# ggsave(loso_eggNOGs_xgb_all.plot, 
#        filename = "data/Machine_Learning_Analysis/XGBoost_eggNOGs_LOSO_barplot_all_prevalencefilter.png",
#        width = 3.5, height = 4)
# 
# 
# 
# 
# write.csv(s2s_eggNOGs_lasso, 
#           file = 'data/Machine_Learning_Analysis/model_stats/eggNOGs/Study2Study_transfer_performance_lasso.csv')
# write.csv(s2s_eggNOGs_rf, 
#           file = 'data/Machine_Learning_Analysis/model_stats/eggNOGs/Study2Study_transfer_performance_rf.csv')
# write.csv(loso_eggNOGs_lasso, file = 'data/Machine_Learning_Analysis/model_stats/eggNOGs/LOSO_performance_lasso.csv')
# write.csv(loso_eggNOGs_rf, file = 'data/Machine_Learning_Analysis/model_stats/eggNOGs/LOSO_performance_rf.csv')
# 
# rm(obj)

























































#______________________________
# 
# obj <- subset_samples(dat.species, donor_id %ni% low_qc[[1]]) %>% 
#   core(detection = 0, prevalence = 0.025) 
# 
# obj_tbc_all <- subset_samples(obj, cohort == "TBC")
# obj_rush_all <- subset_samples(obj, cohort == "Rush")
# 
# cohort_summary_species <-
#   ml_summary(obj_tbc_all = obj_tbc_all, 
#                        obj_rush_all = obj_rush_all,
#                        featSelection =  NULL)
# 
# heatmap <- 
#   cohort_summary_species %>% 
#   mutate(model_perf = coalesce(mean, .estimate)) %>% 
#   filter(.metric == "roc_auc") %>% 
#   pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
#   ggplot(aes(train, test, fill = roc_auc)) +
#   geom_tile() +
#   labs(x = "Training sets", y = "Testing sets", fill = "AUROC") +
#   geom_text(aes(label=round(roc_auc, digits = 2)), size=3, 
#             vjust = 0.77, color = "white") +
#   scale_fill_viridis_c(option = "magma", begin = .9, end = 0, 
#                        na.value = "transparent") +
#   guides(fill = guide_colourbar(barwidth = 1, barheight = 9)) +
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# heatmap
# ggsave(heatmap, filename = "data/Machine_Learning_Analysis/Lasso_Species_Heatmap_Cohort.png",
#        width = 4, height = 3)
# 
# # bars <- 
# #   cohort_summary %>% 
# #   mutate(model_perf = coalesce(mean, .estimate)) %>% 
# #   filter(.metric == "roc_auc") %>% 
# #   dplyr::mutate(model_type = if_else((train == test), "10-fold CV", "Prediction")) %>%
# #   pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
# #   ggplot(aes(train, roc_auc, color = model_type, fill = test)) +
# #   theme_classic() +
# #   labs(x = "Training set", y = "AUROC", fill = "Test") +
# #   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
# #   geom_col(position=position_dodge(), alpha = 0.7) +
# #   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
# #   scale_fill_tableau(palette = "Tableau 10") +
# #   scale_color_manual(values = c("10-fold CV" = "#000000", "Prediction" = "#b3b3b3")) +
# #   theme(panel.grid = element_blank())
# # bars
# 
# barplot_species <- cohort_comparison_bars(cohort_summary_species)
# ggsave(barplot_species, filename = "data/Machine_Learning_Analysis/Lasso_Species_Barplot_Cohort.png",
#        width = 3, height = 3)
# 
# 
# 
# cohort_group_summary_species <- lasso_cohort_x_group_summary(obj_tbc_all, obj_rush_all)
# 
# heatmap2 <- 
#   cohort_group_summary %>% 
#   mutate(model_perf = coalesce(mean, .estimate)) %>% 
#   filter(.metric == "roc_auc") %>% 
#   pivot_wider(names_from = .metric, values_from = "model_perf") %>% 
#   mutate(roc_auc = if_else((train_cohort == test_cohort) & 
#                              (train_group != test_group), 
#                            as.double(NA), roc_auc)) %>% 
#   dplyr::group_by(train) %>% 
#   ggplot(aes(train, test, fill = roc_auc)) +
#   theme_bw() +
#   geom_tile() +
#   labs(x = "Training sets", y = "Testing sets", fill = "AUROC") +
#   geom_text(aes(label=round(roc_auc, digits = 2)), size=3, 
#             vjust = 0.77, color = "white") +
#   scale_fill_viridis_c(option = "magma", begin = .9, end = 0, 
#                        na.value = "transparent") +
#   guides(fill = guide_colourbar(barwidth = 1, barheight = 15)) +
#   theme(panel.border = element_blank(),
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# heatmap2
# ggsave(heatmap2, filename = "data/Machine_Learning_Analysis/Lasso_Species_Heatmap_CohortGroup.png",
#        width = 5, height = 4)
# 
# bars2 <- 
#   cohort_group_summary %>% 
#   dplyr::mutate(model_perf = coalesce(mean, .estimate)) %>% 
#   dplyr::filter(.metric == "roc_auc") %>%
#   dplyr::mutate(traintest = paste(train, test, sep = "_")) %>% 
#   dplyr::mutate(model_type = if_else((train_cohort == test_cohort) & (train_group == test_group), "10-fold CV" , "Prediction")) %>%
#   pivot_wider(names_from = .metric, values_from = "model_perf") %>%
#   dplyr::filter(if_else((train_cohort == test_cohort) & (train_group != test_group), F , T)) %>% 
#   dplyr::group_by(train) %>%
#   dplyr::mutate(ord = row_number()) %>%
#   ggplot(aes(train, roc_auc, group = ord, color = model_type, fill = test)) +
#   theme_classic() +
#   labs(x = "Training set", y = "AUROC", color = "Model", fill = "Test") +
#   scale_y_continuous(expand = expansion(mult = c(0, .1))) +
#   geom_col(position=position_dodge(), alpha = 0.7) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_tableau(palette = "Tableau 20") +
#   scale_color_manual(values = c("10-fold CV" = "#000000", "Prediction" = "#b3b3b3")) +
#   theme(panel.grid = element_blank())
# bars2
# ggsave(bars2, filename = "data/Machine_Learning_Analysis/Lasso_Species_Barplot_CohortGroup.png",
#        width = 5, height = 3)
# 
# 
# 
# 
# 
# #______________________________-
# #                               Pathways
# #______________________________-
# 
# obj <- subset_samples(dat.path.slim, donor_id %ni% low_qc[[1]]) %>% 
#   core(detection = 0, prevalence = 0.025) 
# 
# obj_tbc_all <- subset_samples(obj, cohort == "TBC")
# obj_rush_all <- subset_samples(obj, cohort == "Rush")
# 
# # feature selection via mRMR
# mRMR_pathways <- prep_mRMR_input(obj)
# hits_pathways <- feature_selection(mRMR_pathways, n_features = 100, n_algorithm = 1)
# length(unique(hits_pathways))
# 
# ml_input_tbc <- prep_ml_input(obj_tbc_all) %>% 
#   select(donor_id, PD, cohort, paired, all_of(hits_pathways))
# ml_input_rush <- prep_ml_input(obj_rush_all) %>%
#   select(donor_id, PD, cohort, paired, matches(hits_pathways))
# 
# # Note: some features can't be selected after identification via mRMR
# cohort_summary_pathways <- ml_summary(obj_tbc_all = obj_tbc_all, 
#                                                 obj_rush_all = obj_rush_all,
#                                                 featSelection = NULL)
# 
# cohort_group_summary_pathways <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
#                                                      obj_rush_all = obj_rush_all,
#                                                      featSelection = hits_pathways)
# 
# 
# barplot_pathways <- cohort_comparison_bars(cohort_summary_pathways)
# ggsave(barplot_pathways, filename = "data/Machine_Learning_Analysis/Lasso_Pathways_Barplot_Cohort.png",
#        width = 3, height = 3)
# 
# heatmap_pathways <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_pathways)
# ggsave(heatmap_pathways, filename = "data/Machine_Learning_Analysis/Lasso_Pathways_Heatmap_CohortGroup.png",
#        width = 5, height = 4)
# 
# # # TROUBLE SHOOTING ERROR
# # ml_input_tbc <-
# #   prep_ml_input(obj = obj_tbc_all) %>%
# #   pivot_longer(!c(donor_id, PD, cohort, paired),
# #                names_to = "features", values_to = "nums") %>%
# #   dplyr::select(-c(donor_id, PD, cohort, paired, nums)) %>%
# #   dplyr::distinct() %>%
# #   decode_rfriendly_rows("features")
# # ml_input_rush <- prep_ml_input(obj_rush_all) %>%
# #   pivot_longer(!c(donor_id, PD, cohort, paired),
# #                names_to = "features", values_to = "nums")
# # tst <-
# #   hits_pathways %>%
# #   as.data.frame() %>%
# #   rename("features" = ".") %>%
# #   decode_rfriendly_rows("features") %>%
# #   anti_join(ml_input_tbc, by = "features")
# # 
# # ml_input_tbc <-
# #   prep_ml_input(obj_tbc_all)
# # ml_input_tbc2 <- ml_input_tbc %>% 
# #   dplyr::select(donor_id, PD, cohort, paired, matches(hits_pathways,ignore.case = T))
# # 
# # ml_input_tbc3 <-
# #   prep_ml_input(obj_tbc_all)
# 
# 
# #______________________________-
# #                                  KOs
# #______________________________-
# 
# obj <- subset_samples(dat.KOs.slim, donor_id %ni% low_qc[[1]]) %>% 
#   core(detection = 0, prevalence = 0.15) 
# 
# obj_tbc_all <- subset_samples(obj, cohort == "TBC")
# obj_rush_all <- subset_samples(obj, cohort == "Rush")
# 
# # feature selection via mRMR
# mRMR_KOs <- prep_mRMR_input(obj)
# hits_KOs <- feature_selection(df = mRMR_KOs, n_features = 300, n_algorithm = 1)
# length(unique(hits_KOs))
# 
# # Note: some features can't be selected after identification via mRMR
# cohort_summary_KOs <- ml_summary(obj_tbc_all = obj_tbc_all, 
#                                                 obj_rush_all = obj_rush_all,
#                                                 featSelection = hits_KOs)
# 
# cohort_group_summary_KOs <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
#                                                               obj_rush_all = obj_rush_all,
#                                                               featSelection = hits_KOs)
# 
# 
# barplot_KOs <- cohort_comparison_bars(cohort_summary_KOs)
# ggsave(barplot_KOs, filename = "data/Machine_Learning_Analysis/Lasso_KOs_Barplot_Cohort.png",
#        width = 3, height = 3)
# 
# heatmap_KOs <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_KOs)
# ggsave(heatmap_KOs, filename = "data/Machine_Learning_Analysis/Lasso_KOs_Heatmap_CohortGroup.png",
#        width = 5, height = 4)
# 
# 
# 
# #______________________________-
# #                                  eggNOGs
# #______________________________-
# 
# # add feature selection process
# obj <- subset_samples(dat.EGGNOGs.slim, donor_id %ni% low_qc[[1]]) %>% 
#   core(detection = 0, prevalence = 0.15) 
# 
# obj_tbc_all <- subset_samples(obj, cohort == "TBC")
# obj_rush_all <- subset_samples(obj, cohort == "Rush")
# 
# # feature selection via mRMR
# mRMR_eggNOGs <- prep_mRMR_input(obj)
# hits_eggNOGs <- feature_selection(df = mRMR_eggNOGs, n_features = 100, n_algorithm = 1)
# length(unique(hits_eggNOGs))
# 
# # Note: some features can't be selected after identification via mRMR
# cohort_summary_eggNOGs <- ml_summary(obj_tbc_all = obj_tbc_all, 
#                                            obj_rush_all = obj_rush_all,
#                                            featSelection = hits_eggNOGs)
# 
# cohort_group_summary_eggNOGs <- lasso_cohort_x_group_summary(obj_tbc_all = obj_tbc_all, 
#                                                          obj_rush_all = obj_rush_all,
#                                                          featSelection = hits_eggNOGs)
# 
# 
# barplot_eggNOGs <- cohort_comparison_bars(cohort_summary_eggNOGs)
# ggsave(barplot_eggNOGs, filename = "data/Machine_Learning_Analysis/Lasso_eggNOGs_Barplot_Cohort.png",
#        width = 3, height = 3)
# 
# heatmap_eggNOGs <- heatmap_cohort_groups(cohort_group_summary = cohort_group_summary_eggNOGs)
# ggsave(heatmap_KOs, filename = "data/Machine_Learning_Analysis/Lasso_eggNOGs_Heatmap_CohortGroup.png",
#        width = 5, height = 4)






















  