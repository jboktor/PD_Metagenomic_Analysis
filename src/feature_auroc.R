
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

phylo_objects <- readRDS(file = "files/Phyloseq_all-cohorts_clean_slim.rds")
cohort_list <- c("TBC", "Rush", "Shanghai", "Bonn")

# _______________________________________________________________________________
# Cohort specific ----

aucs <- tibble()
start_time <- Sys.time()
for (obj.name in names(phylo_objects)) {
  for (study in cohort_list) {
    dat <-
      phylo_objects[[obj.name]] %>%
      subset_samples(cohort == study) %>% 
      core(detection = 0, prevalence = 1e-10) %>% 
      microbiome::transform("compositional")

    dat.pd <- subset_samples(dat, donor_group == "PD") %>%
      decoded_abundance()
    dat.cont <- subset_samples(dat, donor_group != "PD") %>%
      decoded_abundance()

    require(foreach)
    require(doParallel)

    cores <- detectCores()
    cl <- makeCluster(cores[1] - 1)
    registerDoParallel(cl)

    roc2add <-
      foreach(
        feat = rownames(dat.pd), .combine = "rbind",
        .packages = c("magrittr", "pROC")
      ) %dopar% {
        x <- dat.pd[feat, ] %>% t()
        y <- dat.cont[feat, ] %>% t()

        # AUROC
        rocdata <- c(roc(controls = y, cases = x, direction = "<", ci = TRUE, auc = TRUE)$ci)
        data.frame(
          "feature" = feat,
          "data_level" = obj.name,
          "study" = study,
          "ci_lower" = rocdata[1],
          "auroc" = rocdata[2],
          "ci_upper" = rocdata[3]
        )
      }
    stopCluster(cl)
    aucs <- rbind(aucs, roc2add)
  }
}

end_time <- Sys.time()
cat(
  "AUROCs calculated in : ",
  end_time - start_time, attr(end_time - start_time, "units"), "\n"
)
saveRDS(aucs, 
        glue("files/Meta_analysis/feature-AUROCs_per-cohort_{Sys.Date()}.rds")
        )

# # _______________________________________________________________________________
# # LODO calculation ----
# 
# aucs <- tibble()
# start_time <- Sys.time()
# for (obj.name in names(datObjs)) {
#   for (study in cohort_list) {
#     dat <-
#       datObjs[[obj.name]] %>%
#       subset_samples(cohort != study) %>%
#       subset_samples(donor_id %ni% low_qc[[1]])
# 
#     dat.pd <- subset_samples(dat, donor_group == "PD") %>%
#       decoded_abundance()
#     dat.cont <- subset_samples(dat, donor_group != "PD") %>%
#       decoded_abundance()
# 
#     require(foreach)
#     require(doParallel)
# 
#     cores <- detectCores()
#     cl <- makeCluster(cores[1] - 1)
#     registerDoParallel(cl)
# 
#     roc2add <-
#       foreach(
#         feat = rownames(dat.pd), .combine = "rbind",
#         .packages = c("magrittr", "pROC")
#       ) %dopar% {
#         x <- dat.pd[feat, ] %>% t()
#         y <- dat.cont[feat, ] %>% t()
# 
#         # AUROC
#         rocdata <- c(roc(controls = y, cases = x, direction = "<", ci = TRUE, auc = TRUE)$ci)
#         data.frame(
#           "feature" = feat,
#           "data_level" = obj.name,
#           "study" = paste0("no_", study),
#           "ci_lower" = rocdata[1],
#           "auroc" = rocdata[2],
#           "ci_upper" = rocdata[3]
#         )
#       }
#     stopCluster(cl)
#     aucs <- rbind(aucs, roc2add)
#   }
# }
# 
# end_time <- Sys.time()
# cat(
#   "AUROCs calculated in : ",
#   end_time - start_time, attr(end_time - start_time, "units"), "\n"
# )
# save(aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")
# 
# 
# # # Save all AUROC data
# # load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")
# # aucs_lodo <- aucs
# # load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.RData")
# # all_aucs <- bind_rows(aucs, aucs_lodo)
# # save(all_aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_all.RData")

