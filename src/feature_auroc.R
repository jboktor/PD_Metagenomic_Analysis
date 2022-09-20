
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
load("files/low_quality_samples.RData")
# load_data("Merged_ML")

phylo_objects <- readRDS("files/Phyloseq_Merged_ML_clean.rds")

cohort_list <- c("TBC", "Rush", "Shanghai", "Bonn")
# datObjs <-
#   list(
#     dat.species,
#     dat.genus,
#     dat.phylum,
#     dat.path.slim,
#     dat.ec.slim,
#     dat.KOs.slim,
#     dat.PFAMs.slim,
#     dat.EGGNOGs.slim
#   )
# names(phylo_objects) <-
#   c(
#     "Species",
#     "Genus",
#     "Phylum",
#     "Pathways.slim",
#     "Enzymes.slim",
#     "KOs.slim",
#     "PFAMs.slim",
#     "eggNOGs.slim"
#   )

# _______________________________________________________________________________
# Cohort specific ----

aucs <- tibble()
start_time <- Sys.time()
for (obj.name in names(phylo_objects)) {
  for (study in cohort_list) {
    dat <-
      phylo_objects[[obj.name]] %>%
      subset_samples(cohort == study)

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
saveRDS(aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.rds")


# _______________________________________________________________________________
# LODO calculation ----

aucs <- tibble()
start_time <- Sys.time()
for (obj.name in names(datObjs)) {
  for (study in cohort_list) {
    dat <-
      datObjs[[obj.name]] %>%
      subset_samples(cohort != study) %>%
      subset_samples(donor_id %ni% low_qc[[1]])

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
          "study" = paste0("no_", study),
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
save(aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")


# # Save all AUROC data
# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO.RData")
# aucs_lodo <- aucs
# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.RData")
# all_aucs <- bind_rows(aucs, aucs_lodo)
# save(all_aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_all.RData")



# _______________________________________________________________________________
# on RAREFIED DATA  ----
# _______________________________________________________________________________

load("files/Phyloseq_Merged_ML_Rarefied.RData")
datObjs <- phyloseq_objs_rare
cohort_list <- c("TBC", "Rush", "Shanghai", "Bonn")

# _______________________________________________________________________________
# Cohort specific ----

aucs <- tibble()
start_time <- Sys.time()
for (obj.name in names(datObjs)) {
  for (study in cohort_list) {
    dat <-
      datObjs[[obj.name]] %>%
      subset_samples(cohort == study) %>%
      subset_samples(donor_id %ni% low_qc[[1]])

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
save(aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_Rarefied.RData")


# _______________________________________________________________________________
# LODO calculation ----

aucs <- tibble()
start_time <- Sys.time()
for (obj.name in names(datObjs)) {
  for (study in cohort_list) {
    dat <-
      datObjs[[obj.name]] %>%
      subset_samples(cohort != study) %>%
      subset_samples(donor_id %ni% low_qc[[1]])

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
          "study" = paste0("no_", study),
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
save(aucs, file = "data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs_LODO_Rarefied.RData")




# # Without parallelization
#
# for (obj.name in names(datObjs)){
#   for (study in cohort_list){
#
#     dat <-
#       datObjs[[obj.name]] %>%
#       subset_samples(cohort == study) %>%
#       subset_samples(donor_id %ni% low_qc[[1]])
#
#     dat.pd <- subset_samples(dat, donor_group == "PD") %>%
#       decoded_abundance()
#     dat.cont <- subset_samples(dat, donor_group != "PD") %>%
#       decoded_abundance()
#
#     for (feat in rownames(dat.pd)){
#       print(feat)
#       x <- dat.pd[feat,] %>% t()
#       y <- dat.cont[feat,] %>% t()
#
#       # AUROC
#       roc <- c(roc(controls=y, cases=x, direction='<', ci=TRUE, auc=TRUE)$ci)
#       roc2add <- cbind(feat, obj.name, study, roc[1], roc[2], roc[3])
#       aucs <- rbind(aucs, roc2add)
#     }
#   }
# }
#
# names(aucs) <- c("feature", "data_level", "study", "ci_lower", "auroc", "ci_upper")
# end_time <- Sys.time()
# cat("Correlations calculated in : ",
#     end_time - start_time, attr(end_time - start_time, "units"), "\n")
