# Differentially Abundant Features (DAFs)

#            Load Data & functions
rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
load("files/low_quality_samples.RData")
wkd <- getwd()

# _______________________________________________________________________________
#                                Cohorts   - TBC & RUSH                      ----
# _______________________________________________________________________________
remove_dats()
# load_data("Merged")
phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")
asdf <- readRDS("files/Phyloseq_Merged_ML_clean.rds")

# ________________________
#          Species ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Species"]])
# PD v PC abundance data
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
# Plot Variance Estimate
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.1)
# PD v HC PAIRED abundance data
dat_pdhc <- subset_samples(dat.object, paired != "No")
# Plot Variance Estimate
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.1)
fit_models(
  dat = dat.object,
  obj.name = "Species",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 1,
  plot_scatter = T
)

# ________________________
#          Genera ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["Genus"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0)
fit_models(
  dat = dat.object,
  obj.name = "Genera",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 1,
  plot_scatter = T
)


# ________________________
#          Phylum  ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["Phylum"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0)
fit_models(
  dat = dat.object,
  obj.name = "Phylum",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 1,
  plot_scatter = T
)

# ________________________
#         Pathways slim ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["Pathways.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)
fit_models(
  dat = dat.object,
  obj.name = "Pathways.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 1,
  plot_scatter = F
)

# ________________________
#         Enzymes slim ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["Enzymes.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)
fit_models(
  dat = dat.object,
  obj.name = "Enzymes.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = F
)

# ________________________
#    Kegg Orthology slim ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["KOs.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(0.2)
fit_models(
  dat = dat.object,
  obj.name = "KOs.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = F
)

# ________________________
#    Gene Ontology ####
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["GOs.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)
fit_models(
  dat = dat.object,
  obj.name = "GOs.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = F
)


# GO Molecular Function
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".MF.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "Merged", tag = "_MF")

# GO Biological Processes
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".BP.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "Merged", tag = "_BP")

# GO Cellular Compartments
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".CC.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "Merged", tag = "_CC")


# ________________________
#          PFAMS slim ----
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["Pfams.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)
fit_models(
  dat = dat.object,
  obj.name = "Pfams.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = F
)

# ________________________
#   EGGNOGS slim
# ________________________

dat.object <- maaslin_prep(phyloseq_objs[["eggNOGs.slim"]])
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.7)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.7)
fit_models(
  dat = dat.object,
  obj.name = "eggNOGs.slim",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = F
)


# _____________________________________________
#              Plot Summaries MERGED    ----
# _____________________________________________

phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")

plot_dafs(obj.name = "Phylum", obj = phyloseq_objs[["Phylum"]], cohort = "Merged")
plot_dafs(obj.name = "Genera", obj = phyloseq_objs[["Genus"]], cohort = "Merged")
plot_dafs(obj.name = "Species", obj = phyloseq_objs[["Species"]], cohort = "Merged")
plot_daf_summary(obj.name = "Species", obj = phyloseq_objs[["Species"]], cohort = "Merged")

plot_dafs(obj.name = "Pathways.slim", obj = phyloseq_objs[["Pathways.slim"]], cohort = "Merged")
plot_daf_summary(obj.name = "Pathways.slim", obj = phyloseq_objs[["Pathways.slim"]], cohort = "Merged", top = 15)

plot_dafs(obj.name = "Enzymes.slim", obj = phyloseq_objs[["Enzymes.slim"]], cohort = "Merged")
plot_daf_summary(obj.name = "Enzymes.slim", obj = phyloseq_objs[["Enzymes.slim"]], cohort = "Merged", top = 20)

plot_dafs(obj.name = "KOs.slim", obj = phyloseq_objs[["KOs.slim"]], cohort = "Merged")
plot_daf_summary(obj.name = "KOs.slim", obj = phyloseq_objs[["KOs.slim"]], cohort = "Merged", top = 15)

plot_dafs(obj.name = "GOs.slim", obj = phyloseq_objs[["GOs.slim"]], cohort = "Merged", tag = "")
plot_daf_summary(obj.name = "GOs.slim", obj = phyloseq_objs[["GOs.slim"]], cohort = "Merged", top = 15)

plot_dafs(obj.name = "eggNOGs.slim", obj = phyloseq_objs[["eggNOGs.slim"]], cohort = "Merged")
plot_daf_summary(obj.name = "eggNOGs.slim", obj = phyloseq_objs[["eggNOGs.slim"]], cohort = "Merged", top = 20)

plot_dafs(obj.name = "Pfams.slim", obj = phyloseq_objs[["Pfams.slim"]], cohort = "Merged")
plot_daf_summary(obj.name = "Pfams.slim", obj = phyloseq_objs[["Pfams.slim"]], cohort = "Merged", top = 15)






# _______________________________________________________________________________
# _______________________________________________________________________________
# _______________________________________________________________________________

# Summary Table for Associations


maaslin2_stat_summary <- tibble()
dat_list <- c(
  "Species",
  "Genera",
  "Phylum",
  "Pathways.slim",
  "Enzymes.slim",
  "KOs.slim",
  "GOs.slim",
  "PFAMs.slim",
  "EGGNOGs.slim"
)

shared_features <- vector(mode = "list")


for (obj.name in dat_list) {
  cohort <- "Merged"

  ### Read-in MaAsLin2 output
  Maas.pd.pc <-
    read_tsv(
      paste0(
        "data/MaAsLin2_Analysis/",
        cohort,
        "/",
        obj.name,
        "_PDvPC_maaslin2_output/all_results.tsv"
      ),
      col_names = T
    ) %>%
    filter(value == "Population Control")
  Maas.pd.pc_assoc <- Maas.pd.pc %>%
    filter(qval < 0.25) %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    dplyr::rename("feature" = "fullnames")

  Maas.pd.hc <-
    read_tsv(
      paste0(
        "data/MaAsLin2_Analysis/",
        cohort, "/", obj.name,
        "_PDvHC_maaslin2_output/all_results.tsv"
      ),
      col_names = T
    ) %>%
    filter(value == "Household Control")
  Maas.pd.hc_assoc <- Maas.pd.hc %>%
    filter(qval < 0.25) %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    dplyr::rename("feature" = "fullnames")

  Maas.pd.pc_sig <- Maas.pd.pc_assoc %>% filter(qval <= 0.1)
  Maas.pd.hc_sig <- Maas.pd.hc_assoc %>% filter(qval <= 0.1)
  ## Split both inputs into enriched and depleted features
  Maas.pd.pc_sig_PD.UP <- Maas.pd.pc_sig %>% filter(coef < 0)
  Maas.pd.pc_sig_PD.DOWN <- Maas.pd.pc_sig %>% filter(coef > 0)
  Maas.pd.hc_sig_PD.UP <- Maas.pd.hc_sig %>% filter(coef < 0)
  Maas.pd.hc_sig_PD.DOWN <- Maas.pd.hc_sig %>% filter(coef > 0)
  Maas.sig.PD.UP.shared <-
    Maas.pd.pc_sig_PD.UP %>% filter(feature %in% Maas.pd.hc_sig_PD.UP$feature)
  Maas.sig.PD.DOWN.shared <-
    Maas.pd.pc_sig_PD.DOWN %>% filter(feature %in% Maas.pd.hc_sig_PD.DOWN$feature)

  Maas.pd.pc_assoc_PD.UP <- Maas.pd.pc_assoc %>% filter(coef < 0)
  Maas.pd.pc_assoc_PD.DOWN <- Maas.pd.pc_assoc %>% filter(coef > 0)
  Maas.pd.hc_assoc_PD.UP <- Maas.pd.hc_assoc %>% filter(coef < 0)
  Maas.pd.hc_assoc_PD.DOWN <- Maas.pd.hc_assoc %>% filter(coef > 0)
  Maas.assoc.PD.UP.shared <-
    Maas.pd.pc_assoc_PD.UP %>% filter(feature %in% Maas.pd.hc_assoc_PD.UP$feature)
  Maas.assoc.PD.DOWN.shared <-
    Maas.pd.pc_assoc_PD.DOWN %>% filter(feature %in% Maas.pd.hc_assoc_PD.DOWN$feature)


  stat_summary <-
    cbind(
      "Feature-level" = obj.name,
      # Significant features
      "PD Enriched vs PC (q<=0.1)" = Maas.pd.pc_sig_PD.UP %>% nrow(),
      "PD Depleted vs PC (q<=0.1)" = Maas.pd.pc_sig_PD.DOWN %>% nrow(),
      "PD Enriched vs HC (q<=0.1)" = Maas.pd.hc_sig_PD.UP %>% nrow(),
      "PD Depleted vs HC (q<=0.1)" = Maas.pd.hc_sig_PD.DOWN %>% nrow(),
      "PD Enriched shared (q<=0.1)" = Maas.sig.PD.UP.shared %>% nrow(),
      "PD Depleted shared (q<=0.1)" = Maas.sig.PD.DOWN.shared %>% nrow(),
      # Percentage of total tested feats
      "PD Enriched vs PC (q<=0.1) % of tested" = (Maas.pd.pc_sig_PD.UP %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Depleted vs PC (q<=0.1) % of tested" = (Maas.pd.pc_sig_PD.DOWN %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Enriched vs HC (q<=0.1) % of tested" = (Maas.pd.hc_sig_PD.UP %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Depleted vs HC (q<=0.1) % of tested" = (Maas.pd.hc_sig_PD.DOWN %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      # Associated features
      "PD Enriched vs PC (q<=0.25)" = Maas.pd.pc_assoc_PD.UP %>% nrow(),
      "PD Depleted vs PC (q<=0.25)" = Maas.pd.pc_assoc_PD.DOWN %>% nrow(),
      "PD Enriched vs HC (q<=0.25)" = Maas.pd.hc_assoc_PD.UP %>% nrow(),
      "PD Depleted vs HC (q<=0.25)" = Maas.pd.hc_assoc_PD.DOWN %>% nrow(),
      "PD Enriched shared (q<=0.25)" = Maas.assoc.PD.UP.shared %>% nrow(),
      "PD Depleted shared (q<=0.25)" = Maas.assoc.PD.DOWN.shared %>% nrow(),
      # Percentage of total tested feats
      "PD Enriched vs PC (q<=0.25) % of tested" = (Maas.pd.pc_assoc_PD.UP %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Depleted vs PC (q<=0.25) % of tested" = (Maas.pd.pc_assoc_PD.DOWN %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Enriched vs HC (q<=0.25) % of tested" = (Maas.pd.hc_assoc_PD.UP %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100,
      "PD Depleted vs HC (q<=0.25) % of tested" = (Maas.pd.hc_assoc_PD.DOWN %>% nrow()) / (Maas.pd.pc %>% nrow()) * 100
    )

  shared_features[[obj.name]][["PD_Enriched_significant"]] <-
    Maas.sig.PD.UP.shared %>% pull(feature)
  shared_features[[obj.name]][["PD_Depleted_significant"]] <-
    Maas.sig.PD.DOWN.shared %>% pull(feature)
  shared_features[[obj.name]][["PD_Enriched_associated"]] <-
    Maas.assoc.PD.UP.shared %>% pull(feature)
  shared_features[[obj.name]][["PD_Depleted_associated"]] <-
    Maas.assoc.PD.DOWN.shared %>% pull(feature)
  maaslin2_stat_summary <- rbind(maaslin2_stat_summary, stat_summary)
}

dt_list <- purrr::map(shared_features, as.data.table)
dt_list$Stats_Summary <- maaslin2_stat_summary
openxlsx::write.xlsx(dt_list, file = "files/Supplementary Tables/Table_S3 MaAsLin2 feature_associations-NEW.xlsx")


#
# #_______________________________________________________________________________
# #                                      TBC ONLY                                ----
# #_______________________________________________________________________________
#
# remove_dats()
# load_tbc()
#
# #________________________
# #          Species ----
# #________________________
#
# ## Filter low QC samples and trim low prevalence features
# dat.object <- maaslin_prep(dat.species)
# # PD v PC abundance data
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# # Plot Variance Estimate
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.1)
# # PD v HC PAIRED abundance data
# dat_pdhc = subset_samples(dat.object, paired !="No")
# # Plot Variance Estimate
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.1)
# fit_models(dat = dat.object,
#            obj.name = "Species",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "TBC")
#
# #________________________
# #          Genera ----
# #________________________
#
# dat.object <- maaslin_prep(dat.genus)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# fit_models(dat = dat.object,
#            obj.name = "Genera",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "TBC")
#
# #________________________
# #          Phylum  ----
# #________________________
#
# dat.object <- maaslin_prep(dat.phylum)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# fit_models(dat = dat.object,
#            obj.name = "Phylum",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "TBC")
#
# #________________________
# #         Pathways slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.path.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# fit_models(dat = dat.object,
#            obj.name = "Pathways.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
# #________________________
# #         Enzymes slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.ec.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# fit_models(dat = dat.object,
#            obj.name = "Enzymes.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
# #________________________
# #    Kegg Orthology slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.KOs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(0.2)
# fit_models(dat = dat.object,
#            obj.name = "KOs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
# #________________________
# #    Gene Ontology ####
# #________________________
#
# dat.object <- maaslin_prep(dat.GOs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# fit_models(dat = dat.object,
#            obj.name = "GOs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
#
# # GO Molecular Function
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".MF.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "TBC", tag = "_MF")
#
# # GO Biological Processes
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".BP.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "TBC", tag = "_BP")
#
# # GO Cellular Compartments
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".CC.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "TBC", tag = "_CC")
#
#
# #________________________
# #   EGGNOGS slim
# #________________________
#
# dat.object <- maaslin_prep(dat.EGGNOGs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.7)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.7)
# fit_models(dat = dat.object,
#            obj.name = "EGGNOGs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
# #________________________
# #          PFAMS slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.PFAMs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# fit_models(dat = dat.object,
#            obj.name = "PFAMs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "TBC")
#
#
# #_______________________________________________________________________________
# #                                      RUSH ONLY                                ----
# #_______________________________________________________________________________
#
# remove_dats()
# load_rush()
#
# ## Filter low QC samples and trim low prevalence features
# dat.object <- maaslin_prep(dat.species)
# # PD v PC abundance data
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# # Plot Variance Estimate
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.1)
# # PD v HC PAIRED abundance data
# dat_pdhc = subset_samples(dat.object, paired !="No")
# # Plot Variance Estimate
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.1)
# fit_models(dat = dat.object,
#            obj.name = "Species",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "RUSH")
#
# #________________________
# #          Genera ----
# #________________________
#
# dat.object <- maaslin_prep(dat.genus)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# fit_models(dat = dat.object,
#            obj.name = "Genera",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "RUSH")
#
# #________________________
# #          Phylum  ----
# #________________________
#
# dat.object <- maaslin_prep(dat.phylum)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0)
# fit_models(dat = dat.object,
#            obj.name = "Phylum",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 1,
#            plot_scatter = T,
#            cohort = "RUSH")
#
# #________________________
# #         Pathways slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.path.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# fit_models(dat = dat.object,
#            obj.name = "Pathways.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")
#
# #________________________
# #         Enzymes slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.ec.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# fit_models(dat = dat.object,
#            obj.name = "Enzymes.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")
#
# #________________________
# #    Kegg Orthology slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.KOs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>% # TROUBLE!
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(0.2)
# fit_models(dat = dat.object,
#            obj.name = "KOs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")
#
# #________________________
# #    Gene Ontology ####
# #________________________
#
# dat.object <- maaslin_prep(dat.GOs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.3)
# fit_models(dat = dat.object,
#            obj.name = "GOs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")
#
#
# # GO Molecular Function
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".MF.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "RUSH", tag = "_MF")
#
# # GO Biological Processes
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".BP.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "RUSH", tag = "_BP")
#
# # GO Cellular Compartments
# abund_rename <-
#   dat.object %>%
#   abundances() %>%
#   as.data.frame() %>%
#   rownames_to_column() %>%
#   filter(grepl(".CC.", rowname) ) %>%
#   column_to_rownames(var = "rowname") %>%
#   otu_table(taxa_are_rows=T)
# my_sample_data <- meta(dat.object) %>% sample_data()
# dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "RUSH", tag = "_CC")
#
#
# #________________________
# #   EGGNOGS slim
# #________________________
#
# dat.object <- maaslin_prep(dat.EGGNOGs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.7)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.7)
# fit_models(dat = dat.object,
#            obj.name = "EGGNOGs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")
#
# #________________________
# #          PFAMS slim ----
# #________________________
#
# dat.object <- maaslin_prep(dat.PFAMs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# fit_models(dat = dat.object,
#            obj.name = "PFAMs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F,
#            cohort = "RUSH")






# _____________________________________________

#
#
#
#
#
#
# #_____________________________________________
# ###              Plot Summaries TBC  ###
# #_____________________________________________
#
# remove_dats()
# load_data("TBC")
#
#
# plot_dafs(obj.name = "Phylum", obj = dat.phylum, cohort = "TBC")
# plot_dafs(obj.name = "Genera", obj = dat.genus, cohort = "TBC")
# plot_dafs(obj.name = "Species", obj = dat.species, cohort = "TBC")
# plot_daf_summary(obj.name = "Species", obj = dat.species, cohort = "TBC")
#
# plot_dafs(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "TBC")
# plot_daf_summary(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "TBC",
#                  repelyup = 0, repelydn = 0)
#
# plot_dafs(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "TBC")
# plot_daf_summary(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "TBC",
#                  repelyup = 0, repelydn = 0)
#
# plot_dafs(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "TBC")
# plot_daf_summary(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "TBC",
#                  repelyup = 0.0025, repelydn = 0.005)
#
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "TBC", tag = "")
# plot_daf_summary(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "TBC",
#                  repelyup = 0.001, repelydn = 0.001)
#
# plot_dafs(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "TBC")
# plot_daf_summary(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "TBC",
#                  repelyup = 0.002, repelydn = 0.002)
#
# plot_dafs(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "TBC")
# plot_daf_summary(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "TBC",
#                  repelyup = 0.001, repelydn = 0.001)
#
#
# #_____________________________________________
# ###                Plot Summaries RUSH    ###
# #_____________________________________________
#
# remove_dats()
# load_data("RUSH")
#
#
# plot_dafs(obj.name = "Phylum", obj = dat.phylum, cohort = "RUSH")
# plot_dafs(obj.name = "Genera", obj = dat.genus, cohort = "RUSH")
# plot_dafs(obj.name = "Species", obj = dat.species, cohort = "RUSH")
# plot_daf_summary(obj.name = "Species", obj = dat.species, cohort = "RUSH")
#
# plot_dafs(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "RUSH")
# plot_daf_summary(obj.name = "Pathways.slim", obj = dat.path.slim, cohort = "RUSH",
#                  repelyup = 0, repelydn = 0)
#
# plot_dafs(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "RUSH")
# plot_daf_summary(obj.name = "Enzymes.slim", obj = dat.ec.slim, cohort = "RUSH",
#                  repelyup = 0, repelydn = 0)
#
# plot_dafs(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "RUSH")
# plot_daf_summary(obj.name = "KOs.slim", obj = dat.KOs.slim, cohort = "RUSH",
#                  repelyup = 0.0025, repelydn = 0.005)
#
# plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "RUSH", tag = "")
# plot_daf_summary(obj.name = "GOs.slim", obj = dat.GOs.slim, cohort = "RUSH",
#                  repelyup = 0.001, repelydn = 0.001)
#
# plot_dafs(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "RUSH")
# plot_daf_summary(obj.name = "EGGNOGs.slim", obj = dat.EGGNOGs.slim, cohort = "RUSH",
#                  repelyup = 0.002, repelydn = 0.002)
#
# plot_dafs(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "RUSH")
# plot_daf_summary(obj.name = "PFAMs.slim", obj = dat.PFAMs.slim, cohort = "RUSH",
#                  repelyup = 0.001, repelydn = 0.001)


# _______________________________________________________________________________
#                META ANALYSIS - TBC/RUSH/Bonn/Shanghai                  ----
# _______________________________________________________________________________
remove_dats()
phyloseq_objs <- readRDS("files/Phyloseq_Merged_ML_clean.rds")

# ________________________
#      Species ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Species"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.1)

fit_models_meta(
  dat = dat.object,
  obj.name = "Species",
  df_input_data = df_input_data_all,
  cores = 1,
  plot_scatter = T
)

# ________________________
#       Genera ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Genus"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0)

fit_models_meta(
  dat = dat.object,
  obj.name = "Genus",
  df_input_data = df_input_data_all,
  cores = 1,
  plot_scatter = T
)

# ________________________
#   Pathways slim ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Pathways.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.1)

fit_models_meta(
  dat = dat.object,
  obj.name = "Pathways.slim",
  df_input_data = df_input_data_all,
  cores = 8,
  plot_scatter = T
)

# ________________________
#     Enzymes slim ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Enzymes.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)

fit_models_meta(
  dat = dat.object,
  obj.name = "Enzymes.slim",
  df_input_data = df_input_data_all,
  cores = 8,
  plot_scatter = T
)

# ________________________
#    Kegg Orthology slim ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["KOs.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)

fit_models_meta(
  dat = dat.object,
  obj.name = "KOs.slim",
  df_input_data = df_input_data_all,
  cores = 8,
  plot_scatter = T
)


# ________________________
#    Gene Ontology ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["GOs.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.3)

fit_models_meta(
  dat = dat.object,
  obj.name = "GOs.slim",
  df_input_data = df_input_data_all,
  cores = 6,
  plot_scatter = T
)

# GO Molecular Function
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".MF.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.MF <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.MF, cohort = "Merged_ML", tag = "_MF")

# GO Biological Processes
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".BP.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.BP <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.BP, cohort = "Merged_ML", tag = "_BP")

# GO Cellular Compartments
abund_rename <-
  dat.object %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(grepl(".CC.", rowname)) %>%
  column_to_rownames(var = "rowname") %>%
  otu_table(taxa_are_rows = T)
my_sample_data <- meta(dat.object) %>% sample_data()
dat.GOs.slim.CC <- phyloseq(abund_rename, my_sample_data)
plot_dafs(obj.name = "GOs.slim", obj = dat.GOs.slim.CC, cohort = "Merged_ML", tag = "_CC")


# ________________________
#   EGGNOGS slim ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["eggNOGs.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.5)

fit_models_meta(
  dat = dat.object,
  obj.name = "eggNOGs.slim",
  df_input_data = df_input_data_all,
  cores = 6,
  plot_scatter = T
)

# ________________________
#     PFAMS slim ----
# ________________________

## Filter low QC samples and trim low prevalence features
dat.object <- maaslin_prep(phyloseq_objs[["Pfams.slim"]])
# Plot Variance Estimate
variance_plot(dat.object)
df_input_data_all <- dat.object %>%
  microbiome::transform("compositional") %>%
  variance_filter(filter.percent = 0.2)

fit_models_meta(
  dat = dat.object,
  obj.name = "PFAMs.slim",
  df_input_data = df_input_data_all,
  cores = 8,
  plot_scatter = T
)

# dat.object <- maaslin_prep(dat.PFAMs.slim)
# dat_pdpc = subset_samples(dat.object, donor_group !="HC")
# variance_plot(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# dat_pdhc = subset_samples(dat.object, paired !="No")
# variance_plot(dat_pdhc)
# df_input_data_pdhc <- dat_pdhc %>%
#   microbiome::transform("compositional") %>%
#   variance_filter(filter.percent = 0.2)
# fit_models(dat = dat.object,
#            obj.name = "PFAMs.slim",
#            dat_pdpc = dat_pdpc,
#            dat_pdhc = dat_pdhc,
#            df_input_data_pdhc = df_input_data_pdhc,
#            df_input_data_pdpc = df_input_data_pdpc,
#            cores = 6,
#            plot_scatter = F)
