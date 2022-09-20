# SGV analysis

######## Load Data & functions
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/daf_functions.R")
source("src/metadata_prep_funcs.R")
load("files/low_quality_samples.RData")

# Load Sample Keys and deletion/variable SGVs

# m <- read.csv(file = "files/metadata_keys.csv", header= TRUE) %>%
#   dplyr::select(c(MBI_Sample_ID, id)) %>%
#   mutate(id = gsub("_", ".", id))

neg <- c(
  "S00A4-ATCC_MSA_1003_S96",
  "S00A4-neg2_S119",
  "S00A4-neg3_S125",
  "S00A4-neg_S118",
  "S00A4NegExt_P00A4_S94",
  "S00A4NegH2O_P00A4_S95",
  "S00A4_stagPos_S117",
  "BLANK"
)

dsgv <-
  read.csv(file = "files/dsgv.csv", row.names = 1, header = TRUE) %>%
  rownames_to_column(var = "MBI_Sample_ID") %>%
  dplyr::filter(MBI_Sample_ID %ni% neg) %>%
  mutate(MBI_Sample_ID = substr(MBI_Sample_ID, 1, 10))
vsgv <-
  read.csv(file = "files/vsgv.csv", row.names = 1, header = TRUE) %>%
  rownames_to_column(var = "MBI_Sample_ID") %>%
  dplyr::filter(MBI_Sample_ID %ni% neg) %>%
  mutate(MBI_Sample_ID = substr(MBI_Sample_ID, 1, 10))

m <- dat.species %>%
  subset_samples(donor_id %ni% low_qc[[1]]) %>%
  meta() %>%
  dplyr::select(donor_id, tube_id) %>%
  dplyr::rename(MBI_Sample_ID = tube_id)

## Map read IDs with group labels
df.dsgv <- dsgv %>%
  dplyr::mutate(MBI_Sample_ID = if_else(
    startsWith(MBI_Sample_ID, "S"), MBI_Sample_ID,
    substr(MBI_Sample_ID, 1, nchar(MBI_Sample_ID) - 4)
  )) %>%
  right_join(m, by = "MBI_Sample_ID") %>%
  column_to_rownames(var = "donor_id") %>%
  dplyr::select(-MBI_Sample_ID)

df.vsgv <- vsgv %>%
  dplyr::mutate(MBI_Sample_ID = if_else(
    startsWith(MBI_Sample_ID, "S"), MBI_Sample_ID,
    substr(MBI_Sample_ID, 1, nchar(MBI_Sample_ID) - 4)
  )) %>%
  right_join(m, by = "MBI_Sample_ID") %>%
  column_to_rownames(var = "donor_id") %>%
  dplyr::select(-MBI_Sample_ID)

# Data clean-up: Select columns with detection in at least 10 individuals
df.dsgv <- df.dsgv[, colSums(is.na(df.dsgv)) < (234 - 10)]
df.vsgv <- df.vsgv[, colSums(is.na(df.vsgv)) < (234 - 10)]

### Metadata
my_sample_data <- meta(dat.species) %>% sample_data()
### Create dsgv Phyloseq Obj
my_dsgv.ab_table <- otu_table(df.dsgv, taxa_are_rows = F)
dat.dsgv <- phyloseq(my_dsgv.ab_table, my_sample_data)
dat.dsgv
### Create dsgv Phyloseq Obj
my_vsgv.ab_table <- otu_table(df.vsgv, taxa_are_rows = F)
dat.vsgv <- phyloseq(my_vsgv.ab_table, my_sample_data)
dat.vsgv
save(dat.dsgv, file = "files/Phyloseq_Merged/Dsgvs_PhyloseqObj.RData")
save(dat.vsgv, file = "files/Phyloseq_Merged/Vsgvs_PhyloseqObj.RData")


#---------------------------------------------
#              DSVG Analysis
#---------------------------------------------

env.merge <- process_meta(dat.species, cohort = "Merged") %>%
  select(donor_id, donor_group, paired)

df.dsgv.stat <- df.dsgv %>%
  rownames_to_column(var = "donor_id") %>%
  left_join(env.merge, by = "donor_id") %>%
  column_to_rownames(var = "donor_id")



dat_pdhc <- subset_samples(dat.dsgv, paired != "No")
dat_pdpc <- subset_samples(dat.dsgv, donor_group != "HC")

pdhc.vals <- dat_pdhc %>%
  abundances() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id")
pdhc.vals <- pdhc.vals %>%
  group_col_from_ids(ids = pdhc.vals$id)
rownames(pdhc.vals) <- NULL
pdhc.vals <- pdhc.vals %>%
  column_to_rownames(var = "id")

pdpc.vals <- dat_pdpc %>%
  abundances() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id")
pdpc.vals <- pdpc.vals %>%
  group_col_from_ids(ids = pdpc.vals$id)
rownames(pdpc.vals) <- NULL
pdpc.vals <- pdpc.vals %>%
  column_to_rownames(var = "id")


pdhc_summary <- nmean_summary(df = pdhc.vals)
pdpc_summary <- nmean_summary(df = pdpc.vals)


#------------------------------------------------
# Fisher's Exact Test - Deletion SVs
#------------------------------------------------

#---------------
# PD v HC
#---------------
pdhc.dsgv.stats <- tibble()
n <- 0
for (i in 1:ncol(pdhc.vals)) {
  var <- colnames(pdhc.vals[i])
  print(var)
  print(n)
  n <- n + 1
  HC.prop <- filter(pdhc_summary, sv == var & group == "HC")$ratio_detected[1]
  HC.n <- filter(pdhc_summary, sv == var & group == "HC")$n[1]
  PD.prop <- filter(pdhc_summary, sv == var & group == "PD")$ratio_detected[1]
  PD.n <- filter(pdhc_summary, sv == var & group == "PD")$n[1]

  if (var != "group") {
    if (!is.na(PD.n) & !is.na(HC.n) > 0 &
      (HC.prop > 0 | PD.prop > 0) &
      (HC.prop < 1 | PD.prop < 1)) {
      tst <- pdhc.vals %>%
        dplyr::count(group, pdhc.vals[[i]]) %>%
        na.omit() %>%
        pivot_wider(
          values_from = n,
          names_from = "pdhc.vals[[i]]"
        ) %>%
        replace(is.na(.), 0) %>%
        dplyr::rename(hit = `1`, nohit = `0`) %>%
        column_to_rownames(var = "group")

      fish <- fisher.test(tst)
      cat("p-value: ", fish$p.value, "\n")

      row2add <-
        cbind(
          var,
          fish$estimate,
          fish$p.value,
          fish$conf.int[[1]],
          fish$conf.int[[2]],
          HC.prop,
          HC.n,
          PD.prop,
          PD.n
        )
      colnames(row2add) <-
        c(
          "sv",
          "estimate",
          "p.value",
          "CI.lower",
          "CI.upper",
          "HC_proportion_detected",
          "HC_n",
          "PD_proportion_detected",
          "PD_n"
        )
      pdhc.dsgv.stats <- rbind(pdhc.dsgv.stats, row2add)
    }
  } else {
    cat("Feature : ", var, " has insufficient data \n")
  }
}
# This correctly converts to numeric
pdhc.dsgv.stats$p.value <- as.numeric(levels(pdhc.dsgv.stats$p.value))[pdhc.dsgv.stats$p.value]
pdhc.dsgv.stats$estimate <- as.numeric(levels(pdhc.dsgv.stats$estimate))[pdhc.dsgv.stats$estimate]
pdhc.dsgv.stats$CI.lower <- as.numeric(levels(pdhc.dsgv.stats$CI.lower))[pdhc.dsgv.stats$CI.lower]
pdhc.dsgv.stats$CI.upper <- as.numeric(levels(pdhc.dsgv.stats$CI.upper))[pdhc.dsgv.stats$CI.upper]
pdhc.dsgv.stats$HC_proportion_detected <- as.numeric(levels(pdhc.dsgv.stats$HC_proportion_detected))[pdhc.dsgv.stats$HC_proportion_detected]
pdhc.dsgv.stats$HC_n <- as.numeric(levels(pdhc.dsgv.stats$HC_n))[pdhc.dsgv.stats$HC_n]
pdhc.dsgv.stats$PD_proportion_detected <- as.numeric(levels(pdhc.dsgv.stats$PD_proportion_detected))[pdhc.dsgv.stats$PD_proportion_detected]
pdhc.dsgv.stats$PD_n <- as.numeric(levels(pdhc.dsgv.stats$PD_n))[pdhc.dsgv.stats$PD_n]

# Save Analysis
# write.csv(pdhc.dsgv.stats, file = 'files/DSGV_PDvHC_FishersExact_02_22_21.csv')






#-----------------------
#####   PD v PC   #####
#-----------------------
pdpc.dsgv.stats <- tibble()
n <- 0
for (i in 1:ncol(pdpc.vals)) {
  var <- colnames(pdpc.vals[i])
  print(var)
  print(n)
  n <- n + 1
  PC.prop <- filter(pdpc_summary, sv == var & group == "PC")$ratio_detected[1]
  PC.n <- filter(pdpc_summary, sv == var & group == "PC")$n[1]
  PD.prop <- filter(pdpc_summary, sv == var & group == "PD")$ratio_detected[1]
  PD.n <- filter(pdpc_summary, sv == var & group == "PD")$n[1]

  if (var != "group") {
    if (!is.na(PD.n) & !is.na(PC.n) > 0 &
      (PC.prop > 0 | PD.prop > 0) &
      (PC.prop < 1 | PD.prop < 1)) {
      tst <- pdpc.vals %>%
        dplyr::count(group, pdpc.vals[[i]]) %>%
        na.omit() %>%
        pivot_wider(
          values_from = n,
          names_from = "pdpc.vals[[i]]"
        ) %>%
        replace(is.na(.), 0) %>%
        dplyr::rename(hit = `1`, nohit = `0`) %>%
        column_to_rownames(var = "group")

      fish <- fisher.test(tst)
      cat("p-value: ", fish$p.value, "\n")

      row2add <-
        cbind(
          var,
          fish$estimate,
          fish$p.value,
          fish$conf.int[[1]],
          fish$conf.int[[2]],
          PC.prop,
          PC.n,
          PD.prop,
          PD.n
        )
      colnames(row2add) <-
        c(
          "sv",
          "estimate",
          "p.value",
          "CI.lower",
          "CI.upper",
          "PC_proportion_detected",
          "PC_n",
          "PD_proportion_detected",
          "PD_n"
        )
      pdpc.dsgv.stats <- rbind(pdpc.dsgv.stats, row2add)
    }
  } else {
    cat("Feature : ", var, " has insufficient data \n")
  }
}
# This correctly converts to numeric
pdpc.dsgv.stats$p.value <- as.numeric(levels(pdpc.dsgv.stats$p.value))[pdpc.dsgv.stats$p.value]
pdpc.dsgv.stats$estimate <- as.numeric(levels(pdpc.dsgv.stats$estimate))[pdpc.dsgv.stats$estimate]
pdpc.dsgv.stats$CI.lower <- as.numeric(levels(pdpc.dsgv.stats$CI.lower))[pdpc.dsgv.stats$CI.lower]
pdpc.dsgv.stats$CI.upper <- as.numeric(levels(pdpc.dsgv.stats$CI.upper))[pdpc.dsgv.stats$CI.upper]
pdpc.dsgv.stats$PC_proportion_detected <- as.numeric(levels(pdpc.dsgv.stats$PC_proportion_detected))[pdpc.dsgv.stats$PC_proportion_detected]
pdpc.dsgv.stats$PC_n <- as.numeric(levels(pdpc.dsgv.stats$PC_n))[pdpc.dsgv.stats$PC_n]
pdpc.dsgv.stats$PD_proportion_detected <- as.numeric(levels(pdpc.dsgv.stats$PD_proportion_detected))[pdpc.dsgv.stats$PD_proportion_detected]
pdpc.dsgv.stats$PD_n <- as.numeric(levels(pdpc.dsgv.stats$PD_n))[pdpc.dsgv.stats$PD_n]

# Save Analysis
# write.csv(pdpc.dsgv.stats, file = 'files/DSGV_PDvPC_FishersExact_02_22_21.csv')





#------------------------------------------------
# Data Visualization
#------------------------------------------------

# # Explore analysis
# ggplot(pdhc.dsgv.stats,
#        aes(x=log2((PD_proportion_detected + 1)/(HC_proportion_detected + 1)),
#            y=log2((PD_n+1)/(HC_n+1)),
#            color=p.value)) +
#   geom_point() +
#   scale_color_viridis_c(direction = -1) +
#   theme_bw()
#
# ggplot(pdhc.dsgv.stats) +
#   geom_point(aes(x=PD_proportion_detected,
#                  y=HC_proportion_detected,
#                  color=p.value)) +
#   scale_color_viridis_c(direction = -1) +
#   theme_bw()

# Ranked SV plot
r1 <-
  pdhc.dsgv.stats %>%
  mutate(PD_deletion_ratio = (1 - PD_proportion_detected + 1) / (1 - HC_proportion_detected + 1)) %>%
  ggplot(aes(
    x = reorder(sv, PD_deletion_ratio),
    y = log2(PD_deletion_ratio),
    color = p.value
  )) +
  labs(
    x = "Ranked Structural Variants",
    y = expression(paste(log[2], "[ PD/HC deletion-ratio ]"))
  ) +
  scale_color_viridis_c(
    option = "magma", direction = -1,
    na.value = "#000000",
    limits = c(0, 0.25)
  ) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )


# Violin Plot
pdhc.dsgv.stats.plot <-
  pdhc.dsgv.stats %>%
  mutate(PD_deletion_ratio = (1 - PD_proportion_detected + 1) / (1 - HC_proportion_detected + 1)) %>%
  mutate(lab = if_else(abs(log2(PD_deletion_ratio)) >= 0.2 & abs(p.value) <= 0.005, "yes", "no"))

v1 <-
  pdhc.dsgv.stats.plot %>%
  ggplot(aes(
    x = log2(PD_deletion_ratio),
    y = -log10(p.value),
    fill = p.value
  )) +
  labs(
    x = expression(paste(log[2], "[ PD/HC detection-ratio ]")),
    y = expression(paste(-log[10], "[ p-value ]")),
    fill = "p-value"
  ) +
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    na.value = "#000000",
    limits = c(0, 0.25)
  ) +
  geom_point(shape = 21) +
  geom_text_repel(
    data = pdhc.dsgv.stats.plot %>%
      filter(lab == "yes"),
    aes(
      x = log2(PD_deletion_ratio),
      y = -log10(p.value), label = sv
    ),
    size = 2, seed = 42, segment.alpha = 0.2, force = 1.5
  ) +
  theme_bw() +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dotted") +
  geom_hline(yintercept = c(-log10(0.05)), linetype = "dotted") +
  theme(panel.grid = element_blank())





# ggplot(pdpc.dsgv.stats) +
#   geom_point(aes(x=log2((PD_proportion_detected + 1)/(PC_proportion_detected + 1)),
#                  y=log2((PD_n+1)/(PC_n+1)),
#                  color=p.value)) +
#   scale_color_viridis_c(direction = -1) +
#   theme_bw()
#
# ggplot(pdpc.dsgv.stats) +
#   geom_point(aes(x=PD_proportion_detected,
#                  y=PC_proportion_detected,
#                  color=p.value)) +
#   scale_color_viridis_c(direction = -1) +
#   theme_bw()

# Ranked SV plot
r2 <-
  pdpc.dsgv.stats %>%
  mutate(PD_deletion_ratio = (1 - PD_proportion_detected + 1) / (1 - PC_proportion_detected + 1)) %>%
  ggplot(aes(
    x = reorder(sv, PD_deletion_ratio),
    y = log2(PD_deletion_ratio),
    color = p.value
  )) +
  labs(
    x = "Ranked Structural Variants",
    y = expression(paste(log[2], "[ PD/PC deletion-ratio ]"))
  ) +
  scale_color_viridis_c(
    option = "magma", direction = -1,
    na.value = "#000000",
    limits = c(0, 0.25)
  ) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )


# Violin Plot
pdpc.dsgv.stats.plot <-
  pdpc.dsgv.stats %>%
  mutate(PD_deletion_ratio = (1 - PD_proportion_detected + 1) / (1 - PC_proportion_detected + 1)) %>%
  mutate(lab = if_else(abs(log2(PD_deletion_ratio)) >= 0.2 & abs(p.value) <= 0.005, "yes", "no"))

v2 <-
  pdpc.dsgv.stats.plot %>%
  ggplot(aes(
    x = log2(PD_deletion_ratio),
    y = -log10(p.value),
    fill = p.value
  )) +
  labs(
    x = expression(paste(log[2], "[ PD/PC detection-ratio ]")),
    y = expression(paste(-log[10], "[ p-value ]")),
    fill = "p-value"
  ) +
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    na.value = "#000000",
    limits = c(0, 0.25)
  ) +
  geom_point(shape = 21) +
  geom_text_repel(
    data = pdpc.dsgv.stats.plot %>%
      filter(lab == "yes"),
    aes(
      x = log2(PD_deletion_ratio),
      y = -log10(p.value), label = sv
    ),
    size = 2, seed = 42, segment.alpha = 0.2, force = 1
  ) +
  theme_bw() +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dotted") +
  geom_hline(yintercept = c(-log10(0.05)), linetype = "dotted") +
  theme(panel.grid = element_blank())


c1 <- cowplot::plot_grid(r1, v1,
  align = "hv", axis = "lrtb",
  ncol = 1, rel_heights = c(1, 2.25)
)
c2 <- cowplot::plot_grid(r2, v2,
  align = "hv", axis = "lrtb",
  ncol = 1, rel_heights = c(1, 2.25)
)
c3 <- cowplot::plot_grid(c1, c2, nrow = 1, align = "hv")


ggsave(c3,
  filename = "data/SGVFinderViz/deletion_sgvs.png",
  height = 8, width = 14, dpi = 1200
)











#------------------------
# VF MaAsLin2 MODELS
#------------------------
# load("files/VFs_PhyloseqObj.RData")

# Distribution Sanity Check
dat.vsgv %>%
  microbiome::abundances() %>%
  distribution_sanity2()

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc <- subset_samples(dat.vsgv, donor_group != "HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  # microbiome::transform("compositional") %>%
  LowVarianceFilter(filter.percent = 0)


# PD v HC paired abundance data
dat_pdhc <- NULL
dat_pdhc <- subset_samples(dat.vsgv, paired != "No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  # microbiome::transform("compositional") %>%
  LowVarianceFilter(filter.percent = 0)


######## Format Metadata
# Run Metadata pre-processing function
process_meta(dat.vsgv)
df_input_metadata <- env %>% column_to_rownames(var = "donor_id")

process_meta(dat_pdpc)
df_input_metadata_pdpc <- env %>% column_to_rownames(var = "donor_id")
df_input_metadata_pdpc$description <- factor(df_input_metadata_pdpc$description,
  levels = c("PD Patient", "Population Control")
)
process_meta(dat_pdhc)
df_input_metadata_pdhc <- env %>% column_to_rownames(var = "donor_id")
df_input_metadata_pdhc$description <- factor(df_input_metadata_pdhc$description,
  levels = c("PD Patient", "Household Control")
)
# set file path
wkd <- getwd()

############  PD v PC - VFs ############

fit_data <- Maaslin2(
  input_data = df_input_data_pdpc,
  input_metadata = df_input_metadata_pdpc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/VSGV_PDvPC_maaslin2_output"),
  fixed_effects = c("description"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  cores = 6
)

############  PD v HC paired - VFs ############

fit_data <- Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/VSGV_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "NONE",
  cores = 6
)







# OLD Chi-Squared Analysis  - Inappropriate for small N

#
# pdhc.dsgv.stats <- tibble()
# for (var in unique(pdhc_summary$sv)) {
#
#   HC.prop <-  filter(pdhc_summary, sv == var & group == "HC")$ratio_detected
#   HC.n <-  filter(pdhc_summary, sv == var & group == "HC")$n
#   PD.prop <- filter(pdhc_summary, sv == var & group == "PD")$ratio_detected
#   PD.n <- filter(pdhc_summary, sv == var & group == "PD")$n
#
#   if (length(PD.n) > 0 & length(HC.n) > 0) {
#
#     cat("testing: ", var, "\n")
#
#     two.pzt <- prop.test(
#       x = c(HC.prop, PD.prop),
#       n = c(HC.n, PD.n),
#       p = NULL, alternative = "two.sided",
#       correct = TRUE)
#
#     row2add <- cbind(var, as.numeric(two.pzt$statistic), as.numeric(two.pzt$p.value),
#                      as.numeric(two.pzt$conf.int[1]), as.numeric(two.pzt$conf.int[2]),
#                      HC.prop, HC.n, PD.prop, PD.n)
#     colnames(row2add) <- c("sv", "statistic", "p.value", "CI.lower", "CI.upper",
#                            "HC_proportion_detected", "HC_n", "PD_proportion_detected", "PD_n")
#     pdhc.dsgv.stats <- rbind(pdhc.dsgv.stats, row2add)
#
#   }
# }
# # This correctly converts to numeric
# pdhc.dsgv.stats$p.value <- as.numeric(levels(pdhc.dsgv.stats$p.value))[pdhc.dsgv.stats$p.value]
# pdhc.dsgv.stats$statistic <- as.numeric(levels(pdhc.dsgv.stats$statistic))[pdhc.dsgv.stats$statistic]
# pdhc.dsgv.stats$CI.lower <- as.numeric(levels(pdhc.dsgv.stats$CI.lower))[pdhc.dsgv.stats$CI.lower]
# pdhc.dsgv.stats$CI.upper <- as.numeric(levels(pdhc.dsgv.stats$CI.upper))[pdhc.dsgv.stats$CI.upper]
# pdhc.dsgv.stats$HC_proportion_detected <- as.numeric(levels(pdhc.dsgv.stats$HC_proportion_detected))[pdhc.dsgv.stats$HC_proportion_detected]
# pdhc.dsgv.stats$HC_n <- as.numeric(levels(pdhc.dsgv.stats$HC_n))[pdhc.dsgv.stats$HC_n]
# pdhc.dsgv.stats$PD_proportion_detected <- as.numeric(levels(pdhc.dsgv.stats$PD_proportion_detected))[pdhc.dsgv.stats$PD_proportion_detected]
# pdhc.dsgv.stats$PD_n <- as.numeric(levels(pdhc.dsgv.stats$PD_n))[pdhc.dsgv.stats$PD_n]


#
# pdpc.dsgv.stats <- tibble()
# for (var in unique(pdpc_summary$sv)) {
#
#   PC.prop <-  filter(pdpc_summary, sv == var & group == "PC")$ratio_detected
#   PC.n <-  filter(pdpc_summary, sv == var & group == "PC")$n
#   PD.prop <- filter(pdpc_summary, sv == var & group == "PD")$ratio_detected
#   PD.n <- filter(pdpc_summary, sv == var & group == "PD")$n
#
#   if (length(PD.n) > 0 & length(PC.n) > 0) {
#
#     cat("testing: ", var, "\n")
#
#     two.pzt <- prop.test(
#       x = c(PC.prop, PD.prop),
#       n = c(PC.n, PD.n),
#       p = NULL, alternative = "two.sided",
#       correct = TRUE)
#
#     row2add <- cbind(var, two.pzt$statistic, two.pzt$p.value, two.pzt$conf.int[1], two.pzt$conf.int[2],
#                      PC.prop, PC.n, PD.prop, PD.n)
#     colnames(row2add) <- c("sv", "statistic", "p.value", "CI.lower", "CI.upper",
#                            "PC_proportion_detected", "PC_n", "PD_proportion_detected", "PD_n")
#     pdpc.dsgv.stats <- rbind(pdpc.dsgv.stats, row2add)
#
#   }
# }
# # This correctly converts to numeric
# pdpc.dsgv.stats$p.value <- as.numeric(levels(pdpc.dsgv.stats$p.value))[pdpc.dsgv.stats$p.value]
# pdpc.dsgv.stats$statistic <- as.numeric(levels(pdpc.dsgv.stats$statistic))[pdpc.dsgv.stats$statistic]
# pdpc.dsgv.stats$CI.lower <- as.numeric(levels(pdpc.dsgv.stats$CI.lower))[pdpc.dsgv.stats$CI.lower]
# pdpc.dsgv.stats$CI.upper <- as.numeric(levels(pdpc.dsgv.stats$CI.upper))[pdpc.dsgv.stats$CI.upper]
# pdpc.dsgv.stats$PC_proportion_detected <- as.numeric(levels(pdpc.dsgv.stats$PC_proportion_detected))[pdpc.dsgv.stats$PC_proportion_detected]
# pdpc.dsgv.stats$PC_n <- as.numeric(levels(pdpc.dsgv.stats$PC_n))[pdpc.dsgv.stats$PC_n]
# pdpc.dsgv.stats$PD_proportion_detected <- as.numeric(levels(pdpc.dsgv.stats$PD_proportion_detected))[pdpc.dsgv.stats$PD_proportion_detected]
# pdpc.dsgv.stats$PD_n <- as.numeric(levels(pdpc.dsgv.stats$PD_n))[pdpc.dsgv.stats$PD_n]
