# GMMxGBM Figures :

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
source("src/daf_functions.R")
# source("src/omixer-rmpR_setup.R")
load("files/low_quality_samples.RData")
load("files/Phyloseq_Merged/GMMs_PhyloseqObj.RData")
wkd <- getwd()


# ArcSinSqrt - Transformation
df.gmm[-1] <- asin(sqrt(df.gmm[-1]))

module.levels <- read_tsv("files/GMM/GMM.hierarchy.v1.07.tsv", col_names = T)

H1.list <- c()
cnt <- 1
# Inspired by https://github.com/zellerlab/crc_meta
for (i in unique(module.levels$HL1)) {

  # Pull Modules belonging to each HL1 class
  label <- gsub(" ", ".", i)
  var <- module.levels %>%
    filter(HL1 == i) %>%
    pull(Module)
  assign(label, var)
  print(label)

  temp_df <- df.gmm %>%
    dplyr::filter(module %in% var)

  temp_df2 <- melt(colSums(temp_df[-1]))
  colnames(temp_df2) <- label
  H1.list[[cnt]] <- temp_df2
  cnt <- cnt + 1
}



df.seed <- H1.list[[1]] %>%
  rownames_to_column(var = "id")
df.seed2 <- H1.list[[2]] %>%
  rownames_to_column(var = "id")
metabolism_df <- left_join(df.seed, df.seed2, by = "id")
# construct data frame with summated
# values for each HL1 level per sample
for (i in 3:length(H1.list)) {
  df2add <- H1.list[[i]] %>%
    rownames_to_column(var = "id")
  metabolism_df <- left_join(metabolism_df, df2add,
    by = "id"
  )
}

metabolism_df <- metabolism_df %>%
  group_col_from_ids(metabolism_df$id) %>%
  column_to_rownames(var = "id")


#---------------------------------------------------------------------------------
#                            Create GMM HL1 Phyloseq Objects
#---------------------------------------------------------------------------------

### Metadata
my_sample_data <- meta(dat.species) %>% sample_data()
my_GMM.HL1.ab_table <- metabolism_df %>%
  dplyr::select(-group) %>%
  otu_table(taxa_are_rows = F)
dat.GMMs.HL1 <- phyloseq(my_GMM.HL1.ab_table, my_sample_data)
print(dat.GMMs.HL1)
save(dat.GMMs.HL1, file = "files/Phyloseq_Merged/GMMs.HL1_PhyloseqObj.RData")


#-------------------------------------------------------------------------------
####                             GMM HL1 MaAsLin2                                ####
#-------------------------------------------------------------------------------

dat.object <- maaslin_prep(dat.GMMs.HL1)
dat_pdpc <- subset_samples(dat.object, donor_group != "HC")
variance_plot(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  variance_filter(0)
dat_pdhc <- subset_samples(dat.object, paired != "No")
variance_plot(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  variance_filter(0)
fit_models(
  dat = dat.object,
  obj.name = "GMMs.HL1",
  dat_pdpc = dat_pdpc,
  dat_pdhc = dat_pdhc,
  df_input_data_pdhc = df_input_data_pdhc,
  df_input_data_pdpc = df_input_data_pdpc,
  cores = 6,
  plot_scatter = T
)


#
# #------------------
# # PD v PC abundance data
# dat_pdpc <- NULL
# dat_pdpc = subset_samples(dat.GMMs.HL1, donor_group !="HC")
# df_input_data_pdpc <- dat_pdpc %>%
#   microbiome::transform("compositional")
# PlotVariance(dat_pdpc)
# df_input_data_pdpc <- dat_pdpc %>%
#   LowVarianceFilter(filter.percent = 0)
#
# #------------------
# # PD v HC PAIRED abundance data
# dat_pdhc <- NULL
# dat_pdhc = subset_samples(dat.GMMs.HL1, Paired !="No")
# df_input_data_pdhc <- dat_pdhc %>%
#   LowVarianceFilter(filter.percent = 0)
#
#
# #------------------Format Metadata
# # Run Metadata pre-processing function
#
# process_meta(dat_pdpc)
# df_input_metadata_pdpc <- env %>% column_to_rownames(var = "donor_id")
# df_input_metadata_pdpc$description <- factor(df_input_metadata_pdpc$description,
#                                              levels = c("PD Patient", "Population Control"))
#
# process_meta(dat_pdhc)
# df_input_metadata_pdhc <- env %>% column_to_rownames(var = "donor_id")
# df_input_metadata_pdhc$description <- factor(df_input_metadata_pdhc$description,
#                                              levels = c("PD Patient", "Household Control"))
#
# #-------------------------------------------------------------------------
# # MaAsLin2 Models
# #--------------------------------------------------------------------------
#
# ############  PD v PC - GMMs HL1  ############
#
# fit_data = Maaslin2(
#   input_data = df_input_data_pdpc,
#   input_metadata = df_input_metadata_pdpc,
#   output = paste0(wkd, "/data/MaAsLin2_Analysis/GMMs.HL1_PDvPC_maaslin2_output"),
#   fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
#   min_prevalence = 0,
#   analysis_method = "LM",
#   normalization = "NONE",
#   transform = "NONE",
#   cores = 1
# )
#
# ############  PD v HC Paired  - GMMs HL1  ############
#
# fit_data = Maaslin2(
#   input_data = df_input_data_pdhc,
#   input_metadata = df_input_metadata_pdhc,
#   output = paste0(wkd, "/data/MaAsLin2_Analysis/GMMs.HL1_PDvHC_maaslin2_output"),
#   min_prevalence = 0,
#   random_effects = c("Paired"),
#   fixed_effects = c("description"),
#   analysis_method = "LM",
#   normalization = "NONE",
#   transform = "NONE",
#   cores = 1
# )

lev <- "GMMs.HL1"
### Read-in Maaslin Files - all features used in significance testing
Maas.pd.pc <- read_tsv(paste0(
  "data/MaAsLin2_Analysis/",
  lev, "_PDvPC_maaslin2_output/all_results.tsv"
),
col_names = T
) %>%
  filter(value == "Population Control")

Maas.pd.hc <- read_tsv(paste0(
  "data/MaAsLin2_Analysis/",
  lev, "_PDvHC_maaslin2_output/all_results.tsv"
),
col_names = T
) %>%
  filter(value == "Household Control")

#--------------------------------------------------------------------------
#                         Plot Metabolic Overview
#--------------------------------------------------------------------------

metabolism_df2 <- metabolism_df %>%
  pivot_longer(-group, names_to = "feature", values_to = "Total_Abundance")

Maas.pd.pc.stats <- Maas.pd.pc %>%
  dplyr::select(c(feature, qval))
Maas.pd.pc.stats$group <- "PC"

Maas.pd.hc.stats <- Maas.pd.hc %>%
  dplyr::select(c(feature, qval))
Maas.pd.hc.stats$group <- "HC"

Maas.stats <- rbind(Maas.pd.pc.stats, Maas.pd.hc.stats)
Maas.stats$qval <- round(Maas.stats$qval, digits = 3)
metabolism_df3 <- left_join(metabolism_df2, Maas.stats, by = c("feature", "group"))
metabolism_df3$qval.print <- sig.symbol.generator(
  Column = metabolism_df3$qval,
  porq = "q", shh = F
)

## Add Max values of each group as a column
feats <- c()
Max.values <- c()
for (process in unique(metabolism_df3$feature)) {
  group.max <- max(filter(metabolism_df3, feature == process)$Total_Abundance)
  feats <- c(feats, process)
  Max.values <- c(Max.values, group.max)
}
featmax <- data.frame("feature" = feats, Max.values)
metabolism_df3 <- left_join(metabolism_df3, featmax, by = c("feature"))

group.cols.fill <- c("PC" = "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
group.cols.rims <- c("PC" = "#c15811", "PD" = "#494949", "HC" = "#2e75b5")
metabolism_df3$group <- factor(metabolism_df3$group, levels = c("PC", "PD", "HC"))
metlabels <- c(
  "amino.acid.degradation" = "Amino acid\n degradation",
  "organic.acid.metabolism" = "Organic acid\n metabolism",
  "glycoprotein.degradation" = "Glycoprotein\n degradation",
  "carbohydrate.degradation" = "Carbohydrate\n degradation",
  "lipid.degradation" = "Lipid\n degradation",
  "central.metabolism" = "Central\n Metabolism",
  "gas.metabolism" = "Gas\n Metabolism",
  "amines.and.polyamines.degradation" = "Amines & Polyamines\n degradation",
  "alcohol.metabolism" = "Alcohol\n Metabolism",
  "inorganic.nutrient.metabolism" = "Inorganic Nutrient\n Metabolism"
)

### Note: No Normalization or Transformation applied -
# Data is Sum of counts for HL1 level of GMMs derived from Relative Abundance counts of
# Kegg Orthologs from HUMANn2

metabolism_overview <-
  ggplot(data = metabolism_df3, aes(x = feature, y = Total_Abundance)) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.2) +
  geom_point(aes(fill = group, color = group),
    position = position_jitterdodge(jitter.width = 0.2),
    shape = 21, size = 1.5, alpha = 0.6
  ) +
  theme_bw() +
  geom_text(
    data = filter(metabolism_df3, group == "PC"),
    aes(
      x = feature, y = (Max.values + Max.values * 0.15),
      label = qval.print
    ), size = 2.9, check_overlap = TRUE, nudge_x = -0.25
  ) +
  geom_text(
    data = filter(metabolism_df3, group == "HC"),
    aes(
      x = feature, y = (Max.values + Max.values * 0.15),
      label = qval.print
    ), size = 2.9, check_overlap = TRUE, nudge_x = 0.25
  ) +
  facet_wrap(. ~ feature,
    scale = "free", nrow = 2,
    labeller = labeller(feature = metlabels)
  ) +
  labs(y = "Summated Normalized Abundance") +
  scale_color_manual(values = group.cols.rims, name = "Group") +
  scale_fill_manual(values = group.cols.fill, name = "Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.length.x = unit(0, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey"),
    legend.background = element_blank()
  )

metabolism_overview


# g <- ggplot_gtable(ggplot_build(metabolism_overview))
# stripr <- which(grepl('strip-t', g$layout$name))
# fills <- c("#c3bc3f", "#bb7693", "#baa094", "#a9b5ae", "#767676",
#   "#6388b4", "#ffae34", "#ef6f6a", "#8cc2ca", "#55ad89")
# k <- 1
# for (i in stripr) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# metabolism_overview2 <- gridExtra::arrangeGrob(g)

ggsave(metabolism_overview,
  filename = "figures/Figure_4/Metabolism_Overview_Boxplot_colored.svg",
  width = 10, height = 5
)
