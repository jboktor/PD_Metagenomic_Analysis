# Feature Correlation Analysis

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
load("files/low_quality_samples.RData")
load_all_cohorts()

#-------------------------------------------------------------------------------
#####                      Correlation Functions                           ##### 
#-------------------------------------------------------------------------------

corr_abund_prep <- function(obj,  obj.name, cohort) {
  
  ### Read-in MaAsLin2 output
  Maas.pd.pc.sig <- read_tsv(
    paste0("data/MaAsLin2_Analysis/", cohort, "/", obj.name, 
           "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Population Control") %>% 
    filter(qval < 0.25)
  
  Maas.pd.hc.sig <- read_tsv(
    paste0("data/MaAsLin2_Analysis/",  cohort, "/", obj.name, 
           "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Household Control") %>% 
    filter(qval < 0.25)
  
  features <- full_join(Maas.pd.pc.sig, Maas.pd.hc.sig, by = "feature") %>% 
    dplyr::select("feature")
  
  abundance_tbl <- 
    obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(rowname %in% features[[1]]) %>% 
    column_to_rownames() %>% 
    t() %>%
    as.data.frame()
  
  return(abundance_tbl)
  
}

corr_loop <- function(metadata, abundance) {
  corr_output <- tibble()
  for (metavar in colnames(metadata)) {
    cat("Calculating correlations for: ", metavar[[1]], "\n")
    for (feature in colnames(abundance)) {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      row2add <-
        cbind(
          "metadata" = metavar,
          "feature" = feature,
          "rho" = spearman$estimate[[1]],
          "S" = spearman$statistic[[1]],
          "n" = length(na.omit(metadata[[metavar]])),
          "p" = spearman$p.value[[1]]
        )
      corr_output <- rbind(corr_output, row2add)
    }
  }
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <- 
    corr_output %>% 
    na.omit() %>% 
    mutate(across(all_of(statvars), as.character),
           across(all_of(statvars), as.numeric)) %>% 
    decode_rfriendly_rows(passed_column = "feature") %>% 
    select(-feature) %>%
    dplyr::rename("feature"="fullnames") %>% 
    dplyr::relocate(feature, .after = metadata)
  corr_output$q <- p.adjust(corr_output$p, method = 'BH')
  return(corr_output)
}


#-------------------------------------------------------------------------------
#####                        Correlation XY Plot                          ##### 
#-------------------------------------------------------------------------------

corr_xy <- function(obj, corr_obj, feature_var, metadata_var){
  
  #' Function creates a scatter plot of a given feature and a metadata column
  abund <- obj %>% 
    abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "feature") %>% 
    decode_rfriendly_rows(passed_column = "feature") %>% 
    select(-feature) %>%
    column_to_rownames(var = "fullnames") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "donor_id")
  df.plot <- obj %>% 
    meta() %>% 
    process_meta(cohort = "Merged") %>% 
    left_join(abund, by = "donor_id")

  stat_col <-
    corr_obj %>%
    dplyr::filter(feature == sym(feature_var)) %>%
    dplyr::filter(metadata == sym(metadata_var))
  stat_title <-
    paste0(
    "Spearman's Rho: ", round(stat_col$rho, digits = 3), "\n",
    "p-value: ", round(stat_col$p, digits = 4),
    "  FDR: ", round(stat_col$q, digits = 4)
    )
  cat("Rho: ", stat_col$rho, "\n")
  cat("P-value: ", stat_col$p, "\n")
  cat("Q-value: ", stat_col$q, "\n")

  df.plot %>%
    drop_na(metadata_var) %>%
    ggplot(aes(x = .data[[feature_var]], y = .data[[metadata_var]])) +
    geom_point(aes(fill = donor_group, color = donor_group), shape=21, alpha = 1) +
    geom_smooth(method = lm, color="darkgrey", linetype="dotted", se = F) +
    theme_bw() +
    labs(x = feature_var, y = metadata_var, title = stat_title) +
    scale_fill_manual(values = cols.pdpchc) +
    scale_color_manual(values = cols.pdpchc.rim) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12))
}

#-------------------------------------------------------------------------------
#####                      Correlation Heatmap                           ##### 
#-------------------------------------------------------------------------------

corr_heatmap <- function(corr.df){
  
  corr.df <- corr.df %>% 
    mutate(metadata = as.character(metadata))
  
  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>% 
    dplyr::filter(p < 0.05) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n)) %>%
    top_n(n = 30, wt = n) %>% 
    slice_head(n = 30)
  
  # Filter correlation df for top 30 features
  corr.df <- corr.df %>% 
    filter(feature %in% corr.df.top$feature)
  
  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <- 
    corr.df %>% 
    dplyr::select(feature, metadata, rho) %>% 
    pivot_wider(names_from = feature, values_from = rho, values_fill = NA) %>% 
    column_to_rownames(var = "metadata")
  
  # Hierarchical clustering of Rho values for features & metadata
  meta.dendro <- 
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  meta.dendro.plot <- ggdendrogram(data = meta.dendro, rotate = TRUE)
  feature.dendro <- 
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature.dendro.plot <- ggdendrogram(data = feature.dendro, rotate = TRUE)
  
  ### Reorder Heatmap axis using order of dendrograms
  feature.order <- order.dendrogram(feature.dendro)
  metadata.order <- order.dendrogram(meta.dendro)
  corr.df$feature <- factor(x = corr.df$feature,
                            levels = corr.df$feature[feature.order],
                            ordered = TRUE)
  corr.df <- corr.df %>% 
    dplyr::arrange(feature)
  corr.df$metadata <- factor(x = corr.df$metadata,
                             levels = corr.df$metadata[metadata.order],
                             ordered = TRUE)
  
  ### Plot Heatmap
  h1 <- 
    corr.df %>% 
    mutate(siglabel = if_else(q < 0.25, "*", "")) %>%
    ggplot(aes(x = metadata, y = feature, fill = rho)) +
    geom_tile() + 
    geom_text(aes(label=siglabel), size=8,vjust = 0.77) +
    labs(fill = "Spearman correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    theme(axis.text.x = element_text(angle=45, hjust =1),
          legend.position = "top",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(3, 3, 3, 6), "cm"))
  
  print(h1)
  return(h1)
}

#-------------------------------------------------------------------------------
#####                      Correlation Analysis                           ##### 
#-------------------------------------------------------------------------------

# Metadata Selection
corr.meta.clinical <-
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  process_meta(cohort = "Merged") %>%
  select(contains(
    c(motor_severity_scores_summary,
      clinical_variables
    )
  ))
corr.meta.diet <-
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  process_meta(cohort = "Merged") %>%
  select(contains(
    c(diet
    )
  ))

#-------------------------------------------------------------------------------
#####                              Species                               ##### 
#-------------------------------------------------------------------------------
# Features with Associations 
# corr.abund.species <-
#   corr_abund_prep(obj.name = "Species",
#                   obj = dat.species,
#                   cohort = "Merged")
corr.abund.species <- 
  dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  abundances() %>% 
  t() %>%
  as.data.frame()

corrs.clinical.species <-
  corr_loop(metadata = corr.meta.clinical, abundance = corr.abund.species)
corrs.diet.species <-
  corr_loop(metadata = corr.meta.diet, abundance = corr.abund.species)

# Heatmaps
h1 <- corr_heatmap(corrs.clinical.species)
h2 <- corr_heatmap(corrs.diet.species)
ggsave(h1, filename = "data/Correlations/Species/Clinical/Heatmap_Species_clinical.svg", 
       height = 6, width = 9)
ggsave(h2, filename = "data/Correlations/Species/Diet/Heatmap_Species_diet.svg", 
       height = 9, width = 9)

# Scatter plots
corrs.clinical.species.sig <- 
  corrs.clinical.species %>% 
  filter(p < 0.05)
corrs.diet.species.sig <- 
  corrs.diet.species %>% 
  filter(p < 0.05)

for (corr in 1:nrow(corrs.clinical.species.sig)) {
  cor_row <- corrs.clinical.species.sig[corr, ]
  p <- corr_xy(obj = dat.species, corrs.clinical.species,
          feature_var = cor_row$feature[[1]],
          metadata_var = as.character(cor_row$metadata[[1]]) )
  print(p)
  plot.name = paste0("data/Correlations/Species/Clinical/", 
                     cor_row$feature[[1]], "_VS_",
                     as.character(cor_row$metadata[[1]]), ".png")
  ggsave(p, filename = plot.name, height = 4, width = 5)
}

for (corr in 1:nrow(corrs.diet.species.sig)) {
  cor_row <- corrs.diet.species.sig[corr, ]
  p <- corr_xy(obj = dat.species, corrs.diet.species,
               feature_var = cor_row$feature[[1]],
               metadata_var = as.character(cor_row$metadata[[1]]) )
  print(p)
  plot.name = paste0("data/Correlations/Species/Diet/", 
                     cor_row$feature[[1]], "_VS_",
                     as.character(cor_row$metadata[[1]]), ".png")
  ggsave(p, filename = plot.name, height = 4, width = 5)
}

#-------------------------------------------------------------------------------
#####                              Pathways                               ##### 
#-------------------------------------------------------------------------------

corr.abund.pathways <-
  corr_abund_prep(obj.name = "Pathways.slim",
                  obj = dat.path.slim,
                  cohort = "Merged")

corrs.clinical.pathways.slim <-
  corr_loop(metadata = corr.meta.clinical, abundance = corr.abund.pathways)
corrs.diet.pathways.slim <-
  corr_loop(metadata = corr.meta.diet, abundance = corr.abund.pathways)

# Heatmaps
h3 <- corr_heatmap(corrs.clinical.pathways.slim)
h4 <- corr_heatmap(corrs.diet.pathways.slim)
ggsave(h3, filename = "data/Correlations/Pathways.slim/Clinical/Heatmap_Pathways.slim_clinical.svg", 
       height = 10, width = 13)
ggsave(h4, filename = "data/Correlations/Pathways.slim/Diet/Heatmap_Pathways.slim_diet.svg", 
       height = 9, width = 11)

# Scatter plots
corrs.clinical.pathways.slim.sig <- 
  corrs.clinical.pathways.slim %>% 
  filter(q < 0.5)
corrs.diet.pathways.slim.sig <- 
  corrs.diet.pathways.slim %>% 
  filter(p < 0.01)

for (corr in 1:nrow(corrs.clinical.pathways.slim.sig)) {
  cor_row <- corrs.clinical.pathways.slim.sig[corr, ]
  p <- corr_xy(obj = dat.path.slim, corrs.clinical.pathways.slim,
               feature_var = cor_row$feature[[1]],
               metadata_var = as.character(cor_row$metadata[[1]]) )
  print(p)
  plot.name = paste0("data/Correlations/Pathways.slim/Clinical/", 
                     cor_row$feature[[1]], "_VS_",
                     as.character(cor_row$metadata[[1]]), ".png")
  ggsave(p, filename = plot.name, height = 4, width = 5)
}

for (corr in 1:nrow(corrs.diet.pathways.slim.sig)) {
  cor_row <- corrs.diet.pathways.slim.sig[corr, ]
  p <- corr_xy(obj = dat.path.slim, corrs.diet.pathways.slim.sig,
               feature_var = cor_row$feature[[1]],
               metadata_var = as.character(cor_row$metadata[[1]]) )
  print(p)
  plot.name = paste0("data/Correlations/Pathways.slim/Diet//", 
                     cor_row$feature[[1]], "_VS_",
                     as.character(cor_row$metadata[[1]]), ".png")
  ggsave(p, filename = plot.name, height = 4, width = 5)
}

#-------------------------------------------------------------------------------
#####                              KOs                               ##### 
#-------------------------------------------------------------------------------

corr.abund.KOs <-
  corr_abund_prep(obj.name = "KOs.slim",
                  obj = dat.KOs.slim,
                  cohort = "Merged")

corrs.clinical.KOs.slim <-
  corr_loop(metadata = corr.meta.clinical, abundance = corr.abund.KOs)
corrs.diet.KOs.slim <-
  corr_loop(metadata = corr.meta.diet, abundance = corr.abund.KOs)

h5 <- corr_heatmap(corrs.clinical.KOs.slim)
h6 <- corr_heatmap(corrs.diet.KOs.slim)
ggsave(h5, filename = "data/Correlations/KOs.slim//Clinical/Heatmap_KOs.slim_clinical.svg", 
       height = 10, width = 12)
ggsave(h6, filename = "data/Correlations/KOs.slim/Diet/Heatmap_KOs.slim_diet.svg", 
       height = 10, width = 13)

# Scatter plots
corrs.clinical.KOs.slim.sig <- 
  corrs.clinical.KOs.slim %>% 
  filter(p < 0.001)
corrs.diet.KOs.slim.sig <- 
  corrs.diet.KOs.slim %>% 
  filter(p < 0.001)

for (corr in 1:nrow(corrs.clinical.KOs.slim.sig)) {
  cor_row <- corrs.clinical.KOs.slim.sig[corr,]
  p <- corr_xy(
    obj = dat.KOs.slim,
    corrs.clinical.KOs.slim,
    feature_var = cor_row$feature[[1]],
    metadata_var = as.character(cor_row$metadata[[1]])
  )
  print(p)
  plot.name = paste0(
    "data/Correlations/KOs.slim/Clinical/",
    cor_row$feature[[1]],
    "_VS_",
    as.character(cor_row$metadata[[1]]),
    ".png"
  )
  ggsave(p,
         filename = plot.name,
         height = 4,
         width = 5)
}

for (corr in 1:nrow(corrs.diet.KOs.slim.sig)) {
  cor_row <- corrs.diet.KOs.slim.sig[corr,]
  p <- corr_xy(
    obj = dat.KOs.slim,
    corrs.diet.KOs.slim,
    feature_var = cor_row$feature[[1]],
    metadata_var = as.character(cor_row$metadata[[1]])
  )
  print(p)
  plot.name = paste0(
    "data/Correlations/KOs.slim/Diet//",
    cor_row$feature[[1]],
    "_VS_",
    as.character(cor_row$metadata[[1]]),
    ".png"
  )
  ggsave(p,
         filename = plot.name,
         height = 4,
         width = 5)
}
#-------------------------------------------------------------------------------
#####                              GOs                               ##### 
#-------------------------------------------------------------------------------

corr.abund.GOs <-
  corr_abund_prep(obj.name = "GOs.slim",
                  obj = dat.GOs.slim,
                  cohort = "Merged")

corrs.clinical.GOs.slim <-
  corr_loop(metadata = corr.meta.clinical, abundance = corr.abund.GOs)
corrs.diet.GOs.slim <-
  corr_loop(metadata = corr.meta.diet, abundance = corr.abund.GOs)

h7 <- corr_heatmap(corrs.clinical.GOs.slim)
h8 <- corr_heatmap(corrs.diet.GOs.slim)
ggsave(h7, filename = "data/Correlations/GOs.slim/Clinical/Heatmap_GOs.slim_clinical.svg", 
       height = 10, width = 13)
ggsave(h8, filename = "data/Correlations/GOs.slim/Diet/Heatmap_GOs.slim_diet.svg", 
       height = 10, width = 13)

#-------------------------------------------------------------------------------
#####                              PFAMs                               ##### 
#-------------------------------------------------------------------------------

corr.abund.PFAMs <-
  corr_abund_prep(obj.name = "PFAMs.slim",
                  obj = dat.PFAMs.slim,
                  cohort = "Merged")

corrs.clinical.PFAMs.slim <-
  corr_loop(metadata = corr.meta.clinical, abundance = corr.abund.PFAMs)
corrs.diet.PFAMs.slim <-
  corr_loop(metadata = corr.meta.diet, abundance = corr.abund.PFAMs)

h9 <- corr_heatmap(corrs.clinical.PFAMs.slim)
h10 <- corr_heatmap(corrs.diet.PFAMs.slim)
ggsave(h9, filename = "data/Correlations/PFAMs.slim/Clinical/Heatmap_PFAMs.slim_clinical.svg", 
       height = 10, width = 12)
ggsave(h10, filename = "data/Correlations/PFAMs.slim/Diet/Heatmap_PFAMs.slim_diet.svg", 
       height = 10, width = 13)

corr_xy(obj = dat.PFAMs.slim, corr_obj = corrs.clinical.PFAMs.slim,
        feature_var = "PF12675: Protein of unknown function (DUF3795)",
        metadata_var = "family_history_pd_degree_relative")

# Scatter plots
corrs.clinical.PFAMs.slim.sig <-
  corrs.clinical.PFAMs.slim %>%
  filter(p < 0.002)
# corrs.diet.PFAMs.slim.sig <-
#   corrs.diet.PFAMs.slim %>%
#   filter(p < 0.001)

for (corr in 1:nrow(corrs.clinical.PFAMs.slim.sig)) {
  cor_row <- corrs.clinical.PFAMs.slim.sig[corr,]
  p <- corr_xy(
    obj = dat.PFAMs.slim,
    corrs.clinical.PFAMs.slim,
    feature_var = cor_row$feature[[1]],
    metadata_var = as.character(cor_row$metadata[[1]])
  )
  print(p)
  plot.name = paste0(
    "data/Correlations/PFAMs.slim/Clinical/",
    cor_row$feature[[1]],
    "_VS_",
    as.character(cor_row$metadata[[1]]),
    ".png"
  )
  ggsave(p,
         filename = plot.name,
         height = 4,
         width = 5)
}

#-------------------------------------------------------------------------------
#####                      Saving as Excel files                           ##### 
#-------------------------------------------------------------------------------
clinical_corrs <-
  list(Species = corrs.clinical.species,
       Pathways = corrs.clinical.pathways.slim,
       Kegg_Orthology = corrs.clinical.KOs.slim,
       Gene_Ontology = corrs.clinical.GOs.slim,
       PFAMs = corrs.clinical.PFAMs.slim)

dietary_corrs <-
  list(Species = corrs.diet.species,
       Pathways = corrs.diet.pathways.slim,
       Kegg_Orthology = corrs.diet.KOs.slim,
       Gene_Ontology = corrs.diet.GOs.slim,
       PFAMs = corrs.diet.PFAMs.slim)

openxlsx::write.xlsx(clinical_corrs, file = 'files/Correlations/Clinical_Correlations.xlsx')
openxlsx::write.xlsx(dietary_corrs, file = 'files/Correlations/Dietary_Correlations.xlsx')



# 
# ggsave(h1, filename = "data/Correlations/Clinical_variables.png",
#        width = 18, height = 7)
# 


















