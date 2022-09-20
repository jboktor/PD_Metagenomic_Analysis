# Meta-Analysis


######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")
source("src/create_phlyoseq_obj_CMG.R")


# Note First run create_phyloseq_obj_CMG

#------------------------------------------------------------------------------------------
# SEQUENCING READS
#------------------------------------------------------------------------------------------

dsranking <- disease_meta %>%
  dplyr::group_by(dataset_name) %>%
  dplyr::summarize(mediandepth = median(number_reads) / 1e6) %>%
  dplyr::mutate(dsorder = rank(mediandepth)) %>%
  dplyr::arrange(dsorder)

meta_reads <- disease_meta %>%
  mutate(ds = factor(disease_meta$dataset_name, levels = dsranking$dataset_name)) %>%
  ggplot(aes(ds, x = dataset_name, y = log10(number_reads / 1e6), fill = study_condition)) +
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition),
    position = position_jitterdodge(jitter.width = 0.4), shape = 21, size = 1, alpha = 0.4
  ) +
  theme_bw() +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  labs(y = "log10 Read Depth (millions)", fill = "Study Condition") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

#------------------------------------------------------------------------------------------
# ALPHA-DIVERSITY
#------------------------------------------------------------------------------------------

## Calculate Alpha Diversity Metrics and add cols to df
disease_meta$Observed <- microbiome::alpha(abundances(physeq), "observed")$observed
disease_meta$Shannon <- microbiome::alpha(abundances(physeq), "shannon")$diversity_shannon
disease_meta$Evenness <- evenness(abundances(physeq), "simpson")$simpson

# Plot histograms to get a sense of data distribution
par(mfrow = c(1, 3))
hist(disease_meta$Observed, main = "observed OTUs", xlab = "", breaks = 10)
hist(disease_meta$Shannon, main = "Shannon diversity", xlab = "", breaks = 10)
hist(disease_meta$Evenness, main = "Simpson's evenness", xlab = "", breaks = 10)


meta_alpha_observed <- ggplot(data = disease_meta, aes(x = dataset_name, y = Observed, fill = study_condition)) +
  theme_bw() +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition),
    position = position_jitterdodge(jitter.width = 0.4), shape = 21, size = 1, alpha = 0.4
  ) +
  labs(y = "Observed Species", fill = "Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

meta_alpha_shannon <- ggplot(data = disease_meta, aes(x = dataset_name, y = Shannon, fill = study_condition)) +
  theme_bw() +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition),
    position = position_jitterdodge(jitter.width = 0.4), shape = 21, size = 1, alpha = 0.4
  ) +
  labs(y = "Shannon Diversity", fill = "Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

meta_alpha_simpsons <- ggplot(data = disease_meta, aes(x = dataset_name, y = Evenness, fill = study_condition)) +
  theme_bw() +
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition),
    position = position_jitterdodge(jitter.width = 0.4), shape = 21, size = 1, alpha = 0.4
  ) +
  labs(y = "Simpson's Evenness", fill = "Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

alp <- cowplot::plot_grid(meta_reads, meta_alpha_observed, meta_alpha_shannon, meta_alpha_simpsons,
  nrow = 2, align = "v", axis = "tblr"
)
ggsave(alp,
  filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_Alpha.png",
  height = 10, width = 14
)



#------------------------------------------------------------------------------------------
# BETA-DIVERSITY
#------------------------------------------------------------------------------------------
# AITCHISONS DISTANCE (EUCLIDIAN DIST ON CLR TRANSFORMED DATA)
#---------------------------------------------------------
obj_clr <- microbiome::transform(physeq, "clr")
iDist <- phyloseq::distance(obj_clr, method = "euclidean")
iMDS <- phyloseq::ordinate(obj_clr, "MDS", distance = iDist)

p <- plot_ordination(obj_clr, iMDS, color = "studyID", axes = c(1, 2))
df12 <- p$data
df12$study_condition <- factor(df12$study_condition,
  levels = c(
    "control", "parkinsons", "ACVD", "BD", "AS", "CDI", "T1D", "hypertension",
    "metabolic_syndrome", "IBD", "T2D", "cirrhosis"
  )
)
df12$dataset_name <- sub(".metaphlan_bugs_list.stool", "", df12$studyID) %>% factor()

df12 <- df12 %>% mutate(disease_category = if_else(study_condition == "control", "Healthy",
  if_else(study_condition == "parkinsons", "Neurodegenerative",
    if_else(study_condition == "AS" | study_condition == "T1D" |
      study_condition == "BD" | study_condition == "IBD", "Inflammatory/Autoimmune",
    if_else(study_condition == "T2D" | study_condition == "metabolic_syndrome" |
      study_condition == "hypertension" | study_condition == "ACVD", "Nutritional",
    if_else(study_condition == "CDI", "Infectious",
      if_else(study_condition == "cirrhosis", "Alcohol related", "ERROR")
    )
    )
    )
  )
))
#---------------------------------------------------------
# Adding metadata columns
#---------------------------------------------------------
sample_data(physeq)$dataset_name <- df12$dataset_name
sample_data(physeq)$disease_category <- df12$disease_category
#---------------------------------------------------------


p1 <-
  ggplot(df12, aes(Axis.1, Axis.2, fill = study_condition, color = study_condition)) +
  geom_point(shape = 21, size = 1.5, alpha = 0.8) +
  geom_point(
    data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
    aes(Axis.1, Axis.2), shape = 21, size = 1.5, color = "black"
  ) +
  theme_bw() +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2]) * 100, digits = 2), "%)")) +
  labs(fill = "Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


p2 <-
  ggplot(df12, aes(Axis.1, Axis.2)) +
  geom_point(shape = 21, size = 1.5, alpha = 0.8, aes(fill = dataset_name, color = dataset_name)) +
  geom_point(
    data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
    aes(Axis.1, Axis.2), shape = 21, size = 1.5, color = "black"
  ) +
  theme_bw() +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2]) * 100, digits = 2), "%)")) +
  labs(fill = "Dataset") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


p3 <-
  ggplot(df12, aes(Axis.1, Axis.2)) +
  geom_point(shape = 21, size = 1.5, alpha = 0.8, aes(fill = country, color = country)) +
  geom_point(
    data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
    aes(Axis.1, Axis.2), shape = 21, size = 1.5, color = "black"
  ) +
  theme_bw() +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2]) * 100, digits = 2), "%)")) +
  labs(fill = "Country of Origin") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


p4 <-
  ggplot(df12, aes(Axis.1, Axis.2)) +
  geom_point(shape = 21, size = 1.5, alpha = 0.8, aes(fill = disease_category, color = disease_category)) +
  geom_point(
    data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
    aes(Axis.1, Axis.2), shape = 21, size = 1.5, color = "black"
  ) +
  theme_bw() +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2]) * 100, digits = 2), "%)")) +
  labs(fill = "Disease Category") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )



ord <- cowplot::plot_grid(p1, p2, p3, p4, nrow = 2, align = "hv")
ggsave(ord,
  filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_Beta_Diversity.png",
  height = 9, width = 14
)





#---------------------------------------------------------------

p <- physeq
m <- "euclidean"

# Make Melted Distance Matrix
wu <- phyloseq::distance(obj_clr, method = "euclidean")
wu.m <- melt(as.matrix(wu))
# Exclude intra-sample distances
wu.m <- wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
# Pull metadata of interest
meta_selected <- meta(p) %>%
  rownames_to_column() %>%
  dplyr::select(c(rowname, dataset_name, study_condition, country, disease_category)) %>%
  mutate_if(is.factor, as.character)

# Add group name for Var1 Column
colnames(meta_selected)[1:2] <- c("Var1", "Type1")
wu.meta_selected <- left_join(wu.m, meta_selected[1:2], by = "Var1")
# Add group name for Var2 Column - add additional metadata for target group
colnames(meta_selected)[1:2] <- c("Var2", "Type2")
wu.meta_selected <- left_join(wu.meta_selected, meta_selected, by = "Var2")
# Select only distances to Population control
wu.meta_selected <- filter(wu.meta_selected, Type1 == "Boktor_Mazmanian")


v.dataset <-
  wu.meta_selected %>%
  mutate(Type2 = fct_reorder(Type2, value)) %>%
  ggplot(aes(x = Type2, y = value)) +
  theme_minimal() +
  geom_violin(aes(fill = Type2), draw_quantiles = c(0.5), trim = T, alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab("Aitchisons Distance to Internal Dataset") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Group") +
  scale_color_manual(values = col_vector) +
  scale_fill_manual(values = col_vector) +
  # scale_colour_tableau(palette = "Classic 20") +
  # scale_fill_tableau(palette = "Classic 20") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
v.dataset


v.study_condition <-
  wu.meta_selected %>%
  mutate(study_condition = fct_reorder(study_condition, value)) %>%
  ggplot(aes(x = study_condition, y = value)) +
  theme_minimal() +
  geom_violin(aes(color = study_condition, fill = study_condition), draw_quantiles = c(0.5), trim = T, alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab("Aitchisons Distance to Internal Dataset") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(fill = "Group") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
v.study_condition

v.disease_category <-
  wu.meta_selected %>%
  dplyr::filter(disease_category != "Healthy") %>%
  mutate(disease_category = fct_reorder(disease_category, value)) %>%
  ggplot(aes(x = disease_category, y = value)) +
  theme_minimal() +
  geom_violin(aes(color = disease_category, fill = disease_category), draw_quantiles = c(0.5), trim = T, alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab("Aitchisons Distance to Internal Dataset") +
  labs(fill = "Group") +
  scale_colour_tableau(palette = "Classic 10") +
  scale_fill_tableau(palette = "Classic 10") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
v.disease_category

v.country <-
  wu.meta_selected %>%
  mutate(country = fct_reorder(country, value)) %>%
  ggplot(aes(x = country, y = value)) +
  theme_minimal() +
  geom_violin(aes(color = country, fill = country), draw_quantiles = c(0.5), trim = T, alpha = 0.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  ) +
  ylab("Aitchisons Distance to Internal Dataset") +
  labs(fill = "Group") +
  # scale_color_manual(values = col_vector) +
  # scale_fill_manual(values = col_vector) +
  scale_colour_tableau(palette = "Classic 10") +
  scale_fill_tableau(palette = "Classic 10") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )
v.country

distance_plots <- cowplot::plot_grid(v.dataset, v.study_condition,
  v.disease_category, v.country,
  nrow = 2, align = "hv"
)
ggsave(distance_plots,
  filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_Distance_Plots.png",
  height = 9, width = 12
)



# library(RColorBrewer)
n <- 60
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=sample(col_vector, n))
