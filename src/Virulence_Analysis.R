#### Virulence/Pathogenicity Analysis


######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")
source("src/Stats_funcs.R")
# source("src/create_phyloseq_obj_ShortBRED.R") # Builds VF Dataframe from individual TSVs

#-------------------------------------------------------------------------


#-------------------------------------------------------------------------

# Temporary shortcut to use previously created VF DF
load("files/shortbred.RData")
# TSS - Normalization
shortbreddata[-1] <- apply(shortbreddata[-1], 2, function(x) x / sum(x))
# colSums(shortbreddata[-1])

merged <- shortbreddata %>% melt()
merged$MBI_Sample_ID <- substr(merged$variable, 1, 10)
# ArcSinSqrt - Transformation
merged$value <- asin(sqrt(merged$value))

# Create DF of Summed VF values
mag <- aggregate(merged$value, by = list(MBI_Sample_ID = merged$MBI_Sample_ID), FUN = sum) %>%
  dplyr::rename(value = x)
m <- read.csv(file = "files/metadata_keys.csv", header = TRUE)
m2 <- dat %>%
  meta() %>%
  dplyr::select(donor_group, Paired, id)
mag2 <- left_join(mag, m, by = "MBI_Sample_ID")
mag2 <- left_join(mag2, m2, by = "id")

#-------------------------------------------------------------------------
# Total VFBD Stats
#-------------------------------------------------------------------------
# Summary statistics
mag2 %>%
  group_by(donor_group) %>%
  get_summary_stats(value, type = "mean_sd")

# Check assumptions - Outliers
mag2 %>%
  group_by(donor_group) %>%
  identify_outliers(value)

# Check assumptions - Normality
model <- lm(value ~ donor_group, data = mag2)
ggqqplot(residuals(model))
# Shapiro-Wilk test of normality
shapiro_test(residuals(model))
ggqqplot(mag2, "value", facet.by = "donor_group")

# Check assumptions - Homogneity of variance
plot(model, 1) # No obvious relationship is good
mag2 %>% levene_test(value ~ donor_group)
# Test reveals there is a significant diff in group variance


process_meta(dat)
mag3 <- mag2 %>%
  dplyr::select(c(id, value)) %>%
  dplyr::rename(donor_id = id, total_VF = value)
mag3$donor_id <- gsub("_", ".", mag3$donor_id)

vfdb_df <- left_join(env, mag3, by = "donor_id", )
# Create new paired column with only household pairs and NAs for rest
vfdb_df$Paired.plot <- as.numeric(levels(vfdb_df$Paired))[vfdb_df$Paired]


#-------------------------------
# PLOT STATS
#-------------------------------
observed.PdPC <- lm.PdPc(metadf = vfdb_df, metric = "total_VF")
observed.PdHC <- lmm.PdHc(metadf = vfdb_df, metric = "total_VF")
## Pull p-values
observed.PdPC.pval <- summary(observed.PdPC)$coefficients["descriptionPopulation Control", "Pr(>|t|)"]
observed.PdHC.pval <- summary(observed.PdHC)$tTable["descriptionHousehold Control", "p-value"]

#-------------------------------
#  TOTAL VFBD PLOT
#-------------------------------

basic_cols <- c("PD Patient" = "#bfbfbf", "Population Control" = "#ed7d31", "Household Control" = "#5b9bd5")
mag2$description <- factor(mag2$description, levels = c("Population Control", "PD Patient", "Household Control"))
cols <- basic_cols
title <- "Total Virulence Factors (VFDB)"
ylabel <- "Summated Normalized Abundance"

set.seed(123)
v1 <- ggplot(data = mag2, aes(x = description, y = value)) +
  geom_boxplot(aes(color = description), outlier.alpha = 0, width = 0.9) +
  geom_point(aes(fill = description),
    position = position_jitterdodge(jitter.width = 1),
    shape = 21, size = 1.5, alpha = 0.8
  ) +
  theme_classic() +
  ggtitle(title) +
  labs(y = ylabel) +
  scale_color_manual(values = cols, name = "Group") +
  scale_fill_manual(values = cols, name = "Group") +
  geom_signif(
    comparisons = list(c("PD Patient", "Household Control")),
    annotations = sig_mapper(observed.PdHC.pval)
  ) +
  geom_signif(
    comparisons = list(c("Population Control", "PD Patient")),
    annotations = sig_mapper(observed.PdPC.pval)
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(colour = "grey"),
    legend.position = "none"
  )
v1

ggsave(v1,
  filename = "figures/Figure_4/Total_Virulence_Factors_VFDB.svg",
  height = 5,
  width = 4
)

#--------------------------------------------------------------------------------------------------
#  TOTAL VFs x Alpha-Diversity : Linear Regression
#--------------------------------------------------------------------------------------------------

# Prep Metadata
meta <- m2 %>% dplyr::select(donor_group, Paired, id)
env <- meta(dat)
env <- left_join(env, mag2, by = "id")
env$description <- factor(env$description, levels = c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels = c("PC", "PD", "HC"))

# Prep abundance table
dat_alpha <- microbiome::transform(dat, "compositional")
## Calculate Alpha Diversity Metrics and add cols to df
vfdb_df$observed <- microbiome::alpha(abundances(dat_alpha), "observed")$observed
vfdb_df$shannon <- microbiome::alpha(abundances(dat_alpha), "shannon")$diversity_shannon
vfdb_df$evenness <- evenness(abundances(dat_alpha), "simpson")$simpson

df <- vfdb_df
x <- vfdb_df$total_VF
y <- vfdb_df$observed
color <- env$donor_group
fill <- env$donor_group
ylabel <- " "
title <- " "

PD.col <- "#bfbfbf"
PC.col <- "#ed7d31"
HC.col <- "#5b9bd5"


p1 <- ggplot(vfdb_df, aes(x = total_VF, y = shannon, color = color, fill = fill)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(fill = fill), shape = 21, size = 2, alpha = 0.9) +
  # labs(x="Filtered Sample Reads", y=ylabel, fill="group") +
  theme_bw() +
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
  #              color = color), label.x = 15e6, label.y.npc=0.25, hjust = 0) +
  ggtitle(title) +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    legend.background = element_blank(),
    legend.position = c(0.9, 0.5)
  )
p1



#--------------------------------------------------------------------------------------------------
#  TOTAL VFs x Beta-Diversity : Linear Regression
#--------------------------------------------------------------------------------------------------

mag4pcoa <- mag2 %>%
  dplyr::select(c(id, value)) %>%
  dplyr::rename(total_VF = value)

obj_clr <- microbiome::transform(dat.KOs.slim, "clr")
iDist <- phyloseq::distance(obj_clr, method = "euclidean")
iMDS <- ordinate(obj_clr, "MDS", distance = iDist)
#  PCoA for Axis 1 and 2

p <- plot_ordination(obj_clr, iMDS, color = "description", axes = c(1, 2))
pcoa.data <- p$data
df.pcoa <- left_join(pcoa.data, mag4pcoa, by = "id")
df.pcoa$donor_group <- factor(df.pcoa$donor_group, levels = c("PC", "PD", "HC"))



df <- df.pcoa
y <- df.pcoa$total_VF
x <- df.pcoa$Axis.1
color <- df.pcoa$donor_group
fill <- df.pcoa$donor_group


p2 <- ggplot(df, aes(x = x, y = y, color = color, fill = fill)) +
  geom_smooth(method = "lm", se = F) +
  geom_point(aes(fill = fill), shape = 21, size = 2, alpha = 0.9) +
  # labs(x="Filtered Sample Reads", y=paste0("PCoA Axis 1: ", feature)) +
  theme_bw() +
  # stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
  #              color = color), label.x = 1.75e7, label.y.npc=0.25, hjust = 0) +
  ggtitle(title) +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
p2





#-------------------------------------------------------------------
# Test individual Features for Significance - MaAsLin2 Models
#-------------------------------------------------------------------


#------------------------
# VF MaAsLin2 MODELS
#------------------------
load("files/VFs_PhyloseqObj.RData")

# Distribtution Sanity Check
dat.VFs %>%
  microbiome::abundances() %>%
  distribution_sanity2()

# PD v PC abundance data
dat_pdpc <- NULL
dat_pdpc <- subset_samples(dat.VFs, donor_group != "HC")
# Plot Variance Estimate
PlotVariance(dat_pdpc)
df_input_data_pdpc <- dat_pdpc %>%
  microbiome::transform("compositional") %>%
  LowVarianceFilter(filter.percent = 0.1)


# PD v HC PAIRED abundance data
dat_pdhc <- NULL
dat_pdhc <- subset_samples(dat.VFs, Paired != "No")
# Plot Variance Estimate
PlotVariance(dat_pdhc)
df_input_data_pdhc <- dat_pdhc %>%
  microbiome::transform("compositional") %>%
  LowVarianceFilter(filter.percent = 0.1)


######## Format Metadata
# Run Metadata pre-processing function
process_meta(dat.VFs)
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
  output = paste0(wkd, "/data/MaAsLin2_Analysis/VF_PDvPC_maaslin2_output"),
  fixed_effects = c("description", "host_age_factor", "sex", "host_body_mass_index"),
  min_prevalence = 0,
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)

############  PD v HC Paired - VFs ############

fit_data <- Maaslin2(
  input_data = df_input_data_pdhc,
  input_metadata = df_input_metadata_pdhc,
  output = paste0(wkd, "/data/MaAsLin2_Analysis/VF_PDvHC_maaslin2_output"),
  min_prevalence = 0,
  random_effects = c("Paired"),
  fixed_effects = c("description"),
  analysis_method = "LM",
  normalization = "NONE",
  transform = "AST",
  cores = 1
)
