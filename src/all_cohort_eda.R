# all cohort exploration

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
load("files/low_quality_samples.RData")

remove_dats()
# load_data("Merged_ML")
load("files/Phyloseq_Merged_ML.RData")
# load("files/Phyloseq_Merged_ML_Rarefied.RData")

# tst <- snm_dat(datObj = phyloseq_objs[["Pathways.slim"]], 
#                ObjName = "Pathways.slim", 
#                distance_metric = "euclidean")

# objInput <- dat.species.SNM # ADJUST AS NEEDED
objInput <- phyloseq_objs_rare[["Species"]] %>% 
  subset_samples(donor_id %ni% low_qc[[1]])

df.meta <- objInput %>%
  meta() %>% 
  mutate_at(c("total_reads", "total_nucleotide_aligned", 
              "total_translated_aligned", "total_species"), as.numeric)

# Sequencing depth vs nucleotide alignment
df.meta %>% 
  ggplot() +
  geom_point(aes(x=total_reads, y=total_translated_aligned, fill = cohort ), shape = 21) +
  scale_fill_d3() +
  scale_x_log10() +
  scale_y_log10()
  
# Sequencing depth distribution
df.meta %>% 
  ggplot() +
  geom_density(aes(x=total_reads, color = cohort), size = 1.3) +
  scale_color_d3() +
  scale_x_log10() +
  my_clean_theme()

# Dimension Reduction
obj <- subset_samples(objInput, donor_id %ni% low_qc[[1]])
iDist <- phyloseq::distance(obj, method="euclidean")
iMDS  <- ordinate(obj, "MDS", distance=iDist)
plot_ordination(obj, iMDS, color="cohort", axes = c(1, 2)) +
  scale_color_d3() +
  my_clean_theme()


# Binning Sample abundance
obj <- phyloseq_objs_rare[["KOs.slim"]] %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) # ALTER DAT AS NEEDED
KO.slim.abun <- obj %>% abundances() %>% as.data.frame()
KO.slim.sums <- as.data.frame(colSums(KO.slim.abun)) %>% 
  dplyr::rename("KOs.slim.sums" = "colSums(KO.slim.abun)") %>% 
  rownames_to_column(var = "donor_id")
df.plot.meta <- left_join(df.meta, KO.slim.sums, by = "donor_id")
df.plot.meta %>%
  ggplot(aes(x=cohort, y=KOs.slim.sums)) +
  geom_violin( draw_quantiles = T) +
  geom_point()

df.quants <- objInput %>% 
  abundances() %>% as.data.frame() %>% rownames_to_column("feature") %>% 
  pivot_longer(!feature, names_to = "donor_id", values_to = "abundance") %>%
  mutate(quants = as.numeric(cut(
    abundance,
    breaks = quantile(abundance, na.rm = T, probs = seq(0.1, 1, by = 0.05)),
    labels = seq(0.1, 0.95, by = 0.05)*100,
    include.lowest = T
  )))

# # Species
# df.quants %>% 
#   left_join(df.meta, by = "donor_id") %>% 
#   ggplot(aes(x=quants, y=abundance), fill = "black") +
#   geom_point(
#     aes(color = quants),
#     position = position_jitterdodge(jitter.width = .4),
#     size = 1,
#     alpha = 0.2) +
#   labs(x = "Percentile of Distribution", y = "CLR Abundance", fill = "quantiles") +
#   facet_grid(PD~cohort) +
#   my_clean_theme()

ecdf_plot <-
  df.quants %>%
  left_join(df.meta, by = "donor_id") %>%
  ggplot(aes(x = abundance, colour = cohort)) +
  stat_ecdf(geom = "step", pad = FALSE) +
  my_clean_theme() +
  labs(y = "ECDF") +
  scale_color_d3() +
  theme(legend.position = c(0.7, 0.4))
ecdf_plot
  




  
#_______________________________________________________________________________

dat.obj <- 
  dat.genus %>% 
  subset_samples(donor_id %ni% low_qc[[1]]) %>%
  core(detection = 0, prevalence = 0.1)

# Create Metadata Column for Cohort x Donor Group
sample_data(dat.obj)$cohort_donor_group <- 
  paste(sample_data(dat.obj)$cohort, sample_data(dat.obj)$donor_group)
# Abundance filter for top 30 Genera
dat.top.30 <- dat.obj %>% 
  get_top_taxa(n=30, relative = TRUE, discard_other = F, other_label = "Other")

dat.top.30 %>% abundance_heatmap(treatment = "cohort_donor_group")

barcols <- c(
  "#386cb0",
  "#7fc97f",
  "#beaed4",
  "#fdc086",
  "#ffff99",
  "#f0027f",
  "#bf5b17",
  "#666666",
  "#7fc97f",
  "#beaed4",
  "#beaed6"
)

# Plot all Samples
barplt1 <- 
  fantaxtic_bar(
    dat.top.30,
    color_by = "Order",
    label_by = "Genus",
    other_label = "Other",
    # facet_type = "grid",
    facet_by = "cohort",
    grid_by = "PD", 
    facet_cols = 4,
    order_alg = "hclust",
    # base_color = "#5b9bd5", 
    palette = barcols
    # color_levels = barcol_ID
  ) +
  labs(y = "Relative Abundance") +
  theme(axis.text.x = element_blank())
barplt1 

ggsave(barplt1,
       filename = "data/Community_Composition/Quality_Control/all_studies_top30Genera.svg",
       width = 14, height = 7)

  
# Pathways


# KOs


# eggNOGs