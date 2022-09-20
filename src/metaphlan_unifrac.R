# Unifrac data

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
load("files/low_quality_samples.RData")

tree_file <- "/Users/josephboktor/Documents/PD_Metagenomic_Analysis/files/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk"
mpa_tree <- ape::read.tree(tree_file)
mpa_tree$tip.label <- gsub(".+\\|s__", "", mpa_tree$tip.label)



#-------------------------------------------------------------------------------
#####                              Merged                                  ##### 
#-------------------------------------------------------------------------------

remove_dats()
load_data("Merged")

bugs.species <- dat.species %>% 
  abundances() %>% 
  as.data.frame()

filt_tree <- ape::keep.tip(mpa_tree, intersect(rownames(bugs.species),mpa_tree$tip.label))
filt_mpa_table <- bugs.species[filt_tree$tip.label,] / 100.0

# Build Phyloseq data with tree
phy_table <- otu_table(filt_mpa_table, taxa_are_rows=T)
phy_meta <- meta(dat.species) %>% sample_data()
phy_tax <- tax_table(as.matrix(as.data.frame(tax_table(dat.species))[filt_tree$tip.label,]))
phy_tree <- phy_tree(filt_tree)
dat.species.tree <- phyloseq(phy_tax, phy_table, phy_meta, phy_tree)
dat.species.tree

save(dat.species.tree, file = "files/Phyloseq_Merged/Species_tree_PhyloseqObj.RData")

#-------------------------------------------------------------------------------
#####                              TBC                                  ##### 
#-------------------------------------------------------------------------------
remove_dats()
load_data("TBC")

bugs.species <- dat.species %>% 
  abundances() %>% 
  as.data.frame()

filt_tree <- ape::keep.tip(mpa_tree, intersect(rownames(bugs.species),mpa_tree$tip.label))
filt_mpa_table <- bugs.species[filt_tree$tip.label,] / 100.0

# Build Phyloseq data with tree
phy_table <- otu_table(filt_mpa_table, taxa_are_rows=T)
phy_meta <- meta(dat.species) %>% sample_data()
phy_tax <- tax_table(as.matrix(as.data.frame(tax_table(dat.species))[filt_tree$tip.label,]))
phy_tree <- phy_tree(filt_tree)
dat.species.tree <- phyloseq(phy_tax, phy_table, phy_meta, phy_tree)
dat.species.tree

save(dat.species.tree, file = "files/Phyloseq_TBC/Species_tree_PhyloseqObj.RData")

#-------------------------------------------------------------------------------
#####                              RUSH                                  ##### 
#-------------------------------------------------------------------------------
remove_dats()
load_data("RUSH")

bugs.species <- dat.species %>% 
  abundances() %>% 
  as.data.frame()

filt_tree <- ape::keep.tip(mpa_tree, intersect(rownames(bugs.species),mpa_tree$tip.label))
filt_mpa_table <- bugs.species[filt_tree$tip.label,] / 100.0

# Build Phyloseq data with tree
phy_table <- otu_table(filt_mpa_table, taxa_are_rows=T)
phy_meta <- meta(dat.species) %>% sample_data()
phy_tax <- tax_table(as.matrix(as.data.frame(tax_table(dat.species))[filt_tree$tip.label,]))
phy_tree <- phy_tree(filt_tree)
dat.species.tree <- phyloseq(phy_tax, phy_table, phy_meta, phy_tree)
dat.species.tree

save(dat.species.tree, file = "files/Phyloseq_RUSH/Species_tree_PhyloseqObj.RData")


#-------------------------------------------------------------------------------
#####                      Unweighted Unifrac                            ##### 
#-------------------------------------------------------------------------------

remove_dats()
load("files/Phyloseq_Merged/Species_tree_PhyloseqObj.RData")

cohort = "Merged"
dist = "unifrac" # unifrac # wunifrac
z <- c("Species")
i <- dat.species.tree
cnt <- 1
load_betadiv_colors()

obj_dist <- dat.species.tree %>%
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  microbiome::transform("compositional")

iDist <- phyloseq::distance(obj_dist, method=dist)
dist_label <- "Unweighted UniFrac"

# Calculate ordination
iMDS  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist)

# ---------------------------------------------------
#  PCoA for Axis 1 and 2
# ---------------------------------------------------

plot_ordination(obj_dist, iMDS, color="cohort", axes = c(1, 2))
p <- plot_ordination(obj_dist, iMDS, color="description", axes = c(1, 2))

df12 = p$data
df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
p <- ggplot(df12, aes(Axis.1, Axis.2, fill = description, color=description))
p <- p + geom_point(shape=21, size=5, alpha=0.7)
ord <- p + 
  theme_bw() + 
  labs(fill="Donor Group") +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Donor Group") +
  scale_fill_manual(values = cols.pdpchc.dark) +
  scale_color_manual(values = cols.pdpchc.rim) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())

# PCoA data to fit ridgeline Plots
my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2

# ---------------------------------------------------
#  Boxplot - Axis 1
# ---------------------------------------------------

r1 <- ggplot(df12, aes(x = Axis.1, y = description)) +
  geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 1), "%)")) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.xrange) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  ggtitle(paste0(dist_label, " Distance PCoA")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")

# ---------------------------------------------------
#  Boxplot - Axis 2
# ---------------------------------------------------

r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
  geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.yrange2) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_blank(),
        legend.position = "none") +
  coord_flip()

# ---------------------------------------------------
#  Beta Diversity Violin Plots
# ---------------------------------------------------
p <- obj_dist
s <- "donor_id"
d <- "description"

# Make Melted Distance Matrix 
wu.m <- melt(as.matrix(iDist))
# Exclude intra-sample distances 
wu.m <- wu.m %>% filter(as.character(Var1) != as.character(Var2)) %>% 
  mutate_if(is.factor, as.character)
# Pull metadata of interest
sd <- meta(p) %>% dplyr::select(all_of(s), all_of(d)) %>% mutate_if(is.factor, as.character)
# Add group name for Var1 Column
colnames(sd) <- c("Var1", "Type1")
wu.sd <- left_join(wu.m, sd, by = "Var1")
# Add group name for Var2 Column
colnames(sd) <- c("Var2", "Type2")
wu.sd <- left_join(wu.sd, sd, by = "Var2")
# Select only distances to Population control
wu.sd <- filter(wu.sd, Type1 == "Population Control")

# Specifying comparisons for analysis
my_comparisons <-
  list(c("Household Control", "PD Patient"),
       c("Population Control", "PD Patient"))
wu.sd$Type2 <-
  factor(wu.sd$Type2,
         levels = c("Population Control", "PD Patient", "Household Control"))

v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
  geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size= 0.2, alpha = 0.5, cex = 0.2) +
  geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha=0) +
  theme_classic() +
  ylab(paste(dist_label, "Distance")) +
  labs(fill="Group") +
  scale_color_manual(values = cols.pdpchc) +
  scale_fill_manual(values = cols.pdpchc) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0) + 
  ggtitle(paste0(dist_label," Distance to Population Control\n", z[cnt], " abundance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")


# ---------------------------------------------------
#  Assemble Plots
# ---------------------------------------------------
cat(paste0("Assembling summary plots for: " , z[cnt], "\n"))
ord.plot <- ord + theme(legend.position = "none")
cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 5), rel_widths = c(5, 1),  align = "vh")
cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 2.75), align = "h")
cow2
ggsave(plot = cow2, 
       filename = paste0("data/Community_Composition/Beta_Diversity_Analysis/",
         z[cnt], "_",
         cohort, "_",
         dist_label,
         ".svg"
       ),
       width = 14, height = 9)


#-------------------------------------------------------------------------------
#####                      Weighted Unifrac                            ##### 
#-------------------------------------------------------------------------------


cohort = "Merged"
dist = "wunifrac" # unifrac # wunifrac
z <- c("Species")
i <- dat.species.tree
cnt <- 1
load_betadiv_colors()

obj_dist <- dat.species.tree %>%
  subset_samples(donor_id %ni% low_qc[[1]]) %>% 
  microbiome::transform("compositional")

iDist <- phyloseq::distance(obj_dist, method=dist)
dist_label <- "Weighted UniFrac"

# Calculate ordination
iMDS  <- phyloseq::ordinate(obj_dist, "MDS", distance=iDist)

# ---------------------------------------------------
#  PCoA for Axis 1 and 2
# ---------------------------------------------------

plot_ordination(obj_dist, iMDS, color="cohort", axes = c(1, 2))
p <- plot_ordination(obj_dist, iMDS, color="description", axes = c(1, 2))

df12 = p$data
df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
p <- ggplot(df12, aes(Axis.1, Axis.2, fill = description, color=description))
p <- p + geom_point(shape=21, size=5, alpha=0.7)
ord <- p + 
  theme_bw() + 
  labs(fill="Donor Group") +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Donor Group") +
  scale_fill_manual(values = cols.pdpchc.dark) +
  scale_color_manual(values = cols.pdpchc.rim) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())

# PCoA data to fit ridgeline Plots
my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2

# ---------------------------------------------------
#  Boxplot - Axis 1
# ---------------------------------------------------

r1 <- ggplot(df12, aes(x = Axis.1, y = description)) +
  geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 1), "%)")) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.xrange) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  ggtitle(paste0(dist_label, " Distance PCoA")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")

# ---------------------------------------------------
#  Boxplot - Axis 2
# ---------------------------------------------------

r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
  geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.yrange2) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_blank(),
        legend.position = "none") +
  coord_flip()

# ---------------------------------------------------
#  Beta Diversity Violin Plots
# ---------------------------------------------------
p <- obj_dist
s <- "donor_id"
d <- "description"

# Make Melted Distance Matrix 
wu.m <- melt(as.matrix(iDist))
# Exclude intra-sample distances 
wu.m <- wu.m %>% filter(as.character(Var1) != as.character(Var2)) %>% 
  mutate_if(is.factor, as.character)
# Pull metadata of interest
sd <- meta(p) %>% dplyr::select(all_of(s), all_of(d)) %>% mutate_if(is.factor, as.character)
# Add group name for Var1 Column
colnames(sd) <- c("Var1", "Type1")
wu.sd <- left_join(wu.m, sd, by = "Var1")
# Add group name for Var2 Column
colnames(sd) <- c("Var2", "Type2")
wu.sd <- left_join(wu.sd, sd, by = "Var2")
# Select only distances to Population control
wu.sd <- filter(wu.sd, Type1 == "Population Control")

# Specifying comparisons for analysis
my_comparisons <-
  list(c("Household Control", "PD Patient"),
       c("Population Control", "PD Patient"))
wu.sd$Type2 <-
  factor(wu.sd$Type2,
         levels = c("Population Control", "PD Patient", "Household Control"))

v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
  geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size= 0.2, alpha = 0.5, cex = 0.3) +
  geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha=0) +
  theme_classic() +
  ylab(paste(dist_label, "Distance")) +
  labs(fill="Group") +
  scale_color_manual(values = cols.pdpchc) +
  scale_fill_manual(values = cols.pdpchc) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0) + 
  ggtitle(paste0(dist_label," Distance to Population Control\n", z[cnt], " abundance")) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none")


# ---------------------------------------------------
#  Assemble Plots
# ---------------------------------------------------
cat(paste0("Assembling summary plots for: " , z[cnt], "\n"))
ord.plot <- ord + theme(legend.position = "none")
cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 5), rel_widths = c(5, 1),  align = "vh")
cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 2.75), align = "h")
cow2
ggsave(plot = cow2, 
       filename = paste0("data/Community_Composition/Beta_Diversity_Analysis/",
                         z[cnt], "_",
                         cohort, "_",
                         dist_label,
                         ".svg"
       ),
       width = 14, height = 9)












# dist_methods <- unlist(distanceMethodList)
# print(dist_methods)
# 
# # These require tree
# dist_methods[(1:3)]
# 
# # Remove them from the vector
# dist_methods <- dist_methods[-(1:3)]
# # This is the user-defined method:
# dist_methods["designdist"]
# 
# # Remove the user-defined distance
# dist_methods = dist_methods[-which(dist_methods=="ANY")]
# 
# plist <- vector("list", length(dist_methods))
# names(plist) = dist_methods
# for( i in dist_methods ){
#   # Calculate distance matrix
#   iDist <- dat.species %>% 
#     subset_samples(donor_id %ni% low_qc[[1]]) %>% 
#     phyloseq::distance(method=i)
#   # Calculate ordination
#   iMDS  <- dat.species %>% 
#     subset_samples(donor_id %ni% low_qc[[1]]) %>% 
#     ordinate("MDS", distance=iDist)
#   ## Make plot
#   # Don't carry over previous plot (if error, p will be blank)
#   p <- NULL
#   # Create plot, store as temp variable, p
#   p <- dat.species %>% 
#     subset_samples(donor_id %ni% low_qc[[1]]) %>% 
#     plot_ordination(iMDS, color="donor_group", shape="donor_group")
#   # Add title to each plot
#   p <- p + ggtitle(paste("MDS using distance method ", i, sep=""))
#   # Save the graphic to file.
#   plist[[i]] = p
# }
# 
# 
# df = ldply(plist, function(x) x$data)
# names(df)[1] <- "distance"
# p = ggplot(df, aes(Axis.1, Axis.2, color=donor_group, shape=donor_group))
# p = p + geom_point(size=3, alpha=0.5)
# p = p + facet_wrap(~distance, scales="free") + theme_minimal()
# p = p + ggtitle("MDS on various distance metrics for Enterotype dataset")
# p
# 


