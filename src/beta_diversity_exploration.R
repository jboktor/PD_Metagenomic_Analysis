# beta diversity exploration

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/community_composition_funcs.R")
load("files/low_quality_samples.RData")

remove_dats()
load_data()

theme_set(theme_bw())
dist_methods <- unlist(distanceMethodList)
print(dist_methods)
# These require tree
dist_methods[(1:3)]
# Remove them from the vector
dist_methods <- dist_methods[-(1:3)]

plist <- vector("list", length(dist_methods))
names(plist) <- dist_methods
dat.obj <- dat.species.shanghai %>%
  subset_samples(donor_id %ni% low_qc[[1]])

for (i in dist_methods) {
  # Calculate distance matrix
  iDist <- phyloseq::distance(dat.obj, method = i)
  # Calculate ordination
  iMDS <- ordinate(dat.obj, "MDS", distance = iDist)
  ## Make plot
  # Don't carry over previous plot (if error, p will be blank)
  p <- NULL
  # Create plot, store as temp variable, p
  p <- plot_ordination(dat.obj, iMDS, color = "donor_group")
  # Add title to each plot
  p <- p + ggtitle(paste("MDS using distance method ", i, sep = ""))
  # Save the graphic to file.
  plist[[i]] <- p
}
iMDS
df.dists <- ldply(plist, function(x) x$data)
names(df.dists)[1] <- "distance"
p <- ggplot(df, aes(Axis.1, Axis.2, color = donor_group))
p <- p + geom_point(size = 3, alpha = 0.5)
p <- p + facet_wrap(~distance, scales = "free")
p <- p + ggtitle("MDS on various distance metrics")
p










iDist <- phyloseq::distance(dat.obj, method = "rlb")
iMDS <- ordinate(dat.obj, "MDS", distance = iDist)
p <- plot_ordination(dat.obj, iMDS, color = "donor_group", axes = c(1, 2))

df12 <- p$data
# df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
p <- ggplot(df12, aes(Axis.1, Axis.2, fill = donor_group, color = donor_group))
p <- p + geom_point(shape = 21, size = 3, alpha = 0.7)
ord <- p +
  theme_bw() +
  labs(fill = "Donor Group") +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2]) * 100, digits = 2), "%)")) +
  labs(fill = "Donor Group") +
  scale_fill_manual(values = cols.pdpchc.dark) +
  scale_color_manual(values = cols.pdpchc.rim) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


# PCoA data to fit ridgeline Plots
my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2

# ---------------------------------------------------
#  Boxplot - Axis 1
# ---------------------------------------------------

r1 <- ggplot(df12, aes(x = Axis.1, y = donor_group)) +
  geom_boxplot(aes(color = donor_group, fill = donor_group), alpha = 0.2) +
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1]) * 100, digits = 1), "%)")) +
  scale_y_discrete(expand = c(0, 0)) + # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.xrange) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  ggtitle(paste0("RLB Distance PCoA")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

# ---------------------------------------------------
#  Boxplot - Axis 2
# ---------------------------------------------------

r2 <- ggplot(df12, aes(x = Axis.2, y = donor_group)) +
  geom_boxplot(aes(color = donor_group, fill = donor_group), alpha = 0.2) +
  scale_y_discrete(expand = c(0, 0)) + # will generally have to set the `expand` option
  scale_x_continuous(limits = my.ggp.yrange2) +
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = cols.pdpchc.dark) +
  scale_fill_manual(values = cols.pdpchc) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
  ) +
  coord_flip()

# ---------------------------------------------------
#  Beta Diversity Violin Plots
# ---------------------------------------------------
p <- dat.obj
s <- "donor_id"
d <- "donor_group"

# Make Melted Distance Matrix
wu.m <- melt(as.matrix(iDist))
# Exclude intra-sample distances
wu.m <- wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
# Pull metadata of interest
sd <- meta(p) %>%
  dplyr::select(all_of(s), all_of(d)) %>%
  mutate_if(is.factor, as.character)
# Add group name for Var1 Column
colnames(sd) <- c("Var1", "Type1")
wu.sd <- left_join(wu.m, sd, by = "Var1")
# Add group name for Var2 Column
colnames(sd) <- c("Var2", "Type2")
wu.sd <- left_join(wu.sd, sd, by = "Var2")
# Select only distances to Population control
wu.sd <- filter(wu.sd, Type1 == "HC")

# Specifying comparisons for analysis
my_comparisons <-
  list(c("HC", "PD"))
# wu.sd$Type2 <-
#   factor(wu.sd$Type2,
#          levels = c("Population Control", "PD Patient", "Household Control"))


# if (cohort == "Merged"){
#   beesize <- 0.2
#   beescale <- 0.2
# } else {
beesize <- 0.3
beescale <- 0.5
# }

v <- ggplot(wu.sd, aes(x = Type2, y = value)) +
  theme_minimal() +
  geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size = beesize, alpha = 0.5, cex = beescale) +
  geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha = 0) +
  theme_classic() +
  ylab(paste("RLB Distance")) +
  labs(fill = "Group") +
  scale_color_manual(values = cols.pdpchc) +
  scale_fill_manual(values = cols.pdpchc) +
  scale_x_discrete(labels = c(
    "Population Control" = "Population \nControl",
    "PD" = "Parkinson's \nDisease",
    "HC" = "Healthy \nControl"
  )) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0) +
  # ggtitle(paste0(" Distance to Population Control\n", z[cnt], " Abundance")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )


# ---------------------------------------------------
#  Assemble Plots
# ---------------------------------------------------
ord.plot <- ord + theme(legend.position = "none")
cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2,
  nrow = 2, ncol = 2,
  rel_heights = c(1, 3.25), rel_widths = c(4.5, 1),
  align = "vh", axis = "tblr"
)
cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 3.5), align = "h")
# print(cow2)

ggsave(
  plot = cow2,
  filename = paste0(
    "data/MISC/Beta_Diversity_Analysis_RLB_distance_Shanghai.svg"
  ),
  width = 7, height = 5.5
)
