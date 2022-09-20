# Differential Abundance of Features Script

library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)
library(grid)
library(devtools)
library(RColorBrewer)
library(ggfortify)
library(vegan)
library(MASS)
library(compositions)
library(zCompositions)
library(phyloseq)
library(Biobase)
library(viridis)
library("foreach")
library("doParallel")
library(ggbeeswarm)
library(FSA)
library(ggpubr)
library(ggsci)
library(microbiome)
library(ggridges)
library(future)
library(cowplot)



### NOTES/UPDATES REQUIRED:
# Merge this script with _functional and pass functions through loop for all desired analyses



rm(list = ls())

######## Load Data & functions
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")

######### SWAP CLADE LEVEL HERE:
LEV <- Phylo_Objects$Phylum
lev <- "Phylum"

### USE if log transformation is to be used for viz
# # find 1/2 min abundance
# abun <- as.data.frame.matrix(abundances(LEV))
# impute <- min(abun[abun > 0])/2
# otu_table(LEV) <- otu_table(LEV) + impute
# #Sanity Check min abundance
# min(abundances(dat_obj))

### Visualization Transformations
taxa_names(LEV) <- gsub("s__", "", taxa_names(LEV))
dat_obj <- microbiome::transform(LEV, "compositional") # %>% transform("log10")
## ArcSinSqrt() Transformation
otu_table(dat_obj) <- asin(sqrt(otu_table(dat_obj)))


############# data prep #############
# PD v PC
dat_pdpc <- subset_samples(dat_obj, donor_group != "HC")
abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
# PD v HC PAIRED
dat_pdhc <- subset_samples(dat_obj, Paired != "No")
abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))


### Read-in Maaslin Files - all features used in significance testing
Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", lev, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>%
  filter(value == "Population Control")
Maas.pd.pc$feature <- gsub("s__", "", Maas.pd.pc$feature)
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/", lev, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>%
  filter(value == "Household Control")
Maas.pd.hc$feature <- gsub("s__", "", Maas.pd.hc$feature)
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25)


### select significant features from abundance tables
# PD v PC
abun.pdpc <- rownames_to_column(abun.pdpc)
abun.pdpc.inpt <- filter(abun.pdpc, rowname %in% Maas.pd.pc.sig$feature) %>%
  column_to_rownames(var = "rowname") %>%
  t() %>%
  melt() %>%
  mutate(group = if_else(grepl(".PC", Var1), "PC", "PD"))
# PD v HC PAIRED
abun.pdhc <- rownames_to_column(abun.pdhc)
abun.pdhc.inpt <- filter(abun.pdhc, rowname %in% Maas.pd.hc.sig$feature) %>%
  column_to_rownames(var = "rowname") %>%
  t() %>%
  melt() %>%
  mutate(group = if_else(grepl(".HC", Var1), "HC", "PD"))

## Color Schemes
cols.pdpc <- c("PD" = "#6495ED", "PC" = "#8D230F")
cols.pdhc <- c("PD" = "#6495ED", "HC" = "#FFD64D")




########################## PD v PC Plots ##########################

# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PD <- subset_samples(dat_pdpc, donor_group == "PD")
dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD)) %>%
  rownames_to_column() %>%
  filter(rowname %in% Maas.pd.pc.sig$feature) %>%
  column_to_rownames()

dat_pdpc.PC <- subset_samples(dat_pdpc, donor_group == "PC")
dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC)) %>%
  rownames_to_column() %>%
  filter(rowname %in% Maas.pd.pc.sig$feature) %>%
  column_to_rownames()

######  Generalized or pseudo-fold change
#' Calculated as geometric mean of the differences between the
#' uantiles for the different classes found in the label

# Initalize counts
i <- 1
probs.fc <- seq(.1, .9, .05)
gfc_data <- c()

for (feature in 1:nrow(dat_pdpc.PCabun)) {
  # Loops through each species row and calculates
  # the Generalized fold change for each feature
  cat("Feature number: ", feature, "\n")
  cat("Testing PD: ", rownames(dat_pdpc.PDabun)[feature], "feature vs ", rownames(dat_pdpc.PCabun)[feature], "\n")
  q.PD <- quantile(dat_pdpc.PDabun[feature, ], probs = probs.fc)
  q.ctl <- quantile(dat_pdpc.PCabun[feature, ], probs = probs.fc)

  gfc <- sum(q.PD - q.ctl) / length(q.ctl)
  print(gfc)

  gfc_data <- c(gfc_data, gfc)
}
PDovrPC <- tibble("feature" = rownames(dat_pdpc.PCabun), "gFC" = gfc_data)


## To get order for Significant feauture plots
PDovrPC.sig <- PDovrPC %>% filter(feature %in% Maas.pd.pc.sig$feature)
PDovrPC.order <- PDovrPC.sig[order(-PDovrPC.sig$gFC), ]



###### Generalized Fold Change (gFC) BarPlot ######
PDovrPC.BP <- PDovrPC.order
PDovrPC.BP$direction <- PDovrPC.BP$gFC > 0
PDovrPC.BP$feature <- factor(PDovrPC.BP$feature, levels = PDovrPC.order$feature)

g0 <- ggplot() +
  geom_col(
    data = PDovrPC.BP, aes(x = gFC, y = feature, fill = direction),
    position = "dodge", width = 0.6, alpha = 0.5
  ) +
  theme_minimal() +
  # labs(x="Generalized Fold Change") +
  ggtitle("Fold Change") +
  scale_fill_manual(values = c("TRUE" = "#6495ED", "FALSE" = "#8D230F")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
g0



###### BoxPlots ######
## PD v PC
set.seed(123)
abun.pdpc.inpt$Var2 <- factor(abun.pdpc.inpt$Var2, levels = PDovrPC.order$feature)
g1 <- ggplot(data = abun.pdpc.inpt, aes(x = value, y = Var2)) +
  geom_boxplot(aes(fill = group), alpha = 0.5, outlier.alpha = 0, width = 0.8) +
  geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = .2), shape = 21, size = .7, alpha = 0.6) +
  theme_minimal() +
  ggtitle(paste0("Differential Abundance: ", lev)) +
  # labs(x = "Abundance") +
  scale_fill_manual(values = cols.pdpc, name = "Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
g1

# ggsave(g1, filename = "PDvPC_DAF_Maaslin_BoxPlots.png",
#        width = 6, height = length(unique(abun.pdpc.inpt$Var2))/2)


###### Significance Plot ######
sigplot.df.pdpc <- dplyr::select(Maas.pd.pc.sig, c("feature", "pval", "qval")) %>% melt()

sigplot.df.pdpc$feature <- factor(sigplot.df.pdpc$feature, levels = PDovrPC.order$feature)
g2 <- ggplot() +
  geom_col(data = sigplot.df.pdpc, aes(x = -log10(value), y = feature, fill = variable), position = "dodge", width = 0.8) +
  theme_minimal() +
  # labs(x=expression('-log'[10]*'(value)')) +
  ggtitle("Significance") +
  scale_fill_manual(values = c("pval" = "#d3d3d3", "qval" = "#676767")) +
  geom_vline(xintercept = -log10(0.1), linetype = "dotted", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
g2

# ggsave(g2, filename = "PDvPC_DAF_Maaslin_Signifiance.png",
#        width = 3, height = length(unique(abun.pdpc.inpt$Var2))/2)


# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PD <- subset_samples(dat_pdpc, donor_group == "PD")
dat_pdpc.PDprev <- tibble::enframe(prevalence(dat_pdpc.PD)) %>% filter(name %in% Maas.pd.pc.sig$feature)
colnames(dat_pdpc.PDprev) <- c("feature", "PD")

dat_pdpc.PC <- subset_samples(dat_pdpc, donor_group == "PC")
dat_pdpc.PCprev <- tibble::enframe(prevalence(dat_pdpc.PC)) %>% filter(name %in% Maas.pd.pc.sig$feature)
colnames(dat_pdpc.PCprev) <- c("feature", "PC")

# Merge DFs
dat_pdpc.PREV <- left_join(dat_pdpc.PDprev, dat_pdpc.PCprev, by = "feature") %>% melt()

###### Prevalence Plot ######
dat_pdpc.PREV$feature <- factor(dat_pdpc.PREV$feature, levels = PDovrPC.order$feature)
dat_pdpc.PREV$variable <- factor(dat_pdpc.PREV$variable, levels = c("PC", "PD"))
g3 <- ggplot() +
  geom_col(
    data = dat_pdpc.PREV, aes(x = value * 100, y = feature, fill = variable),
    position = "dodge", alpha = 0.5, width = 0.8
  ) +
  theme_minimal() +
  labs(x = expression("Prevalence [%]")) +
  # ggtitle("Prevalence") +
  scale_fill_manual(values = cols.pdpc, name = "Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
g3

###### Merge Panels ######
g0a <- g0 + theme(axis.text.y = element_blank(), legend.position = c(.80, .90))
g1a <- g1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
g2a <- g2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .35))
g3a <- g3 + theme(axis.text.y = element_blank(), legend.position = "none")
# cowplot::plot_grid(g1, g2, g0, g3, labels = "A", align = "h", nrow = 1, rel_widths = c(2, 1, 1)) # View to make sure labels are aligned
# fig2a <- cowplot::plot_grid(g1a, g2a, g3a, g0a, labels = "A", align = "h", nrow = 1, rel_widths = c(3, 1, 1, 1))
# fig2a

# ggsave(fig2a, filename = "SignificantFeatures_PDvPC_MaaslinSig.SVG",
# width = 12, height = length(unique(abun.pdpc.inpt$Var2))/2)







########################## PD v HC Plots ##########################

# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PD <- subset_samples(dat_pdhc, donor_group == "PD")
dat_pdhc.PDabun <- as.data.frame.matrix(abundances(dat_pdhc.PD)) %>%
  rownames_to_column() %>%
  filter(rowname %in% Maas.pd.hc.sig$feature) %>%
  column_to_rownames()

dat_pdhc.HC <- subset_samples(dat_pdhc, donor_group == "HC")
dat_pdhc.HCabun <- as.data.frame.matrix(abundances(dat_pdhc.HC)) %>%
  rownames_to_column() %>%
  filter(rowname %in% Maas.pd.hc.sig$feature) %>%
  column_to_rownames()


######  Generalized or pseudo-fold change
#' Calculated as geometric mean of the differences between the
#' uantiles for the different classes found in the label

# Initalize counts
i <- 1
probs.fc <- seq(.1, .9, .05)
gfc_data <- c()

for (feature in 1:nrow(dat_pdhc.PDabun)) {
  # Loops through each species row and calculates
  # the Generalized fold change for each feature
  cat("Feature number: ", feature, "\n")
  cat("Testing PD: ", rownames(dat_pdhc.PDabun)[feature], "feature vs HC", rownames(dat_pdhc.HCabun)[feature], "\n")
  q.PD <- quantile(dat_pdhc.PDabun[feature, ], probs = probs.fc)
  q.ctl <- quantile(dat_pdhc.HCabun[feature, ], probs = probs.fc)

  gfc <- sum(q.PD - q.ctl) / length(q.ctl)
  print(gfc)
  gfc_data <- c(gfc_data, gfc)
}
PDovrHC <- tibble("feature" = rownames(dat_pdhc.PDabun), "gFC" = gfc_data)

## To get order for Significant feauture plots
PDovrHC.sig <- PDovrHC %>% filter(feature %in% Maas.pd.hc.sig$feature)
PDovrHC.order <- PDovrHC.sig[order(-PDovrHC.sig$gFC), ]


###### Generalized Fold Change (gFC) BarPlot ######
PDovrHC.BP <- PDovrHC.order
PDovrHC.BP$direction <- PDovrHC.BP$gFC > 0
PDovrHC.BP$feature <- factor(PDovrHC.BP$feature, levels = PDovrHC.order$feature)

h0 <- ggplot() +
  geom_col(
    data = PDovrHC.BP, aes(x = gFC, y = feature, fill = direction),
    position = "dodge", width = 0.6, alpha = 0.5
  ) +
  theme_minimal() +
  labs(x = "Generalized Fold Change") +
  # ggtitle("Fold Change") +
  scale_fill_manual(values = c("TRUE" = "#6495ED", "FALSE" = "#FFD64D")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
h0



###### BoxPlots ######
abun.pdhc.inpt$Var2 <- factor(abun.pdhc.inpt$Var2, levels = PDovrHC.order$feature)
set.seed(123)
h1 <- ggplot(data = abun.pdhc.inpt, aes(x = value, y = Var2)) +
  geom_boxplot(aes(fill = group), alpha = 0.5, outlier.alpha = 0, width = 0.8) +
  geom_point(aes(fill = group), position = position_jitterdodge(jitter.width = .2), shape = 21, size = .7, alpha = 0.6) +
  theme_minimal() +
  # ggtitle("Differentially Abundant Species") +
  labs(x = "Abundance") +
  scale_fill_manual(values = cols.pdhc, name = "Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
h1

# ggsave(h1, filename = "BoxPlot_PDvHC_MaaslinSig.png",
#        width = 6, height = length(unique(abun.pdhc.inpt$Var2))/2)


###### Significance Plot ######
sigplot.df.pdhc <- dplyr::select(Maas.pd.hc.sig, c("feature", "pval", "qval")) %>% melt()

sigplot.df.pdhc$feature <- factor(sigplot.df.pdhc$feature, levels = PDovrHC.order$feature)
h2 <- ggplot() +
  geom_col(data = sigplot.df.pdhc, aes(x = -log10(value), y = feature, fill = variable), position = "dodge", width = 0.8) +
  theme_minimal() +
  labs(x = expression("-log"[10] * "(value)")) +
  # ggtitle("Significance") +
  scale_fill_manual(values = c("pval" = "#d3d3d3", "qval" = "#676767")) +
  geom_vline(xintercept = -log10(0.1), linetype = "dotted", color = "red") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
h2

# ggsave(h2, filename = "PDvHC_DAF_Maaslin_Signifiance.png",
#        width = 3, height = length(unique(abun.pdpc.inpt$Var2))/2)


# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PDprev <- tibble::enframe(prevalence(dat_pdhc.PD)) %>% filter(name %in% Maas.pd.hc.sig$feature)
colnames(dat_pdhc.PDprev) <- c("feature", "PD")

dat_pdhc.HCprev <- tibble::enframe(prevalence(dat_pdhc.HC)) %>% filter(name %in% Maas.pd.hc.sig$feature)
colnames(dat_pdhc.HCprev) <- c("feature", "HC")

# Merge DFs
dat_pdhc.PREV <- left_join(dat_pdhc.PDprev, dat_pdhc.HCprev, by = "feature") %>% melt()

###### Prevalence Plot ######
dat_pdhc.PREV$feature <- factor(dat_pdhc.PREV$feature, levels = PDovrHC.order$feature)
dat_pdhc.PREV$variable <- factor(dat_pdhc.PREV$variable, levels = c("HC", "PD"))
h3 <- ggplot() +
  geom_col(data = dat_pdhc.PREV, aes(x = value * 100, y = feature, fill = variable), position = "dodge", alpha = 0.5, width = 0.8) +
  theme_minimal() +
  labs(x = expression("Prevalence [%]")) +
  # ggtitle("Prevalence") +
  scale_fill_manual(values = cols.pdhc, name = "Group") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank()
  )
h3


###### Merge Panels ######
h0a <- h0 + theme(axis.text.y = element_blank(), legend.position = c(.80, .85))
h1a <- h1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
h2a <- h2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .85))
h3a <- h3 + theme(axis.text.y = element_blank(), legend.position = "none")
# cowplot::plot_grid(h1, h2, h3, h0, labels = "A", align = "h", nrow = 1, rel_widths = c(3, 1, 1, 1)) # View to make sure labels are aligned
# fig2b <- cowplot::plot_grid(h1a, h2a, h3a, h0a, labels = "B", align = "h", nrow = 1, rel_widths = c(3, 1, 1, 1))
# fig2b

# ggsave(fig2b, filename = "SignificantFeatures_PDvHC_MaaslinSig.SVG",
#        width = 12, height = length(unique(abun.pdhc.inpt$Var2))/2)
#


############ Final View: ############

g0a <- g0a + theme(axis.title.x = element_blank())
g1a <- g1a + theme(axis.title.x = element_blank())
g2a <- g2a + theme(axis.title.x = element_blank())
g3a <- g3a + theme(axis.title.x = element_blank())

# Setting plot length variables - Super extra but it works?
top_len <- length(unique(PDovrPC.BP$feature)) + 1
bottom_len <- length(unique(PDovrHC.BP$feature)) + 1


DAF_part1 <- cowplot::plot_grid(g1a, h1a,
  nrow = 2, align = "hv",
  rel_heights = c(top_len, bottom_len)
)
DAF_part1

DAF_part2 <- cowplot::plot_grid(g2a, g3a, g0a, h2a, h3a, h0a,
  nrow = 2, ncol = 3, align = "h",
  rel_heights = c(top_len, bottom_len)
)
DAF_part2

DAF_final <- cowplot::plot_grid(DAF_part1, DAF_part2, ncol = 2, rel_widths = c(1.5, 1))
# DAF_final

ggsave(DAF_final,
  filename = paste0("data/DAF_Analysis/DAF_", lev, "_MaaslinSig.png"),
  width = 20, height = 10
)






############################  CODE JUNKYARD ############################

## Another option for subsetting instead of making all the extra dat files
# abun[which(!grepl("PC", colnames(abun)) & !grepl("HC", colnames(abun)) )]



###### Volcano Plot
# x_limits <- c(1.5, NA)
# v1 <- ggplot(data=PDovrPC.volc, aes(x=-log10(`pval`), y=`Log2FC(PD/PC)`, color=SigThresh)) + geom_point(alpha = 0.7) +
#   labs(y = expression('log'[2]*'(PD/PC)'), x= expression('log'[10]*'(p-value)')) +
#   scale_colour_manual(values = c("A"= "black", "B"= "#8D230F", "C"= "black")) +
#   ylim(-4, 7) +
#   geom_vline(xintercept = 1.3, linetype = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0.25, linetype = 2, alpha = 0.5) +
#   geom_hline(yintercept = -0.25, linetype = 2, alpha = 0.5) +
#   theme_minimal() +
#   # geom_text_repel(data= subset(PDovrPC.volc, SigThresh != "C"), force = 10,
#   #                 aes(label=feature), size = 3, xlim = x_limits,
#   #                 segment.size = 0.2, segment.alpha = 0.7)+
#   geom_text_repel(data= subset(PDovrPC.volc, SigThresh == "A"), force = 10,
#                   aes(label=feature), size = 3, xlim = x_limits, ylim = c(1, 7),
#                   nudge_x = -0.5, nudge_y = 2,
#                   segment.size = 0.2, segment.alpha = 0.7)+
#   geom_text_repel(data= subset(PDovrPC.volc, SigThresh == "B"), force = 10,
#                   aes(label=feature), size = 3, xlim = x_limits, ylim = c(-1, -4),
#                   nudge_x = -0.5, nudge_y = -2, direction= "y",
#                   segment.size = 0.2, segment.alpha = 0.7)+
#   theme(legend.position="none",
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_blank())
# v1

###### Volcano Plot
# x_limits <- c(NA, 1.5)
# v2 <- ggplot(data=PDovrHC.volc, aes(x=-log10(`pval`), y=`Log2FC(PD/HC)`, color=SigThresh)) + geom_point(alpha = 0.6) +
#   scale_x_reverse() +
#   labs(y = expression('log'[2]*'(PD/HC)'), x= expression('log'[10]*'(p-value)')) +
#   ylim(-4, 7) +
#   scale_colour_manual(values = c("A"= "black", "B"= "#FFD64D", "C"= "black")) +
#   geom_vline(xintercept = 1.3, linetype = 2, alpha = 0.5) +
#   geom_hline(yintercept = 0.25, linetype = 2, alpha = 0.5) +
#   geom_hline(yintercept = -0.25, linetype = 2, alpha = 0.5) +
#   theme_minimal() +
#   geom_text_repel(data= subset(PDovrHC.volc, SigThresh == "A"), force = 10,
#                   aes(label=feature), size = 3, xlim = x_limits, nudge_x = -0.5,
#                   nudge_y = 2, direction= "y",
#                   segment.size = 0.2, segment.alpha = 0.7)+
#   geom_text_repel(data= subset(PDovrHC.volc, SigThresh == "B"), force = 10,
#                   aes(label=feature), size = 3, xlim = x_limits, nudge_x = -0.5,
#                   nudge_y = -2, direction= "y",
#                   segment.size = 0.2, segment.alpha = 0.7)+
#   theme(legend.position="none",
#         panel.grid.minor.y = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_blank())
#
# # v2
#
#### MERGE VOLCANO PLOTS
# v2a <- v2 + theme(axis.text.y = element_blank(), axis.title.y = element_blank())
# fig2d <- cowplot::plot_grid(v1, v2a, labels = "D", align = c("h", "v"), nrow = 1)

# ggsave(fig2d, filename = "SignificantFeatures_VOLCANO.png",
#        width = 20, height = 5)




########## OLD Log2FC sequence
#
# # Subset of phlyoseq obj subset to get samples of interest
# dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
# dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD))
# pdpc.PD.mean <- tibble::enframe(apply(dat_pdpc.PDabun, 1, mean))
#
# dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
# dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC))
# pdpc.PC.mean <- tibble::enframe(apply(dat_pdpc.PCabun, 1, mean))
#
# #  Add .001 to all means - prevents 1/0 situation - may be a more standard method to do this (1/2 min value?)
# PDovrPC <- tibble("feature" = pdpc.PC.mean$name, "Log2FC(PD/PC)" = log2((pdpc.PD.mean$value + .001)/(pdpc.PC.mean$value + .001)) )
# PDovrPC.tested <- PDovrPC %>% filter(feature %in% Maas.pd.pc$feature)
# PDovrPC.volc <- left_join(Maas.pd.pc, PDovrPC.tested, by = "feature")
# PDovrPC.volc <- mutate(PDovrPC.volc, SigThresh = if_else(pval < 0.05 & `Log2FC(PD/PC)` > 0, "A",
#                                                          if_else( pval < 0.05 & `Log2FC(PD/PC)` < 0, "B", "C" )))
#
# ## To get order for Significant feauture plots
# PDovrPC.sig <- PDovrPC %>% filter(feature %in% Maas.pd.pc.sig$feature)
# PDovrPC.order <- PDovrPC.sig[order(-PDovrPC.sig$`Log2FC(PD/PC)`),]
#
#
# ###### LogFC BarPlot
# PDovrPC.BP <- PDovrPC.order
# PDovrPC.BP$direction <- PDovrPC.BP$`Log2FC(PD/PC)` > 0
# PDovrPC.BP$feature <- factor(PDovrPC.BP$feature, levels = PDovrPC.order$feature)
#
# g0 <- ggplot() +
#   geom_col(data=PDovrPC.BP, aes(x=`Log2FC(PD/PC)`, y= feature,  fill = direction),
#            position = "dodge", width = 0.6, alpha = 0.5) +
#   theme_minimal() +
#   labs(x=expression('log'[2]*'(PD/PC)')) +
#   ggtitle("Fold Change") +
#   scale_fill_manual(values = c("TRUE" = "#6495ED", "FALSE" = "#8D230F")) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.y = element_blank(),
#         panel.grid.major.y = element_blank())
# g0
