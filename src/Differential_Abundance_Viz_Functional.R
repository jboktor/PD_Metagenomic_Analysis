# Differential Abundance of Features Script

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot)

rm(list = ls())

######## Load Data & functions
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")
load("files/GBMs_PhyloseqObj.RData")
load("files/GMMs_PhyloseqObj.RData")

################################################################################# 
########################### SWAP FUNCTION LEVEL HERE: ########################### 
################################################################################# 

# LEV <- Phylo_Objects$Enzymes.slim
# lev <- "Enzymes.slim"
LEV <- dat.GMMs
lev <- "GMMs"

################################################################################# 

# ## Color Schemes
cols.pdpc <- c("PD"= "#ed7d31", "PC" = "#bfbfbf")
cols.pdhc <- c("PD"= "#ed7d31", "HC" = "#5b9bd5")
# Rims
cols.pdpc.rim <- c("PD"= "#ed4e31", "PC" = "#999999")
cols.pdhc.rim <- c("PD"= "#ed4e31", "HC" = "#5b7dd5")


# To protect against row/colname errors 
if (lev == "Pathways"| lev == "Pathways.slim") {
  features <- paste0("PATHWAY_", taxa_names(LEV)) 
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  taxa_names(LEV) <- features
} else if (lev == "Enzymes" | lev == "Enzymes.slim") {
  features <- paste0("ENZYME_", taxa_names(LEV)) 
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  taxa_names(LEV) <- features
} else if (lev == "Species") {
  taxa_names(LEV) <- gsub("s__", "", taxa_names(LEV))
} else if (lev == "GMMs" | lev == "GBMs" ) {
  taxa_names(LEV) <- gsub(" ", ".", taxa_names(LEV))
  } else {
  features <- taxa_names(LEV)
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  taxa_names(LEV) <- features
  }

# Check process above
# taxa_names(LEV)

################################################################################# 

### Visualization Transformations
# taxa_names(LEV) <- gsub("s__", "", taxa_names(LEV))
dat_obj <- microbiome::transform(LEV, "compositional")
## ArcSinSqrt() Transformation
otu_table(dat_obj) <- asin(sqrt(otu_table(dat_obj)))


# # Distribution sanity check
# histogram(melt(abundances(dat_obj))$value)
# # ECDF plot slows analysis down 
# abund.melt <- melt(abundances(dat_obj))
# abund.melt <- mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC", if_else(grepl("PC", Var2), "PC", "PD")))
# ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
#   labs(y = "ECDF")



############# data prep ############# 
# PD v PC
dat_pdpc = subset_samples(dat_obj, donor_group !="HC")
abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
# PD v HC PAIRED
dat_pdhc = subset_samples(dat_obj, Paired !="No")
abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))


### Read-in Maaslin Files - all features used in significance testing
Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", lev, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/", lev, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) 

if (lev == "Species") {
  Maas.pd.pc.sig$feature <- gsub("s__", "", Maas.pd.pc.sig$feature)
  Maas.pd.hc.sig$feature <- gsub("s__", "", Maas.pd.hc.sig$feature)
}

########################## Venn Diagram Plots ##########################

v <- VennPlot(Maas.pd.pc.sig, Maas.pd.hc.sig, qval_threshold = 0.1)

# Save Venn Diagrams
if (!is.null(v$venn_depleted)){
  pdf(file = paste0("data/DAF_Analysis/DAF_", lev, "_VennDiagram_PD_depleted.pdf"),
      width = 7, 
      height = 5,
      pointsize = 12)
      # units = "in", pointsize = 12, res=300)
  plot(v$venn_depleted)
  dev.off()
}

if (!is.null(v$venn_enriched)){
  pdf(file = paste0("data/DAF_Analysis/DAF_", lev, "_VennDiagram_PD_enriched.pdf"),
      width = 7, 
      height = 5,
      pointsize = 12)
      # units = "in", pointsize = 12, res=300)
  plot(v$venn_enriched)
  dev.off()
}

#############################################################################

## For plotting select top 50 features by qval rank

if (nrow(Maas.pd.pc.sig) > 25){
  Maas.pd.pc.sig <- Maas.pd.pc.sig[1:25,]
}
if (nrow(Maas.pd.hc.sig) > 25){
  Maas.pd.hc.sig <- Maas.pd.hc.sig[1:25,]
}

### select significant features from abundance tables 
# PD v PC
abun.pdpc <- rownames_to_column(abun.pdpc)
abun.pdpc.inpt <- filter(abun.pdpc, rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames(var="rowname") %>% 
  t() %>% melt() %>% mutate(group = if_else(grepl(".PC", Var1), "PC", "PD"))
# PD v HC PAIRED
abun.pdhc <- rownames_to_column(abun.pdhc)
abun.pdhc.inpt <- filter(abun.pdhc, rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames(var="rowname")  %>% 
  t() %>% melt() %>% mutate(group = if_else(grepl(".HC", Var1), "HC", "PD"))


#############################################################################
##########################     PD v PC Plots       ##########################
#############################################################################

# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
dat_pdpc.PDabun <- as.data.frame.matrix(abundances(dat_pdpc.PD)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()

dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
dat_pdpc.PCabun <- as.data.frame.matrix(abundances(dat_pdpc.PC)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.pc.sig$feature) %>% column_to_rownames()


######  Generalized or pseudo-fold calculation  ######
gfc_data <- generalized_fold_change(dat_pdpc.PDabun, dat_pdpc.PCabun)
PDovrPC <- tibble("feature" = rownames(dat_pdpc.PCabun), "gFC" = gfc_data)
###### To get order for Significant feauture plots
### Stop chunk after next line to order by gFC
PDovrPC.order <- PDovrPC[order(-PDovrPC$gFC),]
### Use for ordering by qval
PDovrPC.order <- left_join(PDovrPC.order, Maas.pd.pc.sig, by = "feature")
PDovrPC.order <- PDovrPC.order[order(-PDovrPC.order$qval),]


###### Generalized Fold Change (gFC) BarPlot ######
PDovrPC.BP <- PDovrPC.order
PDovrPC.BP <- mutate(PDovrPC.BP, direction = if_else(PDovrPC.BP$gFC > 0, "PD",
                                                     if_else(PDovrPC.BP$gFC < 0, "PC",  "error")))
PDovrPC.BP$feature <- factor(PDovrPC.BP$feature, levels = PDovrPC.order$feature)
g0 <- gfc_plot(PDovrPC.BP, cols.pdpc, alfa = 0.8)


###### Significance Plot ######
sigplot.df.pdpc <- dplyr::select(Maas.pd.pc.sig, c("feature", "pval", "qval")) %>%  melt()
sigplot.df.pdpc$feature <- factor(sigplot.df.pdpc$feature, levels = PDovrPC.order$feature)
g2 <- significance_barplot(sigplot.df.pdpc)


###### BoxPlots ######
## Prepping Significance labels
abun.pdpc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdpc, abun.pdpc.inpt)
abun.pdpc.inpt$Var2 <- factor(abun.pdpc.inpt$Var2, levels = PDovrPC.order$feature) 
g1 <- daf_boxplots(abun.pdpc.inpt, cols.pdpc, alfa = 0.2)


###### Prevalence Plot ######
# Subset of phlyoseq obj subset to get samples of interest
dat_pdpc.PD = subset_samples(dat_pdpc, donor_group =="PD")
dat_pdpc.PDprev <- tibble::enframe(prevalence(dat_pdpc.PD)) %>% filter(name %in% Maas.pd.pc.sig$feature) 
colnames(dat_pdpc.PDprev) <- c("feature", "PD")
dat_pdpc.PC = subset_samples(dat_pdpc, donor_group =="PC")
dat_pdpc.PCprev <- tibble::enframe(prevalence(dat_pdpc.PC)) %>% filter(name %in% Maas.pd.pc.sig$feature)
colnames(dat_pdpc.PCprev) <- c("feature", "PC")
dat_pdpc.PREV <- left_join(dat_pdpc.PDprev, dat_pdpc.PCprev, by = "feature") %>% melt()
dat_pdpc.PREV$feature <- factor(dat_pdpc.PREV$feature, levels = PDovrPC.order$feature)
dat_pdpc.PREV$variable <- factor(dat_pdpc.PREV$variable, levels = c("PC", "PD"))
g3 <- prevalence_barplot(dat_pdpc.PREV, cols.pdpc, alfa = 0.2)



#############################################################################
##########################     PD v HC Plots       ##########################
#############################################################################


# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PD = subset_samples(dat_pdhc, donor_group =="PD")
dat_pdhc.PDabun <- as.data.frame.matrix(abundances(dat_pdhc.PD))  %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()

dat_pdhc.HC = subset_samples(dat_pdhc, donor_group =="HC")
dat_pdhc.HCabun <- as.data.frame.matrix(abundances(dat_pdhc.HC)) %>% 
  rownames_to_column() %>%  filter(rowname %in% Maas.pd.hc.sig$feature) %>% column_to_rownames()


######  Generalized or pseudo-fold change 
gfc_data <- generalized_fold_change(pd_abundance=dat_pdhc.PDabun, 
                                    ctrl_abundances=dat_pdhc.HCabun)
PDovrHC <- tibble("feature" = rownames(dat_pdhc.PDabun), "gFC" = gfc_data)


###### To get order for Significant feauture plots
### Stop chunk after next line to order by gFC
PDovrHC.order <- PDovrHC[order(-PDovrHC$gFC),]
### Use for ordering by qval
PDovrHC.order <- left_join(PDovrHC.order, Maas.pd.hc.sig, by = "feature")
PDovrHC.order <- PDovrHC.order[order(-PDovrHC.order$qval),]


###### Generalized Fold Change (gFC) BarPlot ######
PDovrHC.BP <- PDovrHC.order
PDovrHC.BP <- mutate(PDovrHC.BP, direction = if_else(PDovrHC.BP$gFC > 0, "PD",
                                                     if_else(PDovrHC.BP$gFC < 0, "HC",  "error")))
PDovrHC.BP$feature <- factor(PDovrHC.BP$feature, levels = PDovrHC.order$feature) 
h0 <- gfc_plot(PDovrHC.BP, cols.pdhc, alfa = 0.8)


###### Significance Plot ###### 
sigplot.df.pdhc <- dplyr::select(Maas.pd.hc.sig, c("feature", "pval", "qval")) %>%  melt()
sigplot.df.pdhc$feature <- factor(sigplot.df.pdhc$feature, levels = PDovrHC.order$feature) 
h2 <- significance_barplot(sigplot.df.pdhc)


###### BoxPlots ######
## Prepping Significance labels
abun.pdhc.inpt <- daf_boxplot_sigvalues(sigplot.df.pdhc, abun.pdhc.inpt)
abun.pdhc.inpt$Var2 <- factor(abun.pdhc.inpt$Var2, levels = PDovrHC.order$feature)
h1 <- daf_boxplots(abun.pdhc.inpt, cols.pdhc, alfa = 0.2)


###### Prevalence Plot ######
# Subset of phlyoseq obj subset to get samples of interest
dat_pdhc.PDprev <- tibble::enframe(prevalence(dat_pdhc.PD)) %>% filter(name %in% Maas.pd.hc.sig$feature) 
colnames(dat_pdhc.PDprev) <- c("feature", "PD")
dat_pdhc.HCprev <- tibble::enframe(prevalence(dat_pdhc.HC)) %>% filter(name %in% Maas.pd.hc.sig$feature)
colnames(dat_pdhc.HCprev) <- c("feature", "HC")
dat_pdhc.PREV <- left_join(dat_pdhc.PDprev, dat_pdhc.HCprev, by = "feature") %>% melt()
dat_pdhc.PREV$feature <- factor(dat_pdhc.PREV$feature, levels = PDovrHC.order$feature)
dat_pdhc.PREV$variable <- factor(dat_pdhc.PREV$variable, levels = c("HC", "PD"))
h3 <- prevalence_barplot(dat_pdhc.PREV, cols.pdhc, alfa = 0.2)


#############################################################################
#############################################################################


###### Merge Panels ###### 
h0a <- h0 + theme(axis.text.y = element_blank(), legend.position = c(.80, .85))
h1a <- h1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
h2a <- h2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .85))
h3a <- h3 + theme(axis.text.y = element_blank(), legend.position = "none")

g0a <- g0 + theme(axis.text.y = element_blank(), legend.position = c(.80, .80))
g1a <- g1 + theme(legend.position = "none", axis.text.y = element_text(face = "italic"))
g2a <- g2 + theme(axis.text.y = element_blank(), legend.position = c(.90, .35))
g3a <- g3 + theme(axis.text.y = element_blank(), legend.position = "none")

g0a <- g0a + theme(axis.title.x = element_blank())
g1a <- g1a + theme(axis.title.x = element_blank())
g2a <- g2a + theme(axis.title.x = element_blank())
g3a <- g3a + theme(axis.title.x = element_blank())

# Setting plot length variables - 
top_len <- length(unique(PDovrPC.BP$feature)) + 2
bottom_len <- length(unique(PDovrHC.BP$feature)) + 2.5


DAF_part1 <- cowplot::plot_grid(g1a, h1a, nrow = 2, align = "hv", labels = "AUTO",
                       rel_heights = c(top_len, bottom_len))
# DAF_part1

DAF_part2 <- cowplot::plot_grid(g2a, g3a, g0a, h2a, h3a, h0a, nrow = 2, ncol=3, align = "h", 
                       rel_heights = c(top_len, bottom_len) )
# DAF_part2

DAF_final <- cowplot::plot_grid(DAF_part1, DAF_part2, ncol = 2, rel_widths = c(1.5, 1))
# DAF_final

ggsave(DAF_final, filename = paste0("data/DAF_Analysis/DAF_", lev, "_PDvHC_MaaslinSig.svg"),
       width = 20, height = (top_len+bottom_len)/3)


