### QC_Seq_Depth Analysis

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/Community_Composition_Functions.R")
load("files/Pfams.slim_PhyloseqObj.RData")
load("files/Eggnogs.slim_PhyloseqObj.RData")

func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
reads$id <- gsub("_", ".", reads$id)

#################################  Plot Rarefaction Curves per group #################################  

##' Manual Input here: to create the species, genus, phylum, pathways, enzyme, KO, Eggnog, or Pfam 
##' Rarefaction plot input the appropriate __ dat __ frame below

RareFactionPlot(dat, featuretype = "Species", reads)

RareFactionPlot(dat.path.slim, featuretype = "non-stratified Pathways", reads)

RareFactionPlot(dat.KOs.slim, featuretype = "non-stratified KOs", reads)

RareFactionPlot(dat.Pfams.slim, featuretype = "non-stratified Pfams", reads)

# This may crash desktop
# RareFactionPlot(dat.Eggnogs.slim, featuretype = "non-stratified Eggnogs", reads)

################################################################################################
#################################  Plot Seq Depth Distributions by group #################################  

PD.col <- "#bfbfbf"; PC.col <- "#ed7d31"; HC.col <- "#5b9bd5"

df.reads <- group_col_from_ids(reads, reads$id)
df.reads$group <- factor(df.reads$group, levels=c("PC", "PD", "HC"))

histo_plot <- ggplot(df.reads, aes(x=clean_total_reads, color = group), alpha = 0.3) + 
  theme_bw() +
  geom_density() +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col))  +
  theme(axis.title.x = element_blank(),
        legend.position = c(0.9, 0.5),
        panel.grid = element_blank())

ecdf_plot <- ggplot(df.reads, aes(x=clean_total_reads, colour = group)) + 
  stat_ecdf(geom = "step", pad = FALSE) +
  stat_ecdf(geom = "point", pad = FALSE, alpha =0.4) +
  theme_bw() +
  labs(y = "ECDF") +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
  theme(legend.position = c(0.9, 0.5),
        panel.grid = element_blank())

c1 <- cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
ggsave(c1, filename = "data/Quality_Control/Sequencing_Depth_Density_&_ECDF.png",
       width = 6, height = 5, dpi = 600)

################################################################################################
################################# Seq Depth Boxplots by group #################################  

PD.col <- "#FDE725";PD.col2 <- "#fdad19";PD.col3 <- "#d48a02"
PC.col <- "#21908C";PC.col2 <- "#2cc0bb"
HC.col <- "#440154";HC.col2 <- "#73028e"

# Test for normality - Distb is non-normal
shapiro.test(df.reads$clean_total_reads)
# KW AND Wilcoxon test
KW <- kruskal.test(clean_total_reads ~ group, data = df.reads)
DunnB <- dunnTest(clean_total_reads ~ group, data = df.reads, method="bh")
PCvsHC.stat <- DunnB$res$P.adj[1]
PCvsPD.stat <- DunnB$res$P.adj[3]
### Paired Wilcoxon Test between PD and HC 
m <- meta(dat) %>% dplyr::select(Paired) %>% rownames_to_column(var = "id")
paired.donors <- left_join(df.reads, m, var = "id") %>% 
  filter(Paired != "No")
stat.test <- wilcox.test(clean_total_reads ~ group, 
                         paired = F, data = paired.donors)
HCvsPD.stat <- stat.test$p.value


p2 <- ggplot(df.reads, aes(x=group, y=clean_total_reads)) +
  geom_point(aes(fill = group), position = position_jitterdodge(dodge.width=1),shape=21, size=1.25, alpha = 1) +
  geom_boxplot(aes(fill = group), outlier.alpha = 0, alpha = 0.3, width = 0.45) +
  theme_minimal() +
  labs(y="Total Clean Reads per Sample") +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  geom_signif(comparisons = list(c("HC", "PC")), tip_length = 0, 
              annotations = sig_mapper(PCvsHC.stat,  symbols = F, porq = "q")) +
  geom_signif(comparisons = list(c("PC", "PD")), tip_length = 0.01, 
              y_position = 1.95e7, annotations = sig_mapper(PCvsPD.stat,  symbols = F, porq = "q")) +
  geom_signif(comparisons = list(c("HC", "PD")), tip_length = 0.01, 
              y_position = 1.95e7, annotations = sig_mapper(HCvsPD.stat, symbols = F, porq = "p")) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank())

ggsave(p2, filename = "data/Quality_Control/Sequencing_Depth_Boxplt.svg",
       width = 3, height = 6)

################################################################################################
#################################  Plot Abundance Bars by group #################################  

dat2 <- core(dat.genus, detection = 0, prevalence = 0.1)
# Abundance filter for top 20 Genra
dat.top <- get_top_taxa(dat2, n=15, relative = TRUE, discard_other = TRUE,
                              other_label = "Other")
ps1 <- merge_samples(dat.top, "donor_group")

p <- fantaxtic_bar(ps1, color_by = "Order", label_by = "Genus", 
                   palette =  c("#c9ca9c", "#8d5576", "#15476c", "#f2e8d5"))
p <- p + theme(axis.title.x = element_blank()) +
  labs(y = "Relative Abundance")

ggsave(p, filename = "figures/Figure_1/SummaryBarPlot.svg",
       width = 2.5, height = 4)

################################################################################################
#####################  Seq-Depth x Alpha-Diversity : Linear Regression  ######################## 

# Prep Reads 
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")

# Prep Metadata 
env <- meta(dat)
env <- left_join(env, reads, by = "id")
env <- mutate(env, clean_total_reads_factor = as.factor(cut(env$clean_total_reads, 
                              breaks = quantile(env$clean_total_reads), 
                              labels=c(1,2,3,4), include.lowest=T) ))

env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))


# Objects to loop over
obj <- c(dat, dat.path.slim, dat.ec.slim, dat.KOs.slim, dat.Pfams.slim, dat.Eggnogs.slim)
obj.label <- c("Species", "Pathways", "Enzymes", "KOs", "Pfams", "Eggnogs")

# Initiate count
cnt <- 1

for (i in obj) {
  
  cat("\n\n"); cat("Processing", obj.label[cnt], "for Alpha Diversity Analysis \n\n")
  
  # Prep abundnace table
  dat_alpha <-  transform(i, "compositional")
  ## Calculate Alpha Diversity Metrics and add cols to df
  env$observed <- alpha(abundances(dat_alpha), 'observed')$observed
  env$shannon <- alpha(abundances(dat_alpha), 'shannon')$diversity_shannon
  env$evenness <- evenness(abundances(dat_alpha), 'simpson')$simpson
  

  observed.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$observed,
                                                color=env$donor_group, fill=env$donor_group, 
                                                ylabel=paste("Observed", obj.label[cnt]), title= paste("Observed", obj.label[cnt], "by Read Depth"))
  
  shannon.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$shannon,
                                               color=env$donor_group, fill=env$donor_group,
                                               ylabel=paste("Shannon Diversity:", obj.label[cnt]), title="Shannon Index by Read Depth")
  
  evenness.plot <- AlphaLinearRegressionQuantilePlot(df=env, x=env$clean_total_reads, x2=env$clean_total_reads_factor, y=env$evenness,  
                                                color=env$donor_group, fill=env$donor_group,
                                                ylabel= paste("Simpsons Evenness:", obj.label[cnt]), title="Simpsons Evenness by Read Depth")
  
  QC_alphaVsreads <- cowplot::plot_grid(observed.plot, shannon.plot, evenness.plot, 
                                        ncol = 3, align = "h")
  
  ggsave(QC_alphaVsreads, filename = paste0("data/Quality_Control/Read_Depth_vs_Alpha_Diversity/Read_Depth_vs_",  obj.label[cnt], "_AlphaDiversity.png"),
         width = 16, height = 8)
  
  cnt <- cnt + 1
}


################################################################################################
#####################  Seq-Depth x Beta-Diversity : Linear Regression  ######################## 


source("src/Community_Composition_Functions.R")
# Initiate count
cnt <- 1

for (i in obj) {
  
  ## AITCHISONS DISTANCE
  cat("\n\n"); cat("Processing", obj.label[cnt], "for Beta Diversity Analysis \n\n")
  
  obj_clr <- microbiome::transform(i, "clr")
  iDist <- distance(obj_clr, method="euclidean")
  iMDS  <- ordinate(obj_clr, "MDS", distance=iDist)
  #  PCoA for Axis 1 and 2
  
  p <- plot_ordination(obj_clr, iMDS, color="description", axes = c(1, 2))
  pcoa.data = p$data
  df.pcoa <- left_join(pcoa.data, reads, by = "id")
  df.pcoa$donor_group <- factor(df.pcoa$donor_group, levels=c("PC", "PD", "HC"))

  QC_betaVsreads <- BetaLinearRegressionPlot(df=df.pcoa, x=df.pcoa$clean_total_reads, y=df.pcoa$Axis.1, y2=df.pcoa$Axis.2,
                           color=df.pcoa$donor_group, fill=df.pcoa$donor_group,
                           feature=obj.label[cnt], title=paste0(obj.label[cnt], ": Beta Diversity by Read Depth"))

 ggsave(QC_betaVsreads, filename = paste0("data/Quality_Control/Read_Depth_vs_Beta_Diversity/Read_Depth_vs_",  obj.label[cnt], "_BetaDiversity.png"),
         width = 8, height = 8)
  
  cnt <- cnt + 1
}





################################################################################################
#################################  Distrubution Plots by Quantile  ################################# 

# Goal is to see if there are any taxa detected at deeper sequencing at one group and not others
# In other words, is a deeper sequencing responsible for detection of rare taxa that are significant?

############## Species Abundnace by Read Depth : Distributions between donor group
## Pull abundance data - merge with clean_total_reads_factor & donor_group -> melt df

env <- mutate(env, clean_total_reads_factor_tight = 
                as.factor(cut(env$clean_total_reads, 
                              breaks = quantile(env$clean_total_reads, probs=seq(0, 1, 0.2)),
                              labels= c(1:5), include.lowest=TRUE) ))

seqdat <- dplyr::select(env, c("donor_id", "clean_total_reads_factor", "clean_total_reads_factor_tight", "donor_group"))

ab <- dat %>% 
  microbiome::transform("compositional") %>% abundances() %>% 
  t() %>% as.data.frame.matrix() %>% 
  rownames_to_column() %>% cbind(seqdat)
abm <- melt(ab)

# ArcSinSqrt Transform Abundances
abm$value <- asin(sqrt(abm$value))

# Bin Abundance by Quantiles
abm <-  mutate(abm, abundance_factor = 
                 if_else(abm$value > 0.4, 4,
                         if_else(abm$value > 0.1, 3,
                                 if_else(abm$value > 0, 2, 1 
                                         ))))

abm$abundance_factor <- as.factor(abm$abundance_factor)
PD.col <- "#bfbfbf"; PC.col <- "#ed7d31"; HC.col <- "#5b9bd5"

z1 <- abm %>% filter(abundance_factor != 1 ) %>% 
  ggplot(aes(x=value, color = donor_group, fill = donor_group), alpha = 0.4) +
  geom_histogram(position="dodge", boundary = 0) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
  theme_classic() +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "ECDF by Binned Feature Abundance") +
  ggtitle("Sequencing Depth Quantiles") +
  theme(plot.title = element_text(hjust = 0.5))
# z1
ggsave(z1, filename = "data/Quality_Control/SeqDepth_&_Abundance_Quantile_Histogram_Matrix.png",
       width = 9, height = 4.5, dpi = 1200)

z2 <- abm %>% filter(abundance_factor != 1 ) %>% 
  ggplot(aes(x=value, color = donor_group)) +
  stat_ecdf(geom = "point", pad = FALSE, alpha =0.1) +
  stat_ecdf(geom = "step", pad = FALSE) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "Binned Feature Abundance - ECDF") +
  ggtitle("Sequencing Depth Quantiles") +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col))  +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
# z2
ggsave(z2, filename = "data/Quality_Control/SeqDepth_&_Abundance_Quantile_ECDF_matrix.png",
       width = 20, height = 8, dpi = 1200)




z1b <- abm %>% filter(abundance_factor == 2 & clean_total_reads_factor_tight == 1 ) %>% 
  ggplot(aes(x=value, color = donor_group, fill = donor_group), alpha = 0.4) +
  geom_histogram(position="dodge", boundary = 0) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
  theme_classic() +
  ggtitle("Sequencing Depth: First Quantile vs AST transformed Abundance (0 < x < 0.1)") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
# z1b

z2b <- abm %>% filter(abundance_factor == 2 & clean_total_reads_factor_tight == 1 ) %>%
  ggplot(aes(x=value, color = donor_group)) +
  stat_ecdf(geom = "step", pad = FALSE) +
  facet_grid(cols=vars(clean_total_reads_factor_tight), 
             rows=vars(abundance_factor), scales = "free") +
  labs(x = "Arcsin(Sqrt(Abundance))", y = "Binned Feature Abundance - ECDF") +
  scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col))  +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
# z2b

QC_speciesAbundanceVsreads <- cowplot::plot_grid(z1b, z2b, align = "h", nrow = 2)
# QC_speciesAbundanceVsreads

ggsave(QC_speciesAbundanceVsreads, filename = "data/Quality_Control/SeqDepth_Quantile_Distribution.png",
       width = 8, height = 6, dpi = 1200)



