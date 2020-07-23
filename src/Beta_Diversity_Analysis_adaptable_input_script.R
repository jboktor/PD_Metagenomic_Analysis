# Aitchisons Distance Analysis

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot)

rm(list = ls())
####### Load PhyloSeq objects  ####### 
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
load("files/Eggnogs_PhyloseqObj.RData")
load("files/Eggnogs.slim_PhyloseqObj.RData")
load("files/Pfams_PhyloseqObj.RData")
load("files/Pfams.slim_PhyloseqObj.RData")

x <- c(dat, dat.path.slim, dat.KOs.slim, dat.Pfams.slim,  dat.Eggnogs.slim)
z <- c("species", "pathways", "KOs", "Pfams", "Eggnogs")


# ### Distribution Sanity Check - CLR
# for (i in x){
#   print(i)
#   A <- i %>% transform("compositional") %>%
#     transform("clr") %>%
#     abundances() %>% data.frame()
#   plot(distribution_sanity(A))
# }


PD.col <- "#FDE725FF"
PC.col <- "#21908CFF"
HC.col <- "#440154FF"

#################################   Aitchisons PCoA/Ridgeline/Violin Plots Loop (Species, Pathways, Enzymes, Genes-KOs)   ################################# 

cnt <- 1
for (i in x){
  cat("\n\n\n")
  cat("Processing input: ", z[cnt], "\n")
  cat("\n")
  print(i)
  cat("\n")
  
  # print(x[[1]])
  obj <- microbiome::transform(i, "compositional") # transforms species abundance table as (x/sum(x))
  obj_clr <- microbiome::transform(i, "clr") # transforms species abundnce table as (x/geometric-Mean(x)) # Allows for variable independence from Sample total

  ################################# AITCHISONS DISTANCE (EUCLIDIAN DIST ON CLR TRANSFORMED DATA)  ################################# 
  iDist <- distance(obj_clr, method="euclidean")
  iMDS  <- ordinate(obj_clr, "MDS", distance=iDist)
  
  
  #################################  PCoA for Axis 1 and 2
  p <- plot_ordination(obj_clr, iMDS, color="description", axes = c(1, 2))
  
  df12 = p$data
  df12$donor_group <- factor(df12$donor_group, levels=c("PC", "PD", "HC"))
  
  p <- ggplot(df12, aes(Axis.1, Axis.2, fill = description))
  p <- p + geom_point(shape=21, size=5, alpha=0.7, color="black")
  ord <- p + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Donor Group") +
    xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
    ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
    labs(fill="Donor Group") +
    scale_fill_manual(values = c("Household Control" = HC.col, "PD Patient" = PD.col, "Population Control" = PC.col))
  ord
  # ggsave(paste0("data/Beta_Diversity_Analysis/PCoA_Aitchisons_", z[cnt], ".svg"), height = 6, width =9)
  cat(paste0("Complete:  Aitchison Distance PCoA: ", z[cnt], " abundance \n"))
  
  ## PCoA data to fit ridgeline Plots
  my.ggp.xrange <- ggplot_build(ord)$layout$panel_scales_x[[1]]$range$range # For PCoA1
  my.ggp.yrange2 <- ggplot_build(ord)$layout$panel_scales_y[[1]]$range$range # For PCoA2
  
  
  ################################# RidgeLine Plot - Axis 1 #################################

  r1 <- ggplot(df12, aes(x = Axis.1, y = description)) +
    geom_density_ridges(aes(color = description, fill = description),
                        jittered_points = TRUE, quantile_lines = TRUE,
                        position = position_raincloud(adjust_vlines = TRUE),
                        alpha = 0.1, scale = 10) +
    xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 1), "%)")) +
    scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(limits = my.ggp.xrange) +
    coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    scale_color_manual(values = c("Household Control" = HC.col,
                                  "PD Patient" = "#d48a02",
                                  "Population Control" = PC.col)) +
    scale_fill_manual(values = c("Household Control" = HC.col,
                                 "PD Patient" = PD.col,
                                 "Population Control" = PC.col)) +
    theme_classic() +
    ggtitle(paste0("Aitchison Distance PCoA")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x=element_blank(),
          legend.position = "none")
  r1


  # ggsave(paste0("data/Beta_Diversity_Analysis/PCoA1_Aitchisons_RidgePlot_", z[cnt], ".svg"), height = 2, width =9)
  cat(paste0("Complete:  Aitchison Distance PCoA Ridgeline 1 \n"))
  
  ################################# RidgeLine Plot - Axis 2 #################################

  r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
    geom_density_ridges(aes(color = description, fill = description),
                        jittered_points = TRUE, quantile_lines = TRUE,
                        position = position_raincloud(adjust_vlines = TRUE),
                        alpha = 0.1, scale = 10) +
    scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
    scale_x_continuous(limits = my.ggp.yrange2) +
    coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
    scale_color_manual(values = c("Household Control" = HC.col,
                                  "PD Patient" = "#d48a02",
                                  "Population Control" = PC.col)) +
    scale_fill_manual(values = c("Household Control" = HC.col,
                                 "PD Patient" = PD.col,
                                 "Population Control" = PC.col)) +
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
  r2

  # ggsave(paste0("data/Beta_Diversity_Analysis/PCoA2_Aitchisons_RidgePlot_", z[cnt], ".svg"), height = 6, width =2)
  cat(paste0("Complete:  Aitchison Distance PCoA Ridgeline 2 \n"))
  
  ############## ##############  Beta Diversity Violin/Box Plots ############## ############## ##############
  
  p <- obj_clr
  m <- "euclidean"
  s <- "donor_id"
  d <- "description"
  require("phyloseq");require("dplyr");require("reshape2");require("ggplot2")
  
  # Make Melted Distance Matrix 
  wu <- phyloseq::distance(obj_clr, method="euclidean")
  wu.m <- melt(as.matrix(wu))
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
  my_comparisons <- list( c("Household Control", "PD Patient"),c("Population Control", "PD Patient"))
  
  v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
    geom_violin(draw_quantiles = c(0.5), trim = T) +
    geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size=.3, alpha = .5, cex = .5) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab("Aitchison Distance") +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Group") +
    scale_color_manual(values = c("Household Control" = "#340141", 
                                  "PD Patient" = "#d48a02", "Population Control" = "#217890")) +
    scale_fill_manual(values = c("Household Control" = HC.col, 
                                 "PD Patient" = PD.col, "Population Control" = PC.col)) +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0)+ 
    ggtitle(paste0("Aitchison Distance to Population Control\n", z[cnt], " abundance"))
  v
  # ggsave(paste0("data/Beta_Diversity_Analysis/Aitchison_Dissimilarity_ViolinPlot_PC_relative_",  z[cnt], ".svg"), 
  #         height = 10, width =4)
  cat(paste0("Complete:  Aitchison Distance ViolinPlot \n"))
  
  
  cat(paste0("Assembling Cowplot Summary for : " , z[cnt], "\n"))
  ord.plot <- ord + theme(legend.position = "none")
  cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 3), rel_widths = c(9/2, 1),  align = "hv")
  cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 3), align = "h")
  cow2
  ggsave(cow2, filename = paste0("data/Beta_Diversity_Analysis/CowplotSummary_",  z[cnt], ".svg"),
         width = 14, height = 10)
  
  cnt <- cnt + 1
}


