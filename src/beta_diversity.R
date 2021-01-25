# Aitchisons Distance Analysis

####### Load PhyloSeq objects  ####### 
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")

beta_diversity_summary <- function(cohort){
  
  load_betadiv_colors() 
  
  x <- c(dat.species, dat.path, dat.ec, dat.KOs, 
         # dat.EGGNOGs, dat.PFAMs,
         dat.path.slim, dat.ec.slim, dat.KOs.slim, dat.EGGNOGs.slim, dat.PFAMs.slim)

  z <- c("Species", "Pathways", "Enzymes", "KOs", 
         # "Eggnogs", "Pfams",
         "Pathways.slim", "Enzymes.slim", "KOs.slim", "Eggnogs.slim", "Pfams.slim")
  
  
  # --------------------------------------------------- 
  # Aitchisons PCoA/Ridgeline/Violin Plots Loop 
  # --------------------------------------------------- 
  
  cnt <- 1
  for (i in x){
    cat("\n\n\n")
    cat("Processing input: ", z[cnt], "\n")
    cat("\n")
    
    load("files/low_quality_samples.RData")
    i <- i %>% 
      subset_samples(donor_id %ni% low_qc[[1]])
    
    print(i)
    cat("\n")
    
    obj <- microbiome::transform(i, "compositional") # transforms species abundance table as (x/sum(x))
    obj_clr <- microbiome::transform(i, "clr") # transforms species abundnce table as (x/geometric-Mean(x)) # Allows for variable independence from Sample total
    
    # AITCHISONS DISTANCE (EUCLIDIAN DIST ON CLR TRANSFORMED DATA)
    iDist <- phyloseq::distance(obj_clr, method="euclidean")
    iMDS  <- phyloseq::ordinate(obj_clr, "MDS", distance=iDist)
    
    # --------------------------------------------------- 
    #  PCoA for Axis 1 and 2
    # --------------------------------------------------- 
    
    p <- plot_ordination(obj_clr, iMDS, color="description", axes = c(1, 2))
    
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
      scale_fill_manual(values = cols.pdpchc_dark) +
      scale_color_manual(values = cols.pdpchc.rim) +
      theme(plot.title = element_text(hjust = 0.5), 
            panel.grid = element_blank())
    
    cat(paste0("Complete:  Aitchison Distance PCoA: ", z[cnt], " abundance \n"))
    
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
      scale_color_manual(values = cols.pdpchc_dark) +
      scale_fill_manual(values = cols.pdpchc) +
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
    
    cat(paste0("Complete:  Aitchison Distance PCoA Boxplot 1 \n"))
    
    # --------------------------------------------------- 
    #  Boxplot - Axis 2
    # --------------------------------------------------- 
    
    r2 <- ggplot(df12, aes(x = Axis.2, y = description)) +
      geom_boxplot(aes(color = description, fill = description), alpha = 0.2) +
      scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
      scale_x_continuous(limits = my.ggp.yrange2) +
      coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
      scale_color_manual(values = cols.pdpchc_dark) +
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
    
    cat(paste0("Complete:  Aitchison Distance PCoA Boxplot 2 \n"))
    
    # --------------------------------------------------- 
    #  Beta Diversity Violin Plots
    # --------------------------------------------------- 
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
    wu.sd$Type2 <- factor(wu.sd$Type2, levels = c("Population Control", "PD Patient", "Household Control"))
    
    v <- ggplot(wu.sd, aes(x = Type2, y = value)) + theme_minimal() + 
      # geom_beeswarm(aes(color = Type2, fill = Type2), shape = 21, size=.3, alpha = .5, cex = .5) +
      geom_violin(aes(color = Type2), draw_quantiles = c(0.5), trim = T, alpha=0) +
      theme(axis.title.x=element_blank(),
            legend.position = "none") +
      ylab("Aitchison Distance") +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(fill="Group") +
      scale_color_manual(values = cols.pdpchc) +
      scale_fill_manual(values = cols.pdpchc) +
      stat_compare_means(comparisons = my_comparisons, label = "p.signif", tip.length = 0.02, step.increase = 0)+ 
      ggtitle(paste0("Aitchison Distance to Population Control\n", z[cnt], " abundance"))
    
    cat(paste0("Complete:  Aitchison Distance Violin Plot \n"))
    
    # --------------------------------------------------- 
    #  Assemble Plots
    # --------------------------------------------------- 
    cat(paste0("Assembling Cowplot Summary for : " , z[cnt], "\n"))
    ord.plot <- ord + theme(legend.position = "none")
    cow1 <- cowplot::plot_grid(r1, NULL, ord.plot, r2, nrow = 2, ncol = 2, rel_heights = c(1, 5), rel_widths = c(5, 1),  align = "vh")
    cow2 <- cowplot::plot_grid(v, cow1, nrow = 1, rel_widths = c(1, 2.75), align = "h")
    print(cow2)
    cow2
    ggsave(plot = cow2, 
           filename = paste0("data/Community_Composition/Beta_Diversity_Analysis/CowplotSummary_",  z[cnt], "_", cohort, ".svg"),
           width = 13, height = 9)
    
    cnt <- cnt + 1
  }
}



