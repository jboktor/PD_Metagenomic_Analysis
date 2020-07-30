### QC Functions


#--------------------------------------------------------------------------------------------------
# Get Pseduo-Counts from clean reads and abundance table
#--------------------------------------------------------------------------------------------------

PseudoCounts <- function(dat, reads){
  
  psudocnts <- dat %>% microbiome::transform("compositional") %>% 
    microbiome::abundances() %>% as.data.frame()
  cat("TSS \n")
  
  for (i in colnames(psudocnts)){
    donor_reads <- reads[[which(reads$id == i), 2]]
    psudocnts[i] <- psudocnts[i] * donor_reads
  }
  cat("Pseudocount Transformation Complete\n")
  return(psudocnts)
}




#--------------------------------------------------------------------------------------------------
# Rarefaction Analysis 
#--------------------------------------------------------------------------------------------------

RareFactionPlot <- function(dat, featuretype="Species", reads){
  
  cat("\n\n\n"); cat(featuretype, "Rarefaction  \n")

  # Get Pseuo-counts
  psudocnts <- dat %>% transform("compositional") %>% 
    abundances() %>% as.data.frame()
  cat("TSS \n")
  
  for (i in colnames(psudocnts)){
    donor_reads <- reads[[which(reads$id == i), 2]]
    psudocnts[i] <- psudocnts[i] * donor_reads
  }
  cat("Pseudocount Estimation \n")
  
  # Filter 
  psudocnts.HC <- psudocnts %>% dplyr::select(contains("HC")) %>% 
    t() %>% as.data.frame()
  psudocnts.PC <- psudocnts %>% dplyr::select(contains("PC")) %>% 
    t() %>% as.data.frame()
  psudocnts.PD <- psudocnts %>% dplyr::select(!contains(c("HC", "PC"))) %>% 
    t() %>% as.data.frame()
  
  cat("Calculating Rarefaction Estimate for HC : This may take a second -   \n")
  acc.HC <- specaccum(psudocnts.HC, method = "exact")
  cat("Calculating Rarefaction Estimate for PC - Almost there - ༼ つ ◕_◕ ༽つ  \n")
  acc.PC <- specaccum(psudocnts.PC, method = "exact")
  cat("Calculating Rarefaction Estimate for PD - Homestretch -  ༼ つ ಥ_ಥ ༽つ  \n\n")
  acc.PD <- specaccum(psudocnts.PD, method = "exact")
  cat(featuretype, " Rarefaction Calculations Complete:  ヽ༼ຈل͜ຈ༽ﾉ  \n\n")
  
  df.acc.HC <- data.frame(Sites=acc.HC$sites, Richness=acc.HC$richness, SD=acc.HC$sd)
  df.acc.PC <- data.frame(Sites=acc.PC$sites, Richness=acc.PC$richness, SD=acc.PC$sd)
  df.acc.PD <- data.frame(Sites=acc.PD$sites, Richness=acc.PD$richness, SD=acc.PD$sd)
  
  PD.col <- "#FDE725"; PD.col2 <- "#fdad19"; PD.col3 <- "#d48a02"
  PC.col <- "#21908C"; PC.col2 <- "#2cc0bb"
  HC.col <- "#440154"; HC.col2 <- "#73028e"
  
  p1 <- ggplot() +
    theme_bw() +
    # PD data
    geom_point(data=df.acc.PD, aes(x=Sites, y=Richness), alpha=1.5, color = PD.col2) +
    geom_line(data=df.acc.PD, aes(x=Sites, y=Richness), size = 2, alpha=0.6, color = PD.col3) +
    geom_ribbon(data=df.acc.PD, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = PD.col) +
    # PC data
    geom_point(data=df.acc.PC, aes(x=Sites, y=Richness), alpha=1.5, color = PC.col) +
    geom_line(data=df.acc.PC, aes(x=Sites, y=Richness), size = 2, alpha=0.6, color = PC.col) +
    geom_ribbon(data=df.acc.PC, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = PC.col2) +
    # HC data
    geom_point(data=df.acc.HC, aes(x=Sites, y=Richness), alpha=1.5, color = HC.col) +
    geom_line(data=df.acc.HC, aes(x=Sites, y=Richness), size = 2, alpha=0.6, color = HC.col2) +
    geom_ribbon(data=df.acc.HC, aes(x=Sites, ymin=(Richness-2*SD),ymax=(Richness+2*SD)),alpha=0.2, fill = HC.col2) +
    labs(x = "Sample #", y = paste0(featuretype, " Detected")) +
    theme(strip.background = element_blank(),
          panel.grid = element_blank())
  
  ggsave(p1, filename = paste0("data/Quality_Control/Rarefaction_Plots/Rarefaction_Plot_", featuretype, ".svg"),
         width = 8, height = 6)
  return(p1)
  
}

###################################################################################################


AlphaLinearRegressionQuantilePlot <- function(df, x, x2, y, color, fill, ylabel, title){
  
  PD.col <- "#bfbfbf"
  PC.col <- "#ed7d31"
  HC.col <- "#5b9bd5"
  
  p1 <- ggplot(df, aes(x=x, y=y, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=ylabel, fill="group") +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), 
                 color = color), label.x = 15e6, label.y.npc=0.25, hjust = 0) +
    ggtitle(title) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          legend.background = element_blank(),
          legend.position = c(0.9, 0.5))
  p1
  
  p2 <- ggplot(df, aes(x=x2, y=y, color=color)) +
    geom_boxplot() +
    labs(x="Filtered Sample Read Quantiles", y=ylabel) +
    theme_bw() +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())
  p2
  
  c1 <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
  return(c1)
}

###################################################################################################

BetaLinearRegressionPlot <- function(df, x, y, y2, color, fill, feature, title){
  
  PD.col <- "#bfbfbf"
  PC.col <- "#ed7d31"
  HC.col <- "#5b9bd5"
  
  
  p1 <- ggplot(df, aes(x=x, y=y, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=paste0("PCoA Axis 1: ", feature)) +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
                 color = color), label.x = 1.75e7, label.y.npc=0.25, hjust = 0) +
    ggtitle(title) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  
  p2 <- ggplot(df, aes(x=x, y=y2, color=color, fill=fill)) +
    geom_smooth(method="lm", se=F) +
    geom_point(aes(fill=fill),shape=21, size=2, alpha = 0.9) +
    labs(x="Filtered Sample Reads", y=paste0("PCoA Axis 2: ", feature), fill="group") +
    theme_bw() +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"),
                 color = color), label.x = 1.75e7, label.y.npc=0.25, hjust = 0) +
    scale_fill_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col)) +
    scale_color_manual(values = c("HC" = HC.col, "PD" = PD.col, "PC" = PC.col), guide = FALSE)  +
    theme(legend.position = c(0.85, 0.8),
          panel.grid = element_blank())
  
  c1 <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
  return(c1)
}


###################################################################################################


