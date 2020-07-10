# Differentially Abundant Features (DAF) Functions


############################################################################################################

VennPlot <- function(Maaslin_PDvPC, Maaslin_PDvHC, qval_threshold = 0.25){
  
  ###' Function reads in maaslin2 significance results and creates Venn-Diagrams
  ###' for shared and unique features between Household & Population controls in 
  ###' reference to PD patients
  
  ## Note qval_threshold must be 0.25 or less
  require(eulerr)
  
  # Initalize Return variables
  venn_depleted <- NULL
  venn_enriched <- NULL
  
  ## Split both inputs into enriched and depleted features
  Maaslin_PDvPC <- Maaslin_PDvPC %>% filter(qval <= qval_threshold)
  Maaslin_PDvHC <- Maaslin_PDvHC %>% filter(qval <= qval_threshold)
  
  ## If either of the data frames have 0 features after qval-filtering, Stop and print error:
  if (nrow(Maaslin_PDvPC) < 1 | nrow(Maaslin_PDvHC) < 1){
    cat("One or both of comparisons doesn't meet the qval threshold \n")
    return()
  } else {
    cat("### Sufficient features for plotting \n")
  }
  
  PC.model <- mutate(Maaslin_PDvPC, direction = if_else(coef > 0, "Depleted","Enriched"))
  HC.model <- mutate(Maaslin_PDvHC, direction = if_else(coef > 0, "Depleted","Enriched")) 
  a.pc <- dplyr::filter(PC.model, direction == "Depleted")
  b.pc <- dplyr::filter(PC.model, direction == "Enriched")
  a.hc <- dplyr::filter(HC.model, direction == "Depleted")
  b.hc <- dplyr::filter(HC.model, direction == "Enriched")
  
  
  if (nrow(a.pc) > 0 & nrow(a.hc) > 0 ){
    ## Create joined matrix of features 
    down.vs.pc.df <- data.frame("features"=a.pc$feature, "Down_vs_PC" = TRUE)
    down.vs.hc.df <- data.frame("features"=a.hc$feature, "Down_vs_HC" = TRUE)
    depleted_matrix <-full_join(down.vs.pc.df, down.vs.hc.df,  by="features")
    depleted_matrix[is.na(depleted_matrix)] <-  F
    
    venn_depleted <- plot(euler(depleted_matrix[-1], shape = "ellipse"), quantities = TRUE)
    cat("1) PD DEPLETED FEATURES: Plotting shared PD depleted features \n")
    # print(venn_depleted)
    
  } else {
    cat("1) PD DEPLETED FEATURES: At least one comparison without features: no plot \n")
  }
  
  
  if (nrow(b.pc) > 0 & nrow(b.hc) > 0 ){
    
    ## Create joined matrix of features 
    up.vs.pc.df <- data.frame("features"=b.pc$feature, "Up_vs_PC" = T)
    up.vs.hc.df <- data.frame("features"=b.hc$feature, "Up_vs_HC" = T)
    enriched_matrix <-full_join(up.vs.pc.df, up.vs.hc.df,  by="features")
    enriched_matrix[is.na(enriched_matrix)] <-  F
    
    venn_enriched <- plot(euler(enriched_matrix[-1], shape = "ellipse"), quantities = TRUE)
    cat("2) PD ENRICHED FEATURES: Plotting shared PD enriched features \n")
    # print(venn_enriched)
    
  } else {
    cat("2) PD ENRICHED FEATURES: At least one comparison without features: no plot \n")
  }
  
  plots2return <- list( "venn_depleted" = venn_depleted, "venn_enriched" = venn_enriched)

  cat("Venn Diagram Plotting Complete: \n\n")
  return(plots2return)
  
}

############################################################################################################

taxa_genus_phlyum_annotation = function(physeq, selectedTaxa){
  
  ###'  Returns Phylum level of a selected list 
  ###'  of species from Phlyoseq Object
  
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% selectedTaxa)]
  physeq <- prune_taxa(allTaxa, physeq)
  tax.df <- as.data.frame(tax_table(physeq)[,2])
  
  return(tax.df)
}



############################################################################################################

pull.phylo.sig <- function(comparison, phylo.names){
  
  ###' Function takes in a comparion (either HC or PC) and table
  ###' containing Genus and Phylum annotations 
  ###' Returns long format table with features and their qvals from maaslin
  ###' at the Genus and Phylum levels
  
  output <- NULL
  
  if (comparison == "PC"){
    groop <- "Population Control"
  } 
  else if (comparison == "HC"){
    groop <- "Household Control"
  } 
  else {
    return("ERROR: Please enter either PC or HC for comparison: ")
  }
  
  ### Read-in Maaslin Files - PHYLUM LEVEL
  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/Phylum_PDv", comparison, "_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == groop)
  Maas.pd.pc$feature <- gsub("s__", "", Maas.pd.pc$feature)
  Phylum.qvals <- Maas.pd.pc %>% filter(feature %in% unique(phylo.names$Phylum))
  vals <- Phylum.qvals %>% dplyr::select(feature, qval)
  
  colnames(vals)[1] <- "Phylum"
  return(vals)
  
}

############################################################################################################

boxplot_phylobars <- function(inpt.phylo, sigvals, tile.cols){
  
  ##' Function takes in features selected for plotting 
  ##' Returns color bars, legends, and y-axis order separately, 
  ##' mapping each feature's Phylum and Genus
  ##' Adds significance symbol if below threshold in respective
  ##' MaAsLin2 model
  full_plot <- NULL
  Legends <- NULL
  products <- NULL
  
  phylo.plot <- inpt.phylo %>% rownames_to_column(var = "Species")
  phylo.plot$Species <- gsub("s__", "", phylo.plot$Species)
  phylo.plot <-  left_join(phylo.plot, sigvals, by = "Species")
  

  ## Sorting - Phylum then qval
  phylo.plot <- phylo.plot[order(phylo.plot$Phylum, phylo.plot$qval), ]
  phylo.plot$order_col <- 1:nrow(phylo.plot)
  phylo.plot$xaxis <- "Phylum"
  
  ## Select colors from pre-made list
  cut_colors <- tile.cols[names(tile.cols) %in% phylo.plot$Phylum]

  full_plot <- 
    ggplot(phylo.plot, aes(x= xaxis, y=reorder(Species, -order_col), fill= Phylum)) + 
    geom_tile() +
    theme_minimal() +
    scale_fill_manual(values=cut_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none")
  
  Phylum_legend <-
    phylo.plot %>%
    ggplot(aes(x=xaxis, y=reorder(Species, -order_col), fill= Phylum)) +
    geom_tile() +
    scale_fill_manual(values = cut_colors, name = "Phylum")
  
  Legends <- plot_grid(get_legend(Phylum_legend))
  
  products <- list("Bars" = full_plot, "Legends" = Legends, "Axis.order" = unique(phylo.plot$Species))
  
  return(products)
  
}

############################################################################################################
## Adapted from: 


generalized_fold_change <- function(pd_abundance, ctrl_abundances){
  # Initalize counts
  i = 1
  probs.fc <- seq(.1, .9, .05)
  gfc_data <- c()
  
  for (feature in 1:nrow(ctrl_abundances)){
    # Loops through each species row and calculates 
    # the Generalized fold change for each feature
    cat("Feature number: ", feature, "\n")
    cat("Testing PD: ", rownames(pd_abundance)[feature], "feature vs ", rownames(ctrl_abundances)[feature], "\n")
    q.pd <- quantile(pd_abundance[feature,], probs = probs.fc)
    q.ctl <- quantile(ctrl_abundances[feature,], probs = probs.fc)
    
    gfc <- sum(q.pd - q.ctl) / length(q.ctl)
    print(gfc)
    
    gfc_data <- c(gfc_data, gfc)
  }
  return(gfc_data)
}



############################################################################################################

gfc_plot <- function(df, manual_colors, alfa = 0.5){
  ######' Generalized Fold Change (gFC) BarPlot ######
  p <- ggplot() +
    geom_col(data=df, aes(x=gFC, y= feature,  fill = direction), 
             position = "dodge", width = 0.6, alpha = alfa) +
    theme_minimal() +
    labs(x="Generalized Fold Change") +
    ggtitle("Fold Change") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  return(p)
}

############################################################################################################

significance_barplot <- function(df){
  ######' Significance BarPlot ######
  p <- ggplot() +
    geom_col(data=df, aes(x=-log10(value), y= feature,  fill = variable), position = "dodge", width = 0.8) +
    theme_minimal() +
    labs(x=expression('-log'[10]*'(value)')) +
    ggtitle("Significance") +
    scale_fill_manual(values = c("pval" = "#d3d3d3", "qval" = "#676767")) +
    geom_vline(xintercept = -log10(0.1), linetype = "dotted", color = "red") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  return(p)
}

############################################################################################################

daf_boxplot_sigvalues <- function(sigplot.df, abund.input){
  
  sigplot.df.QVAL <- filter(sigplot.df, variable == "qval")
  
  sig.labs <- c()
  for (i in sigplot.df.QVAL$value){
    sig.labs <- c(sig.labs, sig_mapper(i, shh = F, porq = "q"))
  }
  sigplot.df.QVAL$sig.labels <- sig.labs
  sigplot.df.QVAL <- sigplot.df.QVAL %>% dplyr::select(feature, sig.labels) %>% 
    dplyr::rename(Var2 = feature) 
  abund.input2 <- left_join(abund.input, sigplot.df.QVAL, by = "Var2")
  
  return(abund.input2)
  
}

############################################################################################################

daf_boxplots <- function(df, manual_colors, alfa = 0.5){
  
  set.seed(123)
  p <- ggplot(data=df, aes(x=value, y= Var2)) +
    geom_boxplot(aes(fill = group), alpha = alfa, outlier.alpha = 0, width = 0.8) +
    geom_point(aes(fill=group), position = position_jitterdodge(jitter.width = .2), shape=21, size=1, alpha = 0.6) +
    theme_minimal() +
    ggtitle(paste0("Differential Abundance: ", lev)) +
    geom_text(aes(x= max(value) + max(value)*0.15, 
                  y=Var2, label = sig.labels), size = 4, check_overlap = TRUE) +
    labs(x = "Abundance") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  return(p)
}

############################################################################################################

prevalence_barplot <- function(df, manual_colors, alfa = 0.5){
  p <- ggplot() +
    geom_col(data=df, aes(x=value*100, y= feature, fill = variable), position = "dodge",  alpha = alfa, width = 0.8) +
    theme_minimal() +
    labs(x=expression('Prevalence [%]')) +
    ggtitle("Prevalence") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())

  return(p)
}

############################################################################################################
# Top and Bottom Row Mods for Cowplot
# remove x-axis labels from top row and titles from bottom row



############################################################################################################







############################################################################################################
## FUNCTION TO ADD GROUP COLUMNS FOR D3 FlashWeave NETWORK

flashweave_group_colors <- function(features, metadata_added, Maas.pd.pc.sig, Maas.pd.hc.sig){
  
  ##' Loop for identifying significant features and stratifying them
  ##' by enrichment or depletion in PC, HC, or both groups in reference 
  ##' to PD groups 
  ##' Note: a negative coefficient in Maaslin model indicates higher levels in PD
  
  features <- features %>% 
    mutate(group = if_else(Node %in% metadata_added, "Donor_Group",
                           if_else(Node %in% Maas.pd.pc.sig$feature & Node %in% Maas.pd.hc.sig$feature, "BOTH",
                                   if_else(Node %in% Maas.pd.pc.sig$feature, "PC",
                                           if_else(Node %in% Maas.pd.hc.sig$feature, "HC",
                                                   "None")))))
  # Initiate vars for loop
  n = 1
  group2 <- c()
  for (node in features$Node){
    if (features$group[n] == "BOTH") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef  < 0 | 
          filter(Maas.pd.hc.sig, feature == node)$coef  < 0) {
        group2  = c(group2, "Up_PDvBoth")
      }  else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0 & 
                  filter(Maas.pd.hc.sig, feature == node)$coef  > 0) {
        group2  = c(group2, "Down_PDvBoth")
      }
    } else if (features$group[n] == "PC") {
      if (filter(Maas.pd.pc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "Up_PDvPC")
      } else if (filter(Maas.pd.pc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "Down_PDvPC")
      }
    } else if (features$group[n] == "HC") {
      if (filter(Maas.pd.hc.sig, feature == node)$coef < 0) {
        group2  = c(group2, "Up_PDvHC")
      }  else if (filter(Maas.pd.hc.sig, feature == node)$coef > 0) {
        group2  = c(group2, "Down_PDvHC")
      }
    } else if (features$group[n] == "Donor_Group") {
      group2  = c(group2, "Donor_Group")
    } else if (features$group[n] == "None") {
      group2  = c(group2, "None")
    } else {
      group2 = c(group2, "error")
    }
    n = n + 1
  }
  features$group <- group2
  
  ## Complex color codes - DAF PD/PC (up & down) and PD/HC feature (up and down) coloring -
  features2 <- features %>%
    mutate(group_color = if_else(group == "Donor_Group", "#d62728", 
                                 if_else(group == "Up_PDvBoth", "#98df8a",
                                         if_else(group == "Down_PDvBoth", "#2ca02c",
                                                 if_else(group == "Up_PDvPC", "#ffbb78",
                                                         if_else(group == "Down_PDvPC", "#ff7f0e",
                                                                 if_else(group == "Up_PDvHC", "#aec7e8",
                                                                         if_else(group == "Down_PDvHC", "#1f77b4",
                                                                                 if_else(group == "None", "#c7c7c7",
                                                                                         "Error")))))))))
  
  return(features2)
  
}

############################################################################################################

