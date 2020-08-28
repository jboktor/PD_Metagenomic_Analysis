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
  
  Legends <- cowplot::plot_grid(get_legend(Phylum_legend))
  
  products <- list("Bars" = full_plot, "Legends" = Legends, "Axis.order" = unique(phylo.plot$Species))
  
  return(products)
  
}

############################################################################################################


boxplot_hl1.bars <- function(inpt.hl1, tile.cols){
  
  ##' Function takes in features selected for plotting 
  ##' Returns color bars, legends, and y-axis order separately, 
  ##' mapping each feature's Metabolic Module and it's higher level hierarchy
  ##' Adds significance symbol if below threshold in respective
  ##' MaAsLin2 model
  
  full_plot <- NULL
  Legends <- NULL
  products <- NULL
  
  ## Sorting - HL1 then qval
  hl1.plot <- inpt.hl1[order(inpt.hl1$HL1, inpt.hl1$qval), ]
  hl1.plot$order_col <- 1:nrow(hl1.plot)
  hl1.plot$xaxis <- "level_1"
  
  ## Select colors from pre-made list
  cut_colors <- tile.cols[names(tile.cols) %in% hl1.plot$HL1]
  
  full_plot <- 
    ggplot(hl1.plot, aes(x= xaxis, y=reorder(Name, -order_col), fill= HL1)) + 
    geom_tile() +
    theme_minimal() +
    scale_fill_manual(values=cut_colors) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position="none")
  
  Phylum_legend <-
    hl1.plot %>%
    ggplot(aes(x=xaxis, y=reorder(Name, -order_col), fill= HL1)) +
    geom_tile() +
    scale_fill_manual(values = cut_colors, name = "HL1")
  
  Legends <- cowplot::plot_grid(get_legend(Phylum_legend))
  
  products <- list("Bars" = full_plot, "Legends" = Legends, "Axis.order" = unique(hl1.plot$Name))
  
  return(products)
  
}


############################################################################################################
## Adapted from: https://github.com/zellerlab/crc_meta


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
  p <- ggplot(data=df, aes(x=gFC, y= feature, fill = direction)) +
    geom_point(aes(fill=direction), shape=21, size=3, alpha = alfa) +
    geom_segment(aes(x=0, xend=gFC, y=feature, yend=feature, color=direction)) +
    theme_minimal() +
    xlim(-max(abs(df$gFC)), max(abs(df$gFC))) +
    labs(x="Average Difference") +
    ggtitle("Generalized Fold Change") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    scale_color_manual(values = manual_colors, name ="Group") +
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
  sigplot.df.QVAL$value <- round(sigplot.df.QVAL$value, digits = 3)
  
  sig.labs <- c()
  for (i in sigplot.df.QVAL$value){
    sig.labs <- c(sig.labs, sig_mapper(i, porq = "q", symbols = F))
  }
  sigplot.df.QVAL$sig.labels <- sig.labs
  sigplot.df.QVAL <- sigplot.df.QVAL %>% dplyr::select(feature, sig.labels) %>% 
    dplyr::rename(Var2 = feature) 
  abund.input2 <- left_join(abund.input, sigplot.df.QVAL, by = "Var2")
  
  return(abund.input2)
  
}

############################################################################################################

daf_boxplots <- function(df, fill_cols, rim_cols, alfa = 0.5){
  
  set.seed(123)
  p <- ggplot(data=df, aes(x=value, y= Var2)) +
    geom_boxplot(aes(fill = group), alpha = alfa, outlier.alpha = 0, width = 0.8) +
    geom_point(aes(fill=group, color=group), position = position_jitterdodge(jitter.width = .2), shape=21, size=1, alpha = 0.9) +
    theme_minimal() +
    ggtitle(paste0("Differential Abundance: ", lev)) +
    geom_text(aes(x= max(value) + max(value)*0.15, 
                  y=Var2, label = sig.labels), size = 4, check_overlap = TRUE) +
    labs(x = "Abundance") +
    scale_fill_manual(values = fill_cols, name ="Group") +
    scale_color_manual(values = rim_cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          panel.grid.major.y = element_blank())
  
  return(p)
}

############################################################################################################

prevalence_barplot <- function(df, manual_colors, alfa = 0.5){
  p <- ggplot(data=df, aes(x=value*100, y= feature, fill=variable, color=variable)) +
    geom_bar(position = position_dodge(),  alpha = alfa, width = 0.8,
             stat = "identity") +
    theme_minimal() +
    labs(x=expression('Prevalence [%]')) +
    ggtitle("Prevalence") +
    scale_fill_manual(values = manual_colors, name ="Group") +
    scale_color_manual(values = manual_colors, name ="Group") +
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


DAF_Analysis <- function(obj.name, obj){
  
  source("src/miscellaneous_funcs.R")
  
  # ## Color Schemes
  cols.pdpc <- c("PD"= "#ed7d31", "PC" = "#bfbfbf")
  cols.pdhc <- c("PD"= "#ed7d31", "HC" = "#5b9bd5")
  # Rims
  cols.pdpc.rim <- c("PD"= "#ed4e31", "PC" = "#999999")
  cols.pdhc.rim <- c("PD"= "#ed4e31", "HC" = "#5b7dd5")
  
  
  # To protect against row/colname errors 
  if (obj.name == "Pathways"| obj.name == "Pathways.slim") {
    features <- paste0("PATHWAY_", taxa_names(obj)) 
    features <- gsub(":", ".", features)
    features <- gsub("\\|", ".", features)
    features <- gsub(" ", "_", features)
    features <- gsub("-", "_", features)
    taxa_names(obj) <- features
  } else if (obj.name == "Enzymes" | obj.name == "Enzymes.slim") {
    features <- paste0("ENZYME_", taxa_names(obj)) 
    features <- gsub(":", ".", features)
    features <- gsub("\\|", ".", features)
    features <- gsub(" ", "_", features)
    features <- gsub("-", "_", features)
    taxa_names(obj) <- features
  } else if (obj.name == "Species") {
    taxa_names(obj) <- gsub("s__", "", taxa_names(obj))
  } else if (obj.name == "GMMs" | obj.name == "GBMs" ) {
    taxa_names(obj) <- gsub(" ", ".", taxa_names(obj))
  } else {
    features <- taxa_names(obj)
    features <- gsub(":", ".", features)
    features <- gsub("\\|", ".", features)
    features <- gsub(" ", "_", features)
    features <- gsub("-", "_", features)
    taxa_names(obj) <- features
  }
  
  # Check process above
  # taxa_names(obj)
  
  ################################################################################# 
  
  ### Visualization Transformations
  # taxa_names(obj) <- gsub("s__", "", taxa_names(obj))
  dat_obj <- microbiome::transform(obj, "compositional")
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
  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", obj.name, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Population Control")
  Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)
  
  Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/", obj.name, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
    filter(value == "Household Control")
  Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) 
  
  if (obj.name == "Species") {
    Maas.pd.pc.sig$feature <- gsub("s__", "", Maas.pd.pc.sig$feature)
    Maas.pd.hc.sig$feature <- gsub("s__", "", Maas.pd.hc.sig$feature)
  }
  
  ########################## Venn Diagram Plots ##########################
  
  v <- VennPlot(Maas.pd.pc.sig, Maas.pd.hc.sig, qval_threshold = 0.1)
  
  # Save Venn Diagrams
  if (!is.null(v$venn_depleted)){
    pdf(file = paste0("data/DAF_Analysis/DAF_", obj.name, "_VennDiagram_PD_depleted.pdf"),
        width = 7, 
        height = 5,
        pointsize = 12)
    # units = "in", pointsize = 12, res=300)
    plot(v$venn_depleted)
    dev.off()
  }
  
  if (!is.null(v$venn_enriched)){
    pdf(file = paste0("data/DAF_Analysis/DAF_", obj.name, "_VennDiagram_PD_enriched.pdf"),
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
  print(DAF_final)
  
  ggsave(DAF_final, filename = paste0("data/DAF_Analysis/DAF_", obj.name, "_PDvHC_MaaslinSig.svg"),
         width = 20, height = (top_len+bottom_len)/3)
}

