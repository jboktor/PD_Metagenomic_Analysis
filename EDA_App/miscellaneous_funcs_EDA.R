### Miscellaneous Functions

#-----------------------------------------------------------------------------------------------------------
######## p-value significance (integer to symbol function)

sig_mapper <- function(pval, shh = F, porq = "p", symbols = T) {
  ###' Traditional mapping of p-value to symbol 
  ###' prints p-values if below significance
  
  if (symbols == T){
    if (pval <= .001) {
      sigvalue = "***"
    } else if (pval <= .01) {
      sigvalue = "**"
    } else if (pval <= .05) {
      sigvalue = "*"
    } else if (pval > .05 & shh == F) {
      sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
    } else if (pval > .05 & shh == T) {
      sigvalue = ""
    }
  } else if (symbols == F){
    sigvalue = paste0(porq, "=", format.pval(pval, digits=2)) 
  }
  return(sigvalue)
}

#-----------------------------------------------------------------------------------------------------------

sig.symbol.generator <- function(Column){
  sig.symbol <- c()
  for (i in Column){
    sig.symbol <- c(sig.symbol, sig_mapper(i))
  }
  return(sig.symbol)
}

#-----------------------------------------------------------------------------------------------------------

distribution_sanity <- function(df) {
  
  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data 
  ###' colored by donor group
  
  abund.melt <- melt(df)
  abund.melt <- 
    mutate(abund.melt, group = if_else(grepl("HC", variable), "HC",
                                if_else(grepl("PC", variable), "PC","PD")))
  
  histo_plot <- ggplot(abund.melt, aes(x=value, fill = group), alpha = 0.4) + 
    theme_minimal() +
    geom_histogram(color = "black", position="dodge", boundary = 0) +
    theme(axis.title.x = element_blank(),
          legend.position = c(0.9, 0.5))
  
  ecdf_plot <- ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    theme(legend.position = c(0.9, 0.5))
  
  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
}

#--------------------------------- Version 2

distribution_sanity2 <- function(df, binN = 30) {
  
  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data 
  ###' colored by donor group
  
  abund.melt <- melt(df)
  abund.melt <- 
    mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC",
                                       if_else(grepl("PC", Var2), "PC","PD")))
  cols=c("PD"= "#bfbfbf", 
         "PC" = "#ed7d31", 
         "HC" = "#5b9bd5")
  
  histo_plot <- ggplot(abund.melt, aes(x=value, fill = group), alpha = 0.4) + 
    theme_minimal() +
    geom_histogram(color = "black", position="dodge", boundary = 0, bins = binN) +
    scale_fill_manual(values = cols) +
    # scale_size_area(limits = c(0, 1000), max_size = 10, guide = NULL)
    theme(axis.title.x = element_blank(),
          legend.position = c(0.9, 0.5))
  
  ecdf_plot <- ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    scale_colour_manual(values = cols) +
    scale_size_area(limits = c(0, 1000), max_size = 10, guide = NULL)
    theme(legend.position ="none")
  
  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v", axis = "lb")
}


#-----------------------------------------------------------------------------------------------------------
################################  Functions to adjust feature names  ########################################
#-----------------------------------------------------------------------------------------------------------


prep_species_names <- function(phylo_obj){
  taxa_names(phylo_obj) <- gsub("s__", "", taxa_names(phylo_obj))
  return(features)
}

prep_pathway_names <- function(phylo_obj){
  features <- paste0("PATHWAY_", taxa_names(phylo_obj))
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}
prep_enzyme_names <- function(phylo_obj){
  features <- paste0("ENZYME_", taxa_names(phylo_obj))
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}
prep_ko_names <- function(phylo_obj){
  features <- taxa_names(phylo_obj)
  features <- gsub(":", ".", features)
  features <- gsub("\\|", ".", features)
  features <- gsub(" ", "_", features)
  features <- gsub("-", "_", features)
  return(features)
}

#-----------------------------------------------------------------------------------------------------------

group_col_from_ids <- function(df, ids){
  df <- mutate(df, group = if_else(grepl("HC", ids), "HC",
                                   if_else(grepl("PC", ids), "PC","PD")))
  return(df)
}


#-----------------------------------------------------------------------------------------------------------

boxplot_all <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
  ###' Basic all group boxplot function
  
  set.seed(123)
  ggplot(data=df, aes(x=x, y=y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x), position = position_jitterdodge(jitter.width = 1), 
               shape=21, size=1.5, alpha = 0.8) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank())
  
}

#-----------------------------------------------------------------------------------------------------------

alpha_div_boxplots <- function(df, x, y, 
                               df.pairs, df.pairs.x, df.pairs.y, pairs.column, 
                               cols, cols.rim, ylabel,PDvPC.stat, PDvHC.stat){
  
  ###' Function for alpha diveristy
  ###' boxplots with paired lines
  
  set.seed(123)
  p <- ggplot(data = df, aes(x = x, y = y)) + 
    theme_minimal() + 
    geom_point(aes(fill = x, color = x), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 1) +
    geom_boxplot(aes(fill = x, color = x), width=0.3, alpha = 0.1, outlier.alpha = 0) +
    geom_line(data = df.pairs, aes(x = df.pairs.x, y = df.pairs.y, group = pairs.column), 
              linetype = 'solid', color = "grey", alpha = 0.7) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab(ylabel) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    labs(fill="Group") +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols.rim) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(PDvHC.stat)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(PDvPC.stat))
  return(p)
}





#-----------------------------------------------------------------------------------------------------------
# Inspired by MicrobiomeAnalystR: https://github.com/xia-lab/MicrobiomeAnalystR


# Plot (0.1 - 0.9) features by rank : to help decide on percentage cutoff 

PlotVariance <- function(dat) {
  
  int.mat <- abundances(dat) %>% as.data.frame()
  filter.val <- apply(int.mat, 1, function (x) {
    diff(quantile(as.numeric(x), c(0.1, 0.9), na.rm = FALSE, names = FALSE, 
                  type = 7)) })

  var.df <- as.data.frame(filter.val) %>% 
    rownames_to_column(var = "features")
  
  rk <- rank(-filter.val, ties.method='random')
  rws <-  nrow(var.df)

  p <- ggplot(var.df, aes(x= reorder(features, -filter.val), y= filter.val)) +
    geom_point(color="#1170aa") + 
    ggthemes::theme_clean() +
    labs(x = "Ranked Features", y = "[0.1 - 0.9] Quantile Range") +
    geom_vline(xintercept = c(rk[rk == round(rws*.9)], rk[rk == round(rws*.8)], rk[rk == round(rws*.7)], 
                              rk[rk == round(rws*.6)], rk[rk == round(rws*.5)], rk[rk == round(rws*.4)],
                              rk[rk == round(rws*.3)], rk[rk == round(rws*.2)]), linetype = "dashed", alpha = 0.7 ) +
    theme(axis.text.x = element_blank())
  return(p)
}


#-----------------------------------------------------------------------------------------------------------

LowVarianceFilter <- function(dat, filter.percent = 0.1) {
  
  #' This function filters features by their
  #' Inter-quartile range - larger values indicate larger spread
  #' features are ranked by IQR and a specified percentage is trimmed
  
  int.mat <- abundances(dat) %>% as.data.frame()
  filter.val <- apply(int.mat, 1, function (x) {
    diff(quantile(as.numeric(x), c(0.1, 0.9), na.rm = FALSE, names = FALSE, 
                  type = 7)) })
  
  rk <- rank(-filter.val, ties.method='random')
  var.num <- nrow(int.mat);
  remain <- rk <= var.num*(1-filter.percent);
  int.mat <- int.mat[remain,];
  
  cat("A total of", sum(!remain), "low variance features were removed based on the Quantile Range between [0.1 - 0.9]. \n")
  cat("The number of features remaining after filtering is:", nrow(int.mat), "\n")
  
  return(int.mat)

}



#-----------------------------------------------------------------------------------------------------------
#                           Functions to explore select features in dataset
#-----------------------------------------------------------------------------------------------------------

explore_table <- function(obj){
  data.table <- obj %>%
    microbiome::abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() 
  return(data.table)
}

plot_feature <- function(obj, feature){
  
  d <- obj %>% abundances() %>% 
    as.data.frame() %>% 
    rownames_to_column() %>% 
    filter(rowname == feature) %>% 
    melt()
  
  dm <-  group_col_from_ids(d, id=d$variable)
  
  dm$value <- asin(sqrt(dm$value))
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
  
  boxplot_all(dm, x=dm$group, y=dm$value, 
              cols=c("PC"= "#bfbfbf", 
                     "PD" = "#ed7d31", 
                     "HC" = "#5b9bd5"), 
              title=" ", 
              ylabel= paste(unique(dm$rowname), "Abundance"))
}

#-----------------------------------------------------------------------------------------------------------


