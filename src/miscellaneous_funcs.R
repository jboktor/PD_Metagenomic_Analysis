### Miscellaneous Functions

############################################################################################################
######## p-value significance (integer to symbol function)

sig_mapper <- function(pval, shh = F, porq = "p") {
  ###' Traditional mapping of p-value to symbol 
  ###' prints p-values if below significance
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
  return(sigvalue)
}
############################################################################################################
######## p-value significance (integer to symbol function)

sig_mapper2 <- function(pval, shh = F) {
  ###' Traditional mapping of p-value to symbol 
  ###' prints p-values if below significance
  if (pval <= .01) {
    sigvalue = "***"
  } else if (pval <= .05) {
    sigvalue = "**"
  } else if (pval <= .1) {
    sigvalue = "*"
  } else if (pval > .1 & shh == F) {
    sigvalue = paste0("p=", format.pval(pval, digits=2)) 
  } else if (pval > .1 & shh == T) {
    sigvalue = ""
  }
  return(sigvalue)
}

############################################################################################################
sig.symbol.generator <- function(Column){
  sig.symbol <- c()
  for (i in Column){
    sig.symbol <- c(sig.symbol, sig_mapper(i))
  }
  return(sig.symbol)
}

############################################################################################################
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

############################################################################################################

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
  
  plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
}


############################################################################################################
################################  Functions to adjust feature names  ########################################
############################################################################################################


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

############################################################################################################

group_col_from_ids <- function(df, ids){
  df <- mutate(df, group = if_else(grepl("HC", ids), "HC",
                                   if_else(grepl("PC", ids), "PC","PD")))
  return(df)
}


############################################################################################################

boxplot_all <- function(df, x, y, cols, title, ylabel){
  
  ###' Basic all group boxplot function
  
  df$x <- factor(df$x, levels = c("HC", "PD","PC") )
  
  set.seed(123)
  ggplot(data=df, aes(x=x, y=y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x), position = position_jitterdodge(jitter.width = 1), 
               shape=21, size=1.5, alpha = 0.8) +
    theme_minimal() +
    ggtitle(title) +
    labs(y = ylabel) +
    scale_color_manual(values = cols, name ="Group") +
    scale_fill_manual(values = cols, name ="Group") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          panel.grid.major.y = element_blank())
  
}


############################################################################################################