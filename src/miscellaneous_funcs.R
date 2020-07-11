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
  ###' - can specify if both or only one plot is desired
  
  abund.melt <- melt(df)
  abund.melt <- 
    mutate(abund.melt, group = if_else(grepl("HC", variable), "HC",
                                if_else(grepl("PC", variable), "PC","PD")))
  
  histo_plot <- ggplot(abund.melt, aes(x=value, fill = group), alpha = 0.4) + 
    theme_minimal() +
    geom_histogram(color = "black", position="dodge") +
    theme(axis.title.x = element_blank(),
          legend.position = c(0.9, 0.5))
  
  ecdf_plot <- ggplot(abund.melt, aes(x=value, colour = group)) + stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    theme(legend.position = c(0.9, 0.5))
  
  plot_grid(histo_plot, ecdf_plot, ncol = 1, align="v")
}





