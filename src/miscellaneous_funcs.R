### Miscellaneous Functions


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

sig.symbol.generator <- function(Column){
  sig.symbol <- c()
  for (i in Column){
    sig.symbol <- c(sig.symbol, sig_mapper(i))
  }
  return(sig.symbol)
}


stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}

distribution <- function(obj) {
  microbiome::abundances(obj)
  microbiome::meta(obj)
  
}



