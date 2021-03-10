### Miscellaneous Functions


#-------------------------------------------------------------------------------
#####                Functions to adjust feature names    ########
#-------------------------------------------------------------------------------

decode_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df %>% 
    dplyr::mutate(fullnames = gsub("\\.gc.", ":", !!features)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gsqrr.", "\\]", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gsqrl.", "\\[", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gplus.", "\\+", fullnames)) %>%
    dplyr::mutate(fullnames = gsub("\\.gpar.", "\\'", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gpr.", ")", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gpl.", "(", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gm.", ",", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gp.", "\\|", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gs.", " ", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gh.", "-", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("\\.gd.", "/", fullnames)) %>% 
    dplyr::mutate(fullnames = gsub("feat_", "", fullnames))
}

#-------------------------------------------------------------------------------

group_col_from_ids <- function(df, ids){
  df <- dplyr::mutate(df, group = dplyr::if_else(grepl("HC", ids), "HC",
                                   dplyr::if_else(grepl("PC", ids), "PC","PD")))
  return(df)
}


#-------------------------------------------------------------------------------

boxplot_all <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel){
  
  blank.title = " "; blank.ylabel = " "
  group.cols = c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")

  set.seed(123)
  plot <- 
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
cat("functions loaded ..\n")

