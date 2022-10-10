### Miscellaneous Functions


# _______________________________________________________________________________-
# Aesthetic variables

## Color Schemes
cols.pdpc <- c("PD" = "#bfbfbf", "PC" = "#ed7d31")
cols.pdhc <- c("PD" = "#bfbfbf", "HC" = "#5b9bd5")
cols.pdpchc <- c("PD" = "#bfbfbf", "PC" = "#ed7d31", "HC" = "#5b9bd5")
cols.pdpchc.dark <- c("PD" = "#494949", "PC" = "#ed7d31", "HC" = "#5b9bd5")

# Rims
cols.pdpc.rim <- c("PD" = "#494949", "PC" = "#c15811")
cols.pdhc.rim <- c("PD" = "#494949", "HC" = "#2e75b5")
cols.pdpchc.rim <- c("PD" = "#494949", "PC" = "#c15811", "HC" = "#2e75b5")

load_alphadiv_colors <- function() {
  cols.pdpchc <- c(
    "PD" = "#494949",
    "PC" = "#ed7d31",
    "HC" = "#5b9bd5"
  )
  cols.pdpchc.rim <- c(
    "PD" = "#494949",
    "PC" = "#c15811",
    "HC" = "#2e75b5"
  )
  assign("cols.pdpchc", cols.pdpchc, envir = .GlobalEnv)
  assign("cols.pdpchc.rim", cols.pdpchc.rim, envir = .GlobalEnv)
}

load_betadiv_colors <- function() {
  cols.pdpchc <- c(
    "PD Patient" = "#bfbfbf",
    "Population Control" = "#ed7d31",
    "Household Control" = "#5b9bd5"
  )
  cols.pdpchc.dark <- c(
    "PD Patient" = "#494949",
    "Population Control" = "#ed7d31",
    "Household Control" = "#5b9bd5"
  )
  cols.pdpchc.rim <- c(
    "PD Patient" = "#494949",
    "Population Control" = "#c15811",
    "Household Control" = "#2e75b5"
  )

  assign("cols.pdpchc", cols.pdpchc, envir = .GlobalEnv)
  assign("cols.pdpchc.dark", cols.pdpchc.dark, envir = .GlobalEnv)
  assign("cols.pdpchc.rim", cols.pdpchc.rim, envir = .GlobalEnv)
}

colormash <- c(
  # Set 1
  "#e41a1c", "#377eb8", "#4daf4a", "#984ea3",
  "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999",
  # Accent
  "#7fc97f", "#beaed4", "#fdc086", "#ffff99",
  "#386cb0", "#f0027f", "#bf5b17", "#666666",
  # Dark
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666"
)

# _______________________________________________________________________________-

my_clean_theme <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(1, 1, 1, 2), "cm")
    )
  return(th)
}

my_clean_theme2 <- function() {
  th <- ggplot2::theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5),
      plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
    )
  return(th)
}


# _______________________________________________________________________________
# Save legend and plot as separate files to fix width:
save_me_cleanly <- function(ggobj, filename, plot_w, plot_h, leg_w, leg_h,
                            filetype = ".svg") {

  # require(cowplot)
  # require(ggpubr)

  my_legend <- get_legend(ggobj) %>% as_ggplot()
  ggobj_out <- ggobj + theme(legend.position = "none")
  plot_filename <- paste0(filename, filetype)
  legend_filename <- paste0(filename, "__Legend", filetype)

  ggsave(ggobj_out, filename = plot_filename, width = plot_w, height = plot_h)
  ggsave(my_legend, filename = legend_filename, width = leg_w, height = leg_h)
}

# _______________________________________________________________________________
# Remove Objects from environment

remove_dats <- function() {
  obj <-
    c(
      "dat.species",
      "dat.genus",
      "dat.order",
      "dat.class",
      "dat.phylum",
      "dat.kingdom",
      "dat.path",
      "dat.ec",
      "dat.KOs",
      "dat.EGGNOGs",
      "dat.PFAMs",
      "dat.path.slim",
      "dat.ec.slim",
      "dat.KOs.slim",
      "dat.EGGNOGs.slim",
      "dat.PFAMs.slim"
    )
  withCallingHandlers(rm(list = obj), warning = function(w) {
    invokeRestart("muffleWarning")
  })
}


print_line <- function() {
  cat("---------------------------------------------------------------------\n")
}

# _______________________________________________________________________________
# Load number of clean sample reads
# _______________________________________________________________________________

# load_reads <- function(){
#   func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
#   reads <- dplyr::select(func_reads, c("# samples","total reads")) %>%
#     dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
#   return(reads)
# }

load_reads <- function(cohort) {
  negative_controls <- c(
    "S00A4-ATCC_MSA_1003_S96",
    "S00A4-neg2_S119",
    "S00A4-neg3_S125",
    "S00A4-neg_S118",
    "S00A4NegExt_P00A4_S94",
    "S00A4NegH2O_P00A4_S95",
    "S00A4_stagPos_S117",
    "BLANK"
  )

  TBC_keys <- read.csv(file = "files/metadata_keys.csv", header = TRUE) %>%
    dplyr::select(c(MBI_Sample_ID, id)) %>%
    mutate(id = gsub("_", ".", id)) %>%
    mutate(MBI_Sample_ID = as.character(MBI_Sample_ID)) %>%
    mutate(id = as.character(id)) %>%
    dplyr::rename(`# samples` = MBI_Sample_ID)

  RUSH_keys <- read.csv(file = "files/metadata_phyloseq_RUSH.csv", header = TRUE) %>%
    dplyr::filter(study_group == "PD") %>%
    dplyr::select(donor_id, host_subject_id) %>%
    dplyr::mutate(donor_id = as.character(donor_id)) %>%
    dplyr::mutate(host_subject_id = as.character(host_subject_id)) %>%
    dplyr::rename(`# samples` = host_subject_id)

  if (cohort == "TBC") {
    func_reads_TBC <-
      read_tsv(
        "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads <-
      func_reads_TBC %>%
      dplyr::filter(`# samples` %ni% negative_controls) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>%
      left_join(TBC_keys, by = "# samples") %>%
      janitor::clean_names() %>%
      dplyr::select(-number_samples) %>%
      dplyr::rename("donor_id" = "id") %>%
      return(reads)
  } else if (cohort == "RUSH") {
    func_reads_RUSH <-
      read_tsv(
        "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads <-
      func_reads_RUSH %>%
      dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>%
      dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>%
      left_join(RUSH_keys, by = "# samples") %>%
      janitor::clean_names() %>%
      dplyr::select(-number_samples)
    return(reads)
  } else if (cohort == "shanghai") {
    func_reads <-
      read_tsv(
        "files/biobakery_output_SHANGHAI_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads <-
      func_reads %>%
      dplyr::rename(run = `# samples`) %>%
      left_join(shanghai_keys, by = "run") %>%
      dplyr::select(c("donor_id", "total reads")) %>%
      dplyr::rename("clean_total_reads" = "total reads")
    return(reads)
  } else if (cohort == "Merged") {
    func_reads_TBC <-
      read_tsv(
        "files/TBC_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads_TBC <-
      func_reads_TBC %>%
      dplyr::filter(`# samples` %ni% negative_controls) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 10)) %>%
      left_join(TBC_keys, by = "# samples") %>%
      dplyr::select(c("id", "total reads")) %>%
      dplyr::rename("donor_id" = "id", "clean_total_reads" = "total reads") %>%
      dplyr::mutate(cohort = "TBC")
    func_reads_RUSH <-
      read_tsv(
        "files/RUSH_biobakery_output_slim/humann/counts/humann_read_and_species_count_table.tsv",
        col_names = T
      )
    reads_RUSH <-
      func_reads_RUSH %>%
      dplyr::filter(str_detect(`# samples`, "BLANK", negate = TRUE)) %>%
      dplyr::filter(str_detect(`# samples`, "MSA", negate = TRUE)) %>%
      dplyr::mutate(`# samples` = substr(`# samples`, 1, 6)) %>%
      left_join(RUSH_keys, by = "# samples") %>%
      dplyr::select(c("donor_id", "total reads")) %>%
      dplyr::rename("clean_total_reads" = "total reads") %>%
      dplyr::mutate(cohort = "RUSH")
    reads <- rbind(reads_TBC, reads_RUSH)

    return(reads)
  }
}
# _______________________________________________________________________________
######## p-value significance (integer to symbol function)

sig_mapper <- function(pval, shh = F, porq = "p", symbols = T) {
  ###' Traditional mapping of p-value to symbol
  ###' prints p-values if below significance

  if (symbols == T) {
    if (is.na(pval)) {
      sigvalue <- ""
    } else if (pval <= .001) {
      sigvalue <- "***"
    } else if (pval <= .01) {
      sigvalue <- "**"
    } else if (pval <= .05) {
      sigvalue <- "*"
    } else if (pval > .05 & shh == F) {
      sigvalue <- paste0(porq, "=", format.pval(pval, digits = 2))
    } else if (pval > .05 & shh == T) {
      sigvalue <- ""
    }
  } else if (symbols == F) {
    sigvalue <- paste0(porq, "=", format.pval(pval, digits = 2))
  }
  return(sigvalue)
}

# _______________________________________________________________________________

sig.symbol.generator <- function(Column, porq = "p", shh = F) {
  sig.symbol <- c()
  for (i in Column) {
    sig.symbol <- c(sig.symbol, sig_mapper(i, porq = porq, shh = shh))
  }
  return(sig.symbol)
}

# _______________________________________________________________________________

distribution_sanity <- function(df) {

  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data
  ###' colored by donor group

  abund.melt <- melt(df)
  abund.melt <-
    mutate(abund.melt, group = if_else(grepl("HC", variable), "HC",
      if_else(grepl("PC", variable), "PC", "PD")
    ))

  histo_plot <- ggplot(abund.melt, aes(x = value, fill = group), alpha = 0.4) +
    theme_minimal() +
    geom_histogram(color = "black", position = "dodge", boundary = 0) +
    theme(
      axis.title.x = element_blank(),
      legend.position = c(0.9, 0.5)
    )

  ecdf_plot <- ggplot(abund.melt, aes(x = value, colour = group)) +
    stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    theme(legend.position = c(0.9, 0.5))

  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align = "v")
}

#--------------------------------- Version 2

distribution_sanity2 <- function(df, binN = 30) {

  ###' input an abundance table and displays
  ###' a histogram and ECDF distribution of data
  ###' colored by donor group

  abund.melt <- melt(df)
  abund.melt <-
    mutate(abund.melt, group = if_else(grepl("HC", Var2), "HC",
      if_else(grepl("PC", Var2), "PC", "PD")
    ))
  cols <- c(
    "PC" = "#bfbfbf",
    "PD" = "#ed7d31",
    "HC" = "#5b9bd5"
  )

  histo_plot <- ggplot(abund.melt, aes(x = value, fill = group), alpha = 0.4) +
    theme_minimal() +
    geom_histogram(color = "black", position = "dodge", boundary = 0, bins = binN) +
    scale_fill_manual(values = cols) +
    theme(
      axis.title.x = element_blank(),
      legend.position = c(0.9, 0.5)
    )

  ecdf_plot <- ggplot(abund.melt, aes(x = value, colour = group)) +
    stat_ecdf(geom = "step", pad = FALSE) +
    theme_minimal() +
    labs(y = "ECDF") +
    scale_colour_manual(values = cols) +
    theme(legend.position = "none")

  cowplot::plot_grid(histo_plot, ecdf_plot, ncol = 1, align = "v")
}

# _______________________________________________________________________________
#                   Functions to adjust feature names              ----
# _______________________________________________________________________________


make_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df.out <- df %>%
    mutate(simplenames = gsub(":", ".gc.", !!features)) %>%
    mutate(simplenames = gsub("\\|", ".gp.", simplenames)) %>%
    mutate(simplenames = gsub(" ", ".gs.", simplenames)) %>%
    mutate(simplenames = gsub("-", ".gh.", simplenames)) %>%
    mutate(simplenames = gsub("/", ".gd.", simplenames)) %>%
    mutate(simplenames = gsub("\\]", ".gsqrr.", simplenames)) %>%
    mutate(simplenames = gsub("\\[", ".gsqrl.", simplenames)) %>%
    mutate(simplenames = gsub("\\)", ".gpr.", simplenames)) %>%
    mutate(simplenames = gsub("\\(", ".gpl.", simplenames)) %>%
    mutate(simplenames = gsub(",", ".gm.", simplenames)) %>%
    mutate(simplenames = gsub("\\+", ".gplus.", simplenames)) %>%
    mutate(simplenames = gsub("\\'", ".gpar.", simplenames)) %>%
    mutate(simplenames = paste0("feat_", simplenames)) %>%
    tibble::column_to_rownames(var = "simplenames") %>%
    dplyr::select(-!!features)
  return(df.out)
}


decode_rfriendly_rows <- function(df, passed_column) {
  features <- rlang::sym(passed_column)
  df.out <- df %>%
    mutate(fullnames = gsub("\\.gc.", ":", !!features)) %>%
    mutate(fullnames = gsub("\\.gsqrr.", "\\]", fullnames)) %>%
    mutate(fullnames = gsub("\\.gsqrl.", "\\[", fullnames)) %>%
    mutate(fullnames = gsub("\\.gplus.", "\\+", fullnames)) %>%
    mutate(fullnames = gsub("\\.gpar.", "\\'", fullnames)) %>%
    mutate(fullnames = gsub("\\.gpr.", ")", fullnames)) %>%
    mutate(fullnames = gsub("\\.gpl.", "(", fullnames)) %>%
    mutate(fullnames = gsub("\\.gm.", ",", fullnames)) %>%
    mutate(fullnames = gsub("\\.gp.", "\\|", fullnames)) %>%
    mutate(fullnames = gsub("\\.gs.", " ", fullnames)) %>%
    mutate(fullnames = gsub("\\.gh.", "-", fullnames)) %>%
    mutate(fullnames = gsub("\\.gd.", "/", fullnames)) %>%
    mutate(fullnames = gsub("feat_", "", fullnames))
  return(df.out)
}

decoded_abundance <- function(datObj) {
  df <- datObj %>%
    microbiome::abundances() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    decode_rfriendly_rows(passed_column = "rowname") %>%
    column_to_rownames(var = "fullnames") %>%
    dplyr::select(-rowname)
  return(df)
}

# _______________________________________________________________________________

group_col_from_ids <- function(df.in, ids) {
  df.out <- dplyr::mutate(df.in, group = if_else(grepl("HC", ids), "HC",
    if_else(grepl("PC", ids), "PC", "PD")
  ))
  # rownames(df.out) <- rownames(df.in)
  return(df.out)
}

group_col_from_ids2 <- function(df.in) {
  df.out <-
    df.in %>%
    dplyr::mutate(group = if_else(grepl("HC", donor_id), "HC",
      if_else(grepl("PC", donor_id), "PC", "PD")
    ))

  rownames(df.out) <- rownames(df.in)
  return(df.out)
}


# _______________________________________________________________________________

boxplot_all <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel) {
  blank.title <- " "
  blank.ylabel <- " "
  group.cols <- c("PC" = "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
  ###' Basic all group boxplot function

  set.seed(123)
  ggplot(data = df, aes(x = x, y = y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x),
      position = position_jitterdodge(jitter.width = 0.75),
      shape = 21, size = 1.5, alpha = 0.8
    ) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel) +
    scale_color_manual(values = cols, name = "Group") +
    scale_fill_manual(values = cols, name = "Group") +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

boxplot_all_facet <-
  function(df,
           x,
           y,
           cols = c("PC" = "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5"),
           title = blank.title,
           ylabel = blank.ylabel) {
    blank.title <- " "
    blank.ylabel <- " "
    # group.cols = c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5")
    ###' Basic all group boxplot function with a cohort facet wrap

    set.seed(123)
    plot <-
      ggplot(data = df, aes(x = x, y = y)) +
      geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
      geom_point(aes(fill = x),
        position = position_jitterdodge(jitter.width = 1),
        shape = 21, size = 1.5, alpha = 0.8
      ) +
      theme_classic() +
      ggtitle(title) +
      labs(y = ylabel) +
      facet_wrap(~cohort) +
      scale_color_manual(values = cols, name = NULL) +
      scale_fill_manual(values = cols, name = NULL) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        panel.grid.major.y = element_blank()
      )
  }


# _______________________________________________________________________________

boxplot_all_generic <- function(df, x, y, cols = group.cols, title = blank.title, ylabel = blank.ylabel) {
  blank.title <- " "
  blank.ylabel <- " "
  ###' Basic all group boxplot function

  set.seed(123)
  ggplot(data = df, aes(x = x, y = y)) +
    geom_boxplot(aes(color = x), outlier.alpha = 0, width = 0.9) +
    geom_point(aes(fill = x),
      position = position_jitterdodge(jitter.width = 1),
      shape = 21, size = 1.5, alpha = 0.8
    ) +
    theme_classic() +
    ggtitle(title) +
    labs(y = ylabel, color = "Group") +
    guides(fill = FALSE) +
    scale_color_d3() +
    scale_fill_d3() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      panel.grid.major.y = element_blank()
    )
}

# _______________________________________________________________________________

alpha_div_boxplots <- function(df, x, y,
                               df.pairs, df.pairs.x, df.pairs.y, pairs.column,
                               cols, cols.rim, ylabel, PDvPC.stat, PDvHC.stat) {

  ###' Function for alpha diveristy
  ###' boxplots with paired lines

  set.seed(123)
  p <- ggplot(data = df, aes(x = x, y = y)) +
    theme_minimal() +
    geom_point(aes(fill = x, color = x), position = position_jitterdodge(dodge.width = 1), shape = 21, size = 1, alpha = 1) +
    geom_boxplot(aes(fill = x, color = x), width = 0.3, alpha = 0.1, outlier.alpha = 0) +
    geom_line(
      data = df.pairs, aes(x = df.pairs.x, y = df.pairs.y, group = pairs.column),
      linetype = "solid", color = "grey", alpha = 0.7
    ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    ylab(ylabel) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
    ) +
    labs(fill = "Group") +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols.rim) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(PDvHC.stat)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(PDvPC.stat))
  return(p)
}


# _______________________________________________________________________________
#                           Functions to explore select features in dataset
# _______________________________________________________________________________

explore_table <- function(obj) {
  data.table <- obj %>%
    microbiome::abundances() %>%
    as.data.frame() %>%
    rownames_to_column()
  return(data.table)
}

plot_feature <- function(obj, feature) {
  d <- obj %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    filter(rowname == feature) %>%
    melt()

  dm <- group_col_from_ids(d, id = d$variable)

  dm$value <- asin(sqrt(dm$value))
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))

  boxplot_all(dm,
    x = dm$group, y = dm$value,
    cols = c(
      "PC" = "#bfbfbf",
      "PD" = "#ed7d31",
      "HC" = "#5b9bd5"
    ),
    title = " ",
    ylabel = paste(unique(dm$rowname), "Abundance")
  )
}

# _______________________________________________________________________________
#                                            ML FUNCTIONS
# _______________________________________________________________________________
# from: http://jaehyeon-kim.github.io/2015/05/Setup-Random-Seeds-on-Caret-Package.html

setSeeds <- function(method = "repeatedcv", numbers = 1, repeats = 1, tunes = NULL, seed = 42) {

  #' function to set up random seeds for repeated ML cross validation models

  # B is the number of resamples and integer vector of M (numbers + tune length if any)

  B <- if (method == "cv") {
    numbers
  } else if (method == "repeatedcv") {
    numbers * repeats
  } else {
    NULL
  }

  if (is.null(length)) {
    seeds <- NULL
  } else {
    set.seed(seed = seed)
    seeds <- vector(mode = "list", length = B)
    seeds <- lapply(seeds, function(x) sample.int(n = 1000000, size = numbers + ifelse(is.null(tunes), 0, tunes)))
    seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
  }
  # return seeds
  seeds
}

# _______________________________________________________________________________

unregister <- function() {
  #' Unregister a foreach backend
  #' Use when looping functions that run
  #' analyses in parallel

  enviorn <- foreach:::.foreachGlobals
  rm(list = ls(name = enviorn), pos = enviorn)
}


# _______________________________________________________________________________

group_col_from_ids_ML <- function(df.in, ids) {
  df.out <- mutate(df.in, group = if_else(grepl("control", ids), "control", "disease"))
  rownames(df.out) <- rownames(df.in)
  return(df.out)
}

# _______________________________________________________________________________


prep.CMD.Species.ML <- function(study, metafilter = NA) {

  #' Funtion that preps input for ML analysis
  #' Input: Study name from curatedMetagenomicData
  #' Output: dataframe of Species abundnace inluding an
  #' identifying group column

  # alt.disease <- curatedMetagenomicData(paste0(study, ".metaphlan_bugs_list.stool"), dryrun=F)
  #
  # df.spec <- alt.disease[[1]] %>%
  #   ExpressionSet2phyloseq() %>%
  #   subset_taxa(!is.na(Species)) %>%
  #   subset_taxa(is.na(Strain))

  # study <- study

  df.abund <- study %>%
    # microbiome::transform("compositional") %>%
    microbiome::abundances() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column()
  # df.abund[-1] <- asin(sqrt(df.abund[-1]))

  m <- microbiome::meta(study) %>%
    select(study_condition) %>%
    rownames_to_column()

  model.input <- left_join(df.abund, m) %>%
    column_to_rownames()

  if (!is.na(metafilter)) {
    model.input <- filter(model.input, study_condition != metafilter)
  }

  count(model.input$study_condition)
  unique(model.input$study_condition)

  model.input <- model.input %>%
    group_col_from_ids_ML(ids = model.input$study_condition) %>%
    dplyr::select(-study_condition)

  model.input$group <- factor(model.input$group)

  return(model.input)
}

# _______________________________________________________________________________

# helper function for vizualizing Xgboost plots
# https://www.kaggle.com/pelkoja/visual-xgboost-tuning-with-caret/report
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    theme_bw()
}

# _______________________________________________________________________________


# Plot feature of interest (foi)
plot_foi <- function(datObj, foi) {
  asin.input <- datObj %>%
    microbiome::transform("compositional") %>%
    microbiome::abundances()
  df1 <- asin(sqrt(asin.input)) %>%
    as.data.frame() %>%
    rownames_to_column()
  d <- df1 %>%
    filter(rowname == foi) %>%
    melt()
  dm <- group_col_from_ids(d, id = d$variable)
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
  boxplot_all(dm,
    x = dm$group, y = dm$value,
    title = " ",
    ylabel = paste(unique(dm$rowname), "Abundance")
  )
}

# _______________________________________________________________________________

# Plot feature of interest (foi) for any Dataset from CMG
plot_foi_general <- function(datObj, foi, title = "") {
  asin.input <- datObj %>%
    microbiome::transform("compositional") %>%
    microbiome::abundances()
  df1 <- asin(sqrt(asin.input)) %>%
    as.data.frame() %>%
    rownames_to_column()
  met <- meta(df.spec) %>%
    dplyr::select(c(study_condition)) %>%
    rownames_to_column(var = "variable")
  d <- df1 %>%
    filter(rowname == foi) %>%
    melt()
  d2 <- left_join(d, met)
  d2$study_condition <- factor(d2$study_condition)
  boxplot_all_generic(d2,
    x = d2$study_condition, y = d2$value,
    title = "",
    ylabel = paste(unique(d2$rowname), "Abundance")
  )
}

# _______________________________________________________________________________


# 2) CLR Log-odds ratio : facultative anaerobic / obligate anaerobic bacteria

#' Function returns normalized values for feature(s) of interest
#' 1) Adaptive to multiple types of normalization
#' 2) Will sum values if a list of foi's are input

fois <- function(datObj, foi) {
  if (length(foi) > 1) {

    # Init
    cnt <- 1

    for (i in foi) {
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
      if (cnt == 1) {
        summed.vars <- d
      } else {
        summed.vars$value <- summed.vars$value + d$value
      }
      cnt <- cnt + 1
    }
    return(summed.vars)
  } else {
    for (i in foi) {
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
    }
    return(d)
  }
}

# _______________________________________________________________________________
# Negate %in% function

"%ni%" <- Negate("%in%")

# _______________________________________________________________________________

# SGV Stats Function

nmean_summary <- function(df) {
  output.stats <- tibble()
  vars <- colnames(df)[1:ncol(df) - 1]
  for (sv in vars) {
    print(sv)
    stat.col <- df %>%
      dplyr::select(group, sv) %>%
      filter(!is.na(df[[sv]])) %>%
      group_by(group) %>%
      dplyr::summarise(n = n(), across(
        where(is.numeric),
        ~ mean(.x, na.rm = TRUE)
      ))
    colnames(stat.col)[3] <- "ratio_detected"

    row2add <- cbind(stat.col[1:3], sv)
    output.stats <- rbind(output.stats, row2add)
  }
  return(output.stats)
}


# _______________________________________________________________________________

# Functions to clean up column names

clean.cols.tax <- function(x) {
  colnames(x) <- gsub("_taxonomic_profile", "", colnames(x))
  x
}

clean.cols.abund <- function(x) {
  colnames(x) <- gsub("_Abundance", "", colnames(x))
  x
}

clean.cols.abund_RPK <- function(x) {
  colnames(x) <- gsub("_Abundance-RPKs", "", colnames(x))
  x
}

clean.cols.abund_CPM <- function(x) {
  colnames(x) <- gsub("_Abundance-CPM", "", colnames(x))
  x
}


trim_cols <- function(x, cohort) {
  if (cohort == "TBC") {
    colnames(x) <- substr(colnames(x), 1, 10)
    x
  } else if (cohort == "RUSH") {
    colnames(x) <- substr(colnames(x), 1, 6)
    x
  } else if (cohort == "Bonn") {
    colnames(x) <- str_replace(colnames(x), "_1", "")
    x
  }
}


# _______________________________________________________________________________

tss <- function(x) {
  result <- (x / sum(x))
}
hund2relab <- function(x) {
  result <- (x / 100)
}
cpm2relab <- function(x) {
  result <- (x / 1000000)
}
# _______________________________________________________________________________
pseudoCounts_bonn <- function(df, reads) {
  psudocnts <- df
  for (i in colnames(psudocnts)) {
    donor_reads <- reads[[which(reads$run == i), 2]]
    psudocnts[i] <- psudocnts[i] * donor_reads
    # print(donor_reads)
    # print(psudocnts[i])
  }
  cat("Pseudocount Transformation Complete\n")
  return(psudocnts)
}

pseudoCounts <- function(obj) {

  # TROUBLE
  # df <- dat.species %>% abundances() %>%
  #   as.data.frame()
  # reads <- dat.species %>% meta() %>%
  #   dplyr::select(donor_id, total_reads) %>%
  #   as.data.frame()

  df <- obj %>%
    abundances() %>%
    as.data.frame()
  reads <- obj %>%
    meta() %>%
    dplyr::select(donor_id, total_reads) %>%
    as.data.frame()

  psudocnts <- df
  for (i in colnames(psudocnts)) {
    donor_reads <- reads[[which(reads$donor_id == i), "total_reads"]] %>% as.numeric()
    psudocnts[i] <- psudocnts[i] * donor_reads
  }
  cat("Pseudocount Transformation Complete\n")
  return(psudocnts)
}


impute_phyloseq <- function(phy) {
  #' Function replaces all zero values with 
  #' 1/2 the smallest value in a sample
  
  impt <- phy %>%
    abundances() %>%
    as.data.frame() %>%
    mutate_all( ~ replace(., . == 0, min(.[. > 0], na.rm = TRUE) / 2))
  otu_table(phy) <- otu_table(impt, taxa_are_rows = TRUE)
  
  return(phy)
  
}

# _______________________________________________________________________________

calculate_auroc <- function(case_df, control_df) {
  #' This function takes in two data-frames with identical rownames and 
  #' calcuates the AUROC for each feature based on abundance values
  require(foreach)
  require(doParallel)
  aucs <- tibble()
  start_time <- Sys.time()
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1)
  registerDoParallel(cl)
  
  roc2add <-
    foreach(
      feat = rownames(case_df),
      .combine = "rbind",
      .packages = c("magrittr", "pROC")
    ) %dopar% {
      x <- case_df[feat,] %>% t()
      y <- control_df[feat,] %>% t()
      # AUROC calculation
      rocdata <-
        c(roc(
          controls = y,
          cases = x,
          direction = "<",
          ci = TRUE,
          auc = TRUE
        )$ci)
      data.frame(
        "feature" = feat,
        "ci_lower" = rocdata[1],
        "auroc" = rocdata[2],
        "ci_upper" = rocdata[3]
      )
    }
  stopCluster(cl)
  aucs <- rbind(aucs, roc2add)
  end_time <- Sys.time()
  cat(
    "AUROCs calculated in : ",
    end_time - start_time,
    attr(end_time - start_time, "units"),
    "\n"
  )
  return(aucs)
}

# _______________________________________________________________________________

bonn_metadata <- function(dat) {
  meta_stats <- meta(dat) %>% select(contains("total_"))
  meta_df <- cbind(meta_stats, metadata_Bonn.slim)
  sample_data(dat) <- sample_data(meta_df)
  return(dat)
}

# # _______________________________________________________________________________

maaslin_prep <- function(dat) {
  dat <- dat %>%
    subset_samples(donor_id %ni% low_qc[[1]]) %>%
    core(detection = 0, prevalence = 0.1)
  return(dat)
}

# _______________________________________________________________________________-

binarize <- function(x) {
  ifelse(x == 0, 0, 1)
}
# _______________________________________________________________________________-

color_loop_generator <- function(names_col) {
  # Init
  accent_pal <-
    c(
      "#7fc97f", "#beaed4", "#fdc086", "#ffff99",
      "#386cb0", "#f0027f", "#bf5b17", "#666666"
    )
  color_output <- c()
  pal_length <- length(accent_pal)
  feats <- unique(names_col)
  n_cols <- length(feats)
  color_output <- c(
    rep(accent_pal, floor(n_cols / pal_length)),
    accent_pal[1:(n_cols %% pal_length)]
  )
  names(color_output) <- feats

  return(color_output)
}

# _______________________________________________________________________________
#####                      Correlation Functions                           #####
# _______________________________________________________________________________

corr_abund_prep <- function(obj, obj.name, cohort, sigfilter = T) {

  ### Read-in MaAsLin2 output
  Maas.pd.pc.sig <- read_tsv(
    paste0(
      "data/MaAsLin2_Analysis/", cohort, "/", obj.name,
      "_PDvPC_maaslin2_output/all_results.tsv"
    ),
    col_names = T
  ) %>%
    filter(value == "Population Control") %>%
    filter(qval < 0.25)

  Maas.pd.hc.sig <- read_tsv(
    paste0(
      "data/MaAsLin2_Analysis/", cohort, "/", obj.name,
      "_PDvHC_maaslin2_output/all_results.tsv"
    ),
    col_names = T
  ) %>%
    filter(value == "Household Control") %>%
    filter(qval < 0.25)

  features <- full_join(Maas.pd.pc.sig, Maas.pd.hc.sig, by = "feature") %>%
    dplyr::select("feature")

  if (sigfilter) {
    abundance_tbl <-
      obj %>%
      subset_samples(donor_id %ni% low_qc[[1]]) %>%
      core(detection = 0, prevalence = 0.1) %>%
      abundances() %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      filter(rowname %in% features[[1]]) %>%
      column_to_rownames() %>%
      t() %>%
      as.data.frame()
  } else {
    abundance_tbl <-
      obj %>%
      subset_samples(donor_id %ni% low_qc[[1]]) %>%
      core(detection = 0, prevalence = 0.1) %>%
      abundances() %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
  }
  return(abundance_tbl)
}



corr_loop_parallel <- function(metadata, abundance, obj.name) {
  require(foreach)
  require(doParallel)

  # setup parallel back-end to use multiple threads
  start_time <- Sys.time()
  cores <- detectCores()
  cl <- makeCluster(cores[1] - 1) # not to overload your computer
  registerDoParallel(cl)


  corr_output <- tibble()
  looped_df <-
    foreach(metavar = colnames(metadata), .combine = "rbind") %:%
    foreach(feature = colnames(abundance), .combine = "rbind") %dopar% {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )

      data.frame(
        "metadata" = metavar,
        "feature" = feature,
        "object_name" = obj.name,
        "rho" = spearman$estimate[[1]],
        "S" = spearman$statistic[[1]],
        "n" = length(na.omit(metadata[[metavar]])),
        "p" = spearman$p.value[[1]]
      )
    }
  stopCluster(cl)
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    looped_df %>%
    na.omit() %>%
    mutate(
      across(all_of(statvars), as.character),
      across(all_of(statvars), as.numeric)
    ) %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    dplyr::rename("feature" = "fullnames") %>%
    dplyr::relocate(feature, .after = metadata) %>%
    group_by(object_name, metadata) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup()
  end_time <- Sys.time()
  cat(
    "Correlations calculated in : ",
    end_time - start_time, attr(end_time - start_time, "units"), "\n"
  )
  return(corr_output)
}

corr_loop <- function(metadata, abundance, obj.name) {
  corr_output <- tibble()
  for (metavar in colnames(metadata)) {
    cat("Calculating correlations for: ", metavar[[1]], "\n")
    for (feature in colnames(abundance)) {
      # Calculate Spearman's Correlation
      spearman <-
        cor.test(
          x = metadata[[metavar]],
          y = abundance[[feature]],
          method = "spearman",
          na.action = na.exclude,
          alternative = "two.sided"
        )
      row2add <-
        cbind(
          "metadata" = metavar,
          "feature" = feature,
          "object_name" = obj.name,
          "rho" = spearman$estimate[[1]],
          "S" = spearman$statistic[[1]],
          "n" = length(na.omit(metadata[[metavar]])),
          "p" = spearman$p.value[[1]]
        )
      corr_output <- rbind(corr_output, row2add)
    }
  }
  # Remove NAs and add FDR (Benjamini Hochberg)
  statvars <- c("rho", "S", "n", "p")
  corr_output <-
    corr_output %>%
    na.omit() %>%
    mutate(
      across(all_of(statvars), as.character),
      across(all_of(statvars), as.numeric)
    ) %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    dplyr::rename("feature" = "fullnames") %>%
    dplyr::relocate(feature, .after = metadata) %>%
    group_by(object_name, metadata) %>%
    mutate(q = p.adjust(p, method = "BH")) %>%
    ungroup()
  return(corr_output)
}


# _______________________________________________________________________________
#####                        Correlation XY Plot                          #####
# _______________________________________________________________________________


corr_xy <- function(obj, corr_obj, feature_var, metadata_var) {

  #' Function creates a scatter plot of a given feature and a metadata column
  abund <- obj %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column(var = "feature") %>%
    decode_rfriendly_rows(passed_column = "feature") %>%
    dplyr::select(-feature) %>%
    column_to_rownames(var = "fullnames") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "donor_id")
  df.plot <- obj %>%
    meta() %>%
    process_meta(cohort = "Merged") %>%
    left_join(abund, by = "donor_id")

  stat_col <-
    corr_obj %>%
    dplyr::filter(feature == sym(feature_var)) %>%
    dplyr::filter(metadata == sym(metadata_var))
  stat_title <-
    paste0(
      "Spearman's Rho: ", round(stat_col$rho, digits = 3), "\n",
      "P-value: ", format(stat_col$p, digits = 3, scientific = T),
      "  FDR: ", format(stat_col$q, digits = 3, scientific = T)
    )
  cat("Rho: ", stat_col$rho, ", ")
  cat("P-value: ", stat_col$p, ",  ")
  cat("Q-value: ", stat_col$q, "\n")

  df.plot %>%
    drop_na(metadata_var) %>%
    ggplot(aes(x = .data[[feature_var]], y = .data[[metadata_var]])) +
    geom_point(aes(fill = donor_group, color = donor_group), shape = 21, alpha = 1) +
    geom_smooth(method = lm, color = "darkgrey", linetype = "dotted", se = F) +
    theme_bw() +
    labs(x = feature_var, y = metadata_var, title = stat_title) +
    scale_fill_manual(values = cols.pdpchc) +
    scale_color_manual(values = cols.pdpchc.rim) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(size = 12)
    )
}

# _______________________________________________________________________________
#####                      Correlation Heatmap                           #####
# _______________________________________________________________________________

corr_heatmap <- function(corr.df) {

  # Calculate number of significant (P-VALUE < 0.05) associations
  #  and select top 30 features
  corr.df.top <- corr.df %>%
    dplyr::filter(q < 0.25) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(n = n()) %>%
    arrange(desc(n)) %>%
    top_n(n = 30, wt = n) %>%
    slice_head(n = 30)

  # Filter correlation df for top 30 features
  corr.df.trim <- corr.df %>%
    mutate(metadata = as.character(metadata)) %>%
    filter(feature %in% corr.df.top$feature)

  # create dataframe matrix of Rho correlation values for distance functions
  rho.df <-
    corr.df.trim %>%
    dplyr::select(feature, metadata, rho) %>%
    pivot_wider(names_from = feature, values_from = rho, values_fill = NA) %>%
    column_to_rownames(var = "metadata")

  # Hierarchical clustering of Rho values for features & metadata
  meta.dendro <-
    as.dendrogram(hclust(d = dist(x = rho.df), method = "complete"))
  meta.dendro.plot <- ggdendrogram(data = meta.dendro, rotate = TRUE)
  feature.dendro <-
    as.dendrogram(hclust(d = dist(x = as.data.frame(t(rho.df))), method = "complete"))
  feature.dendro.plot <- ggdendrogram(data = feature.dendro, rotate = TRUE)

  ### Reorder Heatmap axis using order of dendrograms
  feature.order <- order.dendrogram(feature.dendro)
  metadata.order <- order.dendrogram(meta.dendro)
  corr.df.trim.ordered <- corr.df.trim %>%
    dplyr::mutate(feature = factor(feature,
      ordered = TRUE,
      levels = unique(corr.df.trim$feature)[feature.order]
    )) %>%
    dplyr::mutate(metadata = factor(metadata,
      ordered = TRUE,
      levels = unique(corr.df.trim$metadata)[metadata.order]
    ))

  ### Plot Heatmap
  h1 <-
    corr.df.trim.ordered %>%
    mutate(siglabel = if_else(q < 0.20, "*", "")) %>%
    ggplot(aes(x = metadata, y = feature, fill = rho)) +
    geom_tile() +
    geom_text(aes(label = siglabel), size = 5, vjust = 0.77, color = "white") +
    labs(fill = "Spearman correlation") +
    scale_y_discrete(position = "right") +
    scale_fill_distiller(palette = "RdBu") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      plot.margin = unit(c(1, 1, 1, 5), "cm")
    )

  print(h1)
  return(h1)
}
# _______________________________________________________________________________
#####                Correlation Plot helper functions                       #####
# _______________________________________________________________________________

trim_sig_helper <- function(df_cors, n = 10) {
  output <-
    df_cors %>%
    filter(q < 0.25) %>%
    slice_min(n = n, order_by = q)
  return(output)
}

top_n_scatterplots <- function(dat, obj.name,
                               df.cors, metaclass,
                               n_plots = 10) {
  df.cors.sig <- trim_sig_helper(df.cors, n = n_plots)

  for (corr in 1:nrow(df.cors.sig)) {
    cor_row <- df.cors.sig[corr, ]
    p <- corr_xy(
      obj = dat, df.cors,
      feature_var = cor_row$feature[[1]],
      metadata_var = as.character(cor_row$metadata[[1]])
    )
    # print(p)
    obj.name.out <- gsub("/", "_", obj.name)
    obj.name.out <- gsub(",", "", obj.name.out)
    plot.name <- paste0(
      "data/Correlations/", obj.name.out, "/", metaclass, "/",
      cor_row$feature[[1]], "_VS_",
      as.character(cor_row$metadata[[1]]), ".svg"
    )
    ggsave(p, filename = plot.name, height = 4, width = 5)
  }
}
# _______________________________________________________________________________


prevalence_filter <- function(datObj, threshold = 0.1) {

  # function conducts a 10% prevalence filter on each individual cohort and joins
  # rows that are shared in all cohorts

  df_Shanghai <- subset_samples(datObj, cohort == "Shanghai") %>%
    core(detection = 0, prevalence = threshold) %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column()
  df_Rush <- subset_samples(datObj, cohort == "Rush") %>%
    core(detection = 0, prevalence = threshold) %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column()
  df_TBC <- subset_samples(datObj, cohort == "TBC") %>%
    core(detection = 0, prevalence = threshold) %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column()
  df_Bonn <- subset_samples(datObj, cohort == "Bonn") %>%
    core(detection = 0, prevalence = threshold) %>%
    abundances() %>%
    as.data.frame() %>%
    rownames_to_column()

  df_abund <-
    inner_join(df_Shanghai, df_Rush, by = "rowname") %>%
    inner_join(df_TBC, by = "rowname") %>%
    inner_join(df_Bonn, by = "rowname") %>%
    column_to_rownames()

  my_dat_table <- otu_table(df_abund, taxa_are_rows = T)
  my_sample_data <- meta(datObj) %>% sample_data()
  dat <- phyloseq(my_dat_table, my_sample_data)
  # print(dat)
  return(dat)
}

# _______________________________________________________________________________

enrichment_formula <- function(N, n, m, k) {

  #' N : the number of total features in both the study and pathway library
  #' n : the number of significant features in both the study and pathway library
  #' m : the number of detected features in the pathway
  #' k : number of significant features in the pathway

  cat("number of total features: ", N, "\n")
  cat("number of total significant features: ", n, "\n")
  cat("number of detected features in pathway: ", m, "\n")
  cat("number of significant features in the pathway: ", k, "\n")

  i <- seq(from = 0, to = k, by = 1)
  binomal_coefs <- gmp::asNumeric(sum(div.bigz((chooseZ(m, i) * chooseZ((N - m), (n - i))), chooseZ(N, n))))
  output <- -log10(1 - binomal_coefs)
  return(output)
}
# _______________________________________________________________________________
