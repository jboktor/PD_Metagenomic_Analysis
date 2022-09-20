# Interactive plot 4 John

######## Load Data & functions
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
load("files/low_quality_samples.RData")
library(plotly)

remove_dats()
load("files/Phyloseq_Merged/GOs_PhyloseqObj.RData")


datObj <- phyloseq_objs[["Enzymes.slim"]] %>% 
  subset_samples(donor_id %ni% low_qc[[1]])


feature <- "1.8.4.12"
df.barplot <-
  datObj %>%
  microbiome::transform("compositional") %>%
  microbiome::abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  decode_rfriendly_rows("rowname") %>% 
  dplyr::select(-rowname) %>% 
  dplyr::rename(rowname = fullnames) %>% 
  filter(grepl(feature, rowname, fixed = TRUE)) %>%
  filter(grepl("\\|", rowname)) %>%
  column_to_rownames() %>%
  t() %>%
  melt()
df.barplot
# all_meta <- process_meta(datObj, cohort = "Merged") %>% 
#   dplyr::select(donor_id, donor_group, cohort)

df.barplot <- datObj %>%
  microbiome::transform("compositional") %>%
  microbiome::abundances() %>%
  as.data.frame() %>%
  rownames_to_column("tempvars") %>%
  decode_rfriendly_rows(passed_column = "tempvars") %>%
  dplyr::select(-tempvars) %>%
  dplyr::rename(features = fullnames) %>%
  dplyr::filter(grepl(feature, features, fixed = TRUE)) %>%
  dplyr::filter(grepl("\\|", features)) %>%
  column_to_rownames(var = "features") %>%
  t() %>% melt()
df.barplot <- df.barplot %>%
  group_col_from_ids(ids = df.barplot$Var1) %>%
  dplyr::mutate(group = factor(group, levels = c("PC", "PD", "HC")))

stratplot.df <- df.barplot %>%
  mutate(Var2 =  sub(".*\\|", "", Var2)) %>%
  group_by(Var2) %>% 
  mutate(mean_size = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  dplyr::mutate(Var2 = fct_reorder(Var2, mean_size)) 
stratplot <- 
  stratplot.df %>% 
  ggplot(aes(
    x = reorder(Var1,-value),
    y = value,
    fill = Var2
  )) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Donor", y = "Relative Abundance", fill = NULL) +
  facet_wrap( ~ group, scales = "free_x", 
              labeller = labeller(group = c(
                "PC" = "Population Controls", 
                "PD" ="Parkinson's Disease",  
                "HC" = "Household Controls"))) +
  scale_fill_manual(values = color_loop_generator(levels(stratplot.df$Var2)), 
                    guide = guide_legend(reverse = TRUE)) +
  theme(
    axis.text.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background =element_rect(fill="white"), 
    plot.margin = unit(c(1, 1, 1, 3), "cm")
  )

# Save Interactive plot
stratPlot_interactive <-
  ggplotly(stratplot, tooltip = c("y", "x", "fill"))
# saveWidget(as_widget(stratPlot), file, selfcontained = TRUE)

p <- ggplotly(stratPlot_interactive, tooltip = c("y", "x", "fill"))
htmlwidgets::saveWidget(as_widget(p), selfcontained = TRUE,
                        paste0("data/Stacked_Barplots/", feature, "stackedbarplot.html"))

# Save Legend
stratplot_legend <- cowplot::plot_grid(get_legend(stratplot))
ggsave(stratplot_legend, filename = paste0("data/Stacked_Barplots/", feature, "stackedbarplot_legend.svg"),
       width = 30, height = 7)
# Save Barplot
stratplot_out <- stratplot + theme(legend.position = "none")
ggsave(stratplot_out, filename = paste0("data/Stacked_Barplots/", feature, "stackedbarplot.svg"),
       width = 12, height = 5)






