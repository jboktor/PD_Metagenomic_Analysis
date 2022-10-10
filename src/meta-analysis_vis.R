
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
library(patchwork)
library(ggpackets)

#_______________________________________________________________________________
# Visualization utility functions
format_cohort_df <- function(mmuphin_results, cohort) {
  mmuphin_results$maaslin_fits[[cohort]] %>% 
    filter(feature %in% plot_df$feature) %>% 
    mutate(feature = factor(feature, levels = plot_df$feature))
}

ggpk_heatbar <- function(...) {
  ggpacket(...) %+%
    geom_tile() %+% 
    coord_flip() %+% 
    theme_void() %+% 
    theme(axis.text.y = element_blank())
}

analysis_levels <- c("Species","Pathways.slim", "Pfams.slim", "KOs.slim", "GOs.slim")
#_______________________________________________________________________________

MMUPHin_models <- readRDS("data/MMUPHin_Meta-Analysis/results_2022-10-07.rds")

# Explore results 
mmuphin_results_meta <- tibble()
for (level in analysis_levels) {
  mmuphin_results_meta %<>% bind_rows(
    MMUPHin_models[[glue("{level}_adj")]][["meta_fits"]] %>% 
      mutate(data_level = level)
  )
}

#_______________________________________________________________________________

# save supplementary table
meta_save <- list()

for (level in analysis_levels) {
  mmuphin_results <- MMUPHin_models[[glue("{level}_adj")]]
  meta_save[[glue("{level}_aggregated")]] <- mmuphin_results$meta_fits
  meta_save[[glue("{level}_TBC")]] <- mmuphin_results$maaslin_fits$TBC
  meta_save[[glue("{level}_RUMC")]] <- mmuphin_results$maaslin_fits$Rush
  meta_save[[glue("{level}_Bonn")]] <- mmuphin_results$maaslin_fits$Bonn
  meta_save[[glue("{level}_Shanghai")]] <-mmuphin_results$maaslin_fits$Shanghai
}

supp_loc <- "files/Supplementary Tables"
openxlsx::write.xlsx(meta_save,
                     file = glue("{supp_loc}/Table_S8 Meta-Anaylsis-statistics_{Sys.Date()}.xlsx"))

#_______________________________________________________________________________


fig <- list()
for (level in analysis_levels){
  mmuphin_results <- MMUPHin_models[[glue("{level}_adj")]]
  mod_df <- mmuphin_results$meta_fits
  
  # meta-summary estimates ----
  plot_df <- 
    mod_df %>%
    filter(qval.fdr <= 0.1) %>% 
    slice_min(order_by = qval.fdr, n = 25, with_ties = F) %>%
    # slice_max(order_by = abs(coef), n = 25, with_ties = F) %>%
    arrange(coef) %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    mutate(qval.fdr.bins = case_when(
      qval.fdr < 0.01 ~ "≤ 0.01",
      qval.fdr < 0.05 ~ "> 0.05 ≤ 0.1",
      qval.fdr < 0.1 ~ "> 0.1 ≤ 0.25"
    ))
  
  indiv_mod <- mod_df %>% bind_rows()
  indiv_mod_trim <- indiv_mod %>% filter(feature %in% plot_df$feature)
  
  fig_meta <- plot_df %>% 
    ggplot(aes(x = feature, y = coef)) +
    geom_pointrange(aes(ymin = coef - stderr, ymax = coef + stderr)) +
    labs(x = NULL, y = "Effect Size") +
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_blank())
  
  fig_meta_q <- plot_df %>%
    ggplot(aes(x = feature, y = 1, fill = qval.fdr.bins)) +
    ggpk_heatbar() +
    scale_fill_brewer(palette = "Greys", na.value="white") +
    labs(title = 'q', fill = 'q-value \n(FDR)') +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.y = element_text(hjust = 1))
  
  fig_meta_tau2 <- plot_df %>%
    ggplot(aes(x = feature, y = 1, fill = tau2)) +
    ggpk_heatbar() +
    scale_fill_distiller(palette = "Greys", na.value="white") +
    labs(title = expression(tau^{2}), fill = expression(tau^{2})) +
    theme(plot.title = element_text(hjust = 0.5))
  fig_meta_I2 <- plot_df %>%
    ggplot(aes(x = feature, y = 1, fill = I2)) +
    ggpk_heatbar() +
    scale_fill_distiller(palette = "Greys", na.value="white") +
    labs(title = expression(I^{2}), fill = expression(I^{2})) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # per-cohort estimates ----
  indp_cohorts_df <-
    bind_rows(
      format_cohort_df(mmuphin_results, "TBC"),
      format_cohort_df(mmuphin_results, "Rush"),
      format_cohort_df(mmuphin_results, "Bonn"),
      format_cohort_df(mmuphin_results, "Shanghai")
    )
  
  # colorscheme for fc heatmap 
  mx <- max(abs(range(indp_cohorts_df$coef, na.rm=TRUE)))
  mx_digits <- floor(-log10(mx)) + 1
  mx <- ifelse(round(mx, digits = mx_digits) < mx, 
               round(mx, digits = mx_digits) + 1*10^(-mx_digits - 1), 
               round(mx, digits = mx_digits))
  col.hm = c(rev(colorRampPalette(brewer.pal(9, 'Blues'))(9)),
             rep('#FFFFFF'),
             colorRampPalette(brewer.pal(9, 'Reds'))(9))
  
  coef.tbc <- format_cohort_df(mmuphin_results, "TBC") %>%
    ggplot(aes(x = feature, y = 1, fill =coef )) +
    ggpk_heatbar() +
    scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), na.value="white") +
    labs(title = 'TBC', fill = "Regression\nCoefficient") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, angle = 45))
  coef.rumc <- format_cohort_df(mmuphin_results, "Rush") %>%
    ggplot(aes(x = feature, y = 1, fill =coef )) +
    ggpk_heatbar() +
    scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), na.value="white") +
    labs(title = 'Rush', fill = "Regression\nCoefficient") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, angle = 45))
  coef.bonn <- format_cohort_df(mmuphin_results, "Bonn") %>%
    ggplot(aes(x = feature, y = 1, fill =coef )) +
    ggpk_heatbar() +
    scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), na.value="white") +
    labs(title = 'Bonn', fill = "Regression\nCoefficient") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, angle = 45))
  coef.shanghai <- format_cohort_df(mmuphin_results, "Shanghai") %>%
    ggplot(aes(x = feature, y = 1, fill =coef )) +
    ggpk_heatbar() +
    scale_fill_gradientn(colours=col.hm, limits=c(-mx, mx), na.value="white") +
    labs(title = 'Shanghai', fill = "Regression\nCoefficient") +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, angle = 45))
  
  fig[[level]] <- fig_meta_q + fig_meta_tau2 + fig_meta_I2 + fig_meta + 
    coef.tbc + coef.rumc + coef.bonn + coef.shanghai + plot_spacer() +
    plot_layout(guides = "collect", widths = c(1, 1, 1, 8, 1, 1, 1, 1, 2))

  ggsave(fig[[level]], filename = glue("data/MMUPHin_Meta-Analysis/figures/{level}_adj.panel.svg"),
         width = 11, height = 6)
}

ggsave(fig$Pfams.slim, filename = glue("data/MMUPHin_Meta-Analysis/figures/Pfams.slim_adj.panel.svg"),
       width = 11, height = 6)
ggsave(fig$Pathways.slim, filename = glue("data/MMUPHin_Meta-Analysis/figures/Pathways.slim_adj.panel.svg"),
       width = 13, height = 4)
ggsave(fig$GOs.slim, filename = glue("data/MMUPHin_Meta-Analysis/figures/GOs.slim_adj.panel.svg"),
       width = 15, height = 6)
ggsave(fig$KOs.slim, filename = glue("data/MMUPHin_Meta-Analysis/figures/KOs.slim_adj.panel.svg"),
       width = 13, height = 6)

