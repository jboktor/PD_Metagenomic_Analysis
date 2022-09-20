#

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
# load("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.RData")
# base::load("files/Phyloseq_Merged_ML.RData") #phyloseq_objs

# level <- "eggNOGs.slim"
### Read-in MaAsLin2 output

# base::load("data/Machine_Learning_Analysis/feature_AUROCs/feature_AUROCs.RData")
phyloseq_objs <- readRDS("files/Phyloseq_Merged_ML_clean.rds")
aucs <- readRDS("data/Machine_Learning_Analysis/feature_AUROCs/feature_slim_AUROCs.rds")


#_______________________________________________________________________________
# Summary Statistics ----

obj_list <- names(phyloseq_objs)[-3]
maaslin_df_sig <- tibble()
for (level in obj_list){

  maaslin_df <- read_tsv(paste0("data/MaAsLin2_Analysis/Meta/", level, 
                                "_maaslin2_output/all_results.tsv"), col_names = T)
  
  sig_df <-  maaslin_df %>% 
    filter(qval <= 0.1) %>% 
    decode_rfriendly_rows(passed_column = "feature") %>% 
    dplyr::select(-feature) %>%
    dplyr::rename("feature"="fullnames") %>% 
    slice_min(qval, n = 25, with_ties = F) %>% 
    mutate(data_level = level)
  
  cat(level, " : ", maaslin_df %>% filter(qval <= 0.1) %>% nrow(), "\n")
  maaslin_df_sig <- bind_rows(maaslin_df_sig, sig_df)
  
}

df_plot <- maaslin_df_sig %>% 
  left_join(aucs, by = c("feature", "data_level")) %>% 
  drop_na(auroc)

meta_plot_main <- 
  df_plot %>% 
  filter(data_level %in% c("Species","Pathways.slim", "GOs.slim", "KOs.slim")) %>% 
  mutate(data_level = factor(data_level, levels = c("Species","Pathways.slim", "GOs.slim", "KOs.slim"))) %>% 
  ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study, 
             ymin=ci_lower, ymax=ci_upper, group = feature)) +
  geom_pointrange(aes(group = feature), 
                  position = position_jitterdodge(jitter.height = 0),
                  shape = 24, stroke = 0.1, colour="grey") +
  # facet_wrap(~data_level, scales = "free", ncol = 1) +
  facet_grid(rows = vars(data_level), scales = "free", space = "free") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
  labs(y = "AUROC", x = NULL) +
  my_clean_theme2() +
  scale_fill_d3() +
  coord_flip()
meta_plot_main

ggsave(meta_plot_main,
       filename = paste0("data/Meta_analysis/AUROC_CI_main.svg"),
       width = 9, height = 12)


meta_plot_supp <- 
  df_plot %>% 
  filter(data_level %nin% c("Species","Pathways.slim", "GOs.slim", "KOs.slim")) %>% 
  mutate(data_level = factor(data_level, levels = c("Genus", "Enzymes.slim", "Pfams.slim", "eggNOGs.slim"))) %>% 
  ggplot(aes(x=fct_reorder(feature, auroc), y = auroc, fill = study, 
             ymin=ci_lower, ymax=ci_upper, group = feature)) +
  geom_pointrange(aes(group = feature), 
                  position = position_jitterdodge(jitter.height = 0),
                  shape = 24, stroke = 0.1, colour="grey") +
  # facet_wrap(~data_level, scales = "free", ncol = 1) +
  facet_grid(rows = vars(data_level), scales = "free", space = "free") +
  geom_hline(yintercept = 0.5, linetype = 2, alpha = 0.8) +
  labs(y = "AUROC", x = NULL) +
  my_clean_theme2() +
  scale_fill_d3() +
  coord_flip()
meta_plot_supp

ggsave(meta_plot_supp,
       filename = paste0("data/Meta_analysis/AUROC_CI_supplement.svg"),
       width = 9, height = 14)




aucs_wide <- aucs %>%
  pivot_wider(names_from = study, values_from = c(ci_lower, auroc , ci_upper))
df_saveme <- maaslin_df_sig %>% 
  left_join(aucs_wide, by = c("feature", "data_level"))
write.csv(df_saveme, file = paste0('files/Meta_analysis/maaslin_auroc_merged.csv'))
openxlsx::write.xlsx(df_saveme, file = 'files/Supplementary Tables/Table_S8_Meta-Analyses.xlsx', overwrite = T)





#______________________________________________________________________________
# List of species of interest for Ana M

level <- "Species"
maaslin_df <- read_tsv(paste0("data/MaAsLin2_Analysis/Meta/", level, 
                              "_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(qval <= 0.1) %>% 
  decode_rfriendly_rows(passed_column = "feature") %>% 
  dplyr::select(-feature) %>%
  dplyr::rename("feature"="fullnames") %>% 
  mutate(data_level = level)

bug.list.4.AnaM <- maaslin_df %>% 
  left_join(aucs, by = c("feature", "data_level")) %>% 
  drop_na(auroc) %>% 
  pivot_wider(names_from = study, values_from = c(ci_lower, auroc , ci_upper))
openxlsx::write.xlsx(bug.list.4.AnaM, file = 'files/Meta-Analysis_significant_species.xlsx', overwrite = T)

