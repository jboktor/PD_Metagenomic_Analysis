# correlation vis

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

load('files/Correlations/Dietary_Correlations.RData') # dietary_corrs
load('files/Correlations/Clinical_Correlations.RData') # clinical_corrs

df_dietary_corrs <- do.call(rbind.data.frame, dietary_corrs)
df_clinical_corrs <- do.call(rbind.data.frame, clinical_corrs)

# number of significant clincial correlations
df_clinical_corrs %>% filter(q<=0.1) %>% nrow()
# number of significant dietary correlations
df_dietary_corrs %>% filter(q<=0.1) %>% nrow()


df_clinical_corrs %>% 
  filter(rho < -0.5) %>%
  ggplot(aes(x=rho, y=-log10(q))) +
  geom_point(
    aes( color = metadata),
    alpha = 0.6) +
  scale_color_futurama() +
  my_clean_theme2() +
  theme(legend.position = "none")

df_clinical_corrs %>% 
  filter(q <= 0.1)
df_dietary_corrs %>% 
  filter(q <= 0.1)


strong_clinical_corrs <- 
  df_clinical_corrs %>% 
  filter(q <= 0.1) %>% 
  filter(rho >= 0.75 | rho <= -0.75) %>% 
  mutate(label_col = paste(metadata, feature, sep = " ")) %>% 
  mutate(facet_col = paste0(metadata, " (n = ", n, ")")) %>% 
  ungroup() %>% 
  dplyr::arrange(metadata, rho) %>% 
  mutate(order = as.factor(row_number()))

strong_clinical_corrs_plot <- 
  strong_clinical_corrs %>% 
  ggplot(aes(x=rho, y= order)) +
  geom_segment(aes(x=0, xend=rho, y=order, yend=order), color="gray") +
  geom_point(aes(fill=object_name), shape=21, size=3) +
  xlim(-1, 1) +
  theme_bw() +
  facet_grid(rows = vars(facet_col), scales = "free", space = "free", switch = "y") +
  labs(x="Spearman's Rho", fill = "") +
  scale_y_discrete(
    position = "right",
    breaks = strong_clinical_corrs$order,
    labels = strong_clinical_corrs$feature) +
  scale_fill_nejm() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major.y = element_blank(), 
        strip.background.y = element_rect(fill="white"),
        legend.position = "top")
strong_clinical_corrs_plot
# ggsave(strong_clinical_corrs_plot, filename = "data/Correlations/Clinical_summary_topfeats_new.svg",
#        width = 10, height  = 4)





strong_dietary_corrs <- 
  df_dietary_corrs %>% 
  filter(rho > 0.4 | rho < -0.4) %>% 
  mutate(label_col = paste(metadata, feature, sep = " ")) %>% 
  mutate(facet_col = paste0(metadata, " (n = ", n, ")")) %>% 
  ungroup() %>% 
  dplyr::arrange(metadata, rho) %>% 
  mutate(order = as.factor(row_number()))

strong_dietary_corrs_plot <- 
  strong_dietary_corrs %>% 
  ggplot(aes(x=rho, y= order)) +
  geom_segment(aes(x=0, xend=rho, y=order, yend=order), color="gray") +
  geom_point(aes(fill=object_name), shape=21, size=3) +
  theme_bw() +
  facet_grid(rows = vars(facet_col), scales = "free", space = "free", switch = "y") +
  labs(x="Spearman's Rho", fill = "") +
  scale_y_discrete(
    position = "right",
    breaks = strong_dietary_corrs$order,
    labels = strong_dietary_corrs$feature) +
  scale_fill_nejm() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        panel.grid.major.y = element_blank(), 
        strip.background.y = element_rect(fill="white"),
        legend.position = "top")
strong_dietary_corrs_plot
# ggsave(strong_dietary_corrs_plot, filename = "data/Correlations/Dietary_summary_topfeats.svg",
#        width = 11, height  = 7)



