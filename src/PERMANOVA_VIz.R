# PERMANOVA PLOTS

# source("src/load_packages.R")
# source("src/miscellaneous_funcs.R")
# source("src/load_phyloseq_obj.R")
# source("src/metadata_prep_funcs.R")
# load("files/low_quality_samples.RData")
# load_all_cohorts()

source("src/PERMANOVA_Analysis.R")


#-------------------------------------------------------------------------------
#####                           Plotting                                 ##### 
#-------------------------------------------------------------------------------

permdf <- permdf %>% filter(metacat != "AA_notmatched")
# All Metadata - pvalue fill
permdf$FDR.symbol <- sig.symbol.generator(permdf$FDR)
permdf.q25 <- permdf %>% filter(FDR <= 0.25)
permdf.q50 <- permdf %>% filter(FDR <= 0.5)
permdf.q75 <- permdf %>% filter(FDR <= 0.75)

p1 <- ggplot() +
  geom_col(data=permdf, aes(x=R2, y=reorder(vars, -FDR), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  # ggtitle("PERMANOVA Covariate Independence") +
  annotate("text", y = unique(permdf$vars), x = (permdf$R2+10), label = permdf$FDR.symbol, size = 2) +
  annotate("text", y = unique(permdf$vars), x = (permdf$R2+3), label = paste0("n =", permdf$n_meta), size = 2) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank(),
        legend.position = c(.90, .50))
p1
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_All_Metdata_Aitchisons.svg", height = 10, width =10)


# All Metadata - METADATA CATEGORY FILL
p2 <- 
  ggplot() +
  geom_col(data=permdf, aes(x=R2, y=reorder(vars, -FDR), fill=metacat)) +
  theme_bw() +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdf$vars), x = (permdf$R2+10), label = permdf$FDR.symbol, size = 2) +
  annotate("text", y = unique(permdf$vars), x = (permdf$R2+3), label = paste0("n =", permdf$n_meta), size = 2) +
  scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.position = c(.80, .50))
p2
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_All_Metdata_Aitchisons_metacatfill.svg", height = 10, width =10)


# Metadata with q value above 0.25- METADATA CATEGORY FILL

p2.2 <- 
  ggplot(data = permdf.q25) +
  geom_col(aes(x=R2, y=reorder(vars, -FDR), fill=metacat), width = 0.6) +
  theme_bw() +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Metadata\nCategory") +
  annotate("text", y = unique(permdf.q25$vars), x = (permdf.q25$R2+6), label = permdf.q25$FDR.symbol, size = 3) +
  annotate("text", y = unique(permdf.q25$vars), x = (permdf.q25$R2+3), label = paste0("n =", permdf.q25$n_meta), size = 3) +
  scale_fill_simpsons() +
  theme(axis.title.y = element_blank(),
        panel.grid.minor  = element_blank(),
        panel.grid.major.y = element_blank())
p2.2
# ggsave(p2.2, filename = "data/PERMANOVA_Analysis/PERMANOVA_Barplt_Q0.25_Aitchisons_metacatfill.svg", height = 5, width = 13)


###########################################   FIGURE 1 Plots #########################################################

############# Metadata Category Summary - Fill by p-Value
p3 <- 
  ggplot(data = permdf.q50) +
  geom_col(aes(x=R2, y=reorder(metacat, R2), fill=FDR), width = 0.6) +
  theme_bw() +
  ggtitle("Cumulative Metadata Explained Variance") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  scale_fill_viridis(direction=1, option="cividis") +
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
p3
ggsave(p3, filename = "data/PERMANOVA_Analysis/PERMANOVA_Barplt_Cumulative_Q0.25_Aitchisons.svg", height = 3, width =5)

#################### FILTERING FOR Essential METADATA

permdfess$FDR.symbol <- sig.symbol.generator(permdfess$FDR)

p6 <- ggplot() +
  geom_col(data=permdfess, aes(x=R2, y=reorder(vars, -FDR), fill=metacat)) +
  theme_bw() +
  ggtitle("Foundational Metadata") +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Category") +
  # annotate("text", y = unique(permdfess$vars), x = (permdfess$R2+4), label = paste0("p =", round(permdfess$FDR, digits = 3)), size = 3) +
  annotate("text", y = unique(permdfess$vars), x = (permdfess$R2+11), label = permdfess$FDR.symbol, size = 3) +
  annotate("text", y = unique(permdfess$vars), x = (permdfess$R2+4), label = paste0("n =", permdfess$n_meta), size = 3) +
  scale_fill_simpsons()+
  theme(axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank())
p6
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_essential_Metdata_Aitchisons.svg", height = 4, width =10)

c1 <- cowplot::plot_grid(p3, p6, ncol=2, align = "hv")
# ggsave(c1, filename = "figures/Figure_1/PERMANOVA_Analysis.svg", height = 3.25, width =20)


####################################################################################################





####################  FILTERING FOR FDR SIGNIFICANT 

permdfsig <- filter(permdf, FDR <= 0.05)
p4 <- ggplot() +
  geom_col(data=permdfsig, aes(x=R2, y=reorder(vars, R2), fill=FDR)) +
  theme_bw() +
  ggtitle("Significant Metadata") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdfsig$vars), x = (permdfsig$R2+3), label = paste0("p =", round(permdfsig$FDR, digits = 3)), size = 2) +
  annotate("text", y = unique(permdfsig$vars), x = (permdfsig$R2+7), label = paste0("n =", permdfsig$n_meta), size = 2) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
p4
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_Significant_Metdata_Aitchisons.svg", height = 4, width =10)





#################### FILTERING FOR CORE METADATA - n > 90 (Metadata present with most samples)

p5 <- ggplot() +
  geom_col(data=permdfcoresig, aes(x=R2, y=reorder(vars, R2), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  ggtitle("PERMANOVA Covariate Independence") +
  annotate("text", y = unique(permdfcoresig$vars), 
           x = (permdfcoresig$R2+0.7), label = paste0("p =", round(permdfcoresig$FDR, digits = 3)), size = 2) +
  annotate("text", y = unique(permdfcoresig$vars), 
           x = (permdfcoresig$R2+1.5), label = paste0("n =", permdfcoresig$n_meta), size = 2) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
p5
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_CoreFeatures_Metdata_Aitchisons.svg", height = 4, width =10)








#################### Extraneous Plots #################### 

#################### FILTERING FOR DIET METADATA

e1 <- ggplot() +
  geom_col(data=permdf_diet, aes(x=R2, y=reorder(vars, -FDR), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdf_diet$vars), x = (permdf_diet$R2+0.5), label = paste0("p =", round(permdf_diet$FDR, digits = 3)), size = 3) +
  annotate("text", y = unique(permdf_diet$vars), x = (permdf_diet$R2+1), label = paste0("n =", permdf_diet$n_meta), size = 3) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
e1
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_Diet_Metdata_Aitchisons.svg", height = 5, width =10)



#################### FILTERING FOR Supplement METADATA

# permdf_supp <- filter(permdf, metacat == "Supplements")

e2 <- ggplot() +
  geom_col(data=permdf_supp, aes(x=R2, y=reorder(vars, -FDR), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdf_supp$vars), x = (permdf_supp$R2+0.5), label = paste0("p =", round(permdf_supp$FDR, digits = 3)), size = 3) +
  annotate("text", y = unique(permdf_supp$vars), x = (permdf_supp$R2+1), label = paste0("n =", permdf_supp$n_meta), size = 3) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
e2
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_Supplements_Metdata_Aitchisons.svg", height = 5, width =10)



#################### FILTERING FOR General Drug use METADATA

# permdf_gendrug <- filter(permdf, metacat == "Other medication")

e3 <- ggplot() +
  geom_col(data=permdf_gendrug, aes(x=R2, y=reorder(vars, -FDR), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdf_gendrug$vars), x = (permdf_gendrug$R2+0.5), label = paste0("p =", round(permdf_gendrug$FDR, digits = 3)), size = 3) +
  annotate("text", y = unique(permdf_gendrug$vars), x = (permdf_gendrug$R2+1), label = paste0("n =", permdf_gendrug$n_meta), size = 3) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_DrugUse_Metdata_Aitchisons.svg", height = 5, width =10)



#################### FILTERING FOR PD Drug use METADATA

# permdf_pddrug <- filter(permdf, metacat == "PD medication")

e4 <- ggplot() +
  geom_col(data=permdf_pddrug, aes(x=R2, y=reorder(vars, -FDR), fill=`FDR`)) +
  theme_bw() +
  # labs(x="Variance Explained - Adjusted R^2 %", fill = "Adjusted\np-value") +
  labs(x=expression("Variance Explained - " ~ R^2~ "(%)"), fill = "Adjusted\np-value") +
  annotate("text", y = unique(permdf_pddrug$vars), x = (permdf_pddrug$R2+0.5), label = paste0("p =", round(permdf_pddrug$FDR, digits = 3)), size = 3) +
  annotate("text", y = unique(permdf_pddrug$vars), x = (permdf_pddrug$R2+1), label = paste0("n =", permdf_pddrug$n_meta), size = 3) +
  scale_fill_viridis(direction=-1, option="cividis")+
  theme(axis.title.y = element_blank())
# ggsave("data/PERMANOVA_Analysis/PERMANOVA_Barplt_PD_DrugUse_Metdata_Aitchisons.svg", height = 5, width =10)


