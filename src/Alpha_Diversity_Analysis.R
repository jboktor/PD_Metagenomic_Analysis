##### Alpha Diversity Boxplots Script

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot);
library(EnvStats);library(sjlabelled);library(sjmisc);library(sjPlot);library(nlme)

rm(list = ls())

######## Load Data & functions
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")


################################################################################# 
########################### SWAP INPUT LEVEL HERE: ########################### 
################################################################################# 

# x <- c(dat, dat.path, dat.ec, dat.KOs.all)
# z <- c("species", "pathways", "enzymes", "genes")


# Phylo_Objects$KOs.all
obj <- "Enzymes"
dat_obj <- Phylo_Objects[[obj]]


################################################################################# 

### WARNING : OTU TABLE CONVERTED TO BINARY - DO NOT USE DOWNSTREAM FOR other analyses
dat_alpha <- dat_obj
# otu_table(dat_alpha) <- ((otu_table(dat_alpha) > 0) + 0)
## Calculate Alpha Diversity Metrics
tab <- microbiome::alpha(dat_alpha, index= c("observed" ,'shannon', "rarity_log_modulo_skewness"))


# Run Metadata pre-processing function
process_meta(dat_alpha)
env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))

# Add specific diversity cols to df
env$observed <- tab$observed
env$Shannon <- tab$diversity_shannon 
env$rarity_log_modulo_skewness <- tab$rarity_log_modulo_skewness 

# Create new paired column with only household pairs and NAs for rest 
env$Paired.plot <- as.numeric(levels(env$Paired))[env$Paired]
env[which(env$Paired.plot > 30),"Paired.plot"] <- NA
env.pairs <- dplyr::filter(env, Paired.plot < 30)

# Plot histograms to get a sense of data distribution
par(mfrow = c(1, 3))
hist(env$observed, main="observed OTUs", xlab="", breaks=10)
hist(env$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(env$rarity_log_modulo_skewness, main="rarity_log_modulo_skewness", xlab="", breaks=10)



######### Test for Normality
# shapiro.test(env$Shannon) # Normal
# shapiro.test(env$diversity_gini_simpson) # Non-normal
# shapiro.test(env$observed) # Normal
# shapiro.test(env$chao1) # Normal



########### Observed Species Plot ########### 
# Linear Mixed Model
lmm.obs <- lme(observed ~ description+host_age+host_body_mass_index+sex, random= ~ 1 | Paired, data=env, na.action = na.omit)
summary(lmm.obs)
# coef(lmm.obs)
# plot_model(lmm.obs, show.values = TRUE, value.offset = .3)
# tab_model(lmm.obs
pc_sig <- summary(lmm.obs)$tTable["descriptionPopulation Control","p-value"]
hc_sub <- summary(lmm.obs)$tTable["descriptionHousehold Control","p-value"]

set.seed(123)
p1 <- ggplot(env, aes(x = donor_group, y = observed)) + theme_minimal() + 
  geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
  geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
  geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
  geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Observed Counts") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # theme(axis.text.x = element_blank())+
  labs(fill="Group") +
  # scale_linetype_manual(values= line_format) +
  scale_fill_manual(values = c("HC" = "#440154", 
                               "PD" = "#FDE725", "PC" = "#21908C")) +
  # stat_pvalue_manual(summary(fit)$tTable, label = "p-value") +
  geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(hc_sub)) +
  geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(pc_sig))
# ggsave("AlphaDiversity_BoxViolinPlot_Observed_species.png", height = 6, width =3)




########### Shannon Diversity Plot ########### 

# Linear Mixed Model
lmm.shannon <- lme(Shannon ~ description+host_age+host_body_mass_index+sex, random= ~ 1 | Paired, data=env, na.action = na.omit)
summary(lmm.shannon)
# plot_model(lmm.shannon, show.values = TRUE, value.offset = .3)
# tab_model(lmm.shannon)
pc_sig <- summary(lmm.shannon)$tTable["descriptionPopulation Control","p-value"]
hc_sub <- summary(lmm.shannon)$tTable["descriptionHousehold Control","p-value"]


set.seed(123)
p2 <- ggplot(env, aes(x = donor_group, y = Shannon)) + theme_minimal() + 
  geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
  geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
  geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
  geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("Shannon's Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # theme(axis.text.x = element_blank())+
  labs(fill="Group") +
  # scale_linetype_manual(values= line_format) +
  scale_fill_manual(values = c("HC" = "#440154", 
                               "PD" = "#FDE725", "PC" = "#21908C")) +
  geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(hc_sub)) +
  geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(pc_sig))
# ggsave("AlphaDiversity_BoxViolinPlot_Shannon.png", height = 6, width =3)




########### Rarity- log-modulo Skewness ########### 

# Linear Mixed Model
lmm.rarity <- lme(rarity_log_modulo_skewness ~ description+host_age+host_body_mass_index+sex, random= ~ 1 | Paired, data=env, na.action = na.omit)
summary(lmm.rarity)
# plot_model(lmm.rarity, show.values = TRUE, value.offset = .3)
# tab_model(lmm.rarity)
pc_sig <- summary(lmm.rarity)$tTable["descriptionPopulation Control","p-value"]
hc_sub <- summary(lmm.rarity)$tTable["descriptionHousehold Control","p-value"]


set.seed(123)
p3 <- ggplot(env, aes(x = donor_group, y = rarity_log_modulo_skewness)) + theme_minimal() + 
  geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
  geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
  geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
  geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
  theme(axis.title.x=element_blank(),
        legend.position = "none") +
  ylab("rarity_log_modulo_skewness's Diversity") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # theme(axis.text.x = element_blank())+
  labs(fill="Group") +
  # scale_linetype_manual(values= line_format) +
  scale_fill_manual(values = c("HC" = "#440154", 
                               "PD" = "#FDE725", "PC" = "#21908C")) +
  geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(hc_sub)) +
  geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(pc_sig))
# ggsave("AlphaDiversity_BoxViolinPlot_rarity_log_modulo_skewness.png", height = 6, width =3)






alpha_cow <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "v")
alpha_cow
# ggsave(alpha_cow, filename = paste0("Alpha_Diversity_Analysis/AlphaDiversity_BoxViolinPlot_", obj,"_Summary.png"),
#        height = 6, width =7)






############ CODE Junkyard ############ 

# Similar option with LMER
# mixed.ranslope <- lmer(Shannon ~ description+host_age+host_body_mass_index+sex + (1 | Paired), data = env) 
# summary(mixed.ranslope)
# plot_model(mixed.ranslope, show.values = TRUE, value.offset = .3)
# tab_model(mixed.ranslope)

# (mm_plot <- ggplot(env, aes(x = bristol_stool_scale, y = Shannon, colour = host_body_mass_index)) +
#     facet_wrap(~description, nrow=1) +   # a panel for each mountain range
#     geom_point(alpha = 0.8) +
#     theme_classic() +
#     # geom_line(data = cbind(env, pred = predict(mixed.ranslope)), aes(y = pred), size = 1) +  # adding predicted line from mixed model
#     theme(legend.position = "bottom",
#           panel.spacing = unit(2, "lines"))  # adding space between panels
# )


# ########### Gini-Simpson Diversity Plot ########### 
# env$description <- factor(env$description, levels=c("PC", "PD", "HC"))
# 
# # Test for normality - Distb is non-normal
# shapiro.test(env$diversity_gini_simpson) # test is Non-Normal!
# # KW AND Wilcoxon test 
# # env$description <- as.factor(env$description)
# kruskal.test(env$diversity_gini_simpson ~ description, data = env)
# sigv <- dunnTest(env$diversity_gini_simpson ~ description, data = env, method="bh")
# sigv
# 
# # Linear Mixed Model
# lmm.gini <- lme(diversity_gini_simpson ~ description+host_age+host_body_mass_index+sex, random= ~ 1 | Paired, data=env, na.action = na.omit)
# summary(lmm.gini)
# plot_model(lmm.gini, show.values = TRUE, value.offset = .3)
# tab_model(lmm.gini)
# 
# set.seed(123)
# ggplot(env, aes(x = description, y = diversity_gini_simpson)) + theme_minimal() + 
#   geom_violin(draw_quantiles = c(0.5), trim = T) +
#   geom_boxplot(aes(fill = donor_group), width=0.1, alpha = 0.6, outlier.alpha = 0) +
#   geom_point(aes(fill=description), position = position_jitterdodge(dodge.width=.3),shape=21, size=2, alpha = 0.6) +
#   theme(axis.title.x=element_blank(),
#         legend.position = "none") +
#   ylab("Gini-Simpson") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   labs(fill="Group") +
#   scale_fill_manual(values = c("HC" = "#440154", 
#                                "PD" = "#FDE725", "PC" = "#21908C")) +
#   geom_signif(comparisons = list(c("PD", "HC")), map_signif_level= sigv$res$P.adj[2]) +
#   geom_signif(comparisons = list(c("PC", "PD")) , annotations="**")
# # ggsave("AlphaDiversity_BoxViolinPlot_Gini_Simpson.png", height = 6, width =3)
# 

######################## GINI-SIMPSION V BSS ######################## 


# # env$bristol_stool_scale[is.na(env$bristol_stool_scale)]
# env[env=="not provided"] <- NA
# env$bristol_stool_scale <- as.numeric(env$bristol_stool_scale)
# 
# 
# ggplot(env, aes(x=bristol_stool_scale, y=diversity_gini_simpson )) +
#   geom_smooth(aes(color = description), method = 'lm', se = F) +
#   theme_minimal() +
#   # geom_point(aes(fill=description), size=2.5, shape=21, alpha=0.6) +
#   geom_point(aes(fill=description), position = position_jitterdodge(dodge.width=.3),shape=21, size=2.5, alpha = 0.6) +
#   scale_fill_manual(values = c("HC" = "#440154", 
#                                "PD" = "#FDE725", "PC" = "#21908C"))+
#   scale_color_manual(values = c("HC" = "#440154", 
#                                "PD" = "#FDE725", "PC" = "#21908C")) +
#   theme(legend.position = "none") +
#   labs(x="Bristol Stool Scale", y="Gini-Simpson")
# ggsave("AlphaDiversity_BSS_v_GiniSimpson.png", height = 6, width =6)


######################## Paired boxplot example with shannon Diversity PD vs HC
# env$Paired
# env$Paired <- as.numeric(levels(env$Paired))[env$Paired]
# env.pairs <- dplyr::filter(env, Paired < 30)
# ggpaired(env.pairs,
#          x = "donor_group", y = "Shannon",
#          line.color = "gray", line.size = 0.4,
#          palette = "jco")+
#   stat_compare_means(paired = TRUE)


# # ANOVA TEST HERE
# aov(env$Shannon~env$description)
# anovaResults <- aov(Shannon~description, data=env)
# summary(anovaResults)
# sigv <- TukeyHSD(anovaResults)
# sigv$description[3,4]

  