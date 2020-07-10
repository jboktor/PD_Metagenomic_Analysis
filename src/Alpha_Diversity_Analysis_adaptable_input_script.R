# Alpha Diversity adaptable

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
###########################  INPUT LEVELS HERE: ########################### 
################################################################################# 

x <- c(dat, dat.path, dat.ec, dat.KOs.all)
z <- c("Species", "Pathways", "Enzymes", "KOs.all")


################################ Begining of Loop################################ 


cnt <- 1
for (i in x){
  
  cat("Processing input: ", z[cnt], "\n")
  print(i)
  cat("\n")
  
  dat_alpha <- i

  ## Calculate Alpha Diversity Metrics
  # tab <- microbiome::alpha(abundances(dat_alpha), index= c("observed" ,'shannon', "rarity_log_modulo_skewness"))
  
  # Run Metadata pre-processing function
  process_meta(dat_alpha)
  env$description <- factor(env$description, levels=c("PD Patient", "Population Control", "Household Control"))
  env$donor_group <- factor(env$donor_group, levels=c("PC", "PD", "HC"))
  
  ## Calculate Alpha Diversity Metrics and add cols to df
  env$Observed <- alpha(abundances(dat_alpha), 'observed')$observed
  env$Shannon <- alpha(abundances(dat_alpha), 'shannon')$diversity_shannon
  env$Evenness <- evenness(abundances(dat_alpha), 'simpson')$simpson
  
  # Create new paired column with only household pairs and NAs for rest 
  env$Paired.plot <- as.numeric(levels(env$Paired))[env$Paired]
  env[which(env$Paired.plot > 30),"Paired.plot"] <- NA
  env.pairs <- dplyr::filter(env, Paired.plot < 30)
  
  # Plot histograms to get a sense of data distribution
  par(mfrow = c(1, 3))
  hist(env$Observed, main="observed OTUs", xlab="", breaks=10)
  hist(env$Shannon, main="Shannon diversity", xlab="", breaks=10)
  hist(env$Evenness, main="Simpson's evenness", xlab="", breaks=10)
  
  ######### Test for Normality
  # shapiro.test(env$Shannon) # Normal
  # shapiro.test(env$rarity_log_modulo_skewness) # Non-normal
  # shapiro.test(env$observed) # Normal
  
  ########################################################################################
  #############################    Functions for Analysis   #############################    
  
  lm.PdPc <- function(metadf, metric){
    ###' Function conducts Linear Model for PD vs PC
    env.PdPc <- filter(metadf, donor_group != "HC")
    formula <- as.formula(
      paste(metric, "~", paste(c("description", "host_age_factor", "host_body_mass_index", "sex"), collapse="+") ) )
    
    linear.model <- lm(formula, data=env.PdPc, na.action = na.omit)
    plot_model(linear.model, show.values = TRUE, value.offset = .3)
    qqnorm(resid(linear.model))
    qqline(resid(linear.model))
    dev.off()
    return(linear.model)
  }

  
  lmm.PdHc <- function(metadf, metric){
    ###' Function conducts Linear Mixed Model for PD vs HC
    env.PdHc <- filter(metadf, Paired.plot < 30)

    formula <- as.formula(paste(metric, "~", paste(c("description"))))
    lmm <- lme(formula, random= ~ 1 | Paired, data=env.PdHc, na.action = na.omit)
    qqnorm(resid(lmm))
    qqline(resid(lmm))
    return(lmm)
  }

  ########################################################################################
  
  
  
  
  
  ########### Observed Species Plot ########### 
  
  
  ### STATS 
  observed.PdPC <- lm.PdPc(metadf=env, metric="Observed")
  observed.PdHC <- lmm.PdHc(metadf=env, metric="Observed")
  ## Pull p-values
  observed.PdPC.pval <- summary(observed.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
  observed.PdHC.pval <- summary(observed.PdHC)$tTable["descriptionHousehold Control","p-value"]
  
  ### PLOT 
  set.seed(123)
  p1 <- ggplot(env, aes(x = donor_group, y = Observed)) + theme_minimal() + 
    geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
    geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
    geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
    geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab(paste0("Observed Counts: ", z[cnt])) +
    theme(plot.title = element_text(hjust = 0.5)) +
    # theme(axis.text.x = element_blank())+
    labs(fill="Group") +
    # scale_linetype_manual(values= line_format) +
    scale_fill_manual(values = c("HC" = "#440154", 
                                 "PD" = "#FDE725", "PC" = "#21908C")) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(observed.PdHC.pval)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(observed.PdPC.pval))

  
  
  ########### Shannon Diversity Plot ########### 
  
  ### STATS 
  observed.PdPC <- lm.PdPc(metadf=env, metric="Shannon")
  observed.PdHC <- lmm.PdHc(metadf=env, metric="Shannon")
  ## Pull p-values
  observed.PdPC.pval <- summary(observed.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
  observed.PdHC.pval <- summary(observed.PdHC)$tTable["descriptionHousehold Control","p-value"]
  
  
  ### PLOT 
  set.seed(123)
  p2 <- ggplot(env, aes(x = donor_group, y = Shannon)) + theme_minimal() + 
    geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
    geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
    geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
    geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab(paste0("Shannon's Diversity: ", z[cnt])) +
    labs(fill="Group") +
    scale_fill_manual(values = c("HC" = "#440154", 
                                 "PD" = "#FDE725", "PC" = "#21908C")) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(observed.PdHC.pval)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(observed.PdPC.pval))

  
  
  ########### EVENESS ########### 
  
  ### STATS 
  observed.PdPC <- lm.PdPc(metadf=env, metric="Evenness")
  observed.PdHC <- lmm.PdHc(metadf=env, metric="Evenness")
  ## Pull p-values
  observed.PdPC.pval <- summary(observed.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
  observed.PdHC.pval <- summary(observed.PdHC)$tTable["descriptionHousehold Control","p-value"]
  
  
  ### PLOT 
  set.seed(123)
  p3 <- ggplot(env, aes(x = donor_group, y = Evenness)) + theme_minimal() + 
    geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
    geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
    geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
    geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
    theme(axis.title.x=element_blank(),
          legend.position = "none") +
    ylab(paste0("Simpson's Evenness: ", z[cnt])) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(fill="Group") +
    scale_fill_manual(values = c("HC" = "#440154", 
                                 "PD" = "#FDE725", "PC" = "#21908C")) +
    geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(observed.PdHC.pval)) +
    geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(observed.PdPC.pval))

  

  alpha_cow <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "v")
  alpha_cow
  ggsave(alpha_cow, filename = paste0("data/Alpha_Diversity_Analysis/AlphaDiversity_BoxViolinPlot_",  z[cnt], "_Summary.svg"),
         height = 6, width =8)

  cnt <- cnt + 1
  
}


# ####################################    JUNKYARD   ####################################################

# #############################    Plotting Functions
# 
# paired_violin_plot <- function(df, metric, PC.pvalue, HC.pvalue){
#   ###' Plotting Alpha Diversity
#   
#   df.pairs <- dplyr::filter(df, Paired.plot < 30)
#   print(df.pairs$Paired.plot)
#   print(df.pairs$donor_group)
#   
#   p1 <- ggplot(df, aes(x = donor_group, y = metric)) + theme_minimal() + 
#     geom_violin(draw_quantiles = c(0.5), trim = T, width = 0.75) +
#     geom_boxplot(aes(fill = donor_group), width=0.15, alpha = 0.6, outlier.alpha = 0) +
#     geom_point(aes(fill = donor_group), position = position_jitterdodge(dodge.width=1),shape=21, size=1, alpha = 0.6) +
#     geom_line(aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
#     
#     # geom_line(data = env.pairs, aes(group = Paired.plot), linetype = 'solid', color = "grey", alpha = 0.7) +
#     
#     labs(fill="Group") +
#     scale_fill_manual(values = c("HC" = "#440154", "PD" = "#FDE725", "PC" = "#21908C")) +
#     geom_signif(comparisons = list(c("PD", "HC")), annotations = sig_mapper(HC.pvalue)) +
#     geom_signif(comparisons = list(c("PC", "PD")), annotations = sig_mapper(PC.pvalue)) +
#     # ylab(paste0("Observed Counts: ", z[cnt])) +
#     theme(axis.title.x=element_blank(),
#           legend.position = "none",
#           plot.title = element_text(hjust = 0.5))
#   p1
#   return(p1)
# }
# 
