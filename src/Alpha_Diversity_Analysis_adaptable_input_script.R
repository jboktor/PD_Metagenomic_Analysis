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
load("files/Eggnogs.slim_PhyloseqObj.RData")
load("files/Eggnogs_PhyloseqObj.RData")
load("files/Pfams.slim_PhyloseqObj.RData")
load("files/Pfams_PhyloseqObj.RData")

################################################################################# 
###########################  INPUT LEVELS HERE: ########################### 
################################################################################# 

x <- c(dat, dat.path, dat.ec, dat.KOs, dat.Eggnogs, dat.Pfams,
       dat.path.slim, dat.ec.slim, dat.KOs.slim, dat.Eggnogs.slim, dat.Pfams.slim)
z <- c("Species", "Pathways", "Enzymes", "KOs", "Eggnogs", "Pfams", 
       "Pathways.slim", "Enzymes.slim", "KOs.slim", "Eggnogs.slim", "Pfams.slim")

color_palette <- c("HC" = "#440154", "PD" = "#FDE725", "PC" = "#21908C")

################################ Begining of Loop################################ 


cnt <- 1
for (i in x){
  
  cat("Processing input: ", z[cnt], "\n")
  print(i)
  cat("\n")
  
  dat_alpha <- i


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
  
  p1 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Observed, 
                           df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Observed,
                           pairs.column=env.pairs$Paired.plot,
                           cols=color_palette, ylabel = paste0("Observed Counts: ", z[cnt]), 
                           PDvPC.stat = observed.PdPC.pval, PDvHC.stat = observed.PdHC.pval)
  
  ########### Shannon Diversity ########### 
  
  ### STATS 
  shannon.PdPC <- lm.PdPc(metadf=env, metric="Shannon")
  shannon.PdHC <- lmm.PdHc(metadf=env, metric="Shannon")
  ## Pull p-values
  shannon.PdPC.pval <- summary(shannon.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
  shannon.PdHC.pval <- summary(shannon.PdHC)$tTable["descriptionHousehold Control","p-value"]
  
  
  p2 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Shannon, 
                           df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Shannon,
                           pairs.column=env.pairs$Paired.plot,
                           cols=color_palette, ylabel = paste0("Shannon's Diversity: ", z[cnt]), 
                           PDvPC.stat = shannon.PdPC.pval, PDvHC.stat = shannon.PdHC.pval)
  
  
  ########### EVENESS ########### 
  
  ### STATS 
  evenness.PdPC <- lm.PdPc(metadf=env, metric="Evenness")
  evenness.PdHC <- lmm.PdHc(metadf=env, metric="Evenness")
  ## Pull p-values
  evenness.PdPC.pval <- summary(evenness.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
  evenness.PdHC.pval <- summary(evenness.PdHC)$tTable["descriptionHousehold Control","p-value"]
  
  
  ### PLOT 
  p3 <- alpha_div_boxplots(df=env, x=env$donor_group, y=env$Evenness, 
                     df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Evenness,
                     pairs.column=env.pairs$Paired.plot,
                     cols=color_palette, ylabel = paste0("Simpson's Evenness: ", z[cnt]), 
                     PDvPC.stat = evenness.PdPC.pval, PDvHC.stat = evenness.PdHC.pval)

  
  ### MERGE PLOTS ### 
  alpha_cow <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = "v")
  alpha_cow
  ggsave(alpha_cow, filename = paste0("data/Alpha_Diversity_Analysis/AlphaDiversity_BoxPlot_",  z[cnt], "_Summary.svg"),
         height = 6, width =8)

  cnt <- cnt + 1
  
}

