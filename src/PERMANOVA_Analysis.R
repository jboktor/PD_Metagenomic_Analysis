## PERMANOVA script

rm(list = ls())
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
load("files/low_quality_samples.RData")
load_data("Merged")

#-------------------------------------------------------------------------------
#####                       PERMANOVA Analysis                             ##### 
#-------------------------------------------------------------------------------

distlist <- c("euclidean", "bray")
sigmetadata <- vector("list", length(distlist))
names(sigmetadata) = distlist

for (j in distlist) {
  sigmetadata[[j]]$siglst <- c()
  sigmetadata[[j]]$almostSiglst <- c()
}
siglst <- c()
env.pa_perm <- tibble()

objs <-
  c(dat.species,
    dat.path,
    dat.ec.slim,
    dat.KOs.slim,
    dat.PFAMs.slim,
    dat.EGGNOGs.slim)
obj_name <-
  c("Species",
    "Pathways",
    "Enzymes.slim",
    "KOs.slim",
    "Pfams.slim",
    "Eggnogs.slim")

#-------------------------------------------------------------------------------
#####                       Aitchisons Distance                            ##### 
#-------------------------------------------------------------------------------

cnt <- 1
for (obj in objs) {
  obj_processed <- obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    microbiome::transform("compositional") %>% 
    microbiome::transform("clr")
  sp <- t(abundances(obj_processed))
  # Run Metadata pre-processing function
  env <- process_meta(obj_processed, cohort = "Merged")
  env <- trim_meta(env, min_ratio = 2/234)
  # remove ID variable and duplicate grouping variable
  env.pa <- dplyr::select(env, -matches(other_metadata))
  ## Adding Alpha Diversity to PERMANOVA
  alpha_diversity <- obj_processed %>% abundances() %>% 
    microbiome::alpha('shannon')
  env.pa$diversity_shannon <- alpha_diversity$diversity_shannon
    
  for (i in 1:length(env.pa)) { 
    a <- env.pa[,i]
    a.narm <- na.omit(a)
    if (any(is.na(a))) {
      sp.narm <- sp[-attr(a.narm, "na.action"), ]
    } else {
      sp.narm <- sp
    }
    cat("Aitchisons", obj_name[cnt], ": ", colnames(env.pa[i]), "\n")
    meta_ano = adonis(vegdist(sp.narm, method = "euclidean") ~ a.narm, permutations = 9999)
    row2add <-
      cbind(colnames(env.pa[i]),
            meta_ano$aov.tab[1, ],
            length(a.narm),
            obj_name[cnt], 
            "Aitchisons")
    env.pa_perm <- rbind(env.pa_perm, row2add)
    sig <- meta_ano$aov.tab$`Pr(>F)`[1]
  }
  cnt = cnt + 1
}
permanova_aitch <- env.pa_perm %>%
  remove_rownames() %>%
  dplyr::rename(
    "p_value" = "Pr(>F)",
    "vars" = "colnames(env.pa[i])",
    "n_meta" = "length(a.narm)",
    "data_type" = "obj_name[cnt]",
    "distance" = '"Aitchisons"'
  ) 

#-------------------------------------------------------------------------------
#####                       Bray-Curtis Distance                            ##### 
#-------------------------------------------------------------------------------

env.pa_perm <- tibble()
cnt <- 1
for (obj in objs) {
  obj_processed <- obj %>% 
    subset_samples(donor_id %ni% low_qc[[1]]) %>% 
    microbiome::transform("compositional")
  sp <- t(abundances(obj_processed))
  # Run Metadata pre-processing function
  env <- process_meta(obj_processed, cohort = "Merged")
  env <- trim_meta(env, min_ratio = 2/234)
  # remove ID variable and duplicate grouping variable
  env.pa <- dplyr::select(env, -matches(other_metadata))
  ## Adding Alpha Diversity to PERMANOVA
  alpha_diversity <- obj_processed %>% abundances() %>% 
    microbiome::alpha('shannon')
  env.pa$diversity_shannon <- alpha_diversity$diversity_shannon
  
  for (i in 1:length(env.pa)) { 
    a <- env.pa[,i]
    a.narm <- na.omit(a)
    if (any(is.na(a))) {
      sp.narm <- sp[-attr(a.narm, "na.action"), ]
    } else {
      sp.narm <- sp
    }
    cat("Bray-Curtis", obj_name[cnt], ": ", colnames(env.pa[i]), "\n")
    meta_ano = adonis(vegdist(sp.narm, method = "bray") ~ a.narm, permutations = 9999)
    row2add <-
      cbind(colnames(env.pa[i]),
            meta_ano$aov.tab[1, ],
            length(a.narm),
            obj_name[cnt],
            "BrayCurtis")
    env.pa_perm <- rbind(env.pa_perm, row2add)
    sig <- meta_ano$aov.tab$`Pr(>F)`[1]
  }
  cnt = cnt + 1
}
permanova_bray <- env.pa_perm %>%
  remove_rownames() %>%
  dplyr::rename(
    "p_value" = "Pr(>F)",
    "vars" = "colnames(env.pa[i])",
    "n_meta" = "length(a.narm)",
    "data_type" = "obj_name[cnt]",
    "distance" = '"BrayCurtis"'
  )  

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

permdf <- 
  bind_rows(permanova_aitch, permanova_bray) %>% 
  dplyr::mutate(vars = as.character(vars)) %>%
  dplyr::mutate(R2 = R2 * 100) %>%
  group_by(data_type, distance) %>%
  mutate(FDR = p.adjust(p_value, method = 'BH')) %>%
  ungroup()


## Add Metadata category column
permdf <-
  mutate(permdf,
         metacat = if_else(
           permdf$vars %in% grouping,
           "Disease-Status/Grouping",
           if_else(
             permdf$vars %in% Anthro,
             "Anthropomorphic",
             if_else(
               permdf$vars %in% GI,
               "GI distress",
               if_else(
                 permdf$vars %in% nonstarters,
                 "Microbiota disruption",
                 if_else(
                   permdf$vars %in% allergy,
                   "Allergies",
                   if_else(
                     permdf$vars %in% smoke,
                     "Smoking",
                     if_else(
                       permdf$vars %in% general_drugs,
                       "Other medication",
                       if_else(
                         permdf$vars %in% pd_drugs,
                         "PD medication",
                         if_else(
                           permdf$vars %in% supplements,
                           "Supplements",
                           if_else(
                             permdf$vars %in% diet,
                             "Diet",
                             if_else(
                               permdf$vars %in% motor_severity_scores_summary,
                               "Motor Symptom Severity",
                               if_else(
                                 permdf$vars %in% clinical_variables,
                                 "Clinical Variables",
                                 if_else(
                                   permdf$vars %in% environmental,
                                   "Environment",
                                   if_else(permdf$vars == "Shannon", "Shannon Diversity",
                                           "AA_notmatched")
                                 )
                               )
                             )
                           )
                         )
                       )
                     )
                   )
                 )
               )
             )
           )
         ))

write.csv(permdf, file = paste0('files/permanova_analysis_', Sys.Date(), '.csv'))

# FILTERING FOR FDR SIGNIFICANT
permdfsig <- filter(permdf, FDR <= 0.05)

# FILTERING FOR CORE METADATA - n > 90 (Metadata present with most samples)
permdfcore <- filter(permdf, n_meta >= 90)
permdfcoresig <- filter(permdfcore, FDR <= 0.5)

# FILTERING FOR Essential METADATA
essential <-
  c(
    "Disease-Status/Grouping",
    "Anthropomorphic",
    "GI distress",
    "Microbiota disruption",
    "Shannon Diversity"
  )
permdfess <- filter(permdf, metacat %in% essential)

# FILTERING FOR DIET METADATA
permdf_diet <- filter(permdf, metacat == "Diet")

# FILTERING FOR Supplement METADATA
permdf_supp <- filter(permdf, metacat == "Supplements")

# FILTERING FOR General Drug use METADATA
permdf_gendrug <- filter(permdf, metacat == "Other medication")

# FILTERING FOR PD Drug use METADATA
permdf_pddrug <- filter(permdf, metacat == "PD medication")



