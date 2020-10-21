# Meta-Analysis  

######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")
# source("src/Meta_Analysis_Community_Composition.R")


abd <- physeq %>% 
  microbiome::abundances() %>% 
  microbiome::transform("compositional")
meta_metadat <- microbiome::meta(physeq)
meta_metadat$studyID <- factor(meta_metadat$studyID)
# meta_metadat$study_condition

# help(adjust_batch)
# The function call indicates for adjust_batch to correct for the effect
# of the batch variable, studyID, while controlling for the effect of the 
# disease variable, study_condition. Many additional options are available 
# through the control parameter, here we specify verbose=FALSE to avoid 
# excessive messages, although they can often be helpful in practice!
fit_adjust_batch <- adjust_batch(feature_abd = abd,
                                 batch = "studyID",
                                 covariates = "study_condition",
                                 data = meta_metadat,
                                 control = list(verbose = T))

# Note that adjust_batch returns a list of more than one components, and 
# feature_abd_adj is the corrected feature abundance matrix. See 
# help(adjust_batch) for the meaning of other components.
abd_adj <- fit_adjust_batch$feature_abd_adj

#-----------------------------------------------------------------
# Replace abundance table with batch corrected abundance
#-----------------------------------------------------------------

abd_adj_input <- otu_table(as.data.frame(abd_adj), taxa_are_rows=T)
my_sample_data <- meta(physeq) %>% sample_data()
physeq_adj <- phyloseq(abd_adj_input, my_sample_data)


#--------------------------------------------
# PERMANOVA QC
#--------------------------------------------

# microbial profiles
D_before <- vegdist(t(abd))
D_after <- vegdist(t(abd_adj))
# fix random seed as adonis runs randomized permutations
set.seed(1)
fit_adonis_before <- adonis(D_before ~ studyID + study_condition + country, data = meta_metadat)
fit_adonis_after <- adonis(D_after ~ studyID + study_condition + country, data = meta_metadat)
print(fit_adonis_before)
print(fit_adonis_after)








#--------------------------------------------------------------------------
# Scrap
#--------------------------------------------------------------------------

# #--------------------------------------------
# # MaAsLin2 with Batch consideration
# #--------------------------------------------
# # lm_meta runs regression and meta-analysis models to identify consistent 
# # effects of the exposure (study_condition, i.e., disease) on feature_abd
# # (microbial feature abundances). Batch variable (studyID) needs to be 
# # specified to identify different studies. Additional covariates to include in 
# # the regression model can be specified via covariates (here set to gender, 
# # age, BMI). Check help(lm_meta) for additional parameter options.
# # Note the warnings: lm_meta can tell if a covariate cannot be meaningfully fit
# # within a batch and will inform the user of such cases through warnings.
# fit_lm_meta <- lm_meta(feature_abd = abd,
#                        batch = "studyID",
#                        exposure = "study_condition",
#                        covariates = c("gender", "age", "BMI"),
#                        data = meta_metadat,
#                        control = list(verbose = T))
# 
# # Again, lm_meta returns a list of more than one components. 
# # meta_fits provides the final meta-analytical testing results. See 
# # help(lm_meta) for the meaning of other components.
# 
# meta_fits <- fit_lm_meta$meta_fits
# 
# meta_fits %>% 
#   filter(qval.fdr < 0.05) %>% 
#   arrange(coef) %>% 
#   mutate(feature = factor(feature, levels = feature)) %>% 
#   ggplot(aes(y = coef, x = feature)) +
#   geom_bar(stat = "identity") +
#   coord_flip()
# 
# #--------------------------------------------------------------
# # Identifying discrete population structures 
# #--------------------------------------------------------------
# 
# # First subset both feature abundance table and metadata to only control samples
# control_meta <- subset(meta_metadat, study_condition == "control")
# control_abd_adj <- abd_adj[, rownames(control_meta)]
# 
# # discrete_discover takes as input sample-by-sample dissimilarity measurements 
# # rather than abundance table. The former can be easily computed from the 
# # latter with existing R packages.
# D_control <- vegdist(t(control_abd_adj))
# fit_discrete <- discrete_discover(D = D_control,
#                                   batch = "studyID",
#                                   data = control_meta,
#                                   control = list(k_max = 8,
#                                                  verbose = FALSE))
# # First subset both feature abundance table and metadata to only control samples
# control_meta <- subset(meta_metadat, study_condition == "control")
# control_abd_adj <- abd_adj[, rownames(control_meta)]
# # discrete_discover takes as input sample-by-sample dissimilarity measurements 
# # rather than abundance table. The former can be easily computed from the 
# # latter with existing R packages.
# D_control <- vegdist(t(control_abd_adj))
# fit_discrete <- discrete_discover(D = D_control,
#                                   batch = "studyID",
#                                   data = control_meta,
#                                   control = list(k_max = 8,
#                                                  verbose = FALSE))
# 
# internal <- data.frame(
#   # By default, fit_discrete evaluates cluster numbers 2-10
#   K = 2:8,
#   statistic = 
#     fit_discrete$internal_mean[, "ZellerG_2014.metaphlan_bugs_list.stool"],
#   se = 
#     fit_discrete$internal_se[, "ZellerG_2014.metaphlan_bugs_list.stool"],
#   type = "internal")
# 
# external <- data.frame(
#   # By default, fit_discrete evaluates cluster numbers 2-10
#   K = 2:8,
#   statistic = 
#     fit_discrete$external_mean[, "ZellerG_2014.metaphlan_bugs_list.stool"],
#   se = 
#     fit_discrete$external_se[, "ZellerG_2014.metaphlan_bugs_list.stool"],
#   type = "external")
# 
# rbind(internal, external) %>% 
#   ggplot(aes(x = K, y = statistic, color = type)) +
#   geom_point(position = position_dodge(width = 0.5)) + 
#   geom_line(position = position_dodge(width = 0.5)) +
#   geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
#                 position = position_dodge(width = 0.5), width = 0.5) +
#   ggtitle("Evaluation of discrete structure in control stool microbiome (ZellerG_2014)")
# 


