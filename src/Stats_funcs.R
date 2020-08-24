## Statistical test functions

# # Total VFBD Stats
# # Summary statistics
# mag2 %>%
#   group_by(donor_group) %>%
#   get_summary_stats(value, type = "mean_sd")
# 
# # Check assumptions - Outliers
# mag2 %>% 
#   group_by(donor_group) %>%
#   identify_outliers(value)
# 
# # Check assumptions - Normality
# model  <- lm(value ~ donor_group, data = mag2)
# ggqqplot(residuals(model))
# # Shapiro-Wilk test of normality
# shapiro_test(residuals(model))
# ggqqplot(mag2, "value", facet.by = "donor_group")
# 
# # Check assumptions - Homogneity of variance
# plot(model, 1) # No obvious relationship is good
# mag2 %>% levene_test(value ~ donor_group) 
# # Test reveals there is a significant diff in group variance
# 

#-------------------------------------------------------------------------
#               LME and LM models 
#-------------------------------------------------------------------------

lm.PdPc <- function(metadf, metric){
  ###' Function conducts Linear Model for PD vs PC
  env.PdPc <- filter(metadf, donor_group != "HC")
  formula <- as.formula(
    paste(metric, "~", paste(c("description", "host_age_factor", "host_body_mass_index", "sex"), collapse="+") ) )
  
  linear.model <- lm(formula, data=env.PdPc, na.action = na.omit)
  # print(plot_model(linear.model, show.values = TRUE, value.offset = .3))
  # plot(linear.model)
  # dev.off()
  return(linear.model)
}

lmm.PdHc <- function(metadf, metric){
  ###' Function conducts Linear Mixed Model for PD vs HC
  env.PdHc <- filter(metadf, Paired.plot < 30)
  env.PdHc$description <- factor(env.PdHc$description, levels = c("PD Patient", "Household Control"))
  # env.PdHc$description <- relevel(env.PdHc$description, ref = "Household Control")
  # formula and model
  formula <- as.formula(paste(metric, "~", paste(c("description"))))
  lmm <- lme(formula, random= ~ 1 | Paired, data=env.PdHc, na.action = na.omit)
  qqnorm(resid(lmm))
  qqline(resid(lmm))
  print(plot(ranef(lmm))) # Random Effect Plot
  print(plot(lmm)) # Plot Model
  return(lmm)
}

#-------------------------------------------------------------------------