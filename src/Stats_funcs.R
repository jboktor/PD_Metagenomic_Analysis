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
#               LMER and LM models 
#-------------------------------------------------------------------------

lm.PdPc <- function(metadf, metric, cohort){
  ###' Function conducts Linear Model for PD vs PC
  env.PdPc <- filter(metadf, donor_group != "HC")
  if (cohort == "Merged"){
    formula <- as.formula(
      paste(metric, "~", paste(c("description", "host_age_factor", "host_body_mass_index", "sex", "(1|cohort)"), collapse = "+")))
    model <- lmer(formula, data = metadf, REML = T, na.action = na.exclude)
  } else {
    formula <- as.formula(
      paste(metric, "~", paste(c("description", "host_age_factor", "host_body_mass_index", "sex"), collapse="+") ) )
    model <- lm(formula, data=metadf, na.action = na.omit)
  }
  # print(plot_model(linear.model, show.values = TRUE, value.offset = .3))
  # print(plot(model))
  # dev.off()
  return(model)
}

lmm.PdHc <- function(metadf, metric, cohort){
  ###' Function conducts Linear Mixed Model for PD vs HC
  metadf$description <- factor(metadf$description, levels = c("PD Patient", "Household Control"))
  
  if (cohort == "Merged"){
    formula <- as.formula(
      paste(metric, "~", paste(c("description", "(1|cohort:paired)"), collapse = "+")))
    model <- lmer(formula, data = metadf, REML = T, na.action = na.exclude)
  } else {  
  formula <- as.formula(
    paste(metric, "~", paste(c("description", "(1|paired)"), collapse = "+")))
  model <- lmer(formula, data = metadf, REML = T, na.action = na.exclude)
  }
  # qqnorm(resid(model))
  # qqline(resid(model))
  # print(plot(ranef(model))) # Random Effect Plot
  # print(plot(model)) # Plot Model
  return(model)
}

#-------------------------------------------------------------------------