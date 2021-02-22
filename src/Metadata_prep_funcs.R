# Metadata_Prep_Functions


#----------------------------------------------------------------------------
##-----           Aggregate testable metadata into genres           ---------
#----------------------------------------------------------------------------


Anthro <- c("sex", "host_age", "host_body_mass_index")
GI <- c("bristol_stool_scale", "bowl_movments_per_day", 
        "GI_symptom_severity")
smoke <-
  c("cigarette_smoked_ever",
    "time_since_last_cigarette_smoke",
    "parents_smoked")
nonstarters <- c("antibiotics", "laxatives")
allergy <- c("non_food_allergy", "food_allergy")
general_drugs <-
  c(
    "antihistamines",
    "suppositories",
    "beta_blockers",
    "statins",
    "proton_pump_inhibitors",
    "amantadine",
    "levothyroxine",
    "tricyclic_antidepressants",
    "ssri_antidepressants",
    "platelet_aggregation_inhibitor",
    "oral_contraceptives",
    "metformin",
    "inhaler",
    "nonsteroidal_anti_inflammatory"
  )

pd_drugs  <-
  c(
    "MAO_B_inhibitors",
    "carbidopa",
    "levodopa",
    "levodopa_daily_doseage",
    "ropinirole",
    "mirapex",
    "dopamine_agonists",
    "rasagiline",
    "selegiline"
  )

supplements  <-
  c(
    "calcium",
    "magnesium",
    "multivitamin",
    "vitamin_C",
    "vitamin_D",
    "vitamin_E",
    "vitamin_B12",
    "vitamin_B9",
    "coenzyme_Q",
    "creatine",
    "biotin",
    "fish_oil",
    "Alpha_lipoc_acid",
    "resveratrol",
    "encapsulated_herbs",
    "encapsulated_spices"
  )
diet  <-
  c(
    "breakfast_frequency",
    "meat_consumption_frequency",
    "fish_consumption_frequency",
    # "egg_consumption_frequency",
    "bread_consumption_frequency",
    "fruit_consumption_frequency",
    "vegetable_consumption_frequency",
    "fermented_product_consumption_frequency",
    "caffeine_consumption_frequency",
    "single_serving_of_milk_cheese_consumption_frequency",
    "dairy_substitute_consumption_frequency",
    "salted_snack_consumption_frequency",
    "sweets_consumption_frequency",
    "chocolate_preferance",
    "chocolate_consumption_frequency",
    "sweetened_softdrink_consumption_frequency",
    "wine_consumption_frequency",
    "hard_alcohol_consumption_frequency",
    "frequency_of_four_cups_water_daily"
  )

grouping <- c("description", "PD")
environmental <- c("cohort", "paired")
other_metadata <- c("donor_id", "donor_group")

mdsupdrs_III <- 
  c("mds_updrs_3_arising_from_chair",
    "mds_updrs_3_constancy_of_rest",
    "mds_updrs_3_facial_expression",
    "mds_updrs_3_finger_tapping_left_hand",
    "mds_updrs_3_finger_tapping_right_hand",
    "mds_updrs_3_freezing_gait",
    "mds_updrs_3_gait",
    "mds_updrs_3_global_spontaneity_of_movement",
    "mds_updrs_3_hand_movements_left_hand",
    "mds_updrs_3_hand_movements_right_hand",
    "mds_updrs_3_kinetic_tremor_left_hand",
    "mds_updrs_3_kinetic_tremor_right_hand",
    "mds_updrs_3_leg_agility_left_leg",
    "mds_updrs_3_leg_agility_right_leg",
    "mds_updrs_3_postural_stability",
    "mds_updrs_3_postural_tremor_left_hand",
    "mds_updrs_3_postural_tremor_right_hand",
    "mds_updrs_3_posture",
    "mds_updrs_3_pronation_supination_movements_left_hand",
    "mds_updrs_3_pronation_supination_movements_right_hand",
    "mds_updrs_3_rest_tremor_amplitude_lipjaw",
    "mds_updrs_3_rest_tremor_amplitude_lle",
    "mds_updrs_3_rest_tremor_amplitude_lue",
    "mds_updrs_3_rest_tremor_amplitude_rle",
    "mds_updrs_3_rest_tremor_amplitude_rue",
    "mds_updrs_3_rigidity_lle",
    "mds_updrs_3_rigidity_lue",
    "mds_updrs_3_rigidity_neck",
    "mds_updrs_3_rigidity_rle",
    "mds_updrs_3_rigidity_rue",
    "mds_updrs_3_speech",
    "mds_updrs_3_toe_tapping_left_foot",
    "mds_updrs_3_toe_tapping_right_foot",
    "mds_updrs_3_total")

updrs <- 
  c("converted_updrs_score",
    "updrs_action_tremor_left_upper_extremity",
    "updrs_action_tremor_right_upper_extremity",
    "updrs_arising_from_chair",
    "updrs_body_bradykinesia",
    "updrs_facial_expression",
    "updrs_finger_taps_left_upper_extremity",
    "updrs_finger_taps_right_upper_extremity",
    "updrs_gait",
    "updrs_hand_movements_left_upper_extremity",
    "updrs_hand_movements_right_upper_extremity",
    "updrs_heel_taps_left_lower_extremity",
    "updrs_heel_taps_right_lower_extremity",
    "updrs_postural_stability",
    "updrs_posture",
    "updrs_prosup_left_upper_extremity",
    "updrs_prosup_right_upper_extremity",
    "updrs_rest_tremor_hn",
    "updrs_rest_tremor_left_lower_extremity",
    "updrs_rest_tremor_left_upper_extremity",
    "updrs_rest_tremor_right_lower_extremity",
    "updrs_rest_tremor_right_upper_extremity",
    "updrs_rigidity_hn",
    "updrs_rigidity_left_lower_extremity",
    "updrs_rigidity_left_upper_extremity",
    "updrs_rigidity_right_lower_extremity",
    "updrs_rigidity_right_upper_extremity",
    "updrs_speech",
    "updrs_total")

mdsupdrs_I.II <- 
  c("mds_updrs_survey_total",
    "sleep_problems",
    "daytime_sleepiness",
    "pain_and_other_sensations",
    "urinary_problems.",
    "constipation_problems",
    "light_headedness_on_standing.",
    "fatigue",
    "speech",
    "saliva_and_drooling.",
    "chewing_and_swallowing",
    "eating_tasks",
    "dressing",
    "hygiene",
    "handwriting",
    "doing_hobbies_and_other_activities",
    "turning_in_bed",
    "tremor",
    "getting_out_of_bed_a_car_or_a_deep_chair",
    "walking_and_balance",
    "freezing")

motor_severity_scores <- 
  c(updrs,
    mdsupdrs_I.II,
    mdsupdrs_III)

motor_severity_scores_summary <- 
  c("updrs_total",
    "mds_updrs_survey_total",
    "mds_updrs_3_total",
    "dyskinesia",
    "hy_stage")

clinical_variables <- 
  c("family_history_pd_degree_relative",
    "olfactory_diagnosis",
    "smell_score_upsit",
    "thought_disorder_scale",
    "year_of_diagnosis",
    "year_of_onset",
    "disease_duration")

summary_metadata <- c(
  other_metadata,
  grouping,
  Anthro,
  GI,
  nonstarters,
  allergy,
  smoke,
  general_drugs,
  pd_drugs,
  motor_severity_scores_summary,
  clinical_variables,
  environmental
)


#----------------------------------------------------------------------------
##--------------           METADATA PREP FUNCTION           --------------
#----------------------------------------------------------------------------

process_meta <- function(dat, cohort, verbose = F){
  
  if (cohort == "TBC"){

    testable_meta <- c(
      other_metadata,
      grouping,
      Anthro,
      GI,
      nonstarters,
      allergy,
      smoke,
      general_drugs,
      pd_drugs,
      motor_severity_scores,
      mdsupdrs_I.II,
      supplements,
      diet,
      environmental
      )
    
    ####### Pull metadata list from Phyloseq Obj
    env <- get_variable(dat, testable_meta)
    # Make all data characters
    env[] <- lapply(env, as.character)
    # Replace "not provided" entries as NA
    env[env == "not provided" ] <- NA
    env[env == "not collected" ] <- NA
    
    #---------------------------------------------
    #-        Anthropomorphic variables  
    #---------------------------------------------
    env$host_age <-  as.numeric(env$host_age)
    env <- env %>%
      mutate(host_age_factor = as.numeric(
        cut(
          env$host_age,
          breaks = quantile(env$host_age, na.rm = T),
          labels = c(1, 2, 3, 4),
          include.lowest = T
        )
      ))
    # Replace # host_age for host_age_factor in anthro
    Anthro <- c("sex", "host_age_factor", "host_body_mass_index")
    # #transform BMI into discrete groups 
    env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
    env$host_body_mass_index[which(env$host_body_mass_index < 18.5)] <- 0
    env$host_body_mass_index[which(env$host_body_mass_index > 35 )] <- 4
    env$host_body_mass_index[which(env$host_body_mass_index < 35 & env$host_body_mass_index >= 30)] <- 3
    env$host_body_mass_index[which(env$host_body_mass_index < 30  & env$host_body_mass_index >= 25)] <- 2
    env$host_body_mass_index[which(env$host_body_mass_index < 25 & env$host_body_mass_index >= 18.5)] <- 1
    
    env <- mutate(env, bowl_movments_per_day = str_replace(bowl_movments_per_day, "Less than one", "0"))
    env$bowl_movments_per_day <- as.numeric(env$bowl_movments_per_day)

    # Cigarette Smoke Qs
    env[env == "not applicable"] <- 0
    env[env == "5 or more years ago"] <- 1
    env[env == "not in the past 6 months but sometime in the past year"] <- 2
    env[env == "not in the past 30 days but sometime in the past 6 months"] <- 3
    env[env == "not in the past 7 days but sometime in the past 30 days"] <- 4
    env[env == "not today but sometime in the past 7 days"] <- 5
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="No" | env$cigarette_smoked_ever=="Never smoked a cigarette before")] <- 0
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever== "Yes, tried at least one cigarette before")] <- 1
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="Yes")] <- 2
    env$time_since_last_cigarette_smoke <- as.numeric(env$time_since_last_cigarette_smoke )
    env$cigarette_smoked_ever <- as.numeric(env$cigarette_smoked_ever)
    
    ## Dietary/Frequency Questions
    env[env == "Never"] <- 0
    env[env == "Rarely"] <- 1
    env[env == "Less than half the time"] <-2
    env[env == "Not in the past 7 days"] <- 3
    env[env == "Almost every day"] <- 4
    env[env == "Every day"] <- 5;   env[env == "every day"] <- 5
    
    env$chocolate_preferance[which(grepl("Milk, Dark", env$chocolate_preferance))] <- 1
    env$chocolate_preferance[which(grepl("Dark", env$chocolate_preferance))] <- 2
    env$chocolate_preferance[which(grepl("Milk", env$chocolate_preferance))] <- 0
  
    # Make Diet/motor score metadata numeric  
    env <- 
      env %>% 
      mutate(across(all_of(diet), as.character),
             across(all_of(diet), as.numeric)) %>% 
      mutate(across(all_of(mdsupdrs_I.II), as.character),
             across(all_of(mdsupdrs_I.II), as.numeric))
    
    #  Make all character variables into Factors
    env <- env %>% mutate_if(is.character, as.factor)
    
    if (verbose == T){
      print(str(env))
    }
    cat ("\n Metadata Processing complete \n\n")
    
  } else if (cohort == "RUSH") {
    
    testable_meta <- c(
      grouping,
      Anthro,
      other_metadata)
    
    ####### Pull metadata list from Phyloseq Obj
    env <- get_variable(dat, testable_meta)
    # Make all data characters
    env[] <- lapply(env, as.character)
    # Replace "not provided" entries as NA
    env[env == "not provided" ] <- NA
    env[env == "not collected" ] <- NA
    
    #---------------------------------------------
    #-        Anthropomorphic variables  
    #---------------------------------------------
    env$host_age <- as.numeric(env$host_age)
    env <- env %>% 
      mutate(host_age_factor = as.numeric(
        cut(
          env$host_age,
          breaks = quantile(env$host_age, na.rm = T),
          labels = c(1, 2, 3, 4),
          include.lowest = T
        )
      ))
    
    # Replace # host_age for host_age_factor in anthro
    Anthro <- c("sex", "host_age_factor", "host_body_mass_index") 
    # #transform BMI into discrete groups 
    env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
    env$host_body_mass_index[which(env$host_body_mass_index < 18.5)] <- 0
    env$host_body_mass_index[which(env$host_body_mass_index > 35 )] <- 4
    env$host_body_mass_index[which(env$host_body_mass_index < 35 & env$host_body_mass_index >= 30)] <- 3
    env$host_body_mass_index[which(env$host_body_mass_index < 30  & env$host_body_mass_index >= 25)] <- 2
    env$host_body_mass_index[which(env$host_body_mass_index < 25 & env$host_body_mass_index >= 18.5)] <- 1
    
    #  Make all character variables into Factors
    env <- env %>% mutate_if(is.character, as.factor)
    if (verbose == T){
      print(str(env))
    }
    cat ("\n Metadata Processing complete \n\n")
    
  } else if (cohort == "Merged") {
    
    # Uses summary of motor scores
    testable_meta <- c(
      other_metadata,
      grouping,
      Anthro,
      GI,
      nonstarters,
      allergy,
      smoke,
      general_drugs,
      pd_drugs,
      supplements,
      diet,
      motor_severity_scores_summary,
      clinical_variables,
      environmental
    )
    
    ####### Pull metadata list from Phyloseq Obj
    env <- get_variable(dat, testable_meta)
    # Make all data characters
    env[] <- lapply(env, as.character)
    # Replace "not provided" entries as NA
    env[env == "not provided" ] <- NA
    env[env == "not collected" ] <- NA
    
    #---------------------------------------------
    #-        Anthropomorphic variables  
    #---------------------------------------------
    env$host_age <- as.numeric(env$host_age)
    env <- env %>%
      mutate(host_age_factor = as.numeric(
        cut(
          env$host_age,
          breaks = quantile(env$host_age, na.rm = T),
          labels = c(1, 2, 3, 4),
          include.lowest = T
        )
      ))
    # Replace # host_age for host_age_factor in anthro
    Anthro <- c("sex", "host_age_factor", "host_body_mass_index")
    # #transform BMI into discrete groups 
    env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
    env$host_body_mass_index[which(env$host_body_mass_index < 18.5)] <- 0
    env$host_body_mass_index[which(env$host_body_mass_index > 35 )] <- 4
    env$host_body_mass_index[which(env$host_body_mass_index < 35 & env$host_body_mass_index >= 30)] <- 3
    env$host_body_mass_index[which(env$host_body_mass_index < 30  & env$host_body_mass_index >= 25)] <- 2
    env$host_body_mass_index[which(env$host_body_mass_index < 25 & env$host_body_mass_index >= 18.5)] <- 1
    
    env <- mutate(env, bowl_movments_per_day = str_replace(bowl_movments_per_day, "Less than one", "0"))
    env$bowl_movments_per_day <- as.numeric(env$bowl_movments_per_day)

    # Cigarette Smoke Qs
    env[env == "not applicable"] <- 0
    env[env == "5 or more years ago"] <- 1
    env[env == "not in the past 6 months but sometime in the past year"] <- 2
    env[env == "not in the past 30 days but sometime in the past 6 months"] <- 3
    env[env == "not in the past 7 days but sometime in the past 30 days"] <- 4
    env[env == "not today but sometime in the past 7 days"] <- 5
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="No" | env$cigarette_smoked_ever=="Never smoked a cigarette before")] <- 0
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever== "Yes, tried at least one cigarette before")] <- 1
    env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="Yes")] <- 2
    env$time_since_last_cigarette_smoke <- as.numeric(env$time_since_last_cigarette_smoke )
    env$cigarette_smoked_ever <- as.numeric(env$cigarette_smoked_ever)
    
    ## Dietary/Frequency Questions
    env[env == "Never"] <- 0
    env[env == "Rarely"] <- 1
    env[env == "Less than half the time"] <-2
    env[env == "Not in the past 7 days"] <- 3
    env[env == "Almost every day"] <- 4
    env[env == "Every day"] <- 5;   env[env == "every day"] <- 5
    
    env$chocolate_preferance[which(grepl("Milk, Dark", env$chocolate_preferance))] <- 1
    env$chocolate_preferance[which(grepl("Dark", env$chocolate_preferance))] <- 2
    env$chocolate_preferance[which(grepl("Milk", env$chocolate_preferance))] <- 0
    
    ## olfactory questions
    env[env == "normosmia"] <- 1
    env[env == "mild_microsmia"] <- 2
    env[env == "moderate_microsmia"] <- 3
    env[env == "severe_microsmia"] <- 4
    env[env == "total_anosmia"] <- 5
    env$smell_score_upsit <- as.numeric(env$smell_score_upsit)
    
    # Metadata to convert to numeric  
    env <- 
      env %>% 
      mutate(across(all_of(diet), as.character),
             across(all_of(diet), as.numeric)) %>% 
      mutate(across(all_of(motor_severity_scores_summary), as.character),
             across(all_of(motor_severity_scores_summary), as.numeric)) %>% 
      mutate(across(all_of(clinical_variables), as.character),
             across(all_of(clinical_variables), as.numeric)) 

    #  Make all character variables into Factors
    env <- env %>% mutate_if(is.character, as.factor)
    if (verbose == T){
      print(str(env))
    }
    cat ("\n Metadata Processing complete \n\n")
  }
  
  # load env into global enviornment
  assign("env",env,envir = .GlobalEnv)
  # assign("Anthro", Anthro, envir = .GlobalEnv)
  # assign("GI", GI, envir = .GlobalEnv)
  # assign("smoke", smoke, envir = .GlobalEnv)
  # assign("nonstarters", nonstarters, envir = .GlobalEnv)
  # assign("allergy", allergy, envir = .GlobalEnv)
  # assign("general_drugs", general_drugs, envir = .GlobalEnv)
  # assign("pd_drugs", pd_drugs, envir = .GlobalEnv)
  # assign("supplements", supplements, envir = .GlobalEnv)
  # assign("diet", diet, envir = .GlobalEnv)
  # assign("grouping", grouping, envir = .GlobalEnv)
  # assign("other_metadata", other_metadata, envir = .GlobalEnv)

}

#############################################################################

##################           TRIM METADATA FUNC           ################## 

#############################################################################


trim_meta <- function(env, min_ratio){
  if (typeof(env) != "list"){
    stop("Please input preprocessed metdata list \n RUN process_meta() function first")
  }
  ###  Remove (Yes/No) questions with less than a certain % (min ratio) of samples as minority
  rmvlst_index <- c()
  rmvlst_names <- c()
  cat("List of Metadata below filtering threshold: \n")
  for (i in 1:length(env)) {
    p <- na.omit(env[, i])
    x <- count(p)
    if (length(x$x) == 2) {
      a <- min(x$freq)
      b <- sum(x$freq)
      if ((a / b) < min_ratio) {
        print(colnames(env)[i])
        rmvlst_names <- c(rmvlst_names, colnames(env)[i])
        rmvlst_index <- c(rmvlst_index, i)
      }
      ###  Remove questions with only a single factor response
    } else if (length(x$x) == 1) {
      print(colnames(env)[i])
      rmvlst_names <- c(rmvlst_names, colnames(env)[i])
      rmvlst_index <- c(rmvlst_index, i)
    }
  }
  if (!is.null(rmvlst_index)) {
    env <- env[-rmvlst_index]
  }
  cat("\nMetadata Filtering Complete: \n\n\n")
  assign("env",env,envir = .GlobalEnv)
}





#############################################################################

##################      STUDY DESIGN PLOT METADTA PREP      ################## 

#############################################################################



process_meta_study_design_plot <- function(dat){
  if (typeof(dat) != "S4"){
    stop("Please input a phyloseq Obj")
  }
  
  ####### Aggregate testable metadata into genres
  Anthro <- c("sex", "host_age", "host_body_mass_index")
  GI <- c("bristol_stool_scale", "bowl_movments_per_day")
  smoke <- c("cigarette_smoked_ever", "time_since_last_cigarette_smoke", "parents_smoked")
  nonstarters <- c("antibiotic_use", "laxatives")
  allergy <- c("non_food_allergy", "food_allergy")
  general_drugs <- c("antihistamines", "suppositories", "beta_blockers", "statins", "proton_pump_inhibitors",
                     "tricyclic_antidepressants", "ssri_antidepressants", "platelet_aggregation_inhibitor", 
                     "oral_contraceptives", "metformin", "inhaler", "nonsteroidal_anti.inflammatory")
  pd_drugs  <- c("carbidopa", "levodopa", "levothyroxine", "ropinrole", "azilect", "mirapex",
                 "rytary", "rasagline")
  supplements  <- c("calcium", "magnesium", "multivitamin", "vitamin_C", "vitamin_D",
                    "vitamin_E", "vitamin_B12", "vitamin_B9", "coenzyme_Q", "creatine",
                    "biotin", "fish_oil", "Alpha_lipoc_acid", "resveratrol", "encapsulated_herbs",
                    "encapsulated_spices")
  diet  <- c("breakfast_frequency", "meat_consumption_frequency", "fish_consumption_frequency",
             "egg_consumption_frequency", "bread_consumption_frequency", "fruit_consumption_frequency", 
             "vegetable_consumption_frequency", "fermented_product_consumption_frequency", "caffeine_consumption_frequency", 
             "single_serving_of_milk_cheese_consumption_frequency", "dairy_substitute_consumption_frequency",
             "salted_snack_consumption_frequency", "sweets_consumption_frequency", "chocolate_preferance", 
             "chocolate_consumption_frequency", "sweetened_softdrink_consumption_frequency", 
             "wine_consumption_frequency", "hard_alcohol_consumption_frequency",
             "frequency_of_four_cups_water_daily")
  grouping <- c("description", "PD", "paired")
  
  other_metadata <- c("donor_id", "donor_group")
  
  # alpha 
  testable_meta <- c(other_metadata, grouping , Anthro, GI, nonstarters, allergy, smoke, 
                     general_drugs, pd_drugs, supplements, diet)
  

  
  ####### Pull metadata list from Phlyoseq Obj
  env <- get_variable(dat, testable_meta)
  # Replace "not provided" entries as NA
  env[env == "not provided" ] <- NA
  
  
  
  #### Manual processing of variables into relevant numeric form 
  # Misc corrections
  # env[which(env$paired == "No"),"paired"] <-NA
  # env$paired <- as.numeric(env$paired)
  # env[which(!is.na(env$paired)),"paired"] <-  paste0("pair_", env[which(!is.na(env$paired)),"paired"])
  
  env$host_age <- as.numeric(env$host_age)
  env <- mutate(env, host_age_factor=as.numeric(cut(env$host_age, breaks = quantile(env$host_age, na.rm = T), labels=c(1,2,3,4), include.lowest=T) ))
  Anthro <- c("sex", "host_age_factor", "host_body_mass_index") # Replace # host_age for host_age_factor in anthro
  env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
  # env[env == "Less than one"] <- 0 # to make "bowel mvmts per day numeric
  env$bowl_movments_per_day <- as.numeric(env$bowl_movments_per_day)
  env$antibiotic_use[ which(grepl("Anti", env$antibiotic_use))] <- c("Yes")
  
  # Cigarette Smoke Qs
  env[env == "not applicable"] <- 0
  env[env == "5 or more years ago"] <- 1
  env[env == "not in the past 6 months but sometime in the past year"] <- 2
  env[env == "not in the past 30 days but sometime in the past 6 months"] <- 3
  env[env == "not in the past 7 days but sometime in the past 30 days"] <- 4
  env[env == "not today but sometime in the past 7 days"] <- 5
  env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="No" | env$cigarette_smoked_ever=="Never smoked a cigarette before")] <- 0
  env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever== "Yes, tried at least one cigarette before")] <- 1
  env$cigarette_smoked_ever[  which(env$cigarette_smoked_ever=="Yes")] <- 2
  env$time_since_last_cigarette_smoke <- as.numeric(env$time_since_last_cigarette_smoke )
  env$cigarette_smoked_ever <- as.numeric(env$cigarette_smoked_ever)
  
  # ## Dietary/Frequency Questions
  # env[env == "Never"] <- 0
  # env[env == "Rarely"] <- 1
  # env[env == "Less than half the time"] <-2
  # env[env == "Not in the past 7 days"] <- 3
  # env[env == "Almost every day"] <- 4
  # env[env == "Every day"] <- 5
  
  env$chocolate_preferance[which(grepl("Milk, Dark", env$chocolate_preferance))] <- 1
  env$chocolate_preferance[which(grepl("Dark", env$chocolate_preferance))] <- 2
  env$chocolate_preferance[which(grepl("Milk", env$chocolate_preferance))] <- 0
  
  #transform BMI into discrete groups  - No need for this actually
  env$host_body_mass_index[which(env$host_body_mass_index < 18.5)] <- 0
  env$host_body_mass_index[which(env$host_body_mass_index < 35 & env$host_body_mass_index >= 30)] <- 3
  env$host_body_mass_index[which(env$host_body_mass_index < 30  & env$host_body_mass_index >= 25)] <- 2
  env$host_body_mass_index[which(env$host_body_mass_index < 25 & env$host_body_mass_index >= 18.5)] <- 1
  
  #  Make all character variables into Factors
  for (i in 1:length(env)){
    if (typeof(env[,i]) == "character") {
      env[,i] <- as.factor(env[,i])
    }
  }
  
  # # Make Diet metadata numeric  
  # j <-  which(colnames(env) == diet[1])
  # for(i in 1:length(diet)) {
  #   # print(j)
  #   # print(typeof(env[,j]))
  #   env[,j] <- as.numeric(env[,j])
  #   # print(typeof(env[,j]))
  #   j = j + 1
  #   # print(j)
  # }
  
  assign("env",env,envir = .GlobalEnv)
  assign("Anthro", Anthro, envir = .GlobalEnv)
  assign("GI", GI, envir = .GlobalEnv)
  assign("smoke", smoke, envir = .GlobalEnv)
  assign("nonstarters", nonstarters, envir = .GlobalEnv)
  assign("allergy", allergy, envir = .GlobalEnv)
  assign("general_drugs", general_drugs, envir = .GlobalEnv)
  assign("pd_drugs", pd_drugs, envir = .GlobalEnv)
  assign("supplements", supplements, envir = .GlobalEnv)
  assign("diet", diet, envir = .GlobalEnv)
  assign("grouping", grouping, envir = .GlobalEnv)
  assign("other_metadata", other_metadata, envir = .GlobalEnv)
  
  cat ("Processing complete : \n\n")
}


#----------------------------------------------------------------------------
##-----                TRIM METADATA BY NAs FUNC               ---------
#----------------------------------------------------------------------------

trim_meta_by_NA <- function(env){
  if (typeof(env) != "list"){
    stop("Please input preprocessed metdata list \n RUN process_meta() function first")
  }
  ###'  Remove questions with less than 70% of samples answered
  rmvlst_index <- c()
  rmvlst_names <- c()
  cat("List of Metadata below filtering threshold: \n")
  for (i in 1:length(env)){
    p <- na.omit(env[,i])
    x <- length(p)
    if (x < 82) {
      print(colnames(env)[i])
      rmvlst_names <- c(rmvlst_names, colnames(env)[i])
      rmvlst_index <- c(rmvlst_index, i)
    }
  }
  cat("\nMetadata Filtering Complete: \n\n\n")
  env <- env[-rmvlst_index]
  assign("env",env,envir = .GlobalEnv)
}



