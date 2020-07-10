# Metadata_Prep_Functions

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot)


# Load Phlyoseq Obj with Metdata
load("files/Species_PhyloseqObj.RData")

#############################################################################

##################         PERMANOVA PREP FUNCTION         ################## 

#############################################################################

process_meta <- function(dat){
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
  grouping <- c("description", "PD", "diagnosed_with_mental_illness", "Paired")
  
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
  env[which(env$Paired == "No"),"Paired"] <- 100:159
  env$host_age <- as.numeric(env$host_age)
  env <- mutate(env, host_age_factor=as.numeric(cut(env$host_age, breaks = quantile(env$host_age, na.rm = T), labels=c(1,2,3,4), include.lowest=T) ))
  Anthro <- c("sex", "host_age_factor", "host_body_mass_index") # Replace # host_age for host_age_factor in anthro
  env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
  env[env == "Less than one"] <- 0 # to make "bowel mvmts per day numeric
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
  
  ## Dietary/Frequency Questions
  env[env == "Never"] <- 0
  env[env == "Rarely"] <- 1
  env[env == "Less than half the time"] <-2
  env[env == "Not in the past 7 days"] <- 3
  env[env == "Almost every day"] <- 4
  env[env == "Every day"] <- 5
  
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
  
  # Make Diet metadata numeric  
  j <-  which(colnames(env) == diet[1])
  for(i in 1:length(diet)) {
    # print(j)
    # print(typeof(env[,j]))
    env[,j] <- as.numeric(env[,j])
    # print(typeof(env[,j]))
    j = j + 1
    # print(j)
  }

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

  cat ("Metadata Processing complete :-D \n\n")
}

#############################################################################

##################           TRIM METADATA FUNC           ################## 

#############################################################################


trim_meta <- function(env){
  if (typeof(env) != "list"){
    stop("Please input preprocessed metdata list \n RUN process_meta() function first")
  }
  ###  Remove (Yes/No) questions with less than 10% of samples as minority
  #initalize vectors
  rmvlst_index <- c()
  rmvlst_names <- c()
  cat("List of Metadata below filtering threshold: \n")
  for (i in 1:length(env)){
    p <- na.omit(env[,i])
    x <- count(p)
    if (length(x$x) == 2) {
      a <- min(x$freq)
      b <- sum(x$freq)
      if ((a/b) < 0.1) {
        print(colnames(env)[i])
        rmvlst_names <- c(rmvlst_names, colnames(env)[i])
        rmvlst_index <- c(rmvlst_index, i)
      }
    }
  }
  cat("\nMetadata Filtering Complete: \n\n\n")
  env <- env[-rmvlst_index]
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
  grouping <- c("description", "PD", "diagnosed_with_mental_illness", "Paired")
  
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
  # env[which(env$Paired == "No"),"Paired"] <-NA
  # env$Paired <- as.numeric(env$Paired)
  # env[which(!is.na(env$Paired)),"Paired"] <-  paste0("pair_", env[which(!is.na(env$Paired)),"Paired"])
  
  env$host_age <- as.numeric(env$host_age)
  env <- mutate(env, host_age_factor=as.numeric(cut(env$host_age, breaks = quantile(env$host_age, na.rm = T), labels=c(1,2,3,4), include.lowest=T) ))
  Anthro <- c("sex", "host_age_factor", "host_body_mass_index") # Replace # host_age for host_age_factor in anthro
  env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
  env[env == "Less than one"] <- 0 # to make "bowel mvmts per day numeric
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


#############################################################################

##################       TRIM METADATA BY NAs FUNC       ################## 

#############################################################################



trim_meta_by_NA <- function(env){
  if (typeof(env) != "list"){
    stop("Please input preprocessed metdata list \n RUN process_meta() function first")
  }
  ###'  Remove questions with less than 70% of samples answered
  #initalize vectors
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

