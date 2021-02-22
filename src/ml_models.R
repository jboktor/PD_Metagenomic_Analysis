# Machine-Learning-Models

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")

# Packages for xgboost
# library(glue)
library(ModelMetrics)
library(OpenMPController) # for Kaggle backend
library(xgboost)
library(parallel)
library(doParallel)
library(mlbench)
library(caret)
library(MLeval)


#----------------------------------------------------------------------------------------
###                       RIDGE, LASSO, ENET LOGISTIC REGRESSION                     ###
#----------------------------------------------------------------------------------------


ridge.lasso.enet.regression.model <- function(obj, comparison = "PDvPC", model.type){
  
  #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
  #' Input: Phyloseq Obj, comparison of interest, and model type
  #' Returns: a list including 
  #' 1) The fitted Model 
  #' 2) A dataframe of values with the optimal parameters (for ROC plot)
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  # Initalize variables
  model <- NULL
  model.input <- NULL
  output.list <- vector(mode="list", length=4)
  names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
  
  if (model.type == "ridge"){
    tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
  } else if (model.type == "lasso"){
    tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
  } else if (model.type == "enet"){
    tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
  }
  
  # Select samples and prep abundance DF 
  if (comparison == "PDvPC"){
    dat_pdpc = subset_samples(obj, donor_group !="HC")
    d <- dat_pdpc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame() 
    d <- asin(sqrt(d))
    pdPC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdPC.df) <- rownames(d)
    pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
    model.input <- pdPC.df
    
  } else if (comparison == "PDvHC"){
    dat_pdhc = subset_samples(obj, Paired !="No")
    d <- dat_pdhc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame() 
    d <- asin(sqrt(d))
    pdHC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdHC.df) <- rownames(d)
    pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
    model.input <- pdHC.df
  }
  
  # Model Parameters
  numbers <- 10
  repeats <- 5  
  set.seed(42)
  seed <- 42
  rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, seed = seed)
  
  mytrainControl <- 
    trainControl(method='repeatedcv',
                 number=numbers, 
                 repeats=repeats,
                 search='grid',
                 seeds = rcvSeeds,
                 savePredictions = TRUE, 
                 classProbs = TRUE, 
                 verboseIter = TRUE)
  
  # Run 10-fold CV Model with 5 Repitions
  model <-train(group ~.,
                data=model.input, 
                method='glmnet', 
                metric='Accuracy', 
                tuneGrid=tune.grid, 
                trControl=mytrainControl,
                verboseIter = T)
  cat("Model Summary")
  print(model)
  cat("\n\n")
  print(plot(model))
  
  # Select model with optimal Lambda
  selectedIndices <- model$pred$lambda == model$bestTune$lambda
  df.output <- model$pred[selectedIndices, ]
  
  # MLeval
  mleval <- evalm(model)
  mleval$roc
  output.list$fitted_model <- model
  output.list$optimal.df <- df.output
  output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
  output.list$MLevaldata <- mleval$stdres
  
  
  return(output.list)
  
}

#----------------------------------------------------------------------------------------
###                              RANDOM FOREST MODELS                                ###
#----------------------------------------------------------------------------------------



random.forest.model <- function(obj, comparison = "PDvPC"){
  
  #' Function to train a Random Forrest Model
  #' Input: Phyloseq Obj and comparison of interest 
  #' Returns: a list including 
  #' 1) The fitted Model 
  #' 2) A dataframe of values with the optimal parameters (for ROC plot)
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  intervalStart <- Sys.time()
  
  # Initalize variables
  model <- NULL
  output.list <- vector(mode="list", length=4)
  names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")

  # Select samples and prep abundance DF 
  if (comparison == "PDvPC"){
    dat_pdpc = subset_samples(obj, donor_group !="HC")
    d <- dat_pdpc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame()
    d <- asin(sqrt(d))
    pdPC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdPC.df) <- rownames(d)
    pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
    model.input <- pdPC.df
    
  } else if (comparison == "PDvHC"){
    dat_pdhc = subset_samples(obj, Paired !="No")
    d <- dat_pdhc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame() 
    d <- asin(sqrt(d))
    pdHC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdHC.df) <- rownames(d)
    pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
    model.input <- pdHC.df
  }
  
  # Model Parameters
  numbers <- 10
  repeats <- 5
  tune.grid = expand.grid(.mtry =   seq(1, 2 * as.integer(sqrt(ncol(model.input) - 1)), by=2)) 
  set.seed(42)
  seed <- 42
  # Repeated cross validation
  rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, 
                       tunes = length(tune.grid$.mtry), seed = seed)
  # c('B + 1' = length(rcvSeeds), M = length(rcvSeeds[[1]]))
  # rcvSeeds[c(1, length(rcvSeeds))]
  mytrainControl <- 
    trainControl(method='repeatedcv',
                 number=numbers, 
                 repeats=repeats,
                 search='grid',
                 savePredictions = TRUE, 
                 classProbs = TRUE, 
                 verboseIter = TRUE,
                 allowParallel = TRUE,
                 seeds = rcvSeeds)
  
  # Run model using multiple threads 
  cluster <- parallel::makeCluster(detectCores() - 1, setup_strategy = "sequential")
  registerDoParallel(cluster)
  
  # Run 10-fold CV Model with 5 Repitions
  system.time(
    model <-train(group ~.,
                data=model.input, 
                method='rf', 
                metric='Accuracy', 
                tuneGrid=tune.grid, 
                ntree = 1000,
                trControl=mytrainControl,
                importance = TRUE,
                verboseIter = T)
  )
  
  stopCluster(cluster)
  unregister()

  cat("Model Summary")
  print(model)
  cat("\n\n")
  print(plot(model))
  
  # Select model with optimal __mtry__ 
  selectedIndices <- model$pred$mtry == model$bestTune$mtry
  df.output <- model$pred[selectedIndices, ]
  
  # MLeval
  mleval <- evalm(model)
  mleval$roc
  output.list$fitted_model <- model
  output.list$optimal.df <- df.output
  output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
  output.list$MLevaldata <- mleval$stdres
  
  intervalEnd <- Sys.time()
  cat("Random Forest analysis took",
        intervalEnd - intervalStart,attr(intervalEnd - intervalStart,"units"))
  cat("\n\n")
  return(output.list)
  
}




#----------------------------------------------------------------------------------------
###                                XBOOST MODELS                                      ###
#----------------------------------------------------------------------------------------


xgboost.model <- function(obj, comparison = "PDvPC"){
  
  
  # ########    DELETE LATER  4 - TESTING ######## 
  # obj = dat
  # comparison = "PDvPC"
  # ########  ########  ########  ########  ######## 
  # 
  
  #' Function to train a Random Forrest Model
  #' Input: Phyloseq Obj and comparison of interest 
  #' Returns: a list including 
  #' 1) A dataframe of values with the optimal parameters (for ROC plot)
  #' 2) MLeval individual model ROC plot
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  intervalStart <- Sys.time()
  
  # Initalize variables
  model <- NULL
  output.list <- vector(mode="list", length=3)
  names(output.list) <- c("optimal.df", "AUCROC", "MLevaldata")

  set.seed(42)
  
  # Select samples and prep abundance DF 
  if (comparison == "PDvPC"){
    dat_pdpc = subset_samples(obj, donor_group !="HC")
    d <- dat_pdpc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame() 
    d <- asin(sqrt(d))
    pdPC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdPC.df) <- rownames(d)
    pdPC.df$group <- factor(pdPC.df$group, levels = c("PC", "PD"))
    model.input <- pdPC.df
    
  } else if (comparison == "PDvHC"){
    dat_pdhc = subset_samples(obj, Paired !="No")
    d <- dat_pdhc %>%
      microbiome::transform("compositional") %>% 
      microbiome::abundances() %>% 
      t() %>% 
      as.data.frame() 
    d <- asin(sqrt(d))
    pdHC.df <- group_col_from_ids(d, id= rownames(d))
    rownames(pdHC.df) <- rownames(d)
    pdHC.df$group <- factor(pdHC.df$group, levels = c("HC", "PD"))
    model.input <- pdHC.df
  }
  
  # Set-up Parallel Processing Threads
  # omp_set_num_threads(3)
  # intervalStart <- Sys.time()
  
  #----------------------------------- 
  # Hyperparameter Tuning 
  #Informed by: https://www.kaggle.com/pelkoja/visual-xgboost-tuning-with-caret/report
  
  #---------------------------------------------------------- 
  # 1) nrounds & eta 
  #---------------------------------------------------------- 

  tune_grid <- expand.grid(
    nrounds = seq(from = 200, to = 1000, by = 50),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 3, 4, 5, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  
  tune_control <-
    trainControl(method='repeatedcv',
                 number=5,
                 repeats=3,
                 search='grid',
                 # savePredictions = TRUE,
                 # classProbs = TRUE,
                 # allowParallel = TRUE,
                 verboseIter = TRUE)
  
  xgb_tune <- caret::train(
    group ~.,
    data=model.input,
    trControl = tune_control,
    tuneGrid = tune_grid,
    method = "xgbTree",
    verbose = TRUE
  )
  
  cat("Tuning Step 1: Maximum Depth, learning rate, and nrounds baseline \n COMPLETE \n")
  print(tuneplot(xgb_tune))
  print(xgb_tune$bestTune)
  
  #---------------------------------------------------------- 
  # 2) Maximum Depth and Minimum Child Weight
  #---------------------------------------------------------- 
  
  if (xgb_tune$bestTune$max_depth == 2) {
    mxdpth <- 
      c(xgb_tune$bestTune$max_depth:4)
  }  else {
    mxdpth <-  
      c((xgb_tune$bestTune$max_depth - 1):(xgb_tune$bestTune$max_depth + 1))
  }

  tune_grid2 <- expand.grid(
    nrounds = seq(from = 50, to = 1000, by = 50),
    eta = xgb_tune$bestTune$eta,
    max_depth = mxdpth,
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = c(1, 2, 3),
    subsample = 1
  )
  
  xgb_tune2 <- caret::train(
    group ~.,
    data=model.input,
    trControl = tune_control,
    tuneGrid = tune_grid2,
    method = "xgbTree",
    verbose = TRUE
  )
  cat("Tuning Step 2: Maximum Depth and Minimum Child Weight \n COMPLETE \n")
  print(tuneplot(xgb_tune2))
  print(xgb_tune2$bestTune)
  
  #---------------------------------------------------------- 
  # 3)  Column and Row Sampling
  #---------------------------------------------------------- 
  
  tune_grid3 <- expand.grid(
    nrounds = seq(from = 50, to = 1000, by = 50),
    eta = xgb_tune$bestTune$eta,
    max_depth = xgb_tune2$bestTune$max_depth,
    gamma = 0,
    colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
    min_child_weight = xgb_tune2$bestTune$min_child_weight,
    subsample = c(0.5, 0.75, 1.0)
  )
  
  xgb_tune3 <- caret::train(
    group ~.,
    data=model.input,
    trControl = tune_control,
    tuneGrid = tune_grid3,
    method = "xgbTree",
    verbose = TRUE
  )
  
  cat("Tuning Step 3: Column and Row Sampling \n COMPLETE \n")
  print(tuneplot(xgb_tune3, probs = .95))
  print(xgb_tune3$bestTune)
  
  #---------------------------------------------------------- 
  # 4)  Gamma
  #---------------------------------------------------------- 
  
  tune_grid4 <- expand.grid(
    nrounds = seq(from = 50, to = 1000, by = 50),
    eta = xgb_tune$bestTune$eta,
    max_depth = xgb_tune2$bestTune$max_depth,
    gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
    colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
    min_child_weight = xgb_tune2$bestTune$min_child_weight,
    subsample = xgb_tune3$bestTune$subsample
  )
  
  xgb_tune4 <- caret::train(
    group ~.,
    data=model.input,
    trControl = tune_control,
    tuneGrid = tune_grid4,
    method = "xgbTree",
    verbose = TRUE
  )
  
  cat("Tuning Step 4: Gamma \n COMPLETE \n")
  print(tuneplot(xgb_tune4))
  print(xgb_tune4$bestTune)
  
  #---------------------------------------------------------- 
  # 5)  Reducing the Learning Rate
  #---------------------------------------------------------- 
  
  tune_grid5 <- expand.grid(
    nrounds = seq(from = 100, to = 1000, by = 100),
    eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
    max_depth = xgb_tune2$bestTune$max_depth,
    gamma = xgb_tune4$bestTune$gamma,
    colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
    min_child_weight = xgb_tune2$bestTune$min_child_weight,
    subsample = xgb_tune3$bestTune$subsample
  )
  
  xgb_tune5 <- caret::train(
    group ~.,
    data=model.input,
    trControl = tune_control,
    tuneGrid = tune_grid5,
    method = "xgbTree",
    verbose = TRUE
  )
  
  cat("Tuning Step 5: Final Learning Rate \n COMPLETE \n")
  print(tuneplot(xgb_tune5))
  print(xgb_tune5$bestTune)
  
  
  #---------------------------------------------------------- 
  # 5)  Fitting the Model
  #---------------------------------------------------------- 
  
  final_grid <- expand.grid(
    nrounds = xgb_tune5$bestTune$nrounds,
    eta = xgb_tune5$bestTune$eta,
    max_depth = xgb_tune5$bestTune$max_depth,
    gamma = xgb_tune5$bestTune$gamma,
    colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
    min_child_weight = xgb_tune5$bestTune$min_child_weight,
    subsample = xgb_tune5$bestTune$subsample
  )
  
  mytrainControl <-
    trainControl(method='repeatedcv',
                 number=5,
                 repeats=5,
                 search='grid',
                 savePredictions = TRUE,
                 classProbs = TRUE,
                 allowParallel = TRUE,
                 verboseIter = TRUE)
  

  # Run 5-fold CV Model with 5 Repetition
  model <-train(group ~.,
                data=model.input, 
                method='xgbTree', 
                metric='Accuracy', 
                tuneGrid=final_grid, 
                trControl=mytrainControl,
                verbose = TRUE)

  
  cat("Model Summary")
  print(model)
  cat("\n\n")
  # print(plot(model))
  
  intervalEnd <- Sys.time()
  cat("XGBoost model tuning completed in: ",
      intervalEnd - intervalStart, attr(intervalEnd - intervalStart, "units"))
  cat("\n\n")
  
  # Select model with optimal hyper-parameters
  df.output <- model$pred 

  # MLeval
  mleval <- evalm(model)
  mleval$roc
  output.list$optimal.df <- df.output
  # output.list$AUCROC <- mleval$roc
  output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
  output.list$MLevaldata <- mleval$stdres


  return(output.list)
  
}



#----------------------------------------------------------------------------------------
###                       ARTIFICAL NEURAL NETWORK MODELS                             ###
#----------------------------------------------------------------------------------------








#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
###         RIDGE, LASSO, ENET LOGISTIC REGRESSION  --  For Curated Datasets            
#----------------------------------------------------------------------------------------


ridge.lasso.enet.regression.model.DS <- function(model.input, model.type){
  
  #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
  #' Input: Phyloseq Obj, comparison of interest, and model type
  #' Returns: a list including 
  #' 1) The fitted Model 
  #' 2) A dataframe of values with the optimal parameters (for ROC plot)
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  # Initialize variables
  model <- NULL
  output.list <- vector(mode="list", length=4)
  names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
  
  if (model.type == "ridge"){
    tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
  } else if (model.type == "lasso"){
    tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
  } else if (model.type == "enet"){
    tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
  }
  
  # AST Transformation on merged data
  groupcol <- factor(model.input$group)
  temp <- dplyr::select(model.input, -group)
  model.input <- asin(sqrt(temp))
  model.input$group <- groupcol
  
  # Model Parameters
  numbers <- 10
  repeats <- 5  
  set.seed(42)
  seed <- 42
  rcvSeeds <- setSeeds(method = "repeatedcv", numbers = numbers, repeats = repeats, seed = seed)
  
  mytrainControl <- 
    trainControl(method='repeatedcv',
                 number=numbers, 
                 repeats=repeats,
                 search='grid',
                 seeds = rcvSeeds,
                 savePredictions = TRUE, 
                 classProbs = TRUE, 
                 verboseIter = TRUE)
  
  # Run 10-fold CV Model with 5 Repitions
  model <-train(group ~.,
                data=model.input, 
                method='glmnet', 
                metric='Accuracy', 
                tuneGrid=tune.grid, 
                trControl=mytrainControl,
                verboseIter = T)
  
  cat("Model Summary")
  print(model)
  cat("\n\n")
  print(plot(model))
  
  # Select model with optimal Lambda
  selectedIndices <- model$pred$lambda == model$bestTune$lambda
  df.output <- model$pred[selectedIndices, ]
  
  # MLeval
  mleval <- evalm(model)
  mleval$roc
  output.list$fitted_model <- model
  output.list$optimal.df <- df.output
  output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
  output.list$MLevaldata <- mleval$stdres
  
  return(output.list)
  
}


#----------------------------------------------------------------------------------------
###   RIDGE, LASSO, ENET LOGISTIC REGRESSION  --  For Curated Datasets & PD COMPARISON        
#----------------------------------------------------------------------------------------


ridge.lasso.enet.regression.model.DSxPD <- function(disease.model.input, model.type, obj = dat){
  
  # # TROUBLESHOOTING
  # disease.model.input = VincentC_2016.model.input;
  # obj = dat;
  # model.type = "enet"
  # 
  #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
  #' Input: Phyloseq Obj, comparison of interest, and model type
  #' Returns: a list including 
  #' 1) The fitted Model 
  #' 2) A dataframe of values with the optimal parameters (for ROC plot)
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  # Initalize variables
  model <- NULL
  output.list <- vector(mode="list", length=4)
  names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
  
  if (model.type == "ridge"){
    tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
  } else if (model.type == "lasso"){
    tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
  } else if (model.type == "enet"){
    tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
  }
  
  # Trim Control Samples from Disease Dataset 
  disease.model.input <- filter(disease.model.input, group != "control")
  # Load Native PD Dataset
  dat_pd = subset_samples(dat, donor_group == "PD")
  d <- dat_pd %>%
    microbiome::transform("compositional") %>% 
    microbiome::abundances() %>% 
    t() %>% 
    as.data.frame() 
  pd.model.input <- group_col_from_ids(d, id= rownames(d))
  rownames(pd.model.input) <- rownames(d)
  
  
# Merge Disease and PD input datasets
  model.input <- full_join(pd.model.input, disease.model.input)
  groupcol <- factor(model.input$group)
  
  # option 1)   Replace all NAs with 0s 
  # model.input[is.na(model.input)] = 0
  # option 2)  Trim all features that aren't shared (Detected in both groups)
  model.input <- model.input[ ,colSums(is.na(model.input)) == 0]
  
  # AST Transformation on merged data
  temp <- dplyr::select(model.input, -group)
  model.input <- asin(sqrt(temp))
  model.input$group <- groupcol
  

  
  # Model Parameters
  numbers <- 10
  repeats <- 5  
  set.seed(42)
  seed <- 42
  rcvSeeds <- setSeeds(method = "repeatedcv", 
                       numbers = numbers, repeats = repeats, seed = seed)
  
  mytrainControl <- 
    trainControl(method='repeatedcv',
                 number=numbers, 
                 repeats=repeats,
                 search='grid',
                 seeds = rcvSeeds,
                 savePredictions = TRUE, 
                 classProbs = TRUE, 
                 verboseIter = TRUE)
  
  # Run 10-fold CV Model with 5 Repitions
  model <-train(group ~.,
                data=model.input, 
                method='glmnet', 
                metric='Accuracy', 
                tuneGrid=tune.grid, 
                trControl=mytrainControl,
                verboseIter = T)
  
  cat("Model Summary")
  print(model)
  cat("\n\n")
  print(plot(model))
  
  # Select model with optimal Lambda
  selectedIndices <- model$pred$lambda == model$bestTune$lambda
  df.output <- model$pred[selectedIndices, ]
  
  # MLeval
  mleval <- evalm(model)
  mleval$roc
  output.list$fitted_model <- model
  output.list$optimal.df <- df.output
  output.list$AUCROC <- mleval$stdres$`Group 1`["AUC-ROC", "Score"]
  output.list$MLevaldata <- mleval$stdres
  
  return(output.list)
  
}



#--------------------------------------------------------------------------------
#                           ML Functions
#--------------------------------------------------------------------------------

prep_ml_input <- function(obj){
  all_abund <- abundances(obj) %>% 
    t() %>% as.data.frame() %>% 
    rownames_to_column("donor_id")
  all_meta <- process_meta(obj, cohort = "Merged") %>% 
    select(donor_id, PD, cohort, paired)
  ml_input <- inner_join(all_meta, all_abund) %>% 
    as_tibble()
  return(ml_input)
}
#---------------------------------------

prep_mRMR_input <- function(obj){
  
  all_abund <- abundances(obj) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column("donor_id")
  all_meta <- process_meta(obj, cohort = "Merged") %>%
    select(donor_id, PD)
  mRMR_input <- left_join(all_meta, all_abund) %>%
    as.data.frame() %>%
    column_to_rownames(var = "donor_id") %>% 
    mutate(PD = factor(PD, levels = c("No", "Yes"))) %>%
    mutate(PD = as.numeric(PD)) 
  
  return(mRMR_input)
}

#---------------------------------------
#             mRMRe wrapper 
#---------------------------------------

feature_selection <- function(df, n_features, n_algorithm = 1){
  mdat <- mRMR.data(data = data.frame(df))
  mRMR <- mRMR.ensemble(data = mdat, target_indices = 1, 
                        feature_count = n_features, solution_count = n_algorithm)
  selections <- mRMR_input[, solutions(mRMR)$`1`] %>% colnames()
  return(selections)
}





#--------------------------------------------------------------------------------
#                   Cross-Validation & Prediction LASSO 
#--------------------------------------------------------------------------------

lasso_model <- function(train, test){
  
  # TROUBLE 
  train = ml_input_tbc
  test = ml_input_rush
  
  combined <- bind_rows(train, test)
  ind <- list(
    analysis = seq(nrow(train)),
    assessment = nrow(train) + seq(nrow(test)))
  splits <- make_splits(ind, combined)
  ml_train <- training(splits)
  ml_test <- testing(splits)
  
  set.seed(42)
  # Generate 10-fold CV model with 5 repetitions, stratified by PD status 
  cv_splits <- rsample::vfold_cv(ml_train, v = 10, repeats = 5, strata = PD)
  # Define tunable Lasso model
  tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
    set_engine("glmnet")
  
  ml_recipe <- recipe(PD ~ ., data = ml_train) %>% 
    update_role(donor_id, new_role = "donor_id") %>% 
    update_role(cohort, new_role = "cohort") %>% 
    update_role(paired, new_role = "paired") %>%
    step_zv(all_numeric(), -all_outcomes()) %>% 
    step_normalize(all_numeric(), -all_outcomes())
  
  wf <- workflow() %>% 
    add_recipe(ml_recipe) %>% 
    add_model(tune_spec)
  
  #-----------------------------------------------------
  #                TUNE LASSO MODEL  
  #-----------------------------------------------------
  # Create a grid of penalty values to test
  lambda_grid <- grid_regular(penalty(), levels = 50)
  ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)
  doParallel::registerDoParallel()
  
  set.seed(42)
  lasso_grid <-
    tune_grid(wf,
              resamples = cv_splits,
              grid = lambda_grid,
              control = ctrl)
  
  # select optimal penalty by filtering largest rocauc
  best_aucroc <- select_best(lasso_grid, "roc_auc")
  # visualize model metrics of grid 
  model_performance <- 
    lasso_grid %>% 
    collect_metrics() %>% 
    ggplot(aes(penalty, mean, color = .metric)) +
    geom_errorbar(aes(ymin = mean - std_err,
                      ymax = mean + std_err),
                  alpha = 0.5) +
    geom_line(size = 1.25, show.legend = F) +
    facet_wrap(~.metric, scales = "free", nrow = 2) +
    theme_bw() + 
    scale_x_log10() +
    scale_color_viridis_d(option = "cividis", begin = .9, end = 0) +
    theme(legend.position = "none")
  print(model_performance)
  
  # Filtering ROCAUC valued from optimal 5x10fold CV model
  train_metrics <- 
    lasso_grid %>% 
    collect_metrics() %>% 
    filter(penalty == best_aucroc$penalty)
  
  # Finalize trained workflow - use on hold out test sets
  train_lasso <-
    wf %>%
    finalize_workflow(best_aucroc) %>%
    fit(ml_train)
  
  # Predictions on test data
  test_lasso <- 
    last_fit(train_lasso, split = splits) 
  test_metrics <- test_lasso %>%
    collect_metrics()
  
  # Test predictions using CV tune model
  lasso_grid %>% 
    collect_predictions()
  lasso_grid %>%
    collect_predictions() %>%
    accuracy(truth = PD, .pred_class)
  
    roc_auc(truth = PD, .pred_Yes)
  
  output <- list(
    "cv_tune_model" = lasso_grid,
    "train_lasso" = train_lasso,
    "test_lasso" = test_lasso,
    "train_metrics" = train_metrics,
    "test_metrics" = test_metrics
  )
  return(output)
}

#--------------------------------------------------------------------------------
#                      LASSO Cohort Cross Testing
#--------------------------------------------------------------------------------

lasso_cohort_summary <- function(obj_tbc_all, obj_rush_all, featSelection = NULL){
  
  # TROUBLE
  obj_tbc_all = obj_tbc_all
  obj_rush_all = obj_rush_all
  featSelection = hits_KOs
  
  ml_input_tbc <- prep_ml_input(obj_tbc_all)
  ml_input_rush <- prep_ml_input(obj_rush_all)
  
  if(!is.null(featSelection)){
    cat("Running mRMR feature selection .. \n")
    ml_input_tbc <- ml_input_tbc %>% 
      dplyr::select(donor_id, PD, cohort, paired, contains(featSelection))
    cat("Features selected from TBC: ", ncol(ml_input_tbc)-4, "\n")
    ml_input_rush <- ml_input_rush %>% 
      dplyr::select(donor_id, PD, cohort, paired, contains(featSelection))
    cat("Features selected from Rush: ", ncol(ml_input_rush)-4, "\n")
  }
  
  tbc_vs_rush <- lasso_model(train = ml_input_tbc,
                             test = ml_input_rush)
  rush_vs_tbc <- lasso_model(train = ml_input_rush,
                             test = ml_input_tbc)

  # PD vs Controls (Rush vs TBC)
  summary_A <-
    bind_rows(
      tbc_vs_rush$train_metrics %>%
        mutate(train = "TBC", test = "TBC"),
      tbc_vs_rush$test_metrics %>%
        mutate(train = "TBC", test = "Rush"),
      rush_vs_tbc$train_metrics %>%
        mutate(train = "Rush", test = "Rush"),
      rush_vs_tbc$test_metrics %>%
        mutate(train = "Rush", test = "TBC")
    )
  
  
  output <- c("tbc_vs_rush" = tbc_vs_rush, "rush_vs_tbc" = rush_vs_tbc, "summary_A" = summary_A)
  return(output)
}

#-----------------------------------------------------

#--------------------------------------------------------------------------------
#                 LASSO Cohort x Donor Group Cross Testing
#--------------------------------------------------------------------------------

lasso_cohort_x_group_summary <- function(obj_tbc_all, obj_rush_all, featSelection = NULL){
  
  ml_input_tbc <- prep_ml_input(obj_tbc_all)
  ml_input_rush <- prep_ml_input(obj_rush_all)
  
  if(!is.null(featSelection)){
    cat("Running mRMR feature selection .. \n")
    ml_input_tbc <- ml_input_tbc %>% 
      dplyr::select(donor_id, PD, cohort, paired, matches(featSelection))
    ml_input_rush <- ml_input_rush %>% 
      dplyr::select(donor_id, PD, cohort, paired, matches(featSelection))
  }
  
  ml_input_tbc_house <- obj_tbc_all %>% 
    subset_samples(paired != "No") %>% prep_ml_input()
  ml_input_tbc_healthy <- obj_tbc_all %>% 
    subset_samples(donor_group != "HC") %>% prep_ml_input() 
  ml_input_rush_house <- obj_rush_all %>% 
    subset_samples(paired != "No") %>% prep_ml_input()
  ml_input_rush_healthy <- obj_rush_all %>% 
    subset_samples(donor_group != "HC") %>% prep_ml_input()
  
  
  rush_PC_vs_rush_HC <- lasso_model(train = ml_input_rush_healthy,
                                    test = ml_input_rush_house)
  rush_PC_vs_tbc_PC <- lasso_model(train = ml_input_rush_healthy,
                                   test = ml_input_tbc_healthy)
  rush_PC_vs_tbc_HC <- lasso_model(train = ml_input_rush_healthy,
                                   test = ml_input_tbc_house)
  
  rush_HC_vs_rush_PC <- lasso_model(train = ml_input_rush_house,
                                    test = ml_input_rush_healthy)
  rush_HC_vs_tbc_PC <- lasso_model(train = ml_input_rush_house,
                                   test = ml_input_tbc_healthy)
  rush_HC_vs_tbc_HC <- lasso_model(train = ml_input_rush_house,
                                   test = ml_input_tbc_house)
  
  tbc_PC_vs_tbc_HC <- lasso_model(train = ml_input_tbc_healthy,
                                  test = ml_input_tbc_house)
  tbc_PC_vs_rush_PC <- lasso_model(train = ml_input_tbc_healthy,
                                   test = ml_input_rush_healthy)
  tbc_PC_vs_rush_HC <- lasso_model(train = ml_input_tbc_healthy,
                                   test = ml_input_rush_house)
  
  tbc_HC_vs_tbc_PC <- lasso_model(train = ml_input_tbc_house,
                                  test = ml_input_tbc_healthy)
  tbc_HC_vs_rush_PC <- lasso_model(train = ml_input_tbc_house,
                                   test = ml_input_rush_healthy)
  tbc_HC_vs_rush_HC <- lasso_model(train = ml_input_tbc_house,
                                   test = ml_input_rush_house)
  
  
  summaryB <-
    bind_rows(
      # RUSH PC vs all
      rush_PC_vs_rush_HC$train_metrics %>%
        mutate(train = "Rush_PC", test = "Rush_PC"),
      rush_PC_vs_rush_HC$test_metrics %>%
        mutate(train = "Rush_PC", test = "Rush_HC"),
      rush_PC_vs_tbc_PC$test_metrics %>%
        mutate(train = "Rush_PC", test = "TBC_PC"),
      rush_PC_vs_tbc_HC$test_metrics %>%
        mutate(train = "Rush_PC", test = "TBC_HC"),
      # RUSH HC vs all
      rush_HC_vs_rush_PC$train_metrics %>%
        mutate(train = "Rush_HC", test = "Rush_HC"),
      rush_HC_vs_rush_PC$test_metrics %>%
        mutate(train = "Rush_HC", test = "Rush_PC"),
      rush_HC_vs_tbc_PC$test_metrics %>%
        mutate(train = "Rush_HC", test = "TBC_PC"),
      rush_HC_vs_tbc_HC$test_metrics %>%
        mutate(train = "Rush_HC", test = "TBC_HC"),
      # TBC PC vs all
      tbc_PC_vs_tbc_HC$train_metrics %>%
        mutate(train = "TBC_PC", test = "TBC_PC"),
      tbc_PC_vs_tbc_HC$test_metrics %>%
        mutate(train = "TBC_PC", test = "TBC_HC"),
      tbc_PC_vs_rush_PC$test_metrics %>%
        mutate(train = "TBC_PC", test = "Rush_PC"),
      tbc_PC_vs_rush_HC$test_metrics %>%
        mutate(train = "TBC_PC", test = "Rush_HC"),
      # TBC PC vs all
      tbc_HC_vs_tbc_PC$train_metrics %>%
        mutate(train = "TBC_HC", test = "TBC_HC"),
      tbc_HC_vs_tbc_PC$test_metrics %>%
        mutate(train = "TBC_HC", test = "TBC_PC"),
      tbc_HC_vs_rush_PC$test_metrics %>%
        mutate(train = "TBC_HC", test = "Rush_PC"),
      tbc_HC_vs_rush_HC$test_metrics %>%
        mutate(train = "TBC_HC", test = "Rush_HC")
    ) %>% 
    mutate(train_cohort = substr(train, 1, nchar(train)-3)) %>% 
    mutate(test_cohort = substr(test, 1, nchar(test)-3)) %>% 
    mutate(train_group = substr(train, nchar(train)-1, nchar(train))) %>% 
    mutate(test_group = substr(test, nchar(test)-1, nchar(test)))
  
  return(summaryB)
}



































# 
# 
# lasso_model <- function(train, test){
#   
#   # TROUBLE
#   ml_input_tbc <- prep_ml_input(obj_tbc_all)
#   ml_input_rush <- prep_ml_input(obj_rush_all)
#   train = ml_input_tbc
#   test = ml_input_rush
# 
#   
# 
#   combined <- bind_rows(train, test)
#   ind <- list(
#     analysis = seq(nrow(train)),
#     assessment = nrow(train) + seq(nrow(test)))
#   splits <- make_splits(ind, combined)
#   ml_train <- training(splits)
#   ml_test <- testing(splits)
#   
#   set.seed(42)
#   # Generate 10-fold CV model with 5 repetitions, stratified by PD status 
#   cv_splits <- rsample::vfold_cv(ml_train, v = 10, repeats = 10, strata = PD)
#   # Define tunable Lasso model
#   tune_spec <- logistic_reg(penalty = tune(), mixture = 1) %>%
#     set_engine("glmnet")
#   
#   ml_recipe <- recipe(PD ~ ., data = ml_train) %>% 
#     update_role(donor_id, new_role = "donor_id") %>% 
#     update_role(cohort, new_role = "cohort") %>% 
#     update_role(paired, new_role = "paired") %>%
#     step_zv(all_numeric(), -all_outcomes()) %>% 
#     step_normalize(all_numeric(), -all_outcomes())
#   
#   wf <- workflow() %>% 
#     add_recipe(ml_recipe) %>% 
#     add_model(tune_spec)
#   
#   #-----------------------------------------------------
#   #                TUNE LASSO MODEL  
#   #-----------------------------------------------------
#   # Create a grid of penalty values to test
#   lambda_grid <- grid_regular(penalty(), levels = 50)
#   ctrl <- control_grid(save_pred = TRUE, verbose = TRUE)
#   doParallel::registerDoParallel()
#   
#   set.seed(42)
#   lasso_grid <-
#     tune_grid(wf,
#               resamples = cv_splits,
#               grid = lambda_grid,
#               control = ctrl)
#   
#   # select optimal penalty by filtering largest rocauc
#   best_aucroc <- select_best(lasso_grid, "roc_auc")
#   # visualize model metrics of grid 
#   model_performance <- 
#     lasso_grid %>% 
#     collect_metrics() %>% 
#     ggplot(aes(penalty, mean, color = .metric)) +
#     geom_errorbar(aes(ymin = mean - std_err,
#                       ymax = mean + std_err),
#                   alpha = 0.5) +
#     geom_line(size = 1.25, show.legend = F) +
#     facet_wrap(~.metric, scales = "free", nrow = 2) +
#     theme_bw() + 
#     scale_x_log10() +
#     scale_color_viridis_d(option = "cividis", begin = .9, end = 0) +
#     theme(legend.position = "none")
#   print(model_performance)
#   
#   train_lasso2 <-
#     wf %>%
#     finalize_workflow(best_aucroc) %>% 
#     fit_resamples(cv_splits, control = control_resamples(save_pred = TRUE))
# 
#   train_lasso2 %>% 
#     collect_metrics() %>% 
#     filter(penalty == best_aucroc$penalty)
#   
#   
#   # Filtering ROCAUC valued from optimal 5x10fold CV model
#   # train_metrics <- 
#   #   lasso_grid %>% 
#   #   collect_metrics() %>% 
#   #   filter(penalty == best_aucroc$penalty)
#   
#   # Finalize trained workflow - use on hold out test sets
#   train_lasso <-
#     wf %>%
#     finalize_workflow(best_aucroc) %>%
#     fit(ml_train)
#   
#   # Predictions on test data
#   test_lasso <- 
#     last_fit(train_lasso2, split = splits) 
#   test_metrics <- test_lasso %>%
#     collect_metrics()
#   
#   output <- list("train_lasso" = train_lasso, 
#                  "test_lasso" = test_lasso, 
#                  "train_metrics" = train_metrics,
#                  "test_metrics" = test_metrics)
#   return(output)
# }