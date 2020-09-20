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
  #' Input: Phlyoseq Obj, comparison of interest, and model type
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
  #' Input: Phlyoseq Obj and comparison of interest 
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
  #' Input: Phlyoseq Obj and comparison of interest 
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
  #' Input: Phlyoseq Obj, comparison of interest, and model type
  #' Returns: a list including 
  #' 1) The fitted Model 
  #' 2) A dataframe of values with the optimal parameters (for ROC plot)
  #' 3) MLeval AUC-ROC values 
  #' 4) Other MLeval analysis parameters
  
  # Initalize variables
  model <- NULL
  # model.input <- NULL
  output.list <- vector(mode="list", length=4)
  names(output.list) <- c("fitted_model", "optimal.df", "AUCROC", "MLevaldata")
  
  if (model.type == "ridge"){
    tune.grid = expand.grid(alpha = 0, lambda=seq(0, 1, 0.1))
  } else if (model.type == "lasso"){
    tune.grid = expand.grid(alpha = 1, lambda=seq(0, 1, 0.1))
  } else if (model.type == "enet"){
    tune.grid = expand.grid(alpha = 0.5, lambda=seq(0, 1, 0.1))
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
###   RIDGE, LASSO, ENET LOGISTIC REGRESSION  --  For Curated Datasets & PD COMPARISON        
#----------------------------------------------------------------------------------------


ridge.lasso.enet.regression.model.DSxPD <- function(disease.model.input, model.type, obj = dat){
  
  #' Function to train a Rigde, lASSO, or ElasticNet Regression Model
  #' Input: Phlyoseq Obj, comparison of interest, and model type
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
  d <- asin(sqrt(d))
  pd.model.input <- group_col_from_ids(d, id= rownames(d))
  rownames(pd.model.input) <- rownames(d)
  
  
# Merge Disease and PD input datasets
  model.input <- full_join(pd.model.input, disease.model.input)
  model.input$group <- factor(model.input$group)
  
  # Replace all NAs with 0s 
  model.input[is.na(model.input)] = 0
  
  # Trim all features that aren't shared (Detected in both groups)
  # model.input <- model.input[ ,colSums(is.na(model.input)) == 0]
  
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




