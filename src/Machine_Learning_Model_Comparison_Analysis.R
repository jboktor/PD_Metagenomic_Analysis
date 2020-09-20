# Machine Learning - Model Comparison Analysis

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Machine_Learning_Models.R")

#---------------------------------------------------------
#                Data Objects for Loop  
#---------------------------------------------------------

x <- c(dat, dat.path.slim, dat.KOs.slim)
z <- c("Species", "Pathways", "KOs")

#---------------------------------------------------------
#           Beginning of ML Comparison Loop
#---------------------------------------------------------

cnt <- 1
for (i in x){
  cat("\n\n\n")
  cat("Processing input: ", z[cnt], "\n")
  cat("\n")
  print(i)
  cat("\n")

  
  # Ridge Regression 
  ridge.PDPC <- ridge.lasso.enet.regression.model(obj = dat.path, comparison = "PDvPC", model.type = "ridge")
  ridge.PDHC <- ridge.lasso.enet.regression.model(obj = i, comparison = "PDvHC", model.type = "ridge")
  # LASSO Regression 
  lasso.PDPC <- ridge.lasso.enet.regression.model(obj = i, comparison = "PDvPC", model.type = "lasso")
  lasso.PDHC <- ridge.lasso.enet.regression.model(obj = i, comparison = "PDvHC", model.type = "lasso")
  # ElasticNet Regression 
  enet.PDPC <- ridge.lasso.enet.regression.model(obj = i, comparison = "PDvPC", model.type = "enet")
  enet.PDHC <- ridge.lasso.enet.regression.model(obj = i, comparison = "PDvHC", model.type = "enet")
  # Random Forrest Models 
  rf.PDPC <- random.forest.model(obj = i, comparison = "PDvPC")
  rf.PDHC <- random.forest.model(obj = i, comparison = "PDvHC")
  # XGBoost Models
  xgboost.PDPC <- xgboost.model(obj = dat, comparison = "PDvPC")
  xgboost.PDHC <- xgboost.model(obj = i, comparison = "PDvHC")
  
  
  # Save Models
  saveRDS(ridge.PDPC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/ridge.PDPC.rds"))
  saveRDS(ridge.PDHC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/ridge.PDHC.rds"))
  saveRDS(lasso.PDPC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/lasso.PDPC.rds"))
  saveRDS(lasso.PDHC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/lasso.PDHC.rds"))
  saveRDS(enet.PDPC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/enet.PDPC.rds"))
  saveRDS(enet.PDHC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/enet.PDHC.rds"))
  saveRDS(rf.PDPC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/RandomForest.PDPC.rds"))
  saveRDS(rf.PDHC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/RandomForest.PDHC.rds"))
  saveRDS(xgboost.PDPC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/XGBoost.PDPC.rds"))
  saveRDS(xgboost.PDHC, file = paste0("files/Machine_Learning_Models/", z[cnt],"/XGBoost.PDHC.rds"))

  
  cols.ml <- brewer.pal(10, "Paired")
  
  rocplot <-  ggplot() + 
    style_roc() +
    theme_bw() +
    geom_abline(intercept = 0, slope = 1, color = "grey") +
    geom_roc(data = ridge.PDHC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
             n.cuts=0, color = cols.ml[1], linealpha = 0.8) +  
    geom_roc(data = ridge.PDPC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "PC"))),
             n.cuts=0, color = cols.ml[2], linealpha = 0.8, labels = TRUE) +  
    geom_roc(data = lasso.PDHC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
             n.cuts=0, color = cols.ml[3], linealpha = 0.8) + 
    geom_roc(data = lasso.PDPC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "PC"))),
             n.cuts=0, color = cols.ml[4], linealpha = 0.8) + 
    geom_roc(data = enet.PDHC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
             n.cuts=0, color = cols.ml[5], linealpha = 0.8) + 
    geom_roc(data = enet.PDPC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "PC"))),
             n.cuts=0, color = cols.ml[6], linealpha = 0.8) + 
    geom_roc(data = rf.PDHC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
             n.cuts=0, color = cols.ml[7], linealpha = 0.8) + 
    geom_roc(data = rf.PDPC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "PC"))),
             n.cuts=0, color = cols.ml[8], linealpha = 0.8) + 
    geom_roc(data = xgboost.PDHC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "HC"))),
             n.cuts=0, color = cols.ml[9], linealpha = 0.8) + 
    geom_roc(data = xgboost.PDPC$optimal.df, aes(m = PD, d = factor(obs, levels = c("PD", "PC"))),
             n.cuts=0, color = cols.ml[10], linealpha = 0.8) + 
    # annotate("text", x = 0.45, y = 0.3, label = paste0("PC: AUC-ROC = ", PC.ROC), color = PC.col, hjust=0) +
    # annotate("text", x = 0.45, y = 0.22, label = paste0("HC: AUC-ROC = ", HC.ROC), color = HC.col, hjust=0) +
    coord_equal() +
    theme(panel.grid = element_blank())
  rocplot
  
  
  ml.models <- list(ridge.PDHC, ridge.PDPC, lasso.PDHC, lasso.PDPC, enet.PDHC, enet.PDPC, 
                    rf.PDHC, rf.PDPC, xgboost.PDHC, xgboost.PDPC)
  names(ml.models) <- c("ridge.PDHC", "ridge.PDPC", "lasso.PDHC", "lasso.PDPC", "enet.PDHC", "enet.PDPC", 
                        "rf.PDHC", "rf.PDPC", "xgboost.PDHC", "xgboost.PDPC")
  names(cols.ml) <- c("ridge.PDHC", "ridge.PDPC", "lasso.PDHC", "lasso.PDPC", "enet.PDHC", "enet.PDPC", 
                        "rf.PDHC", "rf.PDPC", "xgboost.PDHC", "xgboost.PDPC")
  cols.ml.labels <- c("ridge.PDHC" = "Ridge HC", "ridge.PDPC" = "Ridge PC", 
                      "lasso.PDHC" = "Lasso HC", "lasso.PDPC" = "Lasso PC", 
                      "enet.PDHC" = "ElasticNet HC", "enet.PDPC" = "ElasticNet PC", 
                      "rf.PDHC" = "RandomForest HC", "rf.PDPC" = "RandomForest PC", 
                      "xgboost.PDHC" = "XGBoost HC", "xgboost.PDPC" = "XGBoost PC")
  
  
  ml.models.df <- data.frame()
  for(i in 1:length(ml.models)){
    
    model.ID <- names(ml.models)[i]
    AUCROC <- ml.models[[i]]$MLevaldata$`Group 1`["AUC-ROC", "Score"]
    AUCPR <- ml.models[[i]]$MLevaldata$`Group 1`["AUC-PR", "Score"]
    
    temp <- cbind(model.ID, AUCROC, AUCPR)
    ml.models.df <- rbind(ml.models.df, temp)
  }
  
  ml.models.df <- mutate(ml.models.df, comparison = if_else(grepl("HC", model.ID), "PDvsHC", "PDvsPC"))
  ml.models.df <- 
    mutate(ml.models.df, model.type = 
             if_else(grepl("ridge", model.ID), "Ridge", 
                     if_else(grepl("lasso", model.ID), "Lasso",
                             if_else(grepl("ridge", model.ID), "ENet",
                                     if_else(grepl("rf", model.ID), "Random Forest",
                                             if_else(grepl("xgboost", model.ID), "XGBoost",
                                                     "PDvsPC"))))))
  
  ml.models.df$AUCROC <- as.numeric(ml.models.df$AUCROC)
  ml.models.df$AUCPR <- as.numeric(ml.models.df$AUCPR)
  
  barplot.aucroc <- 
    ml.models.df %>% 
    dplyr::group_by(model.type) %>% 
    ggplot(aes(x = model.ID, y=AUCROC, fill = model.ID)) +
    geom_col() +
    theme_bw() +
    coord_cartesian(ylim = c(0.4, max(AUCROC) + 0.1)) +
    labs(y="AUC-ROC", fill = "Models") +
    scale_fill_manual(values = cols.ml, labels = cols.ml.labels) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  barplot.aucpr <- 
    ml.models.df %>% 
    dplyr::group_by(model.type) %>% 
    ggplot(aes(x=model.ID, y=AUCPR, group = model.type, fill = model.ID)) +
    geom_col() +
    theme_bw() +
    coord_cartesian(ylim = c(0.4, max(AUCPR) + 0.1)) +
    labs(y="AUC-PR") +
    scale_fill_manual(values = cols.ml) +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none")
  
  
  legend <- cowplot::plot_grid(get_legend(barplot.aucroc))
  barplot.aucroc2 <- barplot.aucroc + theme(legend.position = "none")
  barplots <- cowplot::plot_grid(barplot.aucroc2, barplot.aucpr, ncol=1, axis = "lrtb", align = "v")
  barplots2 <- cowplot::plot_grid(barplots, legend, ncol=2, rel_widths = c(3,1)) #, axis = "lrtb", align = "v")
  
  
  finalplot <- cowplot::plot_grid(rocplot, barplots2, ncol=2, align = "h")
  finalplot
  
  ggsave(finalplot, filename = paste0("data/Machine_Learning_Analysis/Model_Comparisons_", z[cnt], ".svg"), 
         width = 11.5, height = 5)
  
  cnt <- cnt + 1
}


