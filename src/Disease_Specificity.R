##  Disease Specificity 

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Machine_Learning_Models.R")


# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "devel")
# BiocManager::valid()
# BiocManager::install("curatedMetagenomicData")


library(curatedMetagenomicData)

# ?combined_metadata
# View(combined_metadata)
unique(combined_metadata$body_site)
stool_samples <- combined_metadata %>%
  filter(body_site == "stool")

## Explore studies
# unique(stool_samples$disease)
# testg <- filter(combined_metadata, grepl("metabolic_syndrome", disease))
# testg$dataset_name


#-------------------------------------------------------------------------------------
#                        Explore Dataset Features
#-------------------------------------------------------------------------------------

# Study Summary - CDI
VincentC_2016.metaphlan_bugs_list.stool() %>% 
  experimentData()
# Examine sample conditions
study_vars <- stool_samples %>% 
  filter(dataset_name == "VincentC_2016") %>% 
  select(study_condition)
unique(study_vars$study_condition)
# Prep input data and fit model
VincentC_2016.model.input <- prep.CMD.Species.ML(study = "VincentC_2016")

VincentC_2016.model <- ridge.lasso.enet.regression.model.DS(
  model.input = VincentC_2016.model.input, model.type = "enet")

VincentC_2016.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = VincentC_2016.model.input, obj = dat, model.type = "enet")


# ACVD
JieZ_2017.metaphlan_bugs_list.stool() %>%
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "JieZ_2017") %>% 
  select(study_condition)
unique(study_vars$study_condition)
JieZ_2017.model.input <- prep.CMD.Species.ML(study = "JieZ_2017")
JieZ_2017.model <- ridge.lasso.enet.regression.model.DS(
  model.input = JieZ_2017.model.input, model.type = "enet")
JieZ_2017.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = JieZ_2017.model.input, obj = dat, model.type = "enet")

#  T1D
KosticAD_2015.metaphlan_bugs_list.stool()%>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "KosticAD_2015") %>% 
  select(study_condition)
unique(study_vars$study_condition) # Trim NAs
KosticAD.model.input <- prep.CMD.Species.ML(study = "KosticAD_2015" , metafilter = "NA")
KosticAD.model <- ridge.lasso.enet.regression.model.DS(
  model.input = KosticAD.model.input, model.type = "enet")
KosticAD.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = KosticAD.model.input, obj = dat, model.type = "enet")

# T2D
QinJ_2012.metaphlan_bugs_list.stool() %>%
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "QinJ_2012") %>% 
  select(study_condition)
unique(study_vars$study_condition) # Trim NAs
QinJ_2012.model.input <- prep.CMD.Species.ML(study = "QinJ_2012", metafilter = "NA")
QinJ_2012.model <- ridge.lasso.enet.regression.model.DS(
  model.input = QinJ_2012.model.input, model.type = "enet")
QinJ_2012.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = QinJ_2012.model.input, obj = dat, model.type = "enet")

# cirrhosis
QinN_2014.metaphlan_bugs_list.stool() %>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "QinN_2014") %>% 
  select(study_condition)
unique(study_vars$study_condition)
QinN_2014.model.input <- prep.CMD.Species.ML(study = "QinN_2014")
QinN_2014.model <- ridge.lasso.enet.regression.model.DS(
  model.input = QinN_2014.model.input, model.type = "enet")
QinN_2014.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = QinN_2014.model.input, obj = dat, model.type = "enet")

# IBD
NielsenHB_2014.metaphlan_bugs_list.stool() %>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "NielsenHB_2014") %>% 
  select(study_condition)
unique(study_vars$study_condition)
NielsenHB_2014.model.input <- prep.CMD.Species.ML(study = "NielsenHB_2014")
NielsenHB_2014.model <- ridge.lasso.enet.regression.model.DS(
  model.input = NielsenHB_2014.model.input, model.type = "enet")
NielsenHB_2014.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = NielsenHB_2014.model.input, obj = dat, model.type = "enet")

# Hypertension
LiJ_2017.metaphlan_bugs_list.stool() %>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "LiJ_2017") %>% 
  select(study_condition)
unique(study_vars$study_condition) # Trim pre-hypertension vars
LiJ_2017.model.input <- prep.CMD.Species.ML(study = "LiJ_2017", metafilter = "pre-hypertension")
LiJ_2017.model <- ridge.lasso.enet.regression.model.DS(
  model.input = LiJ_2017.model.input, model.type = "enet")
LiJ_2017.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = LiJ_2017.model.input, obj = dat, model.type = "enet")

# metabolic_syndrome
LiSS_2016.metaphlan_bugs_list.stool() %>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "LiSS_2016") %>% 
  select(study_condition)
unique(study_vars$study_condition) # Trim FMT vars
LiSS_2016.model.input <- prep.CMD.Species.ML(study = "LiSS_2016", metafilter = "FMT")
LiSS_2016.model <- ridge.lasso.enet.regression.model.DS(
  model.input = LiSS_2016.model.input, model.type = "enet")
LiSS_2016.model.PD <- ridge.lasso.enet.regression.model.DSxPD(
  disease.model.input = LiSS_2016.model.input, obj = dat, model.type = "enet")


#-------------------------------------------------------------------------------------


ml.models.0 <- list(VincentC_2016.model, JieZ_2017.model, KosticAD.model, QinJ_2012.model,
                  QinN_2014.model, NielsenHB_2014.model, LiJ_2017.model, LiSS_2016.model)
names(ml.models.0) <- c("CDI", "ACVD", "T1D", "T2D", 
                      "cirrhosis", "IBD", "hypertension", "Metabolic Syndrome")

ml.models.PD <- list(VincentC_2016.model.PD, JieZ_2017.model.PD, KosticAD.model.PD, QinJ_2012.model.PD,
                  QinN_2014.model.PD, NielsenHB_2014.model.PD, LiJ_2017.model.PD, LiSS_2016.model.PD)
names(ml.models.PD) <- c("CDI", "ACVD", "T1D", "T2D", 
                      "cirrhosis", "IBD", "hypertension", "Metabolic Syndrome")


#-------------------------------------------------------------------------------------

ml.models.0.df <- data.frame()
for(i in 1:length(ml.models.0)){
  
  model.ID <- names(ml.models.0)[i]
  AUCROC <- ml.models.0[[i]]$MLevaldata$`Group 1`["AUC-ROC", "Score"]
  AUCPR <- ml.models.0[[i]]$MLevaldata$`Group 1`["AUC-PR", "Score"]
  
  temp <- cbind(model.ID, AUCROC, AUCPR)
  ml.models.0.df <- rbind(ml.models.0.df, temp)
}

ml.models.PD.df <- data.frame()
for(i in 1:length(ml.models.PD)){
  
  model.ID <- names(ml.models.PD)[i]
  AUCROC <- ml.models.PD[[i]]$MLevaldata$`Group 1`["AUC-ROC", "Score"]
  AUCPR <- ml.models.PD[[i]]$MLevaldata$`Group 1`["AUC-PR", "Score"]
  
  temp <- cbind(model.ID, AUCROC, AUCPR)
  ml.models.PD.df <- rbind(ml.models.PD.df, temp)
}



#-------------------------------------------------------------------------------------
#                                Plotting Data
#-------------------------------------------------------------------------------------

library(dendextend)
library(ggdendro)


#-------------------------------------------------------------------------------------
#                              Plotting Function
#-------------------------------------------------------------------------------------

ml.heatmap <- function(ml.models.df){
  
  # Clustering
  ml.matrix <- as.matrix(ml.models.df[, -1])
  rownames(ml.matrix) <- ml.models.df$model.ID 
  ml.dendro <- as.dendrogram(hclust(d = dist(x = ml.matrix), method = "centroid"))
  # Create dendro
  dendro.plot <- ggdendrogram(data = ml.dendro, rotate = TRUE)
  print(dendro.plot)
  # Order the levels according to their position in the cluster
  ml.order <- order.dendrogram(ml.dendro)
  ml.models.df$model.ID <- factor(ml.models.df$model.ID,
                                  levels = ml.models.df$model.ID[ml.order], 
                                  ordered = TRUE)
  ml.models.df$ID <- paste0("ML", 1:nrow(ml.models.df))
  
  # pivot into Long format
  ml.models.df2 <- pivot_longer(ml.models.df, cols = c(AUCROC, AUCPR), names_to = "eval_metric")
  ml.models.df2$value <- as.numeric(ml.models.df2$value)
  ml.models.df2$ID <- paste0("ML", 1:nrow(ml.models.df2))
  
  ## Heatmap  
  hm.plot <- 
    ggplot(data = ml.models.df2, aes(x = model.ID, y =  eval_metric, fill= value)) + 
    geom_tile() +
    scale_fill_distiller(palette = "BuPu", direction = 1,
                         limits = c(min(0.4), max(1)),
                         breaks = c(0.4, 0.6, 0.8, 1)) +
    coord_flip() +
    geom_text(aes(label=value), color = "white") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "top",
          legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
  print(hm.plot)
  
  dendro.plot.merge <- dendro.plot + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  final.plot <- cowplot::plot_grid(hm.plot, dendro.plot.merge, align = "h", axis="lrtb", rel_widths = c(3,1))
  
  return(final.plot)
  
}

HM1 <- ml.heatmap(ml.models.0.df)
HM2 <- ml.heatmap(ml.models.PD.df)


feature.specificity <- 
  cowplot::plot_grid(HM1, NULL, HM2, NULL, ncol = 4, align = "hv", 
                     rel_widths = c(5, 0.5, 5, 0.5), labels = "AUTO")

ggsave(feature.specificity,
       filename = "data/Machine_Learning_Analysis/Feature_Specificity_Heatmap_Species_sharedfeatures.svg",
       width = 9, height = 4.5)




# Neurodegenerative Disorders

# Disorders with Gastrointestinal distress

# Auto Immune Disorders

# Metabolic Disorders

# Infectious Disease / Pathogenic Agents

# Otherwise non-related 






