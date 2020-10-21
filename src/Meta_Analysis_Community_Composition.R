# Meta-Analysis  


######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")
source("src/Community_Composition_Funcs.R")
# source("src/create_phlyoseq_obj_CMG.R")

datasets <- curatedMetagenomicData(
  c("YeZ_2018.metaphlan_bugs_list.stool",
    "ChengpingW_2017.metaphlan_bugs_list.stool",
    "VincentC_2016.metaphlan_bugs_list.stool",
    "JieZ_2017.metaphlan_bugs_list.stool",
    "KosticAD_2015.metaphlan_bugs_list.stool",
    "QinJ_2012.metaphlan_bugs_list.stool",
    "QinN_2014.metaphlan_bugs_list.stool",
    "NielsenHB_2014.metaphlan_bugs_list.stool",
    "LiJ_2017.metaphlan_bugs_list.stool",
    "LiSS_2016.metaphlan_bugs_list.stool"),
  dryrun = FALSE)

# Construct phyloseq object from the five datasets
physeq <-
  # Aggregate the five studies into ExpressionSet
  mergeData(datasets) %>%
  # Convert to phyloseq object
  ExpressionSet2phyloseq()

# Select only disease and control
physeq %>% 
  meta() %>% 
  select(study_condition) %>% 
  unique()
study_groups <- c("control", "ACVD", "BD", "AS", "CDI", "T1D", "hypertension", "metabolic_syndrome", 
  "IBD", "T2D", "cirrhosis")

physeq <- 
  physeq %>%
  # Subset samples to only CRC and controls
  subset_samples(study_condition %in% study_groups) %>% 
  # Subset features to species
  subset_taxa(!is.na(Species) & is.na(Strain)) %>%
  # Normalize abundances to relative abundance scale
  microbiome::transform("compositional") 
  # Filter features to be of at least 1e-5 relative abundance in five samples
  # filter_taxa(kOverA(5, 1e-5), prune = TRUE)


#------------------------------------------
# PREP PD Metagenomics metadata for merge 
#------------------------------------------
# dat.CMG <- load_CMG_formatted_phylo()
physeq <- merge_phyloseq(physeq, dat.CMG)
disease_abd <- microbiome::abundances(physeq) #otu_table(physeq)@.Data
disease_meta <- microbiome::meta(physeq) # data.frame(sample_data(physeq))
disease_meta$studyID <- factor(disease_meta$studyID)

disease_meta$dataset_name <- sub(".metaphlan_bugs_list.stool", "", disease_meta$studyID)
disease_meta$study_condition <- factor(disease_meta$study_condition, 
                                       levels = c("control", "parkinsons", "ACVD", "BD", "AS", "CDI", "T1D", "hypertension", 
                                                "metabolic_syndrome","IBD", "T2D", "cirrhosis"))
    
# unique(disease_meta$study_condition)
# unique(disease_meta$dataset_name)
#------------------------------------------------------------------------------------------
# SEQUENCING READS
#------------------------------------------------------------------------------------------

dsranking <- disease_meta %>%
  dplyr::group_by(dataset_name) %>%
  dplyr::summarize(mediandepth = median(number_reads) / 1e6) %>%
  dplyr::mutate(dsorder = rank(mediandepth)) %>%
  dplyr::arrange(dsorder)

meta_reads <- disease_meta %>%
  mutate(ds = factor(disease_meta$dataset_name, levels=dsranking$dataset_name)) %>%
  ggplot(aes(ds, x = dataset_name, y = log10(number_reads/1e6), fill=study_condition)) + 
  geom_boxplot(alpha = 0.4, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition), 
             position = position_jitterdodge(jitter.width=0.4), shape=21, size=1, alpha = 0.4) +
  theme_bw() + 
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  labs(y="log10 Read Depth (millions)", fill="Study Condition") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())
# ggsave(meta_reads, filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_ReadDepth.svg",
#        height = 5, width =7)

#------------------------------------------------------------------------------------------
# ALPHA-DIVERSITY
#------------------------------------------------------------------------------------------

## Calculate Alpha Diversity Metrics and add cols to df
disease_meta$Observed <- microbiome::alpha(abundances(physeq), 'observed')$observed
disease_meta$Shannon <- microbiome::alpha(abundances(physeq), 'shannon')$diversity_shannon
disease_meta$Evenness <- evenness(abundances(physeq), 'simpson')$simpson

# Plot histograms to get a sense of data distribution
par(mfrow = c(1, 3))
hist(disease_meta$Observed, main="observed OTUs", xlab="", breaks=10)
hist(disease_meta$Shannon, main="Shannon diversity", xlab="", breaks=10)
hist(disease_meta$Evenness, main="Simpson's evenness", xlab="", breaks=10)


meta_alpha_observed <- ggplot(data = disease_meta, aes(x = dataset_name, y = Observed, fill = study_condition)) + 
  theme_bw() + 
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition), 
             position = position_jitterdodge(jitter.width=0.4), shape=21, size=1, alpha = 0.4) +
  labs(y = "Observed Species",  fill="Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

meta_alpha_shannon <- ggplot(data = disease_meta, aes(x = dataset_name, y = Shannon, fill = study_condition)) + 
  theme_bw() + 
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition), 
             position = position_jitterdodge(jitter.width=0.4), shape=21, size=1, alpha = 0.4) +
  labs(y = "Shannon Diversity",  fill="Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

meta_alpha_simpsons <- ggplot(data = disease_meta, aes(x = dataset_name, y = Evenness, fill = study_condition)) + 
  theme_bw() + 
  geom_boxplot(alpha = 0.3, outlier.alpha = 0) +
  geom_point(aes(fill = study_condition, color = study_condition), 
             position = position_jitterdodge(jitter.width=0.4), shape=21, size=1, alpha = 0.4) +
  labs(y = "Simpson's Evenness",  fill="Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

alp <- cowplot::plot_grid(meta_reads, meta_alpha_observed, meta_alpha_shannon, meta_alpha_simpsons, 
                          nrow=2, align = "v", axis = "tblr") 
ggsave(alp, filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_Alpha.svg",
       height = 10, width =14)



#------------------------------------------------------------------------------------------
# BETA-DIVERSITY
#------------------------------------------------------------------------------------------
# AITCHISONS DISTANCE (EUCLIDIAN DIST ON CLR TRANSFORMED DATA)  
#---------------------------------------------------------
obj_clr <- microbiome::transform(physeq, "clr")
iDist <- phyloseq::distance(obj_clr, method="euclidean")
iMDS  <- phyloseq::ordinate(obj_clr, "MDS", distance=iDist)

p <- plot_ordination(obj_clr, iMDS, color="studyID", axes = c(1, 2))
df12 = p$data
df12$study_condition <- factor(df12$study_condition, 
                               levels = c("control", "parkinsons", "ACVD", "BD", "AS", "CDI", "T1D", "hypertension", 
                                          "metabolic_syndrome","IBD", "T2D", "cirrhosis"))
df12$dataset_name <- sub(".metaphlan_bugs_list.stool", "", df12$studyID) %>% factor()


df12 <- df12 %>% mutate(disease_category = if_else(study_condition == "control", "Healthy",
                                                   if_else(study_condition == "parkinsons", "Neurodegenerative",
                                                           if_else(study_condition == "AS" | study_condition == "T1D" | 
                                                                     study_condition == "BD" | study_condition == "IBD", "Inflammatory/Autoimmune",
                                                                   if_else(study_condition == "T2D" | study_condition == "metabolic_syndrome" |
                                                                             study_condition == "hypertension" | study_condition == "ACVD" , "Nutritional",
                                                                           if_else(study_condition == "CDI", "Infectious",
                                                                                   if_else(study_condition == "cirrhosis", "Alcohol related", "ERROR"
                                                                                           )))))))
# unique(df12$disease_category)
# unique(df12$study_condition)

p1 <- 
  ggplot(df12, aes(Axis.1, Axis.2, fill = study_condition, color=study_condition)) + 
  geom_point(shape=21, size=1.5, alpha=0.8) + 
  geom_point(data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
             aes(Axis.1, Axis.2), shape=21, size=1.5, color= "black") +
  theme_bw() + 
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Study Condition") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())


p2 <- 
  ggplot(df12, aes(Axis.1, Axis.2)) +
  geom_point(shape=21, size=1.5, alpha=0.8, aes(fill = dataset_name, color=dataset_name)) + 
  geom_point(data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
             aes(Axis.1, Axis.2), shape=21, size=1.5, color= "black") +
  theme_bw() + 
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Dataset") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())


p3 <- 
  ggplot(df12, aes(Axis.1, Axis.2)) + 
  geom_point(shape=21, size=1.5, alpha=0.8, aes(fill = country, color=country)) + 
  geom_point(data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
             aes(Axis.1, Axis.2), shape=21, size=1.5, color= "black") +
  theme_bw() + 
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Country of Origin") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())


p4 <- 
  ggplot(df12, aes(Axis.1, Axis.2)) + 
  geom_point(shape=21, size=1.5, alpha=0.8, aes(fill = disease_category, color=disease_category)) +
  geom_point(data = dplyr::filter(df12, dataset_name == "Boktor_Mazmanian"),
             aes(Axis.1, Axis.2), shape=21, size=1.5, color= "black") +
  theme_bw() + 
  xlab(paste0("PCoA 1 (", round((iMDS$values$Relative_eig[1])*100, digits = 2), "%)")) +
  ylab(paste0("PCoA 2 (", round((iMDS$values$Relative_eig[2])*100, digits = 2), "%)")) +
  labs(fill="Disease Category") +
  scale_colour_tableau(palette = "Classic 20") +
  scale_fill_tableau(palette = "Classic 20") +
  guides(color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank())



ord <- cowplot::plot_grid(p1, p2, p3, p4, nrow=2, align = "hv") 
ggsave(ord, filename = "data/Disease_Specificity_Meta_Analysis/Meta_Analysis_Beta_Diversity.svg",
       height = 9, width =14)



