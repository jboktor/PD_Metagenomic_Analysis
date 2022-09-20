# Dirichlet Multinomial Mixtures 

rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/community_composition_funcs.R")
load("files/low_quality_samples.RData")


# Load Read Counts
# func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
# reads <- dplyr::select(func_reads, c("# samples", "total reads")) %>% 
#   dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
# reads$id <- gsub("_", ".", reads$id)

# remove_dats()
# load_tbc()
# reads <- load_reads("TBC")

# load_all_cohorts()
# reads <-load_reads("Merged")

# remove_dats()
# load_data("TBC")

DMM_analysis <- function(cohort){
  
  remove_dats()
  load_data(cohort)
  
  #--------------------------------------------------------------------------------------------
  #                                  Species
  #--------------------------------------------------------------------------------------------
  dat <- dat.species %>% 
    subset_samples(donor_id %ni% low_qc[[1]])
  
  DMM.species <- DMM_fit(dat, nmax = 6)
  fit.species <- DMM.species[["fit"]]
  DMM.species[["laplace"]]
  DMM.species[["plot"]]
  
  ggsave(DMM.species[["plot"]],
         filename = paste0("data/Community_Composition/Dirichlet_Multinomial_Mixtures/", 
                           cohort ,"/Species/DMM_cluster_fit_species.svg"),
         width = 6, height = 4)
  
  # Pick the optimal model
  best.species <- fit.species[[which.min(unlist( DMM.species[["laplace"]] ))]]
  # Mixture parameters pi and theta
  mixturewt(best.species)
  # extract summary stats
  sas.stats <- DMM_stats(best.species)[["sas.stats"]]
  # SAVE cluster ID for each sample
  DMM_stats(best.species)[["sas"]] %>% 
    write.csv(file = paste0('data/Community_Composition/Dirichlet_Multinomial_Mixtures/', 
                            cohort ,'/DMM_sample_cluster_mapping/species_cluster_IDs.csv'))
  # Plot cluster distribution
  cluster_distribution_species <- DMM_cluster_plot(sas.stats)
  ggsave(cluster_distribution_species,
         filename = paste0("data/Community_Composition/Dirichlet_Multinomial_Mixtures", 
                           cohort ,"/Species/cluster_distribution_species.svg"),
         width = 5.5, height = 3.5)
  
  #----------------------------------------------------------------
  # Features with greatest Change in Cluster occupancy probability 
  #----------------------------------------------------------------
  melt.species <- melt(fitted(best.species))
  d.k1 <- DMM_select_cluster(df = melt.species, cluster_n = 1)
  d.k2 <- DMM_select_cluster(df = melt.species, cluster_n = 2)
  dfm <- left_join(d.k1, d.k2, by = "feature")
  # Calculate Log2 Fold Change
  dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)
  
  #-----------------------------------------
  # Plot the top 20 most variable features
  #-----------------------------------------
  K2vK1.df <- 
    dfm %>% 
    arrange(desc(abs(K2vK1)))
  K2vK1.df <- K2vK1.df[1:20,]
  K2vK1.df$feature <- gsub("s__", "", K2vK1.df$feature)
  LFC.K2vK1.plot <- DMM_cluster_driver_plot(K2vK1.df, yval =  K2vK1.df$K2vK1, comparison = "2/1")
  
  ggsave(LFC.K2vK1.plot, 
         filename = paste0("data/Community_Composition/Dirichlet_Multinomial_Mixtures/", 
                           cohort ,"/Species/Species_LFC.K2vK1.plot.svg"),
         width = 6, height = 4)
  
}





#--------------------------------------------------------------------------------------------
#                                  Species
#--------------------------------------------------------------------------------------------
dat <- dat.species %>%
  subset_samples(donor_id %ni% low_qc[[1]])

DMM.species <- DMM_fit(dat, nmax = 6)
fit.species <- DMM.species[["fit"]]
DMM.species[["laplace"]]
DMM.species[["plot"]]
# 
# # ggsave(DMM.species[["plot"]],
# #        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/Species/DMM_cluster_fit_species.svg",
# #        width = 6, height = 4)
# 
# # Pick the optimal model
# best.species <- fit.species[[which.min(unlist( DMM.species[["laplace"]] ))]]
# # Mixture parameters pi and theta
# mixturewt(best.species)
# # extract summary stats
# sas.stats <- DMM_stats(best.species)[["sas.stats"]]
# # SAVE cluster ID for each sample
# DMM_stats(best.species)[["sas"]] %>% 
#   write.csv(file = 'data/Community_Composition/Dirichlet_Multinomial_Mixtures/DMM_sample_cluster_mapping/species_cluster_IDs.csv')
# # Plot cluster distribution
# cluster_distribution_species <- DMM_cluster_plot(sas.stats)
# ggsave(cluster_distribution_species,
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/Species/cluster_distribution_species.svg",
#        width = 5.5, height = 3.5)
# 
# #----------------------------------------------------------------
# # Features with greatest Change in Cluster occupancy probability 
# #----------------------------------------------------------------
# melt.species <- melt(fitted(best.species))
# d.k1 <- DMM_select_cluster(df = melt.species, cluster_n = 1)
# d.k2 <- DMM_select_cluster(df = melt.species, cluster_n = 2)
# dfm <- left_join(d.k1, d.k2, by = "feature")
# # Calculate Log2 Fold Change
# dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)
# 
# #-----------------------------------------
# # Plot the top 20 most variable features
# #-----------------------------------------
# K2vK1.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK1)))
# K2vK1.df <- K2vK1.df[1:20,]
# K2vK1.df$feature <- gsub("s__", "", K2vK1.df$feature)
# LFC.K2vK1.plot <- DMM_cluster_driver_plot(K2vK1.df, yval =  K2vK1.df$K2vK1, comparison = "2/1")
# 
# ggsave(LFC.K2vK1.plot, 
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/Species/Species_LFC.K2vK1.plot.svg",
#        width = 6, height = 4)
# 
# 
# #--------------------------------------------------------------------------------------------
# #                                  Pathways
# #--------------------------------------------------------------------------------------------
# 
# DMM.pathways <- DMM_fit(dat.path.slim, nmax = 6)
# fit.pathways <- DMM.pathways[["fit"]]
# DMM.pathways[["laplace"]]
# DMM.pathways[["plot"]]
# 
# ggsave(DMM.pathways[["plot"]],
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/Pathways/Pathways_DMM_cluster_fit.svg",
#        width = 6, height = 4)
# 
# # Pick the optimal model
# best.paths <- fit.pathways[[which.min(unlist(DMM.pathways[["laplace"]]))]]
# # Mixture parameters pi and theta
# mixturewt(best.paths)
# # extract summary stats
# sas.stats <- DMM_stats(best.paths)[["sas.stats"]]
# # SAVE cluster ID for each sample
# DMM_stats(best.paths)[["sas"]] %>% 
#   write.csv(file = 'data/Community_Composition/Dirichlet_Multinomial_Mixtures/DMM_sample_cluster_mapping/pathway_cluster_IDs.csv')
# # Plot cluster distribution
# cluster_distribution_paths <- DMM_cluster_plot(sas.stats)
# ggsave(cluster_distribution_paths,
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/Pathways/Pathways_cluster_distribution.svg",
#        width = 5.5, height = 3.5)
# 
# 
# #----------------------------------------------------------------
# # Features with greatest Change in Cluster occupancy probability 
# #----------------------------------------------------------------
# 
# melt.paths <- melt(fitted(best.paths))
# d.k1 <- DMM_select_cluster(df = melt.paths, cluster_n = 1)
# d.k2 <- DMM_select_cluster(df = melt.paths, cluster_n = 2)
# d.k3 <- DMM_select_cluster(df = melt.paths, cluster_n = 3)
# d.k4 <- DMM_select_cluster(df = melt.paths, cluster_n = 4)
# d.k5 <- DMM_select_cluster(df = melt.paths, cluster_n = 5)
# 
# dfm <- left_join(d.k1, d.k2, by = "feature") %>% 
#   left_join(d.k3, by = "feature") %>% 
#   left_join(d.k4, by = "feature") %>% 
#   left_join(d.k5, by = "feature")
# 
# # Calculate Log2 Fold Change
# dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)
# dfm$K2vK3 <- log2(dfm$value.k2/dfm$value.k3)
# dfm$K2vK4 <- log2(dfm$value.k2/dfm$value.k4)
# dfm$K2vK5 <- log2(dfm$value.k2/dfm$value.k5)
# 
# 
# # Plot the top 20 most variable features
# K2vK1.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK1)))
# K2vK1.df <- K2vK1.df[1:20,]
# 
# K2vK3.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK3)))
# K2vK3.df <- K2vK3.df[1:20,]
# 
# K2vK4.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK4)))
# K2vK4.df <- K2vK4.df[1:20,]
# 
# K2vK5.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK5)))
# K2vK5.df <- K2vK5.df[1:20,]
# 
# 
# # PLOTS
# K2vK1.plot <- DMM_cluster_driver_plot(K2vK1.df, yval =  K2vK1.df$K2vK1, comparison = "2/1")
# K2vK3.plot <- DMM_cluster_driver_plot(K2vK3.df, yval =  K2vK1.df$K2vK1, comparison = "2/3")
# K2vK4.plot <- DMM_cluster_driver_plot(K2vK4.df, yval =  K2vK1.df$K2vK1, comparison = "2/4")
# K2vK5.plot <- DMM_cluster_driver_plot(K2vK5.df, yval =  K2vK1.df$K2vK1, comparison = "2/5")
# 
# plts <- list(K2vK1.plot, K2vK3.plot, K2vK4.plot, K2vK5.plot)
# names(plts) <- c("K2vK1.plot", "K2vK3.plot", "K2vK4.plot", "K2vK5.plot")
# 
# cnt <- 1
# for (i in plts){
#   ggsave(i,
#          filename = paste0("data/Community_Composition/Dirichlet_Multinomial_Mixtures/Pathways/Pathways_LFC.", names(plts[cnt]), ".svg"),
#          width = 9, height = 9)
#   cnt <- cnt + 1
# }
# 
# 
# #--------------------------------------------------------------------------------------------
# #                                  KO Data
# #--------------------------------------------------------------------------------------------
# 
# DMM.kos <- DMM_fit(dat.KOs.slim, nmax = 6)
# fit.kos <- DMM.kos[["fit"]]
# DMM.kos[["laplace"]]
# DMM.kos[["plot"]]
# 
# ggsave(DMM.kos[["plot"]],
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/KOs/KO_DMM_cluster_fit.svg",
#        width = 6, height = 4)
# 
# # Pick the optimal model
# best.kos <- fit.kos[[which.min(unlist(DMM.kos[["laplace"]]))]]
# # Mixture parameters pi and theta
# mixturewt(best.kos)
# # extract summary stats
# sas.stats <- DMM_stats(best.kos)[["sas.stats"]]
# # SAVE cluster ID for each sample
# DMM_stats(best.kos)[["sas"]] %>% 
#   write.csv(file = 'data/Community_Composition/Dirichlet_Multinomial_Mixtures/DMM_sample_cluster_mapping/ko_cluster_IDs.csv')
# # Plot cluster distribution
# cluster_distribution_kos <- DMM_cluster_plot(sas.stats)
# ggsave(cluster_distribution_kos, 
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/KOs/KO_cluster_distribution.svg",
#        width = 5.5, height = 3.5)
# 
# 
# #----------------------------------------------------------------
# # Features with greatest Change in Cluster occupancy probability 
# #----------------------------------------------------------------
# melt.kos <- melt(fitted(best.kos))
# d.k1 <- DMM_select_cluster(df = melt.kos, cluster_n = 1)
# d.k2 <- DMM_select_cluster(df = melt.kos, cluster_n = 2)
# d.k3 <- DMM_select_cluster(df = melt.kos, cluster_n = 3)
# 
# dfm <- left_join(d.k1, d.k2, by = "feature") %>%
#   left_join(d.k3, by = "feature")
# 
# # Calculate Log2 Fold Change
# dfm$K3vK1 <- log2(dfm$value.k3/dfm$value.k1)
# dfm$K3vK2 <- log2(dfm$value.k3/dfm$value.k2)
# dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)
# 
# # Plot the top 20 most variable features
# K3vK1.df <- 
#   dfm %>% 
#   arrange(desc(abs(K3vK1)))
# K3vK1.df <- K3vK1.df[1:20,]
# 
# K3vK2.df <- 
#   dfm %>% 
#   arrange(desc(abs(K3vK2)))
# K3vK2.df <- K3vK2.df[1:20,]
# 
# K2vK1.df <- 
#   dfm %>% 
#   arrange(desc(abs(K2vK1)))
# K2vK1.df <- K2vK1.df[1:20,]
# 
# 
# # PLOTS
# LFC.K3vK1.plot <- DMM_cluster_driver_plot(K3vK1.df, yval =  K3vK1.df$K3vK1, comparison = "3/1")
# LFC.K3vK2.plot <- DMM_cluster_driver_plot(K3vK2.df, yval =  K3vK2.df$K3vK2, comparison = "3/2")
# LFC.K2vK1.plot <- DMM_cluster_driver_plot(K2vK1.df, yval =  K2vK1.df$K2vK1, comparison = "2/1")
# 
# ggsave(LFC.K3vK1.plot, 
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/KOs/KO_LFC.K3vK1.plot.svg",
#        width = 4, height = 4)
# ggsave(LFC.K3vK2.plot, 
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/KOs/KO_LFC.K3vK2.plot.svg",
#        width = 4, height = 4)
# ggsave(LFC.K2vK1.plot, 
#        filename = "data/Community_Composition/Dirichlet_Multinomial_Mixtures/KOs/KO_LFC.K2vK1.plot.svg",
#        width = 4, height = 4)




#---------------------------------------------------------------------------------------------
# Scrap
#---------------------------------------------------------------------------------------------

# # Contribution of each taxonomic group to each component
# for (k in seq(ncol(fitted(best.species)))) {
#   # https://microbiome.github.io/tutorials/DMM.html
#   
#   d <- melt(fitted(best.species))
#   colnames(d) <- c("feature", "cluster", "value")
#   d <- subset(d, cluster == k) %>%
#     # Arrange features by assignment strength
#     arrange(value) %>%
#     mutate(feature = factor(feature, levels = unique(feature))) %>%
#     # Only show the most important drivers
#     filter(abs(value) > quantile(abs(value), 0.9))     
#   
#   d$feature <- gsub("s__", "", d$feature)
#   d$feature <- gsub("_", " ", d$feature)
#   
#   p <- ggplot(d, aes(x = reorder(feature, value), y = value)) +
#     geom_bar(stat = "identity") +
#     coord_flip() +
#     labs(title = paste("Top drivers: community type", k)) +
#     theme_bw() +
#     theme(axis.title.y = element_blank(),
#           axis.text.y = element_text(face = "italic"))
#   print(p)
#   ggsave(p,
#          filename = paste0("data/Community_Composition/Dirichlet_Multinomial_Mixtures/cluster_", k, "_drivers.svg"),
#          width = 5, height = 5)
# }
