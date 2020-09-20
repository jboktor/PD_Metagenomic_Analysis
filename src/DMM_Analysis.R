# Dirichlet Multinomial Mixtures 

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Community_Composition_Funcs.R")

# Load Read Counts
func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
reads$id <- gsub("_", ".", reads$id)

# Add prevalence threshold
dat.DMM <- dat %>% 
  core(detection = 0, prevalence = 0.1)

# Calculate Pseudocounts 
count <- PseudoCounts(dat.DMM, reads) %>% 
  t() %>% as.matrix()

# Fit the DMM model. Set the maximum allowed number of community types to 6 to speed up the analysis
set.seed(42)
fit <- lapply(1:6, dmn, count = count, verbose=TRUE)

# Check model fit with different number of mixture components using standard information criteria
lplc <- sapply(fit, laplace)
aic  <- sapply(fit, AIC) 
bic  <- sapply(fit, BIC) 

dmm.fit <- 
  data.frame(cluster = 1:length(lplc), 
             Laplace = lplc, AIC = aic, BIC = bic) %>% 
  pivot_longer(-cluster, names_to = "model.metrics")

dmm.fit.plot <- 
  ggplot(data = dmm.fit, aes(x = cluster, y = value, color = model.metrics)) +
  geom_point() +
  geom_line() +
  scale_color_d3() +
  labs(x="Number of Dirichlet Components", y="Model Fit", 
       color = "Model Metrics") +
  theme_classic() +
  theme(legend.position = c(0.2, 0.7))
dmm.fit.plot

ggsave(dmm.fit.plot,
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/DMM_cluster_fit.svg",
       width = 6, height = 4)


# Pick the optimal model
best.species <- fit[[which.min(unlist(lplc))]]
# Mixture parameters pi and theta
mixturewt(best.species)
# Sample-component assignments
sas <- apply(mixture(best.species), 1, which.max) %>% 
  as.data.frame() %>% dplyr::rename(cluster = ".")
sas <- group_col_from_ids(sas, rownames(sas))

# Summary Stats by group
sas.stats <- sas %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::count(group) %>% 
  transmute(n, group, Percentage=n/sum(n)*100)
sas.stats



cluster_distribution <- 
  ggplot(data = sas.stats, aes(x = as.factor(cluster), y = n, fill = group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Cluster", y = "Number of Samples") +
  scale_fill_manual(values = cols.pdpchc)
cluster_distribution

cluster_distribution2 <- 
  ggplot(data = sas.stats, aes(x = as.factor(cluster), y = Percentage, fill = group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Cluster", y = "Percentage of Cluster") +
  scale_fill_manual(values = cols.pdpchc)
cluster_distribution2

legend <- cowplot::plot_grid(get_legend(cluster_distribution))
cluster_distribution <- cluster_distribution + theme(legend.position = "none")
cluster_distribution2 <- cluster_distribution2 + theme(legend.position = "none")
cluster_distribution_final <- 
  cowplot::plot_grid(cluster_distribution, cluster_distribution2, legend,
                     ncol = 3, rel_widths = c(4,4,1))

ggsave(cluster_distribution_final,
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/cluster_distribution.svg",
       width = 5.5, height = 3.5)



# Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best.species)))) {
  # https://microbiome.github.io/tutorials/DMM.html
  
  d <- melt(fitted(best.species))
  colnames(d) <- c("feature", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange features by assignment strength
    arrange(value) %>%
    mutate(feature = factor(feature, levels = unique(feature))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.9))     
  
  d$feature <- gsub("s__", "", d$feature)
  d$feature <- gsub("_", " ", d$feature)
  
  p <- ggplot(d, aes(x = reorder(feature, value), y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"))
  print(p)
  ggsave(p,
         filename = paste0("data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/cluster_", k, "_drivers.svg"),
         width = 5, height = 5)
}


#--------------------------------------------------------------------------------
#           Features with greatest Change in Cluster occupancy probabilty 
#--------------------------------------------------------------------------------

d.k1 <- melt(fitted(best.species))
colnames(d.k1) <- c("feature", "cluster", "value.k1")
d.k1 <- subset(d.k1, cluster == "1") %>%
  mutate(feature = factor(feature, levels = unique(feature))) %>% 
  dplyr::select(-cluster)

d.k2 <- melt(fitted(best.species))
colnames(d.k2) <- c("feature", "cluster", "value.k2")
d.k2 <- subset(d.k2, cluster == "2") %>%
  mutate(feature = factor(feature, levels = unique(feature))) %>% 
  dplyr::select(-cluster)

dfm <- left_join(d.k1, d.k2, by = "feature")


#-------------------------------------
# Calculate Log2 Fold Change
dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)

#-------------------------------------
# Plot the top 20 most variable features

K2vK1.df <- 
  dfm %>% 
  arrange(desc(abs(K2vK1)))
K2vK1.df <- K2vK1.df[1:20,]

K2vK1.df$feature <- gsub("s__", "", K2vK1.df$feature)
# K2vK1.df$feature <- gsub("_", " ", K2vK1.df$feature)


LFC.K2vK1.plot <- ggplot(K2vK1.df, aes(x = reorder(feature, K2vK1), y = K2vK1)) +
  geom_bar(stat = "identity") +
  ylim(-8,8)+
  coord_flip() +
  labs(title = paste("Distinguishing features: Clusters 2/1"), y = expression(log[2]*" fold change")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"))
LFC.K2vK1.plot 

ggsave(LFC.K2vK1.plot, 
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/Species_LFC.K2vK1.plot.svg",
       width = 6, height = 4)


#--------------------------------------------------------------------------------------------
#                                  KO Data
#--------------------------------------------------------------------------------------------


# Add prevalence threshold
dat.DMM <- dat.KOs.slim %>% 
  core(detection = 0, prevalence = 0.1)

# Calculate Pseudocounts 
count <- PseudoCounts(dat.DMM, reads) %>% 
  t() %>% as.matrix()

# Fit the DMM model. Set the maximum allowed number of community types to 6 to speed up the analysis
set.seed(42)
fit <- lapply(1:6, dmn, count = count, verbose=TRUE)

# Check model fit with different number of mixture components using standard information criteria
lplc <- sapply(fit, laplace)
aic  <- sapply(fit, AIC) 
bic  <- sapply(fit, BIC) 

dmm.fit <- 
  data.frame(cluster = 1:length(lplc), 
             Laplace = lplc, AIC = aic, BIC = bic) %>% 
  pivot_longer(-cluster, names_to = "model.metrics")

dmm.fit.plot <- 
  ggplot(data = dmm.fit, aes(x = cluster, y = value, color = model.metrics)) +
  geom_point() +
  geom_line() +
  scale_color_d3() +
  labs(x="Number of Dirichlet Components", y="Model Fit", 
       color = "Model Metrics") +
  theme_classic() +
  theme(legend.position = c(0.2, 0.7))
dmm.fit.plot

ggsave(dmm.fit.plot,
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_DMM_cluster_fit.svg",
       width = 6, height = 4)


# Pick the optimal model
best <- fit[[which.min(unlist(lplc))]]
# Mixture parameters pi and theta
mixturewt(best)
# Sample-component assignments
sas <- apply(mixture(best), 1, which.max) %>% 
  as.data.frame() %>% dplyr::rename(cluster = ".")
sas <- group_col_from_ids(sas, rownames(sas))

# Summary Stats by group
sas.stats <- sas %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::count(group) %>% 
  transmute(n, group, Percentage=n/sum(n)*100)
sas.stats



cluster_distribution <- 
  ggplot(data = sas.stats, aes(x = as.factor(cluster), y = n, fill = group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Cluster", y = "Number of Samples") +
  scale_fill_manual(values = cols.pdpchc)
cluster_distribution

cluster_distribution2 <- 
  ggplot(data = sas.stats, aes(x = as.factor(cluster), y = Percentage, fill = group)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Cluster", y = "Percentage of Cluster") +
  scale_fill_manual(values = cols.pdpchc)
cluster_distribution2

legend <- cowplot::plot_grid(get_legend(cluster_distribution))
cluster_distribution <- cluster_distribution + theme(legend.position = "none")
cluster_distribution2 <- cluster_distribution2 + theme(legend.position = "none")
cluster_distribution_final <- 
  cowplot::plot_grid(cluster_distribution, cluster_distribution2, legend,
                     ncol = 3, rel_widths = c(4,4,1))

ggsave(cluster_distribution_final, 
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_cluster_distribution.svg",
       width = 5.5, height = 3.5)



# Contribution of each taxonomic group to each component
for (k in seq(ncol(fitted(best)))) {
  # https://microbiome.github.io/tutorials/DMM.html
  
  d <- melt(fitted(best))
  colnames(d) <- c("feature", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange features by assignment strength
    arrange(value) %>%
    mutate(feature = factor(feature, levels = unique(feature))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.995))     
  
  d$feature <- gsub("s__", "", d$feature)
  d$feature <- gsub("_", " ", d$feature)
  
  p <- ggplot(d, aes(x = reorder(feature, value), y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k)) +
    theme_bw() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(face = "italic"))
  print(p)
  ggsave(p,
         filename = paste0("data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_cluster_", k, "_drivers.svg"),
         width = 5, height = 5)
}


#--------------------------------------------------------------------------------
#           Features with greatest Change in Cluster occupancy probabilty 
#--------------------------------------------------------------------------------

d.k1 <- melt(fitted(best))
colnames(d.k1) <- c("feature", "cluster", "value.k1")
d.k1 <- subset(d.k1, cluster == "1") %>%
  mutate(feature = factor(feature, levels = unique(feature))) %>% 
  dplyr::select(-cluster)

d.k2 <- melt(fitted(best))
colnames(d.k2) <- c("feature", "cluster", "value.k2")
d.k2 <- subset(d.k2, cluster == "2") %>%
  mutate(feature = factor(feature, levels = unique(feature))) %>% 
  dplyr::select(-cluster)

d.k3 <- melt(fitted(best))
colnames(d.k3) <- c("feature", "cluster", "value.k3")
d.k3 <- subset(d.k3, cluster == "3") %>%
  mutate(feature = factor(feature, levels = unique(feature)))  %>% 
  dplyr::select(-cluster)

dfm <- left_join(d.k1, d.k2, by = "feature") %>% 
  left_join(d.k3, by = "feature")

#-------------------------------------
# Calculate Log2 Fold Change

dfm$K3vK1 <- log2(dfm$value.k3/dfm$value.k1)
dfm$K3vK2 <- log2(dfm$value.k3/dfm$value.k2)
dfm$K2vK1 <- log2(dfm$value.k2/dfm$value.k1)

#-------------------------------------
# Plot the top 20 most variable features

K3vK1.df <- 
  dfm %>% 
  arrange(desc(abs(K3vK1)))
K3vK1.df <- K3vK1.df[1:20,]

K3vK2.df <- 
  dfm %>% 
  arrange(desc(abs(K3vK2)))
K3vK2.df <- K3vK2.df[1:20,]

K2vK1.df <- 
  dfm %>% 
  arrange(desc(abs(K2vK1)))
K2vK1.df <- K2vK1.df[1:20,]


LFC.K3vK1.plot <- ggplot(K3vK1.df, aes(x = reorder(feature, K3vK1), y = K3vK1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Distinguishing features: Clusters 3/1"), y = expression(log[2]*" fold change")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"))
LFC.K3vK1.plot

LFC.K3vK2.plot <- ggplot(K3vK2.df, aes(x = reorder(feature, K3vK2), y = K3vK2)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Distinguishing features: Clusters 3/2"), y = expression(log[2]*" fold change")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"))
LFC.K3vK2.plot

LFC.K2vK1.plot <- ggplot(K2vK1.df, aes(x = reorder(feature, K2vK1), y = K2vK1)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = paste("Distinguishing features: Clusters 2/1"), y = expression(log[2]*" fold change")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"))
LFC.K2vK1.plot

ggsave(LFC.K3vK1.plot, 
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_LFC.K3vK1.plot.svg",
       width = 4, height = 4)
ggsave(LFC.K3vK2.plot, 
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_LFC.K3vK2.plot.svg",
       width = 4, height = 4)
ggsave(LFC.K2vK1.plot, 
       filename = "data/Community_Composition/Dirichlet_Mutlinomial_Mixtures/KO_LFC.K2vK1.plot.svg",
       width = 4, height = 4)
