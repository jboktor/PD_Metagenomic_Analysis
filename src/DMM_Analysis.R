# Dirichlet Multinomial Mixtures 

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Community_Composition_Funcs.R")

# load("files/Pfams.slim_PhyloseqObj.RData")
# load("files/Eggnogs.slim_PhyloseqObj.RData")


# Load example data
data(dietswap)
pseq <- dietswap

# To speed up, only consider the core taxa
# that are prevalent at 0.1% relative abundance in 50% of the samples
# (note that this is not strictly correct as information is
# being discarded; one alternative would be to aggregate rare taxa)
pseq.comp <- microbiome::transform(pseq, "compositional")
taxa <- core_members(pseq.comp, detection = 0.1/100, prevalence = 50/100)
pseq <- prune_taxa(taxa, pseq)


func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
reads$id <- gsub("_", ".", reads$id)


source("src/load_phyloseq_obj.R")
# Pick the OTU count matrix and convert it into samples x taxa format
# dat <- core(dat, detection = 0, prevalence = 0.1)
dat.pd = subset_samples(dat.path, donor_group =="PD")

count <- PseudoCounts(dat.pd, reads) %>% 
  t() %>% as.matrix()


# abund <- abundances(dat.pd)
# # dat <- abundances(pseq)
# count <- as.matrix(t(abund))

# Fit the DMM model. Let us set the maximum allowed number of community types to 3 to speed up the example.
fit <- lapply(1:6, dmn, count = count, verbose=TRUE)

# Check model fit with different number of mixture components using standard information criteria
lplc <- sapply(fit, laplace) # AIC / BIC / Laplace
aic  <- sapply(fit, AIC) # AIC / BIC / Laplace
bic  <- sapply(fit, BIC) # AIC / BIC / Laplace

plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
lines(aic, type="b", lty = 2)
lines(bic, type="b", lty = 3)

# Pick the optimal model
best <- fit[[which.min(unlist(lplc))]]
# Mixture parameters pi and theta
mixturewt(best)
# Sample-component assignments
ass <- apply(mixture(best), 1, which.max)


# Contribution of each taxonomic group to each component

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("Feature", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = Feature, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}








