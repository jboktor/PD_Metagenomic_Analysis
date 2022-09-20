## Proudce Input files for FlashWeave Analysis

######## Load Data & functions
rm(list = ls())
source("src/load_phyloseq_obj.R")
source("src/DAF_Functions.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")

######################### Write abundnance and metadata tables #########################

process_meta(dat)

env <- env %>% dplyr::select(donor_group)
write.csv(env, file = "files/FlashWeave_data/metadata_flashweave_friendly.csv", row.names = F)


########## write abundance data for FW input
# Species
dat <- core(dat, detection = 0, prevalence = 0.1)
species.abundance.table <- microbiome::abundances(dat) %>% t()
write.csv(species.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_species_relab.csv")
# Genus
dat.genus <- core(dat.genus, detection = 0, prevalence = 0.1)
genus.abundance.table <- microbiome::abundances(dat.genus) %>% t()
write.csv(genus.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_genus_relab.csv")
# Phylum
dat.phylum <- core(dat.phylum, detection = 0, prevalence = 0.1)
phylum.abundance.table <- microbiome::abundances(dat.phylum) %>% t()
write.csv(phylum.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_phylum_relab.csv")
# Pathways
dat.path.slim <- core(dat.path.slim, detection = 0, prevalence = 0.1)
taxa_names(dat.path.slim) <- prep_pathway_names(dat.path.slim)
pathway.abundance.table <- microbiome::abundances(dat.path.slim) %>% t()
write.csv(pathway.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_pathways_relab.csv")
# Enzymes
dat.ec.slim <- core(dat.ec.slim, detection = 0, prevalence = 0.1)
taxa_names(dat.ec.slim) <- prep_enzyme_names(dat.ec.slim)
enzyme.abundance.table <- microbiome::abundances(dat.ec.slim) %>% t()
write.csv(enzyme.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_enzymes_relab.csv")
# KOs
dat.KOs.slim <- core(dat.KOs.slim, detection = 0, prevalence = 0.1)
taxa_names(dat.KOs.slim) <- prep_ko_names(dat.KOs.slim)
KOs.abundance.table <- microbiome::abundances(dat.KOs.slim) %>% t()
write.csv(KOs.abundance.table, file = "files/FlashWeave_data/flashweave_friendly_kos_relab.csv")



feats <- c("Species", "Genus", "Phylum", "Pathways.slim", "Enzymes.slim", "KOs.slim")
feat_objs <- c(dat, dat.genus, dat.phylum, dat.path.slim, dat.ec.slim, dat.KOs.slim)

i <- 1
for (f in feats) {
  print(f)

  features <- taxa_names(feat_objs[[i]])

  FW.meta <- paste0("donor_group_", unique(env$donor_group))
  features <- c(features, FW.meta)
  features <- data_frame("Node" = features)
  FW.meta <- as.list(FW.meta)

  Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/", f, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>%
    filter(value == "Population Control")
  Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

  Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/", f, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>%
    filter(value == "Household Control")
  Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25)

  network_colors <- flashweave_group_colors(features, FW.meta,
    Maas.pd.pc.sig = Maas.pd.pc.sig,
    Maas.pd.hc.sig = Maas.pd.hc.sig
  )

  write.csv(network_colors, file = paste0("files/FlashWeave_data/", f, "_flashweave_metadata_network_colors.csv"))

  i <- i + 1
}
