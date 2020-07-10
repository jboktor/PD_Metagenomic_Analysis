## Proudce Input files for FlashWeave Analysis

######## Load Data & functions
rm(list = ls())
source("src/load_phyloseq_obj.R")
source("src/DAF_Functions.R")
source("src/Metadata_prep_funcs.R")

######################### Write abundnance and metadata tables ######################### 

process_meta(dat)

env <- env %>% dplyr::select(donor_group)
write.csv(env, file = 'files/metadata_flashweave_friendly.csv', row.names = F)


## write abundance data
dat <- core(dat, detection = 0, prevalence = 0.1)
species.abundance.table <- microbiome::abundances(dat) %>% t()
write.csv(species.abundance.table, file = 'files/flashweave_friendly_species_relab.csv')

dat.genus <- core(dat.genus, detection = 0, prevalence = 0.1)
genus.abundance.table <- microbiome::abundances(dat.genus) %>% t
write.csv(genus.abundance.table, file = 'files/flashweave_friendly_genus_relab.csv')

dat.phylum <- core(dat.phylum, detection = 0, prevalence = 0.1)
phylum.abundance.table <- microbiome::abundances(dat.phylum) %>% t()
write.csv(phylum.abundance.table, file = 'files/flashweave_friendly_phylum_relab.csv')



dat.path.slim <- core(dat.path.slim, detection = 0, prevalence = 0.1)
features <- paste0("PATHWAY_", taxa_names(dat.path.slim))
features <- gsub(":", ".", features)
features <- gsub("\\|", ".", features)
features <- gsub(" ", "_", features)
features <- gsub("-", "_", features)
taxa_names(dat.path.slim) <- features
pathway.abundance.table <- microbiome::abundances(dat.path.slim) %>% t()
write.csv(pathway.abundance.table, file = 'files/flashweave_friendly_pathways_relab.csv')

dat.ec.slim <- core(dat.ec.slim, detection = 0, prevalence = 0.1)
enzyme.abundance.table <- microbiome::abundances(dat.ec.slim) %>% t()
write.csv(enzyme.abundance.table, file = 'files/flashweave_friendly_enzymes_relab.csv')

dat.KOs.slim <- core(dat.KOs.slim, detection = 0, prevalence = 0.1)
KOs.abundance.table <- microbiome::abundances(dat.KOs.slim) %>% t()
write.csv(KOs.abundance.table, file = 'files/flashweave_friendly_kos_relab.csv')



######################## Write df to bind to dictionary & color nodes ######################### 

features <- taxa_names(dat)
FW.meta <- paste0("donor_group_", unique(env$donor_group))
features <- c(features, FW.meta)
features <- data_frame("Node" = features)
FW.meta <- as.list(FW.meta)

############# Read-in Maaslin Files - all features used in significance testing ############# 

Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/Species_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/Species_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) 


## Complex color labels - DAF PD/PC (up & down) and PD/HC feature (up and down) coloring -

species_colors <- flashweave_group_colors(features, FW.meta, 
                        Maas.pd.pc.sig = Maas.pd.pc.sig, 
                        Maas.pd.hc.sig = Maas.pd.hc.sig)

# write.csv(species_colors, file = 'files/flashweave_metadata_network_colors.csv')





######################## Write df to bind to dictionary & color nodes ######################### 

features <- taxa_names(dat.path.slim)

FW.meta <- paste0("donor_group_", unique(env$donor_group))
features <- c(features, FW.meta)
features <- data_frame("Node" = features)
FW.meta <- as.list(FW.meta)

############# Read-in Maaslin Files - all features used in significance testing ############# 

Maas.pd.pc <- read_tsv(paste0("data/MaAsLin2_Analysis/Pathways.slim_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv(paste0("data/MaAsLin2_Analysis/Pathways.slim_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25) 


## Complex color labels - DAF PD/PC (up & down) and PD/HC feature (up and down) coloring -

pathway_colors <- flashweave_group_colors(features, FW.meta, 
                             Maas.pd.pc.sig = Maas.pd.pc.sig, 
                             Maas.pd.hc.sig = Maas.pd.hc.sig)

write.csv(pathway_colors, file = 'files/flashweave_metadata_pathways_network_colors.csv')

