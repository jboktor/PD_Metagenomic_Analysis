## Julia based

using FlashWeave

# Load  metadata files
meta_data_path = "files/FlashWeave_data/metadata_flashweave_friendly.csv"

# Load SPECIES data files
species_data_path = "files/FlashWeave_data/flashweave_friendly_species_relab.csv"
####### Species Model #######
species_results = learn_network(species_data_path, meta_data_path, sensitive=true, heterogeneous=false, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/species_network_output.edgelist", species_results)
save_network("data/Co_Abundance_Analysis/species_network_output.gml", species_results)


# Load GENUS data files
genus_data_path = "files/FlashWeave_data/flashweave_friendly_genus_relab.csv"
####### Species Model #######
genus_results = learn_network(genus_data_path, meta_data_path, sensitive=true, heterogeneous=false, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/genus_network_output.edgelist", genus_results)
save_network("data/Co_Abundance_Analysis/genus_network_output.gml", genus_results)


# Load PHYLUM data files
phylum_data_path = "files/FlashWeave_data/flashweave_friendly_phylum_relab.csv"
####### Species Model #######
phylum_results = learn_network(phylum_data_path, meta_data_path, sensitive=true, heterogeneous=false, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/phylum_network_output.edgelist", phylum_results)
save_network("data/Co_Abundance_Analysis/phylum_network_output.gml", phylum_results)


# Load PATHWAY data files
pathways_data_path = "files/FlashWeave_data/flashweave_friendly_pathways_relab.csv"
####### PATHWAY Model #######
pathways_results = learn_network(pathways_data_path, meta_data_path, sensitive=true, heterogeneous=false, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/pathways_network_output.edgelist", pathways_results)
save_network("data/Co_Abundance_Analysis/pathways_network_output.gml", pathways_results)


# Load ENZYME data files
enzymes_data_path = "files/FlashWeave_data/flashweave_friendly_enzymes_relab.csv"
####### ENZYME Model #######
enzymes_results = learn_network(enzymes_data_path, meta_data_path, sensitive=true, heterogeneous=true, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/enzymes_network_output.edgelist", enzymes_results)
save_network("data/Co_Abundance_Analysis/enzymes_network_output.gml", enzymes_results)


# Load Kegg Ortholog data files
kos_data_path = "files/FlashWeave_data/flashweave_friendly_kos_relab.csv"
####### KO Model #######
kos_results = learn_network(kos_data_path, meta_data_path, sensitive=true, heterogeneous=true, verbose=true, feed_forward= true, alpha=0.05, FDR=true)
save_network("data/Co_Abundance_Analysis/kos_network_output.edgelist", kos_results)
save_network("data/Co_Abundance_Analysis/kos_network_output.gml", kos_results)