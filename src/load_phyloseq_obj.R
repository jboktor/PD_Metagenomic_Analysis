# Load Phyloseq Objects

load_data <- function(folder){
  #' folders include: TBC, RUSH, SHANGHAI, Merged, Merged_ML, Merged_ML_Rarefied

  #### Taxa
  load(paste0("files/Phyloseq_", folder, "/Species_PhyloseqObj.RData"))
  cat("Loading: Species \n"); assign("dat.species", dat.species, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Genus_PhyloseqObj.RData"))
  cat("Loading: Genera \n"); assign("dat.genus", dat.genus, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Family_PhyloseqObj.RData"))
  cat("Loading: Families \n"); assign("dat.family", dat.family, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Order_PhyloseqObj.RData"))
  cat("Loading: Orders \n"); assign("dat.order", dat.order, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Class_PhyloseqObj.RData"))
  cat("Loading: Classes \n"); assign("dat.class", dat.class, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Phylum_PhyloseqObj.RData"))
  cat("Loading: Phyla \n"); assign("dat.phylum", dat.phylum, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Kingdom_PhyloseqObj.RData"))
  cat("Loading: Kindoms \n"); assign("dat.kingdom", dat.kingdom, envir = .GlobalEnv)
  ## Function - Pathways/Enzymes/Genes
  load(paste0("files/Phyloseq_", folder, "/Pathways_PhyloseqObj.RData"))
  cat("Loading: MetaCyc Pathways: Stratified \n"); assign("dat.path", dat.path, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Pathways.slim_PhyloseqObj.RData"))
  cat("Loading: MetaCyc Pathways \n"); assign("dat.path.slim", dat.path.slim, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Enzymes_PhyloseqObj.RData"))
  cat("Loading: Enzymes: Stratified \n"); assign("dat.ec", dat.ec, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/Enzymes.slim_PhyloseqObj.RData"))
  cat("Loading: Enzymes \n"); assign("dat.ec.slim", dat.ec.slim, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/KOs_PhyloseqObj.RData"))
  cat("Loading: Kegg Orthology: Stratified \n"); assign("dat.KOs", dat.KOs, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/KOs.slim_PhyloseqObj.RData"))
  cat("Loading: Kegg Orthology \n"); assign("dat.KOs.slim", dat.KOs.slim, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/GOs_PhyloseqObj.RData"))
  cat("Loading: Gene Ontology: Stratified \n"); assign("dat.GOs", dat.GOs, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/GOs.slim_PhyloseqObj.RData"))
  cat("Loading: Gene Ontology \n"); assign("dat.GOs.slim", dat.GOs.slim, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/PFAMs_PhyloseqObj.RData"))
  cat("Loading: PFAMs: Stratified \n"); assign("dat.PFAMs", dat.PFAMs, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/PFAMs.slim_PhyloseqObj.RData"))
  cat("Loading: PFAMs \n"); assign("dat.PFAMs.slim", dat.PFAMs.slim, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/EGGNOGs_PhyloseqObj.RData"))
  cat("Loading: eggNOGs: Stratified \n"); assign("dat.EGGNOGs", dat.EGGNOGs, envir = .GlobalEnv)
  load(paste0("files/Phyloseq_", folder, "/EGGNOGs.slim_PhyloseqObj.RData"))
  cat("Loading: eggNOGs \n"); assign("dat.EGGNOGs.slim", dat.EGGNOGs.slim, envir = .GlobalEnv)
}

