## Phyloseq Obj Creation - HUMANn2 REGROUP OBJECTS

############  load packages  ############
source("src/ load_packages.R")
source("src/load_phyloseq_obj.R")

### Metadata 
my_sample_data <- meta(dat) %>% sample_data()


build.pfams <- function(strat = F){
  #################### PFAMS #################### 
  # Regrouped Gene Data
  Pfams.abund.all <- read_tsv("files/genefamilies_relab_rescaled_Pfam.tsv", col_names = T)       #All genes (Pfams and ungrouped)
  Pfams.abund <- filter(Pfams.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only Pfams 
  Pfams.abund.slim <- filter(Pfams.abund, !grepl("g__", `# Gene Family`))      
  Pfams.abund.slim <- filter(Pfams.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only Pfams with no stratification
  # Pfams with Strat 
  colnames(Pfams.abund) <- gsub("_Abundance-RPKs", "", colnames(Pfams.abund))
  colnames(Pfams.abund) <- gsub("_", ".", colnames(Pfams.abund))
  Pfams.abund <- Pfams.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.ab_table <- otu_table(Pfams.abund, taxa_are_rows=T)
  dat.Pfams <- phyloseq(my_Pfam.ab_table, my_sample_data)
  print(dat.Pfams)
  save(dat.Pfams, file = "files/Pfams_PhyloseqObj.RData")
  #Pfams no stratification
  colnames(Pfams.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(Pfams.abund.slim))
  colnames(Pfams.abund.slim) <- gsub("_", ".", colnames(Pfams.abund.slim))
  Pfams.abund.slim <- Pfams.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.slim.ab_table <- otu_table(Pfams.abund.slim, taxa_are_rows=T)
  dat.Pfams.slim <- phyloseq(my_Pfam.slim.ab_table, my_sample_data)
  print(dat.Pfams.slim)
  save(dat.Pfams.slim, file = "files/Pfams.slim_PhyloseqObj.RData")
  
  if (strat == F){
    assign("dat.Pfams.slim", dat.Pfams.slim, envir = .GlobalEnv)
  } else if (strat == T){
    assign("dat.Pfams", dat.Pfams, envir = .GlobalEnv)
  }
  
}

build.eggnogs <- function(strat = F) {
  #################### EggNogs #################### 
  # Regrouped Gene Data
  Eggnogs.abund.all <- read_tsv("files/genefamilies_relab_rescaled_Eggnog.tsv", col_names = T)       #All genes (Eggnogs and ungrouped)
  Eggnogs.abund <- filter(Eggnogs.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only Eggnogs 
  Eggnogs.abund.slim <- filter(Eggnogs.abund, !grepl("g__", `# Gene Family`))      
  Eggnogs.abund.slim <- filter(Eggnogs.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only Eggnogs with no stratification
  # Eggnogs with Strat 
  colnames(Eggnogs.abund) <- gsub("_Abundance-RPKs", "", colnames(Eggnogs.abund))
  colnames(Eggnogs.abund) <- gsub("_", ".", colnames(Eggnogs.abund))
  Eggnogs.abund <- Eggnogs.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.ab_table <- otu_table(Eggnogs.abund, taxa_are_rows=T)
  dat.Eggnogs <- phyloseq(my_Pfam.ab_table, my_sample_data)
  print(dat.Eggnogs)
  save(dat.Eggnogs, file = "files/Eggnogs_PhyloseqObj.RData")
  #Eggnogs no stratification
  colnames(Eggnogs.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(Eggnogs.abund.slim))
  colnames(Eggnogs.abund.slim) <- gsub("_", ".", colnames(Eggnogs.abund.slim))
  Eggnogs.abund.slim <- Eggnogs.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.slim.ab_table <- otu_table(Eggnogs.abund.slim, taxa_are_rows=T)
  dat.Eggnogs.slim <- phyloseq(my_Pfam.slim.ab_table, my_sample_data)
  print(dat.Eggnogs.slim)
  save(dat.Eggnogs.slim, file = "files/Eggnogs.slim_PhyloseqObj.RData")

  if (strat == F){
    assign("dat.Eggnogs.slim", dat.Eggnogs.slim, envir = .GlobalEnv)
  } else if (strat == T){
    assign("dat.Eggnogs", dat.Eggnogs, envir = .GlobalEnv)
  }
}

build.GOs <- function(strat = F) {
  #################### GOs #################### 
  # Regrouped Gene Data
  GOs.abund.all <- read_tsv("files/genefamilies_relab_rescaled_GO.tsv", col_names = T)       #All genes (GOs and ungrouped)
  GOs.abund <- filter(GOs.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only GOs 
  GOs.abund.slim <- filter(GOs.abund, !grepl("g__", `# Gene Family`))      
  GOs.abund.slim <- filter(GOs.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only GOs with no stratification
  # GOs with Strat 
  colnames(GOs.abund) <- gsub("_Abundance-RPKs", "", colnames(GOs.abund))
  colnames(GOs.abund) <- gsub("_", ".", colnames(GOs.abund))
  GOs.abund <- GOs.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.ab_table <- otu_table(GOs.abund, taxa_are_rows=T)
  dat.GOs <- phyloseq(my_Pfam.ab_table, my_sample_data)
  print(dat.GOs)
  save(dat.GOs, file = "files/GOs_PhyloseqObj.RData")
  #GOs no stratification
  colnames(GOs.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(GOs.abund.slim))
  colnames(GOs.abund.slim) <- gsub("_", ".", colnames(GOs.abund.slim))
  GOs.abund.slim <- GOs.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.slim.ab_table <- otu_table(GOs.abund.slim, taxa_are_rows=T)
  dat.GOs.slim <- phyloseq(my_Pfam.slim.ab_table, my_sample_data)
  print(dat.GOs.slim)
  save(dat.GOs.slim, file = "files/GOs.slim_PhyloseqObj.RData")

  if (strat == F){
    assign("dat.GOs.slim", dat.GOs.slim, envir = .GlobalEnv)
  } else if (strat == T){
    assign("dat.GOs", dat.GOs, envir = .GlobalEnv)
  }
}

build.Rxns <- function(strat = F) {
  #################### Rxns #################### 
  # Regrouped Gene Data
  Rxns.abund.all <- read_tsv("files/genefamilies_relab_rescaled_Rxn.tsv", col_names = T)       #All genes (Rxns and ungrouped)
  Rxns.abund <- filter(Rxns.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only Rxns 
  Rxns.abund.slim <- filter(Rxns.abund, !grepl("g__", `# Gene Family`))      
  Rxns.abund.slim <- filter(Rxns.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only Rxns with no stratification
  # Rxns with Strat 
  colnames(Rxns.abund) <- gsub("_Abundance-RPKs", "", colnames(Rxns.abund))
  colnames(Rxns.abund) <- gsub("_", ".", colnames(Rxns.abund))
  Rxns.abund <- Rxns.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.ab_table <- otu_table(Rxns.abund, taxa_are_rows=T)
  dat.Rxns <- phyloseq(my_Pfam.ab_table, my_sample_data)
  print(dat.Rxns)
  save(dat.Rxns, file = "files/Rxns_PhyloseqObj.RData")
  #Rxns no stratification
  colnames(Rxns.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(Rxns.abund.slim))
  colnames(Rxns.abund.slim) <- gsub("_", ".", colnames(Rxns.abund.slim))
  Rxns.abund.slim <- Rxns.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.slim.ab_table <- otu_table(Rxns.abund.slim, taxa_are_rows=T)
  dat.Rxns.slim <- phyloseq(my_Pfam.slim.ab_table, my_sample_data)
  print(dat.Rxns.slim)
  save(dat.Rxns.slim, file = "files/Rxns.slim_PhyloseqObj.RData")

  if (strat == F){
    assign("dat.Rxns.slim", dat.Rxns.slim, envir = .GlobalEnv)
  } else if (strat == T){
    assign("dat.Rxns", dat.Rxns, envir = .GlobalEnv)
  }
  
}

build.InfoGOs <- function(strat = F) {
  #################### InfoGOs #################### 
  # Regrouped Gene Data
  InfoGOs.abund.all <- read_tsv("files/genefamilies_relab_rescaled_Infogo.tsv", col_names = T)       #All genes (InfoGOs and ungrouped)
  InfoGOs.abund <- filter(InfoGOs.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only InfoGOs 
  InfoGOs.abund.slim <- filter(InfoGOs.abund, !grepl("g__", `# Gene Family`))      
  InfoGOs.abund.slim <- filter(InfoGOs.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only InfoGOs with no stratification
  # InfoGOs with Strat 
  colnames(InfoGOs.abund) <- gsub("_Abundance-RPKs", "", colnames(InfoGOs.abund))
  colnames(InfoGOs.abund) <- gsub("_", ".", colnames(InfoGOs.abund))
  InfoGOs.abund <- InfoGOs.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.ab_table <- otu_table(InfoGOs.abund, taxa_are_rows=T)
  dat.InfoGOs <- phyloseq(my_Pfam.ab_table, my_sample_data)
  print(dat.InfoGOs)
  save(dat.InfoGOs, file = "files/InfoGOs_PhyloseqObj.RData")
  #InfoGOs no stratification
  colnames(InfoGOs.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(InfoGOs.abund.slim))
  colnames(InfoGOs.abund.slim) <- gsub("_", ".", colnames(InfoGOs.abund.slim))
  InfoGOs.abund.slim <- InfoGOs.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
  my_Pfam.slim.ab_table <- otu_table(InfoGOs.abund.slim, taxa_are_rows=T)
  dat.InfoGOs.slim <- phyloseq(my_Pfam.slim.ab_table, my_sample_data)
  print(dat.InfoGOs.slim)
  save(dat.InfoGOs.slim, file = "files/InfoGOs.slim_PhyloseqObj.RData")

  if (strat == F){
    assign("dat.InfoGOs.slim", dat.InfoGOs.slim, envir = .GlobalEnv)
  } else if (strat == T){
    assign("dat.InfoGOs", dat.InfoGOs, envir = .GlobalEnv)
  }
}



