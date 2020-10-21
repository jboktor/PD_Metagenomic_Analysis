# create_phyloseq_obj_CMG

source("src/load_packages.R")
source("src/Community_Composition_Funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")


  metadata <- read.csv('files/metadata_phyloseq_CMG.csv', header= TRUE) %>% 
    as.data.frame()
  metadata <- left_join(metadata, reads, by = "id")
  rownames(metadata) <- metadata$donor_id
  metadata[is.na(metadata)] <- "not provided"
  reads <- load_reads() %>% 
    dplyr::rename(number_reads = clean_total_reads)
  
  # Prep Abundance
  abd <- read.csv(file = "files/metaphlan2_taxonomic_profiles.csv", 
                  row.names = 1,  header= TRUE)
  # Filter Species
  abd <- abd %>% 
    tibble::rownames_to_column() %>% 
    filter(grepl("s__", rowname)) %>% 
    filter(!grepl("t__", rowname)) %>% 
    tibble::column_to_rownames()
  
  # source("https://raw.githubusercontent.com/waldronlab/presentations/master/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R")
  dat.CMG <- metaphlanToPhyloseq_Waldron(tax = abd, metadat = metadata)
  

  
# Clear workspace
# rm(abd, metadata)

