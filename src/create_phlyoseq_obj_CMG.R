# create_phyloseq_obj_CMG

source("src/load_packages.R")
source("src/Community_Composition_Funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Community_Composition_Funcs.R")


reads <- load_reads() %>% 
  dplyr::rename(number_reads = clean_total_reads)
metadata <- read.csv('files/metadata_phyloseq_CMG.csv', header= TRUE) %>% 
  as.data.frame()
metadata <- left_join(metadata, reads, by = "id")
rownames(metadata) <- metadata$donor_id
metadata[is.na(metadata)] <- "not provided"

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


#----------------------------------------------------------
# Build CuratedMetagenomicData Phyloseq Obj
#----------------------------------------------------------

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

#----------------------------------------------------------
# Merge PD Metagenomics data
#----------------------------------------------------------
physeq <- merge_phyloseq(physeq, dat.CMG)
disease_abd <- microbiome::abundances(physeq)
disease_meta <- microbiome::meta(physeq)
disease_meta$studyID <- factor(disease_meta$studyID)

disease_meta$dataset_name <- sub(".metaphlan_bugs_list.stool", "", disease_meta$studyID)
disease_meta$study_condition <- factor(disease_meta$study_condition, 
                                       levels = c("control", "parkinsons", "ACVD", "BD", "AS", "CDI", "T1D", "hypertension", 
                                                  "metabolic_syndrome","IBD", "T2D", "cirrhosis"))


# Clear workspace
# rm(abd, metadata)

