# Scratch code

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")

negative_controls <- c(
  "S00A4-ATCC_MSA_1003_S96",
  "S00A4-neg2_S119",
  "S00A4-neg3_S125",
  "S00A4-neg_S118",
  "S00A4NegExt_P00A4_S94",
  "S00A4NegH2O_P00A4_S95",
  "S00A4_stagPos_S117",
  "BLANK"
  )

# #   SCRAP
# # ID Key Maps
# ppmi <-
#   read.csv(file = "Whole_Blood_Sample_Collection_PPMI.csv", header = TRUE) 
# # %>% 
# #   janitor::get_dupes(PATNO)
# length(unique(ppmi$PATNO))




# Load UniRef90 Mapping to Entity Names


# Load UniprotKB hits from Hidden Markov Model
csgA_hits <-
  read_tsv(
    "files/Bacterial_Amyloids/0FA1E7F0-6F0B-11EB-AE5D-187D53F04F9B.1.tsv",
    col_names = T)

# ID Key Maps
TBC_keys <-
  read.csv(file = "files/metadata_keys.csv", header = TRUE) %>%
  dplyr::select(c(MBI_Sample_ID, id)) %>%
  mutate(id = gsub("_", ".", id)) %>%
  mutate(MBI_Sample_ID = as.character(MBI_Sample_ID)) %>%
  mutate(id = as.character(id))
TBC_keymap <- TBC_keys$MBI_Sample_ID
names(TBC_keymap) <- TBC_keys$id

metadata_RUSH <-
  read.csv(file = "files/metadata_phyloseq_RUSH.csv", header = TRUE)
RUSH_keys <-
  metadata_RUSH %>%
  select(donor_id, host_subject_id) %>%
  mutate(donor_id = as.character(donor_id)) %>%
  mutate(host_subject_id = as.character(host_subject_id))
RUSH_keymap <- RUSH_keys$host_subject_id
names(RUSH_keymap) <- RUSH_keys$donor_id

# Load UniRef90 Gene hits
uniref90_TBC <-
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/genefamilies_relab.tsv",
    col_names = T)  %>% 
  dplyr::rename(clusterID = `# Gene Family`) %>% 
  dplyr::mutate(clusterID = gsub("UniRef90_", "", clusterID)) %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_RPK() %>% 
  column_to_rownames(var = "clusterID") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  rownames_to_column(var = "clusterID")

uniref90_Rush <-
  read_tsv(
    "files/RUSH_biobakery_output_slim/humann/merged/genefamilies_relab.tsv",
    col_names = T) %>% 
  dplyr::rename(clusterID = `# Gene Family`) %>% 
  dplyr::mutate(clusterID = gsub("UniRef90_", "", clusterID)) %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_RPK() %>% 
  column_to_rownames(var = "clusterID") %>% 
  trim_cols("RUSH") %>%
  dplyr::rename(RUSH_keymap) %>% 
  dplyr::select(-contains("MSA")) %>% 
  rownames_to_column(var = "clusterID")


# ------------------------------------------------------------------------------
#                      Data Wrangling Functions 
# ------------------------------------------------------------------------------

hmmsearch_output_preprocessing <- function(df.hits){
  df.temp <-
    as.data.frame(lapply(df.hits, as.character), stringsAsFactors = FALSE)
  for (row in 1:nrow(df.temp)) {
    if (row %% 2 == 0) {
      cat("Processing: ", row, "\n")
      df.temp[row-1, 6:25] <- df.temp[row, 3:22]
    }
  }
  hit_output <- 
    df.temp %>% 
    dplyr::filter(row_number() %% 2 == 1) %>% 
    dplyr::mutate(clusterID = substr(Target.Accession, 1, nchar(Target.Accession)-6)) %>% 
    dplyr::relocate(clusterID, .after = Target.Accession)
  return(hit_output)
}




# ------------------------------------------------------------------------------

csgA_hmm_hits <- hmmsearch_output_preprocessing(csgA_hits)

csgA_uniref90_TBC <-
  uniref90_TBC %>%
  filter(clusterID %in% csgA_hmm_hits$clusterID) %>% 
  column_to_rownames(var = "clusterID") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "donor_id")

csgA_uniref90_Rush <-
  uniref90_Rush %>%
  filter(clusterID %in% csgA_hmm_hits$clusterID) %>% 
  column_to_rownames(var = "clusterID") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column(var = "donor_id") 

csga_total <- 
  dplyr::bind_rows(csgA_uniref90_TBC, csgA_uniref90_Rush) %>%
  rowwise() %>%
  mutate(csgA_total = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
  dplyr::mutate(group = if_else(grepl("HC", donor_id), "HC",
                                if_else(grepl("PC", donor_id), "PC", "PD")))

boxplot_all(df = csga_total,
            x = csga_total$group,
            y = csga_total$csgA_total)


load_all_cohorts()

KOs <- 
  dat.KOs.slim %>% 
  abundances() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "feature_temp") %>%
  # decode_rfriendly_rows(passed_column = "feature_temp") %>%
  # dplyr::select(-feature_temp) %>% 
  # dplyr::rename(feature = fullnames) %>% 
  pivot_longer(!feature_temp, names_to = "donor_id", values_to = "Abundance") %>% 
  group_col_from_ids2()

KOs %>% 
  filter(gsub())

plt <- plot_foi(dat.KOs.slim, foi = "feat_K04334.gc..gs.major.gs.curlin.gs.subunit")
  



dat.species

# # SCRAP
# remove_dats()
# load_all_cohorts()
# 
# # tst <- abundances(dat.EGGNOGs)
# # # head(tst)
# # colSums(tst)
# tst <- 
#   read_tsv(
#     "files/TBC_biobakery_output_slim/humann/merged/pfam-cpm-named.tsv", 
#     col_names = T)
# # head(ec.abund)
# colSums(tst[-1])


# csgA_hits_temp <- as.data.frame(lapply(csgA_hits, as.character), stringsAsFactors=FALSE)
# for (row in 1:nrow(csgA_hits_temp)) {
#   if (row %% 2 == 0) {
#     print(row)
#     csgA_hits_temp[row-1, 6:25] <- csgA_hits_temp[row, 3:22]
#   }
# }
# 
# # Select only odd rows & create cluster ID column by trimming Target.Accession 
# csgA_hmm_hits <- 
#   csgA_hits_temp %>% 
#   filter(row_number() %% 2 == 1) %>% 
#   mutate(clusterID = substr(Target.Accession, 1, nchar(Target.Accession)-6)) %>% 
#   dplyr::relocate(clusterID, .after = Target.Accession)

