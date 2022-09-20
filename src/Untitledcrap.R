# scrap

# rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")
load("files/low_quality_samples.RData")


#---------------------------------------------------------- 
#-                 Metadata Prep - Bonn   
#---------------------------------------------------------- 
metadata_Bonn <- 
  read.csv(file = "files/metadata_phyloseq_Bonn.csv", header = T, stringsAsFactors = F) %>%
  janitor::clean_names() %>% 
  dplyr::rename("PD" = "pd") %>% 
  dplyr::select(run:cohort, library_layout) 

BonnPE_reads <-
  read_tsv("files/biobakery_output_BONN_PE_slim/humann/counts/humann_read_and_species_count_table.tsv",
           col_names = T)
BonnSE_reads <-
  read_tsv("files/biobakery_output_BONN_SE_slim/humann/counts/humann_read_and_species_count_table.tsv",
           col_names = T)
Bonn_reads <- 
  bind_rows(BonnPE_reads, BonnSE_reads) %>% 
  janitor::clean_names() %>% 
  dplyr::rename("run" = "number_samples") %>% 
  dplyr::mutate(run = str_replace(run, "_1", ""))

metadata_Bonn <- full_join(metadata_Bonn, Bonn_reads, by = "run")
rownames(metadata_Bonn) <- metadata_Bonn$run

metadata_Bonn_cleanreads <- metadata_Bonn %>% 
  dplyr::select(run, total_reads)

metadata_Bonn.slim <- metadata_Bonn %>% 
  distinct(id, donor_id, description, donor_group, PD, paired, cohort)
rownames(metadata_Bonn.slim) <- metadata_Bonn.slim$donor_id

#---------------------------------------------------------- 
#-                Bonn -  Taxonomy 
#---------------------------------------------------------- 

met.table.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T)
met.table.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T)
met.table <- 
  full_join(met.table.PE, met.table.SE, by = "# taxonomy") %>% 
  replace(is.na(.), 0 )

# Select only species rows from 
bugs.species <- 
  met.table %>% 
  dplyr::rename("taxonomy" = `# taxonomy`) %>% 
  filter(grepl("s__", taxonomy)) %>% 
  filter(!grepl("t__", taxonomy)) %>% 
  mutate(taxonomy = gsub("s__", "", taxonomy)) %>% 
  column_to_rownames(var = "taxonomy") %>% 
  clean.cols.tax() %>%
  trim_cols("Bonn") %>% 
  dplyr::mutate_if(is.numeric, ~ hund2relab(.))

dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species.pseudo, metadat = metadata_Bonn) %>% 
  merge_samples("donor_id", fun = sum) %>% 
  bonn_metadata()

# Re-normalize 
species.abund <- dat.species %>% abundances(transform = "compositional")
my_species_table <- otu_table(species.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% 
  dplyr::select(-c(run, alias, library_layout)) %>% sample_data()
dat.species <- phyloseq(my_species_table, my_sample_data, tax_table(dat.species))


#-------- Species Level Object --------
save(dat.species, file = "files/Phyloseq_Bonn/Species_PhyloseqObj.RData")
#-------- Genus Level Object --------
dat.genus = tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[,6]
save(dat.genus, file = "files/Phyloseq_Bonn/Genus_PhyloseqObj.RData")
#-------- Family Level Object --------
dat.family = tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[,5]
save(dat.family, file = "files/Phyloseq_Bonn/Family_PhyloseqObj.RData")
#-------- Order Level Object --------
dat.order = tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[,4]
save(dat.order, file = "files/Phyloseq_Bonn/Order_PhyloseqObj.RData")
#-------- Class Level Object --------
dat.class = tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[,3]
save(dat.class, file = "files/Phyloseq_Bonn/Class_PhyloseqObj.RData")
#-------- Phylum Level Object --------
dat.phylum = tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[,2]
save(dat.phylum, file = "files/Phyloseq_Bonn/Phylum_PhyloseqObj.RData")
#-------- Kingdom Level Object --------
dat.kingdom = tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[,1]
save(dat.kingdom, file = "files/Phyloseq_Bonn/Kingdom_PhyloseqObj.RData")


#---------------------------------------------------------- 
#                   TBC -  Pathways 
#---------------------------------------------------------- 

path.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/pathabundance_relab.tsv",
  col_names = T)
path.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/pathabundance_relab.tsv",
  col_names = T)
path.abund <- 
  full_join(path.abund.PE, path.abund.SE, by = "# Pathway") %>% 
  replace(is.na(.), 0 ) %>% 
  filter(!grepl("UNMAPPED", `# Pathway`)) %>% 
  filter(!grepl("UNINTEGRATED", `# Pathway`)) %>% 
  filter(!grepl("g__", `# Pathway`)) %>% 
  filter(!grepl("unclassified", `# Pathway`)) %>%
  column_to_rownames(var = "# Pathway") %>% 
  colSums()

path.abund <- 
  path.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund() %>% 
  filter(!grepl("UNMAPPED", `# Pathway`)) %>% 
  filter(!grepl("UNINTEGRATED", `# Pathway`))
path.abund.slim <- path.abund %>% 
  filter(!grepl("g__", `# Pathway`)) %>% 
  filter(!grepl("unclassified", `# Pathway`)) %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
path.abund <- path.abund %>% 
  make_rfriendly_rows(passed_column = "# Pathway") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
dat.path
save(dat.path, file = "files/Phyloseq_Bonn/Pathways_PhyloseqObj.RData")
# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data)
dat.path.slim
save(dat.path.slim, file = "files/Phyloseq_Bonn/Pathways.slim_PhyloseqObj.RData")


#---------------------------------------------------------- 
#                   TBC -  Enzymes 
#---------------------------------------------------------- 

ec.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/ecs_relab.tsv", 
    col_names = T)
ec.abund <- 
  ec.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_RPK()
ec.abund.slim <- 
  ec.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
ec.abund <- 
  ec.abund %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec
save(dat.ec, file = "files/Phyloseq_Bonn/Enzymes_PhyloseqObj.RData")
# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec.slim
# Save Phyloseq obj as .RData file 
save(dat.ec.slim, file = "files/Phyloseq_Bonn/Enzymes.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Kegg Orthology
#---------------------------------------------------------- 

KOs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/ko-cpm-named.tsv", 
    col_names = T)
KOs.abund <- 
  KOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
KOs.abund.slim <- 
  KOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>% 
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
KOs.abund <- 
  KOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All KOs 
my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
dat.KOs
save(dat.KOs, file = "files/Phyloseq_Bonn/KOs_PhyloseqObj.RData")
# Slim KOs - no stratification
my_KOs.ab_table.slim <- otu_table(KOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs.slim <- phyloseq(my_KOs.ab_table.slim, my_sample_data)
dat.KOs.slim
# Save Phyloseq obj as .RData file 
save(dat.KOs.slim, file = "files/Phyloseq_Bonn/KOs.slim_PhyloseqObj.RData")


#---------------------------------------------------------- 
#                   TBC -  Gene Ontology
#---------------------------------------------------------- 

GOs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/go-cpm-named.tsv", 
    col_names = T)
GOs.abund <- 
  GOs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>% 
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
GOs.abund.slim <- GOs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
GOs.abund <- 
  GOs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
# All GOs 
my_GOs.ab_table <- otu_table(GOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs <- phyloseq(my_GOs.ab_table, my_sample_data)
dat.GOs
save(dat.GOs, file = "files/Phyloseq_Bonn/GOs_PhyloseqObj.RData")
# Slim GOs - no stratification
my_GOs.ab_table.slim <- otu_table(GOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs.slim <- phyloseq(my_GOs.ab_table.slim, my_sample_data)
dat.GOs.slim
# Save Phyloseq obj as .RData file 
save(dat.GOs.slim, file = "files/Phyloseq_Bonn/GOs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Pfam
#---------------------------------------------------------- 

PFAMs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/pfam-cpm-named.tsv", 
    col_names = T)
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
PFAMs.abund.slim <- PFAMs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
PFAMs.abund <- 
  PFAMs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap)) %>% 
  as.data.frame.matrix()
# All PFAMs 
my_PFAMs.ab_table <- otu_table(PFAMs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs <- phyloseq(my_PFAMs.ab_table, my_sample_data)
dat.PFAMs
save(dat.PFAMs, file = "files/Phyloseq_Bonn/PFAMs_PhyloseqObj.RData")
# Slim PFAMs - no stratification
my_PFAMs.ab_table.slim <- otu_table(PFAMs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs.slim <- phyloseq(my_PFAMs.ab_table.slim, my_sample_data)
dat.PFAMs.slim
# Save Phyloseq obj as .RData file 
save(dat.PFAMs.slim, file = "files/Phyloseq_Bonn/PFAMs.slim_PhyloseqObj.RData")

#---------------------------------------------------------- 
#                   TBC -  Eggnog
#---------------------------------------------------------- 

EGGNOGs.abund <- 
  read_tsv(
    "files/TBC_biobakery_output_slim/humann/merged/eggnog-cpm.tsv", 
    col_names = T)
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::select(-contains(negative_controls)) %>%
  clean.cols.abund_CPM() %>% 
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>% 
  filter(!grepl("UNGROUPED", `# Gene Family`))
EGGNOGs.abund.slim <- EGGNOGs.abund %>%  
  filter(!grepl("g__", `# Gene Family`)) %>% 
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
EGGNOGs.abund <- 
  EGGNOGs.abund %>% 
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>% 
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("TBC") %>%
  dplyr::rename(all_of(TBC_keymap))
# All EGGNOGs 
my_EGGNOGs.ab_table <- otu_table(EGGNOGs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs <- phyloseq(my_EGGNOGs.ab_table, my_sample_data)
dat.EGGNOGs
save(dat.EGGNOGs, file = "files/Phyloseq_Bonn/EGGNOGs_PhyloseqObj.RData")
# Slim EGGNOGs - no stratification
my_EGGNOGs.ab_table.slim <- otu_table(EGGNOGs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs.ab_table.slim, my_sample_data)
dat.EGGNOGs.slim
# Save Phyloseq obj as .RData file 
save(dat.EGGNOGs.slim, file = "files/Phyloseq_Bonn/EGGNOGs.slim_PhyloseqObj.RData")








