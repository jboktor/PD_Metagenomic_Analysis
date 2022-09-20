# Create phyloseq objects

rm(list = ls())
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/metaphlanToPhyloseq_Waldron.R")
load("files/low_quality_samples.RData")

#-------------------------------------------------------------------------------
#######                         shanghai - Cohort                        #######
#-------------------------------------------------------------------------------

metadata_shanghai <-
  read.csv(file = "files/metadata_phyloseq_shanghai.csv", header = TRUE, stringsAsFactors = F) %>%
  dplyr::rename(donor_id = Sample.Name, run = Run) %>%
  dplyr::select(donor_id, donor_group, description, PD, cohort, paired, run) %>%
  dplyr::mutate(paired = as.character(paired)) %>%
  filter(run %ni% c("SRR10983032", "SRR10983014")) %>%
  as.data.frame()

shanghai_keys <-
  metadata_shanghai %>%
  dplyr::select(run, donor_id) %>%
  mutate(donor_id = as.character(donor_id))

reads.shanghai <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  ) %>%
  janitor::clean_names() %>%
  dplyr::rename(run = number_samples) %>%
  left_join(shanghai_keys, by = "run")

metadata_shanghai <-
  left_join(metadata_shanghai, reads.shanghai, by = "donor_id")

shanghai_keymap <- shanghai_keys$run
names(shanghai_keymap) <- shanghai_keys$donor_id
metadata_shanghai[is.na(metadata_shanghai)] <- "not provided"
rownames(metadata_shanghai) <- metadata_shanghai$donor_id
metadata_shanghai <- as.data.frame(metadata_shanghai)

#----------------------------------------------------------
#-                shanghai -  Taxonomy
#----------------------------------------------------------
met.table <- read_tsv(
  file = "files/biobakery_output_shanghai_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T
)

# Select only species rows from
bugs.species <- met.table %>%
  dplyr::rename("taxonomy" = `# taxonomy`) %>%
  filter(grepl("s__", taxonomy)) %>%
  filter(!grepl("t__", taxonomy)) %>%
  mutate(taxonomy = gsub("s__", "", taxonomy)) %>%
  column_to_rownames(var = "taxonomy") %>%
  clean.cols.tax() %>%
  dplyr::rename(shanghai_keymap)

dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species, metadat = metadata_shanghai
)

tst <- meta(dat.species)

#-------- Species Level Object --------
save(dat.species, file = "files/Phyloseq_SHANGHAI/Species_PhyloseqObj.RData")
#-------- Genus Level Object --------
dat.genus <- tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[, 6]
save(dat.genus, file = "files/Phyloseq_SHANGHAI/Genus_PhyloseqObj.RData")
#-------- Family Level Object --------
dat.family <- tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[, 5]
save(dat.family, file = "files/Phyloseq_SHANGHAI/Family_PhyloseqObj.RData")
#-------- Order Level Object --------
dat.order <- tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[, 4]
save(dat.order, file = "files/Phyloseq_SHANGHAI/Order_PhyloseqObj.RData")
#-------- Class Level Object --------
dat.class <- tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[, 3]
save(dat.class, file = "files/Phyloseq_SHANGHAI/Class_PhyloseqObj.RData")
#-------- Phylum Level Object --------
dat.phylum <- tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[, 2]
save(dat.phylum, file = "files/Phyloseq_SHANGHAI/Phylum_PhyloseqObj.RData")
#-------- Kingdom Level Object --------
dat.kingdom <- tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[, 1]
save(dat.kingdom, file = "files/Phyloseq_SHANGHAI/Kingdom_PhyloseqObj.RData")

#-------------------- --------------------------------------
#                   shanghai -  Pathways
#----------------------------------------------------------

path.abund <-
  read_tsv(file = "files/biobakery_output_SHANGHAI_slim/humann/merged/pathabundance_relab.tsv", col_names = T)
path.abund <-
  path.abund %>%
  clean.cols.abund() %>%
  filter(!grepl("UNMAPPED", `# Pathway`)) %>%
  filter(!grepl("UNINTEGRATED", `# Pathway`))
path.abund.slim <- path.abund %>%
  filter(!grepl("g__", `# Pathway`)) %>%
  filter(!grepl("unclassified", `# Pathway`)) %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>%
  dplyr::rename(shanghai_keymap)
path.abund <- path.abund %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>%
  dplyr::rename(shanghai_keymap)
# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
dat.path
save(dat.path, file = "files/Phyloseq_SHANGHAI/Pathways_PhyloseqObj.RData")
# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data)
dat.path.slim
save(dat.path.slim, file = "files/Phyloseq_SHANGHAI/Pathways.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   shanghai -  Enzymes
#----------------------------------------------------------

ec.abund <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/merged/ecs_relab.tsv",
    col_names = T
  )
ec.abund <-
  ec.abund %>%
  clean.cols.abund_RPK()
ec.abund.slim <-
  ec.abund %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
ec.abund <-
  ec.abund %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec
save(dat.ec, file = "files/Phyloseq_SHANGHAI/Enzymes_PhyloseqObj.RData")
# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data)
dat.ec.slim
# Save Phyloseq obj as .RData file
save(dat.ec.slim, file = "files/Phyloseq_SHANGHAI/Enzymes.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   shanghai -  Kegg Orthology
#----------------------------------------------------------

KOs.abund <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/merged/ko-cpm-named.tsv",
    col_names = T
  )
KOs.abund <-
  KOs.abund %>%
  clean.cols.abund_CPM() %>%
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>%
  filter(!grepl("UNGROUPED", `# Gene Family`))
KOs.abund.slim <-
  KOs.abund %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
KOs.abund <-
  KOs.abund %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
# All KOs
my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
dat.KOs
save(dat.KOs, file = "files/Phyloseq_SHANGHAI/KOs_PhyloseqObj.RData")
# Slim KOs - no stratification
my_KOs.ab_table.slim <- otu_table(KOs.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.KOs.slim <- phyloseq(my_KOs.ab_table.slim, my_sample_data)
dat.KOs.slim
# Save Phyloseq obj as .RData file
save(dat.KOs.slim, file = "files/Phyloseq_SHANGHAI/KOs.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   shanghai -  Gene Ontology
#----------------------------------------------------------

GOs.abund <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/merged/go-cpm-named.tsv",
    col_names = T
  )
GOs.abund <-
  GOs.abund %>%
  clean.cols.abund_CPM() %>%
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>%
  filter(!grepl("UNGROUPED", `# Gene Family`))
GOs.abund.slim <- GOs.abund %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap) %>%
  as.data.frame.matrix()
GOs.abund <-
  GOs.abund %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap) %>%
  as.data.frame.matrix()
# All GOs
my_GOs.ab_table <- otu_table(GOs.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs <- phyloseq(my_GOs.ab_table, my_sample_data)
dat.GOs
save(dat.GOs, file = "files/Phyloseq_SHANGHAI/GOs_PhyloseqObj.RData")
# Slim GOs - no stratification
my_GOs.ab_table.slim <- otu_table(GOs.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.GOs.slim <- phyloseq(my_GOs.ab_table.slim, my_sample_data)
dat.GOs.slim
# Save Phyloseq obj as .RData file
save(dat.GOs.slim, file = "files/Phyloseq_SHANGHAI/GOs.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   shanghai -  Pfam
#----------------------------------------------------------

PFAMs.abund <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/merged/pfam-cpm-named.tsv",
    col_names = T
  )
PFAMs.abund <-
  PFAMs.abund %>%
  clean.cols.abund_CPM() %>%
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>%
  filter(!grepl("UNGROUPED", `# Gene Family`))
PFAMs.abund.slim <- PFAMs.abund %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap) %>%
  as.data.frame.matrix()
PFAMs.abund <-
  PFAMs.abund %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap) %>%
  as.data.frame.matrix()
# All PFAMs
my_PFAMs.ab_table <- otu_table(PFAMs.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs <- phyloseq(my_PFAMs.ab_table, my_sample_data)
dat.PFAMs
save(dat.PFAMs, file = "files/Phyloseq_SHANGHAI/PFAMs_PhyloseqObj.RData")
# Slim PFAMs - no stratification
my_PFAMs.ab_table.slim <- otu_table(PFAMs.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.PFAMs.slim <- phyloseq(my_PFAMs.ab_table.slim, my_sample_data)
dat.PFAMs.slim
# Save Phyloseq obj as .RData file
save(dat.PFAMs.slim, file = "files/Phyloseq_SHANGHAI/PFAMs.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   shanghai -  Eggnog
#----------------------------------------------------------

EGGNOGs.abund <-
  read_tsv(
    "files/biobakery_output_SHANGHAI_slim/humann/merged/eggnog-cpm.tsv",
    col_names = T
  )
EGGNOGs.abund <-
  EGGNOGs.abund %>%
  clean.cols.abund_CPM() %>%
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>%
  filter(!grepl("UNGROUPED", `# Gene Family`))
EGGNOGs.abund.slim <- EGGNOGs.abund %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
EGGNOGs.abund <-
  EGGNOGs.abund %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  dplyr::rename(shanghai_keymap)
# All EGGNOGs
my_EGGNOGs.ab_table <- otu_table(EGGNOGs.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs <- phyloseq(my_EGGNOGs.ab_table, my_sample_data)
dat.EGGNOGs
save(dat.EGGNOGs, file = "files/Phyloseq_SHANGHAI/EGGNOGs_PhyloseqObj.RData")
# Slim EGGNOGs - no stratification
my_EGGNOGs.ab_table.slim <- otu_table(EGGNOGs.abund.slim, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>% sample_data()
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs.ab_table.slim, my_sample_data)
dat.EGGNOGs.slim
# Save Phyloseq obj as .RData file
save(dat.EGGNOGs.slim, file = "files/Phyloseq_SHANGHAI/EGGNOGs.slim_PhyloseqObj.RData")


# temporary rename
dat.kingdom.shanghai <- dat.kingdom
dat.phylum.shanghai <- dat.phylum
dat.class.shanghai <- dat.class
dat.order.shanghai <- dat.order
dat.family.shanghai <- dat.family
dat.genus.shanghai <- dat.genus
dat.species.shanghai <- dat.species

dat.path.shanghai <- dat.path
dat.path.slim.shanghai <- dat.path.slim
dat.ec.shanghai <- dat.ec
dat.ec.slim.shanghai <- dat.ec.slim
dat.KOs.shanghai <- dat.KOs
dat.KOs.slim.shanghai <- dat.KOs.slim
dat.GOs.shanghai <- dat.GOs
dat.GOs.slim.shanghai <- dat.GOs.slim
dat.PFAMs.shanghai <- dat.PFAMs
dat.PFAMs.slim.shanghai <- dat.PFAMs.slim
dat.EGGNOGs.shanghai <- dat.EGGNOGs
dat.EGGNOGs.slim.shanghai <- dat.EGGNOGs.slim


#-------------------------------------------------------------------------------
#######                         Bonn - Cohort                        #######
#-------------------------------------------------------------------------------

metadata_Bonn <-
  read.csv(file = "files/metadata_phyloseq_Bonn.csv", header = T, stringsAsFactors = F) %>%
  janitor::clean_names() %>%
  dplyr::rename("PD" = "pd") %>%
  dplyr::select(run:cohort, library_layout)

BonnPE_reads <-
  read_tsv("files/biobakery_output_BONN_PE_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
BonnSE_reads <-
  read_tsv("files/biobakery_output_BONN_SE_slim/humann/counts/humann_read_and_species_count_table.tsv",
    col_names = T
  )
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
  col_names = T
)
met.table.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/metaphlan/merged/metaphlan_taxonomic_profiles.tsv",
  col_names = T
)
met.table <-
  full_join(met.table.PE, met.table.SE, by = "# taxonomy") %>%
  replace(is.na(.), 0)

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
bugs.species.pseudo <-
  pseudoCounts_bonn(df = bugs.species, reads = metadata_Bonn_cleanreads)

dat.species <- metaphlanToPhyloseq_Waldron(
  tax = bugs.species.pseudo, metadat = metadata_Bonn
) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re-normalize
species.abund <- dat.species %>% abundances(transform = "compositional") * 100
my_species_table <- otu_table(species.abund, taxa_are_rows = T)
my_sample_data <- meta(dat.species) %>%
  dplyr::select(-c(run, alias, library_layout)) %>%
  sample_data()
dat.species <- phyloseq(my_species_table, my_sample_data, tax_table(dat.species))


#-------- Species Level Object --------
save(dat.species, file = "files/Phyloseq_Bonn/Species_PhyloseqObj.RData")
#-------- Genus Level Object --------
dat.genus <- tax_glom(dat.species, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[, 6]
save(dat.genus, file = "files/Phyloseq_Bonn/Genus_PhyloseqObj.RData")
#-------- Family Level Object --------
dat.family <- tax_glom(dat.species, taxrank = "Family", NArm = F)
taxa_names(dat.family) <- tax_table(dat.family)[, 5]
save(dat.family, file = "files/Phyloseq_Bonn/Family_PhyloseqObj.RData")
#-------- Order Level Object --------
dat.order <- tax_glom(dat.species, taxrank = "Order", NArm = F)
taxa_names(dat.order) <- tax_table(dat.order)[, 4]
save(dat.order, file = "files/Phyloseq_Bonn/Order_PhyloseqObj.RData")
#-------- Class Level Object --------
dat.class <- tax_glom(dat.species, taxrank = "Class", NArm = F)
taxa_names(dat.class) <- tax_table(dat.class)[, 3]
save(dat.class, file = "files/Phyloseq_Bonn/Class_PhyloseqObj.RData")
#-------- Phylum Level Object --------
dat.phylum <- tax_glom(dat.species, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[, 2]
save(dat.phylum, file = "files/Phyloseq_Bonn/Phylum_PhyloseqObj.RData")
#-------- Kingdom Level Object --------
dat.kingdom <- tax_glom(dat.species, taxrank = "Kingdom", NArm = F)
taxa_names(dat.kingdom) <- tax_table(dat.kingdom)[, 1]
save(dat.kingdom, file = "files/Phyloseq_Bonn/Kingdom_PhyloseqObj.RData")


#----------------------------------------------------------
#                   Bonn -  Pathways
#----------------------------------------------------------

path.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/pathabundance_relab.tsv",
  col_names = T
)
path.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/pathabundance_relab.tsv",
  col_names = T
)
path.abund.slim <-
  full_join(path.abund.PE, path.abund.SE, by = "# Pathway") %>%
  replace(is.na(.), 0) %>%
  filter(!grepl("UNMAPPED", `# Pathway`)) %>%
  filter(!grepl("UNINTEGRATED", `# Pathway`)) %>%
  filter(!grepl("g__", `# Pathway`)) %>%
  filter(!grepl("unclassified", `# Pathway`)) %>%
  clean.cols.abund() %>%
  make_rfriendly_rows(passed_column = "# Pathway") %>%
  trim_cols("Bonn")

path.abund.pseudo.slim <-
  pseudoCounts_bonn(df = path.abund.slim, reads = metadata_Bonn_cleanreads)

# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.pseudo.slim, taxa_are_rows = T)
my_sample_data.all <- metadata_Bonn %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re group and re-normalize
path.abund <- dat.path.slim %>% abundances(transform = "compositional")
my_path_table <- otu_table(path.abund, taxa_are_rows = T)
dat.path.slim <- phyloseq(my_path_table, my_sample_data)

save(dat.path.slim, file = "files/Phyloseq_Bonn/Pathways_PhyloseqObj.RData")


#----------------------------------------------------------
#                   Bonn -  Enzymes
#----------------------------------------------------------

ec.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/ecs_relab.tsv",
  col_names = T
)
ec.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/ecs_relab.tsv",
  col_names = T
)
ec.abund.slim <-
  full_join(ec.abund.PE, ec.abund.SE, by = "# Gene Family") %>%
  replace(is.na(.), 0) %>%
  filter(!grepl("UNMAPPED", `# Gene Family`)) %>%
  filter(!grepl("UNINTEGRATED", `# Gene Family`)) %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  clean.cols.abund_RPK() %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("Bonn")

ec.abund.pseudo.slim <-
  pseudoCounts_bonn(df = ec.abund.slim, reads = metadata_Bonn_cleanreads)

# Slim Pathway Data
my_ec_table <- otu_table(ec.abund.pseudo.slim, taxa_are_rows = T)
dat.ec.slim <- phyloseq(my_ec_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re group and re-normalize
ec.abund <- dat.ec.slim %>% abundances(transform = "compositional")
my_ec_table <- otu_table(ec.abund, taxa_are_rows = T)
dat.ec.slim <- phyloseq(my_ec_table, my_sample_data)
# Save Phyloseq obj as .RData file
save(dat.ec.slim, file = "files/Phyloseq_Bonn/Enzymes.slim_PhyloseqObj.RData")

#----------------------------------------------------------
#                   Bonn -  Kegg Orthology
#----------------------------------------------------------

KOs.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/ko-cpm-named.tsv",
  col_names = T
)
KOs.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/ko-cpm-named.tsv",
  col_names = T
)
KOs.abund.slim <-
  full_join(KOs.abund.PE, KOs.abund.SE, by = "# Gene Family") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  clean.cols.abund_CPM() %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("Bonn")

KOs.abund.pseudo.slim <-
  pseudoCounts_bonn(df = KOs.abund.slim, reads = metadata_Bonn_cleanreads) %>%
  mutate_all(~ tss(.))

# Slim Data
my_KOs_table <- otu_table(KOs.abund.pseudo.slim, taxa_are_rows = T)
dat.KOs.slim <- phyloseq(my_KOs_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Generate Total Sum Scale on non-stratified data
# filter out un-mapped and un-grouped data
KOs.abund <- dat.KOs.slim %>%
  abundances() %>%
  as.data.frame() %>%
  mutate_all(~ tss(.)) %>%
  rownames_to_column() %>%
  filter(!grepl("UNMAPPED", rowname)) %>%
  filter(!grepl("UNGROUPED", rowname)) %>%
  column_to_rownames()
my_KOs_table <- otu_table(KOs.abund, taxa_are_rows = T)
dat.KOs.slim <- phyloseq(my_KOs_table, my_sample_data)
# Save Phyloseq obj as .RData file
save(dat.KOs.slim, file = "files/Phyloseq_Bonn/KOs.slim_PhyloseqObj.RData")


#----------------------------------------------------------
#                   Bonn -  Gene Ontology
#----------------------------------------------------------

GOs.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/GO-cpm-named.tsv",
  col_names = T
)
GOs.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/GO-cpm-named.tsv",
  col_names = T
)
GOs.abund.slim <-
  full_join(GOs.abund.PE, GOs.abund.SE, by = "# Gene Family") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  clean.cols.abund_CPM() %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("Bonn")

GOs.abund.pseudo.slim <-
  pseudoCounts_bonn(df = GOs.abund.slim, reads = metadata_Bonn_cleanreads) %>%
  mutate_all(~ tss(.))

# Slim Data
my_GOs_table <- otu_table(GOs.abund.pseudo.slim, taxa_are_rows = T)
dat.GOs.slim <- phyloseq(my_GOs_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re group and re-normalize
GOs.abund <- dat.GOs.slim %>%
  abundances() %>%
  as.data.frame() %>%
  mutate_all(~ tss(.)) %>%
  rownames_to_column() %>%
  filter(!grepl("UNMAPPED", rowname)) %>%
  filter(!grepl("UNGROUPED", rowname)) %>%
  column_to_rownames()
my_GOs_table <- otu_table(GOs.abund, taxa_are_rows = T)
dat.GOs.slim <- phyloseq(my_GOs_table, my_sample_data)
# Save Phyloseq obj as .RData file
save(dat.GOs.slim, file = "files/Phyloseq_Bonn/GOs.slim_PhyloseqObj.RData")


#----------------------------------------------------------
#                   Bonn -  Pfam
#----------------------------------------------------------

PFAMs.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/pfam-cpm-named.tsv",
  col_names = T
)
PFAMs.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/pfam-cpm-named.tsv",
  col_names = T
)
PFAMs.abund.slim <-
  full_join(PFAMs.abund.PE, PFAMs.abund.SE, by = "# Gene Family") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  clean.cols.abund_CPM() %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("Bonn")

PFAMs.abund.pseudo.slim <-
  pseudoCounts_bonn(df = PFAMs.abund.slim, reads = metadata_Bonn_cleanreads) %>%
  mutate_all(~ tss(.))

# Slim Pathway Data
my_PFAMs_table <- otu_table(PFAMs.abund.pseudo.slim, taxa_are_rows = T)
dat.PFAMs.slim <- phyloseq(my_PFAMs_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re group and re-normalize
PFAMs.abund <- dat.PFAMs.slim %>%
  abundances() %>%
  as.data.frame() %>%
  mutate_all(~ tss(.)) %>%
  rownames_to_column() %>%
  filter(!grepl("UNMAPPED", rowname)) %>%
  filter(!grepl("UNGROUPED", rowname)) %>%
  column_to_rownames()
my_PFAMs_table <- otu_table(PFAMs.abund, taxa_are_rows = T)
dat.PFAMs.slim <- phyloseq(my_PFAMs_table, my_sample_data)
# Save Phyloseq obj as .RData file
save(dat.PFAMs.slim, file = "files/Phyloseq_Bonn/PFAMs.slim_PhyloseqObj.RData")


#----------------------------------------------------------
#                   Bonn -  Eggnog
#----------------------------------------------------------

EGGNOGs.abund.PE <- read_tsv(
  file = "files/biobakery_output_BONN_PE_slim/humann/merged/eggnog-cpm.tsv",
  col_names = T
)
EGGNOGs.abund.SE <- read_tsv(
  file = "files/biobakery_output_BONN_SE_slim/humann/merged/eggnog-cpm.tsv",
  col_names = T
)
EGGNOGs.abund.slim <-
  full_join(EGGNOGs.abund.PE, EGGNOGs.abund.SE, by = "# Gene Family") %>%
  replace(is.na(.), 0) %>%
  dplyr::mutate_if(is.numeric, ~ (. / 1000000)) %>%
  filter(!grepl("g__", `# Gene Family`)) %>%
  filter(!grepl("unclassified", `# Gene Family`)) %>%
  clean.cols.abund_CPM() %>%
  make_rfriendly_rows(passed_column = "# Gene Family") %>%
  trim_cols("Bonn")

EGGNOGs.abund.pseudo.slim <-
  pseudoCounts_bonn(df = EGGNOGs.abund.slim, reads = metadata_Bonn_cleanreads) %>%
  mutate_all(~ tss(.))

# Slim Data
my_EGGNOGs_table <- otu_table(EGGNOGs.abund.pseudo.slim, taxa_are_rows = T)
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs_table, my_sample_data.all) %>%
  merge_samples("donor_id", fun = sum) %>%
  bonn_metadata()

# Re group and re-normalize
EGGNOGs.abund <- dat.EGGNOGs.slim %>%
  abundances() %>%
  as.data.frame() %>%
  mutate_all(~ tss(.)) %>%
  rownames_to_column() %>%
  filter(!grepl("UNMAPPED", rowname)) %>%
  filter(!grepl("UNGROUPED", rowname)) %>%
  column_to_rownames()
my_EGGNOGs_table <- otu_table(EGGNOGs.abund, taxa_are_rows = T)
dat.EGGNOGs.slim <- phyloseq(my_EGGNOGs_table, my_sample_data)
# Save Phyloseq obj as .RData file
save(dat.EGGNOGs.slim, file = "files/Phyloseq_Bonn/EGGNOGs.slim_PhyloseqObj.RData")


# temporary rename
dat.kingdom.bonn <- dat.kingdom
dat.phylum.bonn <- dat.phylum
dat.class.bonn <- dat.class
dat.order.bonn <- dat.order
dat.family.bonn <- dat.family
dat.genus.bonn <- dat.genus
dat.species.bonn <- dat.species

dat.path.slim.bonn <- dat.path.slim
dat.ec.slim.bonn <- dat.ec.slim
dat.KOs.slim.bonn <- dat.KOs.slim
dat.GOs.slim.bonn <- dat.GOs.slim
dat.PFAMs.slim.bonn <- dat.PFAMs.slim
dat.EGGNOGs.slim.bonn <- dat.EGGNOGs.slim

#-------------------------------------------------------------------------------
## --------------------------------   Merge - Cohorts  --------------------------
#-------------------------------------------------------------------------------


phylo_objects <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")

dat.kingdom <- merge_phyloseq(phylo_objects[["Kingdom"]], dat.kingdom.shanghai, dat.kingdom.bonn)
dat.phylum <- merge_phyloseq(phylo_objects[["Phylum"]], dat.phylum.shanghai, dat.phylum.bonn)
dat.class <- merge_phyloseq(phylo_objects[["Class"]], dat.class.shanghai, dat.class.bonn)
dat.order <- merge_phyloseq(phylo_objects[["Order"]], dat.order.shanghai, dat.order.bonn)
dat.family <- merge_phyloseq(phylo_objects[["Family"]], dat.family.shanghai, dat.family.bonn)
dat.genus <- merge_phyloseq(phylo_objects[["Genus"]], dat.genus.shanghai, dat.genus.bonn)
dat.species <- merge_phyloseq(phylo_objects[["Species"]], dat.species.shanghai, dat.species.bonn)

dat.path <- merge_phyloseq(phylo_objects[["Pathways"]], dat.path.shanghai)
dat.path.slim <- merge_phyloseq(phylo_objects[["Pathways.slim"]], dat.path.slim.shanghai, dat.path.slim.bonn)
dat.ec <- merge_phyloseq(phylo_objects[["Enzymes"]], dat.ec.shanghai)
dat.ec.slim <- merge_phyloseq(phylo_objects[["Enzymes.slim"]], dat.ec.slim.shanghai, dat.ec.slim.bonn)
dat.KOs <- merge_phyloseq(phylo_objects[["KOs"]], dat.KOs.shanghai)
dat.KOs.slim <- merge_phyloseq(phylo_objects[["KOs.slim"]], dat.KOs.slim.shanghai, dat.KOs.slim.bonn)
dat.GOs <- merge_phyloseq(phylo_objects[["GOs"]], dat.GOs.shanghai)
dat.GOs.slim <- merge_phyloseq(phylo_objects[["GOs.slim"]], dat.GOs.slim.shanghai, dat.GOs.slim.bonn)
dat.PFAMs <- merge_phyloseq(phylo_objects[["Pfams"]], dat.PFAMs.shanghai)
dat.PFAMs.slim <- merge_phyloseq(phylo_objects[["Pfams.slim"]], dat.PFAMs.slim.shanghai, dat.PFAMs.slim.bonn)
dat.EGGNOGs <- merge_phyloseq(phylo_objects[["eggNOGs"]], dat.EGGNOGs.shanghai)
dat.EGGNOGs.slim <- merge_phyloseq(phylo_objects[["eggNOGs.slim"]], dat.EGGNOGs.slim.shanghai, dat.EGGNOGs.slim.bonn)

save(dat.species, file = "files/Phyloseq_Merged_ML/Species_PhyloseqObj.RData")
save(dat.genus, file = "files/Phyloseq_Merged_ML/Genus_PhyloseqObj.RData")
save(dat.family, file = "files/Phyloseq_Merged_ML/Family_PhyloseqObj.RData")
save(dat.order, file = "files/Phyloseq_Merged_ML/Order_PhyloseqObj.RData")
save(dat.class, file = "files/Phyloseq_Merged_ML/Class_PhyloseqObj.RData")
save(dat.phylum, file = "files/Phyloseq_Merged_ML/Phylum_PhyloseqObj.RData")
save(dat.kingdom, file = "files/Phyloseq_Merged_ML/Kingdom_PhyloseqObj.RData")
save(dat.path, file = "files/Phyloseq_Merged_ML/Pathways_PhyloseqObj.RData")
save(dat.path.slim, file = "files/Phyloseq_Merged_ML/Pathways.slim_PhyloseqObj.RData")
save(dat.ec, file = "files/Phyloseq_Merged_ML/Enzymes_PhyloseqObj.RData")
save(dat.ec.slim, file = "files/Phyloseq_Merged_ML/Enzymes.slim_PhyloseqObj.RData")
save(dat.KOs, file = "files/Phyloseq_Merged_ML/KOs_PhyloseqObj.RData")
save(dat.KOs.slim, file = "files/Phyloseq_Merged_ML/KOs.slim_PhyloseqObj.RData")
save(dat.GOs, file = "files/Phyloseq_Merged_ML/GOs_PhyloseqObj.RData")
save(dat.GOs.slim, file = "files/Phyloseq_Merged_ML/GOs.slim_PhyloseqObj.RData")
save(dat.PFAMs, file = "files/Phyloseq_Merged_ML/PFAMs_PhyloseqObj.RData")
save(dat.PFAMs.slim, file = "files/Phyloseq_Merged_ML/PFAMs.slim_PhyloseqObj.RData")
save(dat.EGGNOGs, file = "files/Phyloseq_Merged_ML/EGGNOGs_PhyloseqObj.RData")
save(dat.EGGNOGs.slim, file = "files/Phyloseq_Merged_ML/EGGNOGs.slim_PhyloseqObj.RData")


phyloseq_objs <- list(
  dat.species,
  dat.genus,
  dat.phylum,
  dat.path.slim,
  dat.ec.slim,
  dat.KOs.slim,
  dat.GOs.slim,
  dat.PFAMs.slim,
  dat.EGGNOGs.slim
)
names(phyloseq_objs) <- c(
  "Species",
  "Genus",
  "Phylum",
  "Pathways.slim",
  "Enzymes.slim",
  "KOs.slim",
  "GOs.slim",
  "Pfams.slim",
  "eggNOGs.slim"
)
saveRDS(phyloseq_objs, file = "files/Phyloseq_Merged_ML_clean.rds")
