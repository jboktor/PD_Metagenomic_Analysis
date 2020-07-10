## Phyloseq Obj Creation


library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot);
library(EnvStats)



##### New Method
# reads in data from merged metaphlan2_taxonomic_profiles output - contains all species

# Prep Metadata
metadata <- read_xlsx('files/metadata_phyloseq.xlsx', sheet = 'Sheet1')
trans1 <- gsub("_taxonomic_profile.tsv", "", metadata$metaphlan2_ID)
rownames(metadata) <- gsub("_", ".", trans1)
metadata <- as.data.frame(metadata)
metadata[is.na(metadata)] <- "not provided"

# Get NEW (trimmed) metaphlanToPhyloseq.R function 
# source("https://raw.githubusercontent.com/waldronlab/presentations/master/Waldron_2016-06-07_EPIC/metaphlanToPhyloseq.R")
source("src/metaphlanToPhyloseq_Waldron.R")
# Load merged metaphlan2 taxonomic data
met.table <-read.csv(file = "files/metaphlan2_taxonomic_profiles.csv", row.names = 1,  header= TRUE)
# Select only species rows from 
met.table$temp <- rownames(met.table)
t <- filter(met.table, grepl("s__", met.table$temp))
t2 <- filter(t, !grepl("t__", t$temp))
rownames(t2) <- t2$temp
t2 <- t2[1:118]

# Run metaphlanToPhyloseq.R function
dat <- metaphlanToPhyloseq_Waldron(tax = t2, metadat = metadata)
print(dat)

######### Species Level Object ######### 
# Save Phyloseq obj as .RData file 
save(dat, file = "files/Species_PhyloseqObj.RData")

######### Genus Level Object ######### 
dat.genus = tax_glom(dat, taxrank = "Genus", NArm = F)
taxa_names(dat.genus) <- tax_table(dat.genus)[,6]
# Save Phylum level Phyloseq obj as .RData file 
save(dat.genus, file = "files/Genus_PhyloseqObj.RData")

######### Phylum Level Object ######### 
dat.phylum = tax_glom(dat, taxrank = "Phylum", NArm = F)
taxa_names(dat.phylum) <- tax_table(dat.phylum)[,2]
# Save Phylum level Phyloseq obj as .RData file 
save(dat.phylum, file = "files/Phylum_PhyloseqObj.RData")





############# Phyloseq Obj for Pathways -  ############# 

path.abund <- read_tsv("files/pathabundance_relab.tsv", col_names = T) 
colnames(path.abund) <- gsub("_Abundance", "", colnames(path.abund))
colnames(path.abund) <- gsub("_", ".", colnames(path.abund))
path.abund.slim <- path.abund %>%  filter(!grepl("g__", `# Pathway`))  %>%  filter(!grepl("unclassified", `# Pathway`))

path.abund <- path.abund %>% column_to_rownames(var = "# Pathway") %>% as.data.frame.matrix()
path.abund.slim <- path.abund.slim %>% column_to_rownames(var = "# Pathway") %>% as.data.frame.matrix()

# All Pathway Data
my_pathab_table <- otu_table(path.abund, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()
dat.path <- phyloseq(my_pathab_table, my_sample_data)
print(dat.path)
# Save Phyloseq obj as .RData file 
save(dat.path, file = "files/Pathways_PhyloseqObj.RData")


# Slim Pathway Data
my_pathab_table <- otu_table(path.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()
dat.path.slim <- phyloseq(my_pathab_table, my_sample_data)
print(dat.path.slim)
# Save Phyloseq obj as .RData file 
save(dat.path.slim, file = "files/Pathways.slim_PhyloseqObj.RData")





############# Phyloseq Obj for Enzymes -#############

ec.abund <- read_tsv("files/ecs_relab.tsv", col_names = T)
colnames(ec.abund) <- gsub("_Abundance-RPKs", "", colnames(ec.abund))
colnames(ec.abund) <- gsub("_", ".", colnames(ec.abund))
ec.abund.slim <- ec.abund %>%  filter(!grepl("g__", `# Gene Family`))  %>%  filter(!grepl("unclassified", `# Gene Family`))

ec.abund <- ec.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()
ec.abund.slim <- ec.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()

# All Enzyme Data
my_EC.ab_table <- otu_table(ec.abund, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()
dat.ec <- phyloseq(my_EC.ab_table, my_sample_data)
print(dat.ec)
# Save Phyloseq obj as .RData file 
save(dat.ec, file = "files/Enzymes_PhyloseqObj.RData")

# Slim Enzyme Data - no stratification
my_EC.ab_table <- otu_table(ec.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()
dat.ec.slim <- phyloseq(my_EC.ab_table, my_sample_data)
print(dat.ec.slim)
# Save Phyloseq obj as .RData file 
save(dat.ec.slim, file = "files/Enzymes.slim_PhyloseqObj.RData")



############# Phyloseq Obj for KOs - Humann2 regroupped #############
### NOTE: gene abundance (relab) was multiplied by 1,000,000 to avoid rounding bug in humann2_regroup function
# Setting up groups - All genes (KO and ungrouped), all genes with no stratification, just KOs, and just KOs with no stratification

KOs.abund.all <- read_tsv("files/genefamilies_relab_rescaled_KO.tsv", col_names = T)            #All genes (KO and ungrouped)
KOs.abund.all.slim <- filter(KOs.abund.all, !grepl("g__", `# Gene Family`))
KOs.abund.all.slim <- filter(KOs.abund.all.slim, !grepl("unclassified", `# Gene Family`))  #All genes (KO and ungrouped) with no stratification
KOs.abund <- filter(KOs.abund.all, !grepl("UNGROUPED", `# Gene Family`))                   #Only KOs 
KOs.abund.slim <- filter(KOs.abund, !grepl("g__", `# Gene Family`))      
KOs.abund.slim <- filter(KOs.abund.slim, !grepl("unclassified", `# Gene Family`))          #Only KOs with no stratification


#All genes (KO and ungrouped)
colnames(KOs.abund.all) <- gsub("_Abundance-RPKs", "", colnames(KOs.abund.all))
colnames(KOs.abund.all) <- gsub("_", ".", colnames(KOs.abund.all))
KOs.abund.all <- KOs.abund.all %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()

my_KOs.ab_table <- otu_table(KOs.abund.all, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()

dat.KOs.all <- phyloseq(my_KOs.ab_table, my_sample_data)
print(dat.KOs.all)
save(dat.KOs.all, file = "files/KOs.all_PhyloseqObj.RData")


#All genes (KO and ungrouped) with no stratification
colnames(KOs.abund.all.slim) <- gsub("_Abundance-RPKs", "", colnames(KOs.abund.all.slim))
colnames(KOs.abund.all.slim) <- gsub("_", ".", colnames(KOs.abund.all.slim))
KOs.abund.all.slim <- KOs.abund.all.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()

my_KOs.ab_table <- otu_table(KOs.abund.all.slim, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()

dat.KOs.all.slim <- phyloseq(my_KOs.ab_table, my_sample_data)
print(dat.KOs.all.slim)
save(dat.KOs.all.slim, file = "files/KOs.all.slim_PhyloseqObj.RData")


#Only KOs 
colnames(KOs.abund) <- gsub("_Abundance-RPKs", "", colnames(KOs.abund))
colnames(KOs.abund) <- gsub("_", ".", colnames(KOs.abund))
KOs.abund <- KOs.abund %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()

my_KOs.ab_table <- otu_table(KOs.abund, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()

dat.KOs <- phyloseq(my_KOs.ab_table, my_sample_data)
print(dat.KOs)
save(dat.KOs, file = "files/KOs_PhyloseqObj.RData")


#Only KOs with no stratification
colnames(KOs.abund.slim) <- gsub("_Abundance-RPKs", "", colnames(KOs.abund.slim))
colnames(KOs.abund.slim) <- gsub("_", ".", colnames(KOs.abund.slim))
KOs.abund.slim <- KOs.abund.slim %>% column_to_rownames(var = "# Gene Family") %>% as.data.frame.matrix()

my_KOs.ab_table <- otu_table(KOs.abund.slim, taxa_are_rows=T)
my_sample_data <- meta(dat) %>% sample_data()

dat.KOs.slim <- phyloseq(my_KOs.ab_table, my_sample_data)
print(dat.KOs.slim)
save(dat.KOs.slim, file = "files/KOs.slim_PhyloseqObj.RData")







##### Old Method 
# Reads data directly from Metaphlan output folder - creates tree but filters features not matching in Rxml database (unclassified species)

#  Run local version of :  metaphlanToPhyloseq.R

# metadata <- read_xlsx('metadata_phyloseq.xlsx', sheet = 'Sheet1')
# rownames(metadata) <- metadata$metaphlan2_ID
# metadata <- as.data.frame(metadata)
# metadata[is.na(metadata)] <- "not provided"
# dat <- metaphlanToPhyloseq("/Users/joeboktor/Documents/PD_Metagenomic_Analysis/FA_neph/metaphlan2/main", metadat = metadata)



