#################################### Fig 1 Table stats #################################### 

##### Alpha Diversity Boxplots Script

source("src/load_packages.R")
source("src/Metadata_prep_funcs.R")
load("files/Species_PhyloseqObj.RData")

#### Summary Stats ####
## Prep Metadata
# process_meta(dat)

# Create Function for the Standard Error of the Mean 
sem <- function(x) sd(x)/sqrt(length(x))

env <- meta(dat) 
env[env == "not provided" ] <- NA

## host-age
env$host_age <- as.numeric(env$host_age)
env %>% 
  group_by(donor_group) %>% 
  dplyr::summarise(count = n(), 
                   mean = mean(host_age, na.rm = TRUE),
                   sd = sd(host_age, na.rm = TRUE),
                   se = sem(na.omit(host_age)))

## BMI
env$host_body_mass_index <- as.numeric(env$host_body_mass_index)
env %>% 
  group_by(donor_group) %>% 
  dplyr::summarise(count = n(), 
                   mean = mean(host_body_mass_index, na.rm = TRUE),
                   sd = sd(host_body_mass_index, na.rm = TRUE),
                   se = sem(na.omit(host_body_mass_index)))
# BSS
env$bristol_stool_scale <- as.numeric(env$bristol_stool_scale)
env %>% 
  group_by(donor_group) %>% 
  dplyr::summarise(count = n(), 
                   mean = mean(bristol_stool_scale, na.rm = TRUE),
                   sd = sd(bristol_stool_scale, na.rm = TRUE),
                   se = sem(na.omit(bristol_stool_scale)))





pd_dat = subset_samples(dat, donor_group == "PD")
hc_dat = subset_samples(dat, donor_group == "HC")
pc_dat = subset_samples(dat, donor_group == "PC")

# pull metadata
pd_dat_meta <- meta(pd_dat)
hc_dat_meta <- meta(hc_dat)
pc_dat_meta <- meta(pc_dat)
pd_dat_meta[pd_dat_meta=="not provided"] <-  NA
hc_dat_meta[hc_dat_meta=="not provided"] <-  NA
pc_dat_meta[pc_dat_meta=="not provided"] <-  NA

# Sex Count
count(na.omit(pd_dat_meta$sex))
count(na.omit(hc_dat_meta$sex))
count(na.omit(pc_dat_meta$sex))

## READS 
func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
reads <- dplyr::select(func_reads, c("# samples","total reads")) %>% 
  dplyr::rename( "id" = "# samples", "clean_total_reads" = "total reads")
reads$id <- gsub("_", ".", reads$id)
df.reads <- group_col_from_ids(reads, reads$id)

df.reads %>% 
  group_by(group) %>% 
  dplyr::summarise(count = n(), 
                   mean = mean(clean_total_reads, na.rm = TRUE),
                   sd = sd(clean_total_reads, na.rm = TRUE),
                   se = sem(na.omit(clean_total_reads)))

# to show estimate by millions
mili <- 1000000 
df.reads %>% 
  group_by(group) %>% 
  dplyr::summarise(count = n(), 
                   mean = mean(clean_total_reads/mili, na.rm = TRUE),
                   sd = sd(clean_total_reads/mili, na.rm = TRUE),
                   se = sem(na.omit(clean_total_reads/mili)))





# ## READS 
# QC_data <- read.table(file = "qc_counts_pairs_table.tsv")
# redz <- mutate(QC_data, group = if_else(grepl("HC", Helper), "HC", 
#                                         if_else(grepl("PC", Helper), "PC",
#                                                 "PD")))
# redz_pd <- filter(redz, group=="PD")
# redz_pc<- filter(redz, group=="PC")
# redz_hc<- filter(redz, group=="HC")
# 
# # Mean Reads per group
# mean(redz_pd$Raw)/1000000
# mean(redz_pc$Raw)/1000000
# mean(redz_hc$Raw)/1000000
# 
# # Mean Reads per group
# sd(redz_pd$Raw)/1000000
# sd(redz_pc$Raw)/1000000
# sd(redz_hc$Raw)/1000000

