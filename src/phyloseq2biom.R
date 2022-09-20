# Write out to biom

# Run the below system commands once to download phyloseq to biom function (requires wget)
# system(command = "wget https://raw.githubusercontent.com/itsmisterbrown/microfiltR/master/microfiltR_source_code.R")
# system(command = "mv microfiltR_source_code.R src/microfiltR_source_code.R")

source("src/load_packages.R")
source("src/microfiltR_source_code.R")
load("files/low_quality_samples.RData")

load_data("Merged")
dat.species2biom <- dat.species %>% 
  subset_samples(donor_id %ni% low_qc[[1]])

sample_data(dat.species2biom)$latitude <- "not provided"
sample_data(dat.species2biom)$longitude <- "not provided"
sample_data(dat.species2biom)$datetime <- "not provided"

countdata <- as.data.frame(abundances(dat.species2biom))
metadata4biom <- as.data.frame(meta(dat.species2biom)) %>% 
  dplyr::select(donor_group, PD, latitude, longitude, datetime)

# countdata2 <- as.matrix(phyloseq::otu_table(dat.species2biom))
# write.table(countdata, file = "files/merged_abund_data.txt", append = FALSE, sep = " ", dec = ".",
#             row.names = TRUE, col.names = TRUE)
write.csv(countdata, file = paste0('files/merged_abund_data.csv'))
write.csv(metadata4biom, file = paste0('files/metadata4biom.csv'))

# countdata %>% 
#   rownames_to_column(var = "donor_id") %>% 
#   write_tsv(file = paste0('files/merged_abund_data2.tsv'))
metadata4biom %>% 
  rownames_to_column(var = "sample") %>% 
  write_tsv(file = paste0('files/metadata4biom2.tsv'))



# write.dataset.biom(
#   ps = dat.species2biom,
#   filePATH = "files/",
#   filePREFIX = "merged_species",
#   writeFASTA = F
# )

# devtools::install_github("jladau/Nestedness")
# githubinstall("Nestedness")