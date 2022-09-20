# create ShortBRED / Phyloseq Objects

######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")


#### Virulence/Pathogenicity Analysis

#-------------------------------------------------------------------------
# Build a dataframe from ShortBRED VFDB individual TSVs
#-------------------------------------------------------------------------

l <- list.files(
  path = "files/shortBRED_VFDB/", pattern = NULL, all.files = F,
  full.names = FALSE, recursive = FALSE,
  ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE
)
nwname <- substr(l[1], 19, nchar(l[1]) - 4)
dfseed <- read_tsv(paste0("files/shortBRED_VFDB/", l[1]), col_names = T) %>%
  dplyr::select(c(Family, Count)) %>%
  dplyr::rename(!!nwname := Count)
for (i in l[2:length(l)]) {
  print(i)
  nwname <- print(substr(i, 19, nchar(i) - 4))
  print(nwname)

  dfadd <- read_tsv(paste0("files/shortBRED_VFDB/", i), col_names = T) %>%
    dplyr::select(c(Family, Count)) %>%
    dplyr::rename(!!nwname := Count)

  dfseed <- dplyr::left_join(dfseed, dfadd, by = "Family")
}

shortbreddata.untrimmed <- dfseed
shortbreddata <- shortbreddata.untrimmed[which(rowSums(shortbreddata.untrimmed[-1]) > 0), ]
save(shortbreddata, file = "files/PhyloseqObjects/shortbred.RData")

#-------------------------------------------------------------------------
# Create Phlyoseq Obj for VF data
#-------------------------------------------------------------------------


# shortbreddata
merged <- shortbreddata %>% melt()
merged$MBI_Sample_ID <- substr(merged$variable, 1, 10)

# Create Data Table of VF with 0 sum columns removed
family_merge <- aggregate(merged$value, by = list(MBI_Sample_ID = merged$MBI_Sample_ID, Family = merged$Family), FUN = sum) %>%
  dplyr::rename(value = x)
family_merge$MBI_Sample_ID <- factor(family_merge$MBI_Sample_ID)
family_merge_wide <- family_merge %>%
  spread(Family, value)

m <- read.csv(file = "files/metadata_keys.csv", header = TRUE) %>%
  dplyr::select(c(MBI_Sample_ID, id))
m$id <- gsub("_", ".", m$id)
vf.data <- left_join(family_merge_wide, m, by = "MBI_Sample_ID") %>%
  column_to_rownames(var = "id") %>%
  dplyr::select(-MBI_Sample_ID)

## Save cleaned VF data with trimmed non-counted features as csv
write.csv(vf.data, file = "files/ShortBRED_regrouped.csv")
vf.data.t <- vf.data %>%
  t() %>%
  as.data.frame()
rs <- rowSums(vf.data.t) %>%
  as.data.frame()
colnames(rs) <- "sum"
vf.data.t$sum <- rs$sum
rank(vf.data.t$sum)
vf.data.t.ranked <- vf.data.t %>%
  rownames_to_column() %>%
  arrange(-sum) %>%
  column_to_rownames()
write.csv(vf.data.t.ranked, file = "files/ShortBRED_regrouped_ranked.csv", row.names = TRUE)



### Metadata
my_sample_data <- meta(dat) %>% sample_data()

### Create new Phyloseq Obj
my_VFDB.ab_table <- otu_table(vf.data, taxa_are_rows = F)
dat.VFs <- phyloseq(my_VFDB.ab_table, my_sample_data)
print(dat.VFs)
cat(" \n\nTrimming all VF's not detected \n")
dat.VFs <- core(dat.VFs, detection = 0, prevalence = 1 / 118)
save(dat.VFs, file = "files/PhyloseqObjects/VFs_PhyloseqObj.RData")



rm("l", "nwname", "dfadd", "dfseed")
