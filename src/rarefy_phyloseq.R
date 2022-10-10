# Caltech - Mazmanian Lab
# Joe Boktor
# October 2021

source("src/miscellaneous_funcs.R")

# _______________________________________________________________________________

rarefy_biobakery_phyloseq <- function(obj, pseudo_counts = FALSE, scaling_factor) {

  # Returns Rarefied Phyloseq Object
  df <- obj %>%
    core(detection = 0, prevalence = 1 / length(sample_names(obj))) %>%
    abundances() %>%
    as.data.frame()
  reads <- obj %>%
    meta() %>%
    dplyr::select(donor_id, total_reads) %>%
    as.data.frame() %>%
    mutate(total_reads = as.numeric(total_reads))
  
  if (pseudo_counts) {
    print(paste0("Column sum max is :", max(colSums(df))))
    
    for (i in colnames(df)) {
      donor_reads <-
        reads[[which(reads$donor_id == i), "total_reads"]] %>% as.numeric()
      if (max(colSums(df)) > 1) {
        df[i] <- df[i] / 100 * as.integer(donor_reads)
      } else {
        df[i] <- df[i] * as.integer(donor_reads)
      }
    }
  }
  
  df <- df * scaling_factor
  df.ints <- df %>%
    mutate_if(is.numeric, round) %>%
    t()
  
  print(paste0("Min (scaled): ", min(rowSums(df.ints))))
  print(paste0("Max (scaled): ", max(rowSums(df.ints))))
  psudocnts.rarefied <- rrarefy(df.ints, sample = min(rowSums(df.ints))) %>% t()
  print(paste0("Min rarefied (scaled): ", min(colSums(psudocnts.rarefied))))
  print(paste0("Max rarefied (scaled): ", max(colSums(psudocnts.rarefied))))
  cat("Rarefaction Complete\n")
  
  # Descaling
  psudocnts.rarefied.desc <- psudocnts.rarefied / scaling_factor
  # Recreate phyloseq obj
  dat_table <- otu_table(psudocnts.rarefied.desc, taxa_are_rows = T)
  my_sample_data <- meta(obj) %>% sample_data()
  datObjOut <- phyloseq(dat_table, my_sample_data)
  return(datObjOut)
}

# _______________________________________________________________________________

##### Rarefy data ---- TBC & RUSH ONLY

phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_slim_clean.rds")

psuedocount_levels <-
  c("Species",
    "Genus",
    "Family",
    "Order",
    "Class",
    "Phylum",
    "Kingdom")

phyloseq_objs_rare <- list()

for (obj in names(phyloseq_objs)) {
  
  if (obj %in% psuedocount_levels) {
    pseudo <- TRUE
    sfactor <- 1
  } else {
    pseudo <- FALSE
    sfactor <- 1
  }
  
  cat("\nRarify: ", obj, 
      "\nPseduo:", pseudo, 
      "\nScaling Factor:", sfactor, 
      "\n")
  
  phyloseq_objs_rare[[obj]] <-
    rarefy_biobakery_phyloseq(obj = phyloseq_objs[[obj]],
                              pseudo_counts = pseudo,
                              scaling_factor = sfactor)
}

names(phyloseq_objs_rare) <- names(phyloseq_objs)
saveRDS(phyloseq_objs_rare, file = "files/Phyloseq_Merged/PhyloseqObj_slim_clean_rarefied.rds")


# ##### Rarefy data ---- (FOR all four cohorts)
# load(file = "files/Phyloseq_Merged_ML.RData")
# 
# phyloseq_objs_rare <- list()
# 
# for (obj in names(phyloseq_objs)) {
#   
#   if (obj %in% psuedocount_levels) {
#     pseudo <- TRUE
#     sfactor <- 1
#   } else {
#     pseudo <- FALSE
#     sfactor <- 1e6
#   }
#   
#   cat("\nRarify: ", obj, 
#       "\nPseduo:", pseudo, 
#       "\nScaling Factor:", sfactor, 
#       "\n")
#   
#   phyloseq_objs_rare[[obj]] <-
#     rarefy_biobakery_phyloseq(obj = phyloseq_objs[[obj]],
#                               pseudo_counts = pseudo,
#                               scaling_factor = sfactor)
# }
# 
# names(phyloseq_objs_rare) <- names(phyloseq_objs)
# saveRDS(phyloseq_objs_rare, file = "files/Phyloseq_Merged_ML_Rarefied.rds")
# save(phyloseq_objs_rare, file = "files/Phyloseq_Merged_ML_Rarefied.RData")


# start_time <- Sys.time()
# cores <- detectCores()
# cl <- makeCluster(cores[1] - 1)
# registerDoParallel(cl)
# 
# phyloseq_objs_rare <-
#   foreach(
#     obj = names(phyloseq_objs),
#     .packages = c(
#       "phyloseq", "magrittr", "microbiome", "Hmisc",
#       "dplyr", "vegan"
#     )
#   ) %dopar% {
#     output <-
#       phyloseq_objs[[obj]] %>%
#       subset_samples(donor_id %nin%
#         c(
#           "PC.79", "PD.59", "PC.87", "PC.69", "PC.89",
#           "PC.72", "PD.58", "PC.59", "PD.56", "PC.60"
#         )) %>%
#       pseudoCountsRarefy()
#   }
# 
# stopCluster(cl)
# end_time <- Sys.time()
# cat(
#   "Rarefaction calculated in : ",
#   end_time - start_time,
#   attr(end_time - start_time, "units"), "\n"
# )
# names(phyloseq_objs_rare) <- names(phyloseq_objs)
# save(phyloseq_objs_rare, file = "files/Phyloseq_Merged_ML_Rarefied.RData")
