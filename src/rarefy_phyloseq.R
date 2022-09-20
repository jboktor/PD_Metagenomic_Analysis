# Caltech - Mazmanian Lab
# Joe Boktor
# October 2021

source("src/miscellaneous_funcs.R")

#_______________________________________________________________________________

pseudoCountsRarefy <- function(obj){
  
  # Returns Rarefied Phyloseq Object
  
  # # TROUBLE
  # obj <- phyloseq_objs[["KOs.slim"]] %>%
  #   subset_samples(donor_id %ni% low_qc[[1]])
  # df <- obj %>%
  #   core(detection = 0, prevalence = 1/length(sample_names(obj))) %>%
  #   abundances() %>%
  #   as.data.frame()
  # reads <- obj %>% meta() %>%
  #   dplyr::select(donor_id, total_reads) %>%
  #   as.data.frame()
  
  load("files/low_quality_samples.RData")
  
  df <- obj %>% 
    core(detection = 0, prevalence = 1/length(sample_names(obj))) %>% 
    abundances() %>%
    as.data.frame()
  reads <- obj %>% meta() %>%
    dplyr::select(donor_id, total_reads) %>%
    as.data.frame() %>% 
    mutate(total_reads = as.numeric(total_reads))
  
  psudocnts <- df
  print(paste0("Column sum max is :", max(colSums(df))))
  
  for (i in colnames(psudocnts)){
    donor_reads <- reads[[which(reads$donor_id == i), "total_reads"]] %>% as.numeric()
    if (max(colSums(df)) > 1){
      psudocnts[i] <- psudocnts[i]/100 * as.integer(donor_reads)
    } else{
      psudocnts[i] <- psudocnts[i] * as.integer(donor_reads)
    }
  }
  psudocnts.ints <- psudocnts %>% mutate_if(is.numeric, round) %>% t()
  print(paste0("Minimum non-rarefied: ", min(rowSums(psudocnts.ints))))
  print(paste0("Maximum non-rarefied: ", max(rowSums(psudocnts.ints))))
  psudocnts.rarefied <- rrarefy(psudocnts.ints, sample = min(rowSums(psudocnts.ints))) %>% t() 
  print(paste0("Minimum rarefied: ", min(colSums(psudocnts.rarefied))))
  print(paste0("Maximum rarefied: ", max(colSums(psudocnts.rarefied))))
  
  
  cat("Pseudocount Transformation Complete\n")
  dat_table <- otu_table(psudocnts.rarefied, taxa_are_rows=T)
  my_sample_data <- meta(obj) %>% sample_data()
  datObjOut <- phyloseq(dat_table, my_sample_data)
  return(datObjOut)
}

#_______________________________________________________________________________

##### Rarefy data ---- TBC & RUSH ONLY

phyloseq_objs <- readRDS("files/Phyloseq_Merged/PhyloseqObj_clean.rds")

start_time <- Sys.time()
cores = detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

phyloseq_objs_rare <-
  foreach(
    obj = names(phyloseq_objs),
    .packages = c("phyloseq", "magrittr", "microbiome", "Hmisc", 
                  "dplyr", "vegan")
  ) %dopar% {
    
    output <-
      phyloseq_objs[[obj]] %>%
      pseudoCountsRarefy()
  }

stopCluster(cl)
end_time <- Sys.time()
cat("Rarefaction calculated in : ",
    end_time - start_time,
    attr(end_time - start_time, "units"), "\n")

names(phyloseq_objs_rare) <- names(phyloseq_objs)
# save(phyloseq_objs_rare, file = "files/Phyloseq_Merged/PhyloseqObj_clean_rarefied.rds")
saveRDS(phyloseq_objs_rare, file = "files/Phyloseq_Merged/PhyloseqObj_clean_rarefied.rds")



##### Rarefy data ----

load(file = "files/Phyloseq_Merged_ML.RData")

start_time <- Sys.time()
cores = detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

phyloseq_objs_rare <-
  foreach(
    obj = names(phyloseq_objs),
    .packages = c("phyloseq", "magrittr", "microbiome", "Hmisc", 
                  "dplyr", "vegan")
  ) %dopar% {
    
    output <-
      phyloseq_objs[[obj]] %>%
      subset_samples(donor_id %nin% 
                       c("PC.79", "PD.59", "PC.87", "PC.69", "PC.89", 
                         "PC.72", "PD.58", "PC.59", "PD.56", "PC.60")) %>%
      pseudoCountsRarefy()
  }

stopCluster(cl)
end_time <- Sys.time()
cat("Rarefaction calculated in : ",
    end_time - start_time,
    attr(end_time - start_time, "units"), "\n")
names(phyloseq_objs_rare) <- names(phyloseq_objs)
save(phyloseq_objs_rare, file = "files/Phyloseq_Merged_ML_Rarefied.RData")





