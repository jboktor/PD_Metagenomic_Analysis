# SIAMCAT analysis

## ----start, message=FALSE, warning=FALSE--------------------------------------
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
load("files/low_quality_samples.RData")
remove_dats()
load_data("Merged_ML")
library(SIAMCAT)

# this is data from Zeller et al., Mol. Syst. Biol. 2014
fn.feat.fr  <-
  'https://www.embl.de/download/zeller/FR-CRC/FR-CRC-N141_tax-ab-specI.tsv'
fn.meta.fr  <-
  'https://www.embl.de/download/zeller/FR-CRC/FR-CRC-N141_metadata.tsv'

# this is the external dataset from Yu et al., Gut 2017
fn.feat.cn  <-
  'https://www.embl.de/download/zeller/CN-CRC/CN-CRC-N128_tax-ab-specI.tsv'
fn.meta.cn  <-
  'https://www.embl.de/download/zeller/CN-CRC/CN-CRC-N128_metadata.tsv'

## ----siamcat_fr---------------------------------------------------------------
# # features
# # be vary of the defaults in R!!!
# feat.fr  <- read.table(fn.feat.fr, sep='\t', quote="",
#                        check.names = FALSE, stringsAsFactors = FALSE)
# # the features are counts, but we want to work with relative abundances
# feat.fr.rel <- prop.table(as.matrix(feat.fr), 2)
# 
# # metadata
# meta.fr  <- read.table(fn.meta.fr, sep='\t', quote="",
#                        check.names=FALSE, stringsAsFactors=FALSE)
# 
# # create SIAMCAT object
# siamcat.fr <- siamcat(feat=feat.fr.rel, meta=meta.fr,
#                       label='Group', case='CRC')

# # features
# # be vary of the defaults in R!!!
# feat.fr  <- read.table(fn.feat.fr, sep='\t', quote="",
#                        check.names = FALSE, stringsAsFactors = FALSE)
# # the features are counts, but we want to work with relative abundances
# feat.fr.rel <- prop.table(as.matrix(feat.fr), 2)
# 
# # metadata
# meta.fr  <- read.table(fn.meta.fr, sep='\t', quote="",
#                        check.names=FALSE, stringsAsFactors=FALSE)

# create SIAMCAT object
obj <- subset_samples(dat.KOs.slim, donor_id %ni% low_qc[[1]]) #%>% 
  # microbiome::transform('compositional')
obj_noShanghai <- subset_samples(obj, cohort != "Shanghai")
obj_Shanghai <- subset_samples(obj, cohort == "Shanghai")
# obj_noRush <- subset_samples(obj, cohort != "Rush")
# obj_Rush <- subset_samples(obj, cohort == "Rush")
# obj_noTBC <- subset_samples(obj, cohort != "TBC")
# obj_TBC <- subset_samples(obj, cohort == "TBC")
# obj_noBonn <- subset_samples(obj, cohort != "Bonn")
# obj_Bonn <- subset_samples(obj, cohort == "Bonn")

abundances_noShanghai <- abundances(obj_noShanghai) 
meta_noShanghai <- meta(obj_noShanghai)
abundances_Shanghai <- abundances(obj_Shanghai)
meta_Shanghai <- meta(obj_Shanghai)

siamcat.noShanghai <- siamcat(feat=abundances_noShanghai, meta=meta_noShanghai,
                      label='PD', case='Yes')

siamcat.Shanghai <- siamcat(feat=abundances_Shanghai, meta=meta_Shanghai,
                       label='PD', case='Yes')

## ----preprocessing_fr---------------------------------------------------------
siamcat.noShanghai <- filter.features(
  siamcat.noShanghai,
  filter.method = 'abundance',
  cutoff = 1e-08,
  rm.unmapped = TRUE,
  verbose=2
)

siamcat.noShanghai <- normalize.features(
  siamcat.noShanghai,
  norm.method = "log.std",
  norm.param = list(log.n0 = 1e-06, sd.min.q = 0.1),
  verbose = 2
)

## ----build_model_fr, results='hide'-------------------------------------------
siamcat.noShanghai <-  create.data.split(
  siamcat.noShanghai,
  num.folds = 5,
  num.resample = 2
)

siamcat.noShanghai <- train.model(
  siamcat.noShanghai,
  method = "lasso"
)

## ----predict_evaluate_fr, results='hide'--------------------------------------
siamcat.noShanghai <- make.predictions(siamcat.noShanghai)

siamcat.noShanghai <-  evaluate.predictions(siamcat.noShanghai)

## ----normalize_cn-------------------------------------------------------------

siamcat.Shanghai <- normalize.features(
  siamcat.Shanghai,
  norm.param = norm_params(siamcat.noShanghai),
  feature.type = 'original',
  verbose = 2
)


## ----predict_cn, results='hide'-----------------------------------------------
siamcat.Shanghai <- make.predictions(
  siamcat = siamcat.noShanghai,
  siamcat.holdout = siamcat.Shanghai,
  normalize.holdout = FALSE)

## ----alternative_pipeline_cn, eval=FALSE--------------------------------------
 ## Alternative Code, not run here
 # siamcat.cn <- siamcat(feat=feat.cn.rel, meta=meta.cn,
 #     label='Group', case='CRC')
 # siamcat.cn <- make.predictions(siamcat = siamcat.fr,
 #     siamcat.holdout = siamcat.cn,
 #     normalize.holdout = TRUE)

## ----eval_cn, message=FALSE---------------------------------------------------
siamcat.Shanghai <- evaluate.predictions(siamcat.Shanghai)

## ----eval_plot, eval=FALSE----------------------------------------------------
 model.evaluation.plot('noShanghai'=siamcat.noShanghai,
     'Shanghai'=siamcat.Shanghai,
     colours=c('dimgrey', 'orange'))

## ----session_info-------------------------------------------------------------
sessionInfo()
