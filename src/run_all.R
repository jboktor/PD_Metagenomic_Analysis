#' This script will run all analyses downstream of the biobakery pipeline to
#' generate all data-tables and figures used in our publication

#' run once initialize environment  and install necessary packages
#' if errors persist, run the configure_enviornment.R script or install 
#' missing packages individually
renv::init()

# Analysis Prep
source("src/create_phyloseq_obj.R")
source("src/quality_control.R")
source("src/rarefy_phyloseq.R")
source("src/create_phyloseq_obj_external_cohorts.R")

# Data exploration and quality assurance
source("src/sample_summary.R")
source("src/tbc-survey-heatmap.R") # Figure 1S:

# Community Composition Analyses
source("src/permanova_analysis.R") 
source("src/permanova_vis.R") # Figure 1B
source("src/community_composition_overview.R") # Figure 1/1S

# Figure 2/3: Differential Abundance Analysis
source("src/dafs.R")
source("src/dafs_vis_Figure-2.R")
source("src/dafs_control_analyses.R")
source("src/dafs_PC-vs-HC.R")
source("src/enrichment_analysis.R")
source("src/enrichment_vis.R")

# Figure 4: Metadata Correlations
source("src/metadata_correlations.R")

# Figure 5/5S: Meta-Analyses 
source("src/feature_auroc.R")
source("src/meta-analysis-MMUPHin.R")
source("src/meta-analysis_vis.R")