## Configure enviornment

#### Install Bioconductor Packages
install.packages("BiocManager")
BiocManager::install(
  c("phyloseq", "microbiome", "Biobase", "Maaslin2"))

#### Install CRAN Packages
install.packages(
  c("ggplot2", "tidyverse", "readxl", "dplyr", "ggrepel", "gridExtra", "reshape2", "plyr", 
    "devtools", "RColorBrewer", "ggfortify", "vegan", "MASS", "compositions", "zCompositions",
    "gplots", "viridis", "lme4", "phangorn", "plotly", "VennDiagram", "viridis","foreach", "doParallel" , 
    "ggbeeswarm", "FSA", "ggpubr", "ggsci", "ggridges", "future", "svglite", "cowplot", "coin", "EnvStats", 
    "sjlabelled", "sjmisc", "sjPlot", "nlme", "eulerr", "ggthemes", "ggforce", "huge", "Matrix"))

## Load Packages into enviornment
library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot);
library(EnvStats);library(sjlabelled);library(sjmisc);library(sjPlot);library(nlme);library(eulerr);library(ggthemes);
library(ggforce);library(huge);library(Matrix)



