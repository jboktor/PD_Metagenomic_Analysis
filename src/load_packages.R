# Load Packages

## Load Packages into enviornment
library(ggplot2)
library(tidyverse)
library(readxl)
library(openxlsx)
library(data.table)
# library(xlsx)
library(plyr)
library(dplyr)
library(ggrepel)
library(grid)

library(gridExtra)
library(reshape2)
library(grid)
library(devtools)
library(RColorBrewer)

library(ggfortify)
library(vegan)
library(MASS)
library(compositions)
library(zCompositions)
library(phyloseq)

library(Biobase)
library(viridis)
library(ggbeeswarm)

library(FSA)
library(ggpubr)
library(ggsci)
library(microbiome)
library(ggridges)
library(future)
library(cowplot)
library(ggdendro)
library(ggdist)
library(GGally)

library(EnvStats)
library(sjlabelled)
library(sjmisc)
library(sjPlot)
library(nlme)
library(lme4)
library(jtools)
library(eulerr)
library(ggthemes)

library(ggforce)
library(huge)
library(Matrix)
library(Maaslin2)
library(MMUPHin)
library(fantaxtic)
library(igraph)
library(slam)
# library(curatedMetagenomicData)

library(scales)
library(DirichletMultinomial)
library(magrittr)
library(randomForest)
library(pROC)
library(plotROC)
library(SGL)

library(rstatix)
library(ppcor)
library(parallel)
library(doParallel)
library(foreach)
library(mlbench)

library(caret)
library(MLeval)
library(tidymodels)
library(vip)
library(Rcpp)
library(mRMRe)
library(ranger)
# library(lightgbm)
library(snm)

library(rsconnect)
# library(phylosmith)

library(lazyeval)

# library(rbiom)
library(ape)
# library(wesanderson)
# library(ClusterR)
library(mclust)
library(gmp)
library(jsonlite)
library(conflicted)

conflict_prefer("print", "base")
conflict_prefer("matrix", "base")
conflict_prefer("apply", "base")

conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("summarize", "dplyr")
conflict_prefer("summarise", "dplyr")

conflict_prefer("desc", "plyr")

conflict_prefer("melt", "reshape2")

conflict_prefer("replace_na", "tidyr")

conflict_prefer("map", "purrr")

conflict_prefer("%nin%", "sjmisc")

conflict_prefer("get_legend", "cowplot")

conflict_prefer("dist", "stats")

conflict_prefer("sd", "BiocGenerics")

conflict_prefer("roc", "pROC")

conflict_prefer("box", "shinydashboard")
