# Plot Select Features of interest

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")
source("src/Community_Composition_Funcs.R")
base::load("files/Phyloseq_Merged/PhyloseqObj.RData")

#-----------------------------------------------------------------------------------------------------
#           This script provides a scaffold to explore selected features of the dataset



#-------------------------------------------------------------------------------
# Step 1) Select an object type

# See github table
obj <- Phylo_Objects[["KOs.slim"]]

#-------------------------------------------------------------------------------
# Step 2) Browse data table
# (command click on data.table)

data_table <- explore_table(obj)

#-------------------------------------------------------------------------------
# Step 3) Select a feature of interest -
# Make sure spelling is correct and feature is found in data table (Step 2)
# Only one feature can match, Use unstratified objs if looking for general
# pathays, enzymes, or genes

feature <- "feat_K01186.gc..gs.sialidase.gh.1.gs..gsqrl.EC.gc.3.2.1.18.gsqrr."
# 1.8.4.12

#-------------------------------------------------------------------------------
# Step 4) Plot a feature of interest
# Note: uses ArcSin Sqrt Transformation

p1 <- plot_feature(obj, feature)
p1


# alpha_div_boxplots(df=env, x=env$donor_group, y=env$Evenness,
#                    df.pairs=env.pairs, df.pairs.x = env.pairs$donor_group, df.pairs.y=env.pairs$Evenness,
#                    pairs.column=env.pairs$Paired.plot,
#                    cols=color_palette, ylabel = paste0("Simpson's Evenness: ", z[cnt]),
#                    PDvPC.stat = evenness.PdPC.pval, PDvHC.stat = evenness.PdHC.pval)

library(plotly)
ggplotly(p1, tooltip = "all")




#-------------------------------------------------------------------------------
# Random metadata plotting

d <- obj %>%
  abundances() %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  filter(rowname == feature) %>%
  column_to_rownames() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column()

metadat <- meta(obj) %>%
  rownames_to_column()

plotme <- d %>%
  full_join(metadat)

plotme %>%
  filter(PD == "Yes") %>%
  ggplot(aes(x = amantadine, y = feat_K01186.gc..gs.sialidase.gh.1.gs..gsqrl.EC.gc.3.2.1.18.gsqrr.)) +
  geom_boxplot() +
  stat_compare_means(method = "anova") +
  geom_point()
