### SCFA Analysis ### 

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/daf_functions.R")
load("files/Phyloseq_Merged/GBMs_PhyloseqObj.RData")
load("files/Phyloseq_Merged/GMMs_PhyloseqObj.RData")
load("files/Phyloseq_Merged/Species_PhyloseqObj.RData")

# Script to print out boxplots with all three groups for 
# all associated features (via MaAsLin2 analysis)
# Note: this is only designed for merged analysis

meta_join <- dat.species %>% 
  meta() %>% 
  select(donor_id, donor_group, cohort)

dat <- dat.GMMs ## Dat object swap as needed
datName <- "GMMs"

df.input <- 
  dat %>% 
  abundances() %>% 
  sqrt() %>% 
  asin() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "feature") %>% 
  pivot_longer(!feature, names_to = "donor_id") %>% 
  left_join(meta_join, by = "donor_id") %>% 
  mutate(donor_group = factor(donor_group, levels = c("PC", "PD","HC")))

# Pull stats from MaAsLin2
Maas.pd.pc <- read_tsv(
  paste0("data/MaAsLin2_Analysis/Merged/", datName, "_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Population Control")
Maas.pd.hc <- read_tsv(
  paste0("data/MaAsLin2_Analysis/Merged/", datName, "_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
  filter(value == "Household Control")

#####################################################################################
PlotGBMs <- function(df.plt, df.targets, Maas.pd.pc, Maas.pd.hc) {
  
  for (i in unique(df.targets$feature)) {
    
    cat("Plotting Feature: ", i, "\n")
    Maas.pd.pc.Q <- Maas.pd.pc %>% filter(feature == i)
    cat("MaAsLin2 PDvPC GLM Q-value: ", Maas.pd.pc.Q$qval, "\n")
    Maas.pd.hc.Q <- Maas.pd.hc %>% filter(feature == i)
    cat("MaAsLin2 PDvHC GLM Q-value: ", Maas.pd.hc.Q$qval, "\n")
    
    df.plt1 <- df.plt %>% 
      filter(feature == i )
    
    p <- boxplot_all(df.plt1, x=df.plt1$donor_group, y=df.plt1$value,
                     cols=c("PC"= "#ed7d31", "PD" = "#bfbfbf", "HC" = "#5b9bd5"), 
                     ylabel="Normalized Abundance", 
                     title= gsub("\\.", " ",  df.plt1$feature))
    p <- p + geom_signif(comparisons = list(c("PD", "HC")), #tip_length = 0.02,
                         annotations = sig_mapper(Maas.pd.hc.Q$qval, porq = "q", symbols = F)) +
      geom_signif(comparisons = list(c("PC", "PD")), #tip_length = 0.02,
                  annotations = sig_mapper(Maas.pd.pc.Q$qval, porq = "q", symbols = F)) +
      theme(legend.position = "none",
            panel.grid.major.x = element_blank(),
            plot.title = element_text(size = 14),
            axis.title.y = element_text(size = 12),
            axis.line.x = element_line(colour="grey"),
            axis.line.y = element_line(colour="grey"))
    
    plot(p)
    ggsave(p, filename = paste0("data/boxplots/", datName, "_",i ,".svg"), height = 4, width =3)
  }
}

#####################################################################################


##### Filter all significant plots from GBM Model 

Maas.pd.hc.sig <- Maas.pd.hc %>% 
  dplyr::filter(qval < 0.25) %>% 
  dplyr::select(feature)
Maas.pd.pc.sig <- Maas.pd.pc %>% 
  dplyr::filter(qval < 0.25) %>% 
  dplyr::select(feature)

df.pdhc.sig <- df.input %>% 
  filter(feature %in% Maas.pd.hc.sig$feature | 
           feature %in% Maas.pd.pc.sig$feature)

PlotGBMs(df.plt = df.input,
         df.targets = df.pdhc.sig,
         Maas.pd.pc = Maas.pd.pc, 
         Maas.pd.hc = Maas.pd.hc)



#__________________
#  Scrap

# ##### Acetate Plots ##### 
# ## Explore acetate variables
# df.acetate <- df.plt %>% 
#   filter(grepl("Acetate", feature, ignore.case = T))
# unique(df.acetate$feature)
# 
# PlotGBMs(df.plt = df.plt,
#          df.targets = df.acetate,
#          Maas.pd.pc = Maas.pd.pc, 
#          Maas.pd.hc = Maas.pd.hc)
# 
# ##### Propionate 
# ## Explore Propionate variables
# df.propionate <- df.plt %>% 
#   filter(grepl("Propionate", feature, ignore.case = T))
# unique(df.propionate$feature)
# 
# PlotGBMs(df.plt = df.plt,
#          df.targets = df.propionate,
#          Maas.pd.pc = Maas.pd.pc, 
#          Maas.pd.hc = Maas.pd.hc)
# 
# ##### Butyrate
# ## Explore Butyrate variables
# df.butyrate <- df.plt %>% 
#   filter(grepl("Butyrate", feature, ignore.case = T))
# unique(df.butyrate$feature)
# 
# PlotGBMs(df.plt = df.plt,
#          df.targets = df.butyrate,
#          Maas.pd.pc = Maas.pd.pc, 
#          Maas.pd.hc = Maas.pd.hc)
# 
# ##### Valerate
# ## Explore Butyrate variables
# df.valerate <- df.plt %>% 
#   filter(grepl("Valer", feature, ignore.case = T))
# unique(df.valerate$feature)
# 
# PlotGBMs(df.plt = df.plt,
#          df.targets = df.valerate,
#          Maas.pd.pc = Maas.pd.pc, 
#          Maas.pd.hc = Maas.pd.hc)
