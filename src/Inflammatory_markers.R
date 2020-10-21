# Inflammatory Markers

######## Load Data & functions
source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/DAF_Functions.R")
source("src/Stats_funcs.R")

#-------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------

# Plot feature of interest (foi)
plot.foi <- function(datObj, foi, title="" ){
  
  asin.input <- datObj %>%
    microbiome::transform("compositional") %>% 
    microbiome::abundances()
  df1 <- asin(sqrt(asin.input)) %>% 
    as.data.frame() %>% 
    rownames_to_column() 
  d <- df1 %>% 
    filter(rowname == foi) %>% 
    melt()
  dm <-  group_col_from_ids(d, id=d$variable)
  dm$group <- factor(dm$group, levels = c("PC", "PD", "HC"))
  boxplot_all(dm, x=dm$group, y=dm$value, 
              title=title, 
              ylabel= paste(unique(dm$rowname), "Abundance"))
  
}

#-------------------------------------------------------------------


# 2) CLR Log-odds ratio : facultative anaerobic / obligate anaerobic bacteria

#' Function returns normalized values for feature(s) of interest
#' 1) Adaptive to multiple types of normalization
#' 2) Will sum values if a list of foi's are input

fois <- function(datObj, foi) {
  
  if (length(foi) > 1) {
    
    # Init
    cnt <- 1
    
    for (i in foi) {
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
      if (cnt == 1){
        summed.vars <- d
      } else {
        summed.vars$value <- summed.vars$value + d$value
      }
      cnt <- cnt + 1
    }
    return(summed.vars)
    
  } else {
    for (i in foi){
      cat("selecting feature: ", i, "\n")
      d <- datObj %>%
        microbiome::transform("compositional") %>%
        microbiome::abundances() %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        filter(rowname == i) %>%
        melt() %>%
        dplyr::select(c(variable, value))
    }
    return(d)
  }
}

#-------------------------------------------------------------------
#-------------------------------------------------------------------


# Taxon level analysis of Inflammatory Markers
# Dominant obligate anaerobic bacteria Class (Clostridia, Bacteroidia ) 
# Dominant facultative anaerobic bacteria Phylum (Proteobacteria), Family (Enterobacteriaceae)


# obligate anaerobic bacteria
Clostridia <-
  plot.foi(dat.class, foi = "Clostridia")
ggsave(Clostridia,filename = "data/Inflammation_Markers/Clostridia.png", 
       height = 4, width =3)

Bacteroidia <-
  plot.foi(dat.class, foi = "Bacteroidia")
ggsave(Bacteroidia, filename = "data/Inflammation_Markers/Bacteroidia.png", 
       height = 4, width =3)


# facultative anaerobic bacteria
Proteobacteria <-
  plot.foi(dat.phylum, foi = "Proteobacteria")
ggsave(Proteobacteria, filename = "data/Inflammation_Markers/Proteobacteria.png", 
       height = 4, width =3)

Enterobacteriaceae <-
  plot.foi(dat.family, foi = "Enterobacteriaceae")
ggsave(Enterobacteriaceae, filename = "data/Inflammation_Markers/Enterobacteriaceae.png", 
       height = 4, width =3)


#------------------------------------------------------
# Ratio of Facultative to Obligate Anaerobes
#------------------------------------------------------

FA <- fois(datObj = dat.family, foi = "Enterobacteriaceae")
OA <- fois(datObj = dat.class, foi = c("Clostridia", "Bacteroidia"))
ratio.plot <- FA
ratio.plot$ratio <- FA$value/OA$value
ratio.plot <- group_col_from_ids(ratio.plot, id = ratio.plot$variable) %>% as.data.frame()
ratio.plot$group <- factor(ratio.plot$group, levels = c("PC", "PD", "HC"))

FAvOA <- boxplot_all(df=ratio.plot, x=ratio.plot$group, y=(ratio.plot$ratio),
            title = "Facultative / Obligate Anerobe Ratio",
            ylabel = "Normalized Abundance Ratio")
ggsave(FAvOA, filename = "data/Inflammation_Markers/Facultative_vs_Obligate_Anaerobe_Abundance.png", 
       height = 4, width =4)
FAvOA_LR <- boxplot_all(df=ratio.plot, x=ratio.plot$group, y=log(ratio.plot$ratio),
            title = "Facultative / Obligate Anerobe Ratio",
            ylabel = "Natural Log-Ratio")
ggsave(FAvOA_LR, filename = "data/Inflammation_Markers/Facultative_vs_Obligate_Anaerobe_LogRatio.png", 
       height = 4, width =4)

   




#-------------------------------------------------------------------------------------
# AEROBIC RESPIRATION
#-------------------------------------------------------------------------------------

# Genes involved in Facultative Anerobes distinction (02 and nitrate consumption)
 
# ENZYMES INVOLVED IN AEROBIC RESPIRATION
# The cytochrome oxidases
# Cytochrome bo oxidase (CyoABCD):   1.10.3.10
# Cytochrome bd-I/bd-II oxidase (CydAB/AppBC).  1.10.3.14      
   
# The dehydrogenases:
# NADH dehydrogenase I (NuoABCDEFGHIJKLMN) : EC:1.6.5.3
# NADH dehydrogenase II (Ndh-2) :  EC 1.6.5.11 
# Succinate dehydrogenase (SdhCDAB) : EC 1.3.5.1  
# L-lactate dehydrogenase: [EC:1.1.2.3]
# glycerol-phosphate dehydrogenase:  
# proline dehydrogenase


# Cytochrome bo oxidase (CyoABCD):   1.10.3.10
CyoABCD <- plot.foi(dat.ec.slim, foi = "1.10.3.10" )
CyoABCD

# Cytochrome bd-I/bd-II oxidase (CydAB/AppBC).  1.10.3.14      
CydAB <- plot.foi(dat.ec.slim, foi = "1.10.3.14")
CydAB

# NADH dehydrogenase I (NuoABCDEFGHIJKLMN) : EC:1.6.5.3
NuoABCDEFGHIJKLMN <- plot.foi(dat.ec.slim, foi = "1.6.5.3")
NuoABCDEFGHIJKLMN
# NADH dehydrogenase II (Ndh-2) :  EC 1.6.5.11 
Ndh2 <- plot.foi(dat.ec.slim, foi = "1.6.5.11")
Ndh2
# Succinate dehydrogenase (SdhCDAB) : EC 1.3.5.1  
SdhCDAB <- plot.foi(dat.ec.slim, foi = "1.3.5.1")
SdhCDAB
# L-lactate dehydrogenase: [EC:1.1.2.3]
ldh <- plot.foi(dat.ec.slim, foi = "1.1.2.3")
ldh


ggsave(CyoABCD, filename = "data/Inflammation_Markers/aerobic_respiration_Cytochrome_oxidase_CyoABCD.png", 
       height = 4, width =4)
ggsave(CydAB, filename = "data/Inflammation_Markers/aerobic_respiration_Cytochrome_oxidase_CydAB.png", 
       height = 4, width =4)
ggsave(NuoABCDEFGHIJKLMN, filename = "data/Inflammation_Markers/aerobic_respiration_NuoABCDEFGHIJKLMN.png", 
       height = 4, width =4)
ggsave(Ndh2, filename = "data/Inflammation_Markers/aerobic_respiration_dehydrogenase_Ndh2.png", 
       height = 4, width =4)
ggsave(SdhCDAB, filename = "data/Inflammation_Markers/aerobic_respiration_dehydrogenase_SdhCDAB.png", 
       height = 4, width =4)
ggsave(ldh, filename = "data/Inflammation_Markers/aerobic_respiration_dehydrogenase_ldh.png", 
       height = 4, width =4)




# Cytochrome oxidases Aerobic respiration
anerobic.respiration.cytochrome_oxidases <- c("1.10.3.10", "1.10.3.14")
air.ecs <- fois(datObj = dat.ec.slim, foi = anerobic.respiration.cytochrome_oxidases)

ratio.plot <- air.ecs
ratio.plot <- group_col_from_ids(ratio.plot, id=ratio.plot$variable) %>% as.data.frame()
ratio.plot$group <- factor(ratio.plot$group, levels = c("PC", "PD", "HC"))

boxplot_all(df=ratio.plot, x=ratio.plot$group, y=(ratio.plot$value),
            title = "Aerobic Respiration Cytochrome Oxidases",
            ylabel = "Normalized Abundance Ratio")
boxplot_all(df=ratio.plot, x=ratio.plot$group, y=log(ratio.plot$value),
            title = "Aerobic Respiration Cytochrome Oxidases",
            ylabel = "Natural Log-Ratio")




#-------------------------------------------------------------------------------------
# ANAEROBIC RESPIRATION
#-------------------------------------------------------------------------------------
### nitrate reductases
#---------------------------------------------------

# narGHJI/narZYWV nitrate reductase (quinone);  EC 1.7.5.1   
# OLD ENZYME ID 1.7.99.4 
nitrate_reductase <- plot.foi(dat.ec.slim, foi = "1.7.99.4" )
nitrate_reductase

# napFDAGHBC KO: K02019
#' The molybdate-responsive Escherichia coli ModE transcriptional 
#' regulator coordinates periplasmic nitrate reductase (napFDAGHBC) 
#' operon expression with nitrate and molybdate availability.
nitrate_reductase_reg <- plot.foi(dat.KOs.slim, foi = "K02019" )
nitrate_reductase_reg


ggsave(nitrate_reductase, filename = "data/Inflammation_Markers/anaerobic_respiration_nitrate_reductase.png", 
       height = 4, width =4)
ggsave(nitrate_reductase_reg, filename = "data/Inflammation_Markers/anaerobic_respiration_nitrate_reductase_reg.png", 
       height = 4, width =4)

#---------------------------------------------------
# S-oxide reductases
#---------------------------------------------------

# dmsABC, ynfDEFGH
# anaerobic dimethyl sulfoxide reductase subunit A [EC:1.8.5.3]

sulfoxide_reductase <- plot.foi(dat.ec.slim, foi = "1.8.5.3" )
sulfoxide_reductase

ggsave(sulfoxide_reductase, filename = "data/Inflammation_Markers/anaerobic_respiration_sulfoxide_reductase.png", 
       height = 4, width =4)
#---------------------------------------------------
# N-oxide reductases
#---------------------------------------------------

# torCAD - K07772
# two-component system, OmpR family, torCAD operon response regulator TorR
torR <- plot.foi(dat.KOs.slim, foi = "K07772" )
torR

# torYZ - K07821
# trimethylamine-N-oxide reductase (cytochrome c), cytochrome c-type subunit TorY
torY <- plot.foi(dat.KOs.slim, foi = "K07821" )
torY

# yedYZ - K07147
# methionine sulfoxide reductase catalytic subunit [EC:1.8.-.-]
yedY <- plot.foi(dat.KOs.slim, foi = "K07147" )
yedY


ggsave(torR, filename = "data/Inflammation_Markers/anaerobic_respiration_N-oxide_reductases_torR.png", 
       height = 4, width =4)
ggsave(torY, filename = "data/Inflammation_Markers/anaerobic_respiration_N-oxide_reductases_torY.png", 
       height = 4, width =4)
ggsave(yedY, filename = "data/Inflammation_Markers/anaerobic_respiration_N-oxide_reductases_yedY.png", 
       height = 4, width =4)



#------------------------------------------------------
# Anaerobic Process
#------------------------------------------------------

aerobic.respiration.ecs <- c("1.7.99.4", "1.8.5.3")
aerobic.respiration.kos <- c("K02019", "K07772", "K07821", "K07147")
anaer.ecs <- fois(datObj = dat.ec.slim, foi = aerobic.respiration.ecs)
anaer.kos <- fois(datObj = dat.KOs.slim, foi = aerobic.respiration.kos)


ratio.plot <- anaer.ecs
ratio.plot <- group_col_from_ids(ratio.plot, id=ratio.plot$variable) %>% as.data.frame()
ratio.plot$group <- factor(ratio.plot$group, levels = c("PC", "PD", "HC"))

boxplot_all(df=ratio.plot, x=ratio.plot$group, y=(ratio.plot$value),
            title = "Anaerobic Respiration - Enzymes",
            ylabel = "Normalized Abundance Ratio")
boxplot_all(df=ratio.plot, x=ratio.plot$group, y=log(ratio.plot$value),
            title = "Anaerobic Respiration - Enzymes",
            ylabel = "Natural Log")

ratio.plot <- anaer.kos
ratio.plot <- group_col_from_ids(ratio.plot, id=ratio.plot$variable) %>% as.data.frame()
ratio.plot$group <- factor(ratio.plot$group, levels = c("PC", "PD", "HC"))

boxplot_all(df=ratio.plot, x=ratio.plot$group, y=(ratio.plot$value),
            title = "Anaerobic Respiration - KOs",
            ylabel = "Normalized Abundance Ratio")
boxplot_all(df=ratio.plot, x=ratio.plot$group, y=log(ratio.plot$value),
            title = "Anaerobic Respiration - KOs",
            ylabel = "Natural Log")


















#-------------------------------------------------------------------------------------
# Inflammatory Diseases
#-------------------------------------------------------------------------------------


# # Obligate Anerobes also encode a  broad spectrum of enzymes for hydrolyzing different complex carbohydrates
## Enrichment 



# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "devel")
# BiocManager::valid()
# BiocManager::install("curatedMetagenomicData")


library(curatedMetagenomicData)

# ?combined_metadata
# View(combined_metadata)
unique(combined_metadata$body_site)
stool_samples <- combined_metadata %>%
  filter(body_site == "stool")

## Explore studies
unique(stool_samples$disease)
# testg <- filter(combined_metadata, grepl("RA", disease))
# testg$dataset_name


#------------------------------------------------------------
# IBD
#------------------------------------------------------------

NielsenHB_2014.metaphlan_bugs_list.stool() %>% 
  experimentData()
study_vars <- stool_samples %>% 
  filter(dataset_name == "NielsenHB_2014") %>% 
  select(study_condition)
unique(study_vars$study_condition)

study <- "NielsenHB_2014"
dat.IBD.bugs <- curatedMetagenomicData(paste0(study, ".metaphlan_bugs_list.stool"), dryrun=F)
dat.IBD.paths <- curatedMetagenomicData(paste0(study, ".pathabundance_relab.stool"), dryrun=F)

df.spec <- dat.IBD.bugs[[1]] %>%
  ExpressionSet2phyloseq() %>% 
  subset_taxa(!is.na(Species)) %>% 
  subset_taxa(is.na(Strain)) 


# obligate anaerobic bacteria
Clostridia <-
  plot.foi(df.spec, foi = "Clostridia")
# ggsave(Clostridia,filename = "data/Inflammation_Markers/Clostridia.png", 
#        height = 4, width =3)

Bacteroidia <-
  plot.foi(dat.class, foi = "Bacteroidia")
# ggsave(Bacteroidia, filename = "data/Inflammation_Markers/Bacteroidia.png", 
#        height = 4, width =3)


# facultative anaerobic bacteria
Proteobacteria <-
  plot.foi(dat.phylum, foi = "Proteobacteria")
# ggsave(Proteobacteria, filename = "data/Inflammation_Markers/Proteobacteria.png", 
#        height = 4, width =3)

Enterobacteriaceae <-
  plot.foi(dat.family, foi = "Enterobacteriaceae")
# ggsave(Enterobacteriaceae, filename = "data/Inflammation_Markers/Enterobacteriaceae.png", 
#        height = 4, width =3)



# genefamilies_relab.stool
# pathabundance_relab.stool




