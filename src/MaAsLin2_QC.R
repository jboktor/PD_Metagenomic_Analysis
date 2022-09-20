# Maaslin_QC

library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(reshape2)
library(ggdendro)
library(plyr)
library(fastcluster)
library(dendextend)
library(grid)
library(devtools)
library(RColorBrewer)
library(ggfortify)
library(vegan)
library(MASS)
library(compositions)
library(zCompositions)
library(phyloseq)
library(gplots)
library(ape)
library(lme4)
library(phangorn)
library(plotly)
library(VennDiagram)
library(ggvegan)
library(Biobase)
library(BiocInstaller)
library(viridis)
library("foreach")
packageVersion("foreach")
library("doParallel")
packageVersion("doParallel")
library(jakR)
library(ggbeeswarm)
library(FSA)
library(ggpubr)
library(ggsci)
library(microbiome)
library(ggridges)
library(future)
library(cowplot)
library(scales)


# Get Data
load("Metaphlan2_PhyloseqObj.RData")
taxa_names(dat) <- gsub("s__", "", taxa_names(dat))

# Implement new FUNC:
# A <- dat %>% transform("compositional") %>%
#   abundances()
# asin(sqrt(A)) %>% data.frame() %>%
#   distribution_sanity()


############# data prep #############
# PD v PC
dat_pdpc <- subset_samples(dat, donor_group != "HC")
abun.pdpc <- as.data.frame.matrix(abundances(dat_pdpc))
# PD v HC PAIRED
dat_pdhc <- subset_samples(dat, Paired != "No")
abun.pdhc <- as.data.frame.matrix(abundances(dat_pdhc))


Maas.pd.pc <- read_tsv("maaslin2-export2/taxa_PDvPC_maaslin2_output/all_results.tsv", col_names = T) %>% filter(value == "Population Control")
Maas.pd.pc$feature <- gsub("s__", "", Maas.pd.pc$feature)
Maas.pd.pc.sig <- Maas.pd.pc %>% filter(qval < 0.25)

Maas.pd.hc <- read_tsv("maaslin2-export2/taxa_PDvHC_maaslin2_output/all_results.tsv", col_names = T) %>% filter(value == "Household Control")
Maas.pd.hc$feature <- gsub("s__", "", Maas.pd.hc$feature)
Maas.pd.hc.sig <- Maas.pd.hc %>% filter(qval < 0.25)


# Pull Species Name / P-value /Q-value (from Maaslin files) - Join with (Species Name & Prevalence from phyloseq Obj)

m.pc <- select(Maas.pd.pc, c("feature", "coef", "pval", "qval"))
colnames(m.pc) <- c("feature", paste0(c("coef", "pval", "qval"), ".PC"))

m.hc <- select(Maas.pd.hc, c("feature", "coef", "pval", "qval"))
colnames(m.hc) <- c("feature", paste0(c("coef", "pval", "qval"), ".HC"))
## Join Maaslin Tables
mas.QC <- left_join(m.pc, m.hc, by = "feature")
## Join prevelance col to maaslin data
prev <- data.frame(prevalence = prevalence(dat)) %>% rownames_to_column()
colnames(prev)[1] <- "feature"
mas.QC <- left_join(mas.QC, prev)

############### Create variables specifying unique or shared p & q value significance  ###############
############### P-values
mas.QC <- mutate(mas.QC, pcol = if_else(pval.PC <= 0.05 & pval.HC > 0.05, "PC_sig",
  if_else(pval.HC <= 0.05 & pval.PC > 0.05, "HC_sig",
    if_else(pval.PC <= 0.05 & pval.HC <= 0.05, "Both_Sig",
      if_else(pval.PC > 0.05 & pval.HC > 0.05, "Not_Sig",
        "error"
      )
    )
  )
))

############### Q-values
mas.QC <- mutate(mas.QC, qcol = if_else(qval.PC <= 0.25 & qval.HC > 0.25, "PC_sig",
  if_else(qval.HC <= 0.25 & qval.PC > 0.25, "HC_sig",
    if_else(qval.PC <= 0.25 & qval.HC <= 0.25, "Both_Sig",
      if_else(qval.PC > 0.25 & qval.HC > 0.25, "Not_Sig",
        "error"
      )
    )
  )
))



# cols <- c()

############### Plots ###############

#### P-val x P-val
p1 <- ggplot(mas.QC, aes(x = pval.PC, y = pval.HC)) +
  geom_point(aes(color = pcol)) +
  theme_minimal() +
  geom_vline(xintercept = 0.05, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0.05, linetype = 2, alpha = 0.5)
p1
#### Q-val x Q-val
p2 <- ggplot(mas.QC, aes(x = qval.PC, y = qval.HC)) +
  geom_point(aes(color = qcol)) +
  theme_minimal() +
  geom_vline(xintercept = 0.25, linetype = 2, alpha = 0.5) +
  geom_hline(yintercept = 0.25, linetype = 2, alpha = 0.5)
p2
####### Filter significant features in any comp
mas.QC.sig <- filter(mas.QC, pcol != "Not_Sig")

p3 <- ggplot(mas.QC.sig, aes(x = pval.PC, y = prevalence)) +
  geom_point(aes(color = pcol), size = 2, alpha = 0.8) +
  theme_minimal()
p3

p4 <- ggplot(mas.QC.sig, aes(x = pval.HC, y = prevalence)) +
  geom_point(aes(color = pcol), size = 2, alpha = 0.8) +
  theme_minimal()
p4

#### Cowplot Merge

merge1 <- plot_grid(p1, p2, align = c("h", "v"), nrow = 2)
ggsave(merge1,
  filename = "MaaslinQC_pqvalue_Plots.png",
  width = 14, height = 12
)

merge3 <- plot_grid(p3, p4, align = c("h"), nrow = 1)
ggsave(merge3,
  filename = "MaaslinQC_prevalence-pvalue_Plots.png",
  width = 14, height = 12
)



#### Coef x Coef
p5 <- ggplot(mas.QC, aes(x = coef.PC, y = coef.HC)) +
  geom_point(aes(color = pcol), size = 2, alpha = 0.8) +
  theme_minimal() +
  geom_text_repel(data = subset(mas.QC.sig, pcol == "Both_Sig"), aes(label = feature), force = 10)

#### P-val x Coef - PC
p6 <- ggplot(mas.QC, aes(x = pval.PC, y = coef.PC)) +
  geom_point(aes(color = pcol), size = 2, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none")

#### P-val x Coef - HC
p7 <- ggplot(mas.QC, aes(x = pval.HC, y = coef.HC)) +
  geom_point(aes(color = pcol), size = 2, alpha = 0.8) +
  theme_minimal() +
  theme(legend.position = "none")

merge2.bottom <- plot_grid(p6, p7, align = c("h", "v"), nrow = 1)
merge2 <- plot_grid(p5, merge2.bottom, align = c("h", "v"), ncol = 1)
merge2

ggsave(merge2,
  filename = "MaaslinQC_Coef-pvalue_Plots.png",
  width = 12, height = 12
)



# PD
dat_pd <- subset_samples(dat, donor_group == "PD")
# HC
dat_hc <- subset_samples(dat, donor_group == "HC")
# PC
dat_pc <- subset_samples(dat, donor_group == "PC")

# PrevPD <- data.fraprevalence(dat_pd)
# PrevHC <- prevalence(dat_hc)
# PrevPC <- prevalence(dat_pc)
prev.dat <- data.frame(PD.prev = prevalence(dat_pd), HC.prev = prevalence(dat_hc), PC.prev = prevalence(dat_pc))




PD.prev.under.1counts <- sum(prevalence(dat_pd) < 0.1)
HC.prev.under.1counts <- sum(prevalence(dat_hc) < 0.1)
PC.prev.under.1counts <- sum(prevalence(dat_pc) < 0.1)

# PD.prev.under.1 <- sum(filter(prev.dat, PD.prev < 0.1)$PD.prev)
# HC.prev.under.1 <- sum(filter(prev.dat, HC.prev < 0.1)$HC.prev)
# PC.prev.under.1 <- sum(filter(prev.dat, PC.prev < 0.1)$PC.prev)

grep(PrevPD, prevalence(dat_pd) < 0.1)

prev.dat <- data.frame(PD.prev = prevalence(dat_pd), HC.prev = prevalence(dat_hc), PC.prev = prevalence(dat_pc))


prev.datm <- melt(prev.dat)
ggplot() +
  geom_bar(data = prev.datm, aes(x = value)) +
  theme_minimal() +
  facet_wrap(~variable, ncol = 1, scales = "free") +
  geom_vline(xintercept = 0.1, linetype = "solid", color = "red") +
  ggtitle("Prevalence of species distribution by Group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "prevalence")

# Code Junkyard.


# pal_npg(palette = ("nrc"), alpha = 1)
# show_col(pal_npg("nrc")(4))
#
# pal_simpsons(palette = c("springfield"), alpha = 1)
# show_col(pal_simpsons("springfield")(4))
# c("#FED439FF", "#709AE1FF", "#8A9197FF", "#D2AF81FF")
