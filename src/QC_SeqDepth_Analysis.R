### QC_Seq_Depth Analysis

library(ggplot2); library(tidyverse); library(readxl);library(dplyr); library(ggrepel);library(grid);
library(gridExtra);library(reshape2);library(plyr);library(grid);library(devtools);library(RColorBrewer);
library(ggfortify);library(vegan);library(MASS);library(compositions);library(zCompositions);library(phyloseq);
library(Biobase);library(viridis);library("foreach");library("doParallel");library(ggbeeswarm);
library(FSA);library(ggpubr);library(ggsci);library(microbiome);library(ggridges);library(future);library(cowplot)


# Load Phyloseq Obj
load("Species_PhyloseqObj.RData")
dat <- transform(dat, "compositional")

func_reads <- read_tsv("files/humann2_read_and_species_count_table.tsv", col_names = T)
reads <- dplyr::select(func_reads, c("# samples","total reads"))
colnames(reads) <- c("id", "clean_total_reads")


########### quick spin-off : get average clean reads for each group ########### 

reads2 <- reads %>% tibble() %>%  mutate(group = if_else(grepl("HC", reads$id), "HC",
                                            if_else(grepl("PC", reads$id), "PC",
                                                    if_else(grepl("PD", reads$id), "PD",
                                                            "error"))))

reads2 %>%
  group_by(group) %>%
  summarize(mean(clean_total_reads, na.rm = TRUE))



################################################################## 

env <- meta(dat)
env <- left_join(env, reads, by = "id")
env <- mutate(env, clean_total_reads_factor= as.factor(cut(env$clean_total_reads, breaks = quantile(env$clean_total_reads), labels=c(1,2,3,4), 
                                                      include.lowest=T) ))

histogram(env$clean_total_reads)


p <- ggplot(env, aes(x = clean_total_reads, y = description)) +
  geom_density_ridges(aes(color = description, fill = description),
                      jittered_points = TRUE, quantile_lines = TRUE,
                      position = position_raincloud(adjust_vlines = TRUE),
                      alpha = 0.1, scale = 10) + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  scale_color_manual(values = c("Household Control" = "#440154FF", 
                                "PD Patient" = "#d48a02", 
                                "Population Control" = "#21908CFF")) +
  scale_fill_manual(values = c("Household Control" = "#440154FF", 
                               "PD Patient" = "#FDE725FF", 
                               "Population Control" = "#21908CFF")) +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        # axis.text.y=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x = element_blank(),
        plot.title=element_blank(),
        legend.position = "none")
p + coord_flip()




# Make Alpha Div specific obj so binary OTUs aren't accidentally used elsewhere 
dat_alpha <- dat
### WARNING : OTU TABLE CONVERTED TO BINARY - DO NOT USE DOWNSTREAM FOR other analyses
otu_table(dat_alpha) <- ((otu_table(dat_alpha) > 0) + 0)
# Removes any taxa present in no samples 
dat_alpha <- prune_taxa(taxa_sums(dat_alpha) > 0, dat_alpha) # none pruned
## Calculate Alpha Diversity Metrics
tab <- global(dat_alpha, index = "all")


env$observed <- tab$observed
env$shannon <- tab$diversities_shannon
env$evenness_simpson <- tab$evenness_simpson

# Test for normality - Distb is non-normal
shapiro.test(env$observed)
shapiro.test(env$shannon)
shapiro.test(env$evenness_simpson)


############### Observed species ############### 

# ANOVA & Tukey Post-Hoc
anovaResults <- aov(observed~clean_total_reads_factor, data=env)
summary(anovaResults)
sigv <- TukeyHSD(anovaResults)
sigv

o1 <- ggplot(env, aes(x=(clean_total_reads), y=observed)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group, color = donor_group),shape=21, size=2, alpha = 0.9) +
  xlab("Filtered Sample Reads") +
  ylab("Observed Species") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 15e6, label.y = c(90, 85, 80)) +
  stat_regline_equation(label.x = 15e6, label.y = c(88, 83, 78)) +
  ggtitle("Species Detected by Read Depth - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
o1

o2 <- ggplot(env, aes(x=(clean_total_reads), y=observed, color = donor_group)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group),shape=21, size=2, alpha = 0.9) +
  xlab("Filtered Sample Reads") +
  ylab("Observed Species") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color = donor_group), label.x = 15e6, label.y = c(90, 85, 80)) +
  stat_regline_equation(label.x = 15e6, label.y = c(88, 83, 78)) +
  ggtitle("Species Detected by Read Depth (Group Specific) - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
o2

o3 <- ggplot(env, aes(x=(clean_total_reads_factor), y=observed, color = donor_group)) +
  geom_boxplot() +
  xlab("Filtered Sample Reads") +
  ylab("Observed Species") +
  theme_minimal() +
  # ggtitle("Total csga Relative Abundance \nLinear Regression qPCR vs Metagenomics \ncsga-Positive samples") +
  theme(plot.title = element_text(hjust = 0.5))
o3


############### Shannon Diversity Index ############### 

# ANOVA & Tukey Post-Hoc
anovaResults <- aov(shannon~clean_total_reads_factor, data=env)
summary(anovaResults)
sigv <- TukeyHSD(anovaResults)
sigv


s1 <- ggplot(env, aes(x=(clean_total_reads), y=shannon)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group, color = donor_group),shape=21, size=2, alpha = 0.9) +
  ylab("Shannon Index") +
  xlab("Filtered Sample Reads") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 15e6, label.y = c(3.8)) +
  stat_regline_equation(label.x = 15e6, label.y = c(3.75)) +
  ggtitle("Shannon Index by Read Depth - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
s1

s2 <- ggplot(env, aes(x=(clean_total_reads), y=shannon, color = donor_group)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group),shape=21, size=2, alpha = 0.9) +
  ylab("Shannon Index") +
  xlab("Filtered Sample Reads") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color = donor_group), label.x = 15e6, label.y = c(3.95, 3.88, 3.81)) +
  stat_regline_equation(label.x = 15e6, label.y = c(3.92, 3.85, 3.78)) +
  ggtitle("Shannon Index by Read Depth (Group Specific) - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
s2

s3 <- ggplot(env, aes(x=(clean_total_reads_factor), y=shannon, color = donor_group)) +
  geom_boxplot() +
  ylab("Shannon Index") +
  xlab("Filtered Sample Read Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
s3



############### Simpsons Evenness Index ############### 

# ANOVA & Tukey Post-Hoc
anovaResults <- aov(evenness_simpson~clean_total_reads_factor, data=env)
summary(anovaResults)
sigv <- TukeyHSD(anovaResults)
sigv

e1 <- ggplot(env, aes(x=(clean_total_reads), y=evenness_simpson)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group, color = donor_group),shape=21, size=2, alpha = 0.9) +
  ylab("Simpsons Evenness") +
  xlab("Filtered Sample Reads") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 15e6, label.y = c(0.16)) +
  stat_regline_equation(label.x = 15e6, label.y = c(0.155)) +
  ggtitle("Simpsons Evenness by Read Depth - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
e1
  
e2 <-ggplot(env, aes(x=(clean_total_reads), y=evenness_simpson, color = donor_group)) +
  geom_smooth(method="lm", se=F) +
  geom_point(aes(fill=donor_group),shape=21, size=2, alpha = 0.9) +
  ylab("Simpsons Evenness") +
  xlab("Filtered Sample Reads") +
  theme_minimal() +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"), color = donor_group), label.x = 15e6, label.y = c(0.175, 0.16, 0.145)) +
  stat_regline_equation(label.x = 15e6, label.y = c(0.17, 0.155, 0.14)) +
  ggtitle("Simpsons Evenness by Read Depth (Group Specific) - Linear Regression") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
e2

e3 <-ggplot(env, aes(x=(clean_total_reads_factor), y=evenness_simpson, color = donor_group)) +
  geom_boxplot() +
  ylab("Simpsons Evenness") +
  xlab("Filtered Sample Read Quantiles") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
e3



QC_speciesVsreads <- cowplot::plot_grid(o1, o2, o3, 
                               s1, s2, s3,
                               e1, e2, e3, align = "h", nrow = 3)


# ggsave(QC_speciesVsreads, filename = "Read_Depth_vs_AlphaDiv.png",
#        width = 18, height = 18)




################## Density Plots by Quantile ################## 

z <-ggplot(env, aes(x=(clean_total_reads), color = donor_group)) +
  geom_density() +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
z
# ggsave(z, filename = "Clean_Read_Depth_DensityPlot_by_Group.png",
#        width = 7, height = 5)

############## Species Abundnace by Read Depth : Distributions between donor group
## Pull abundance data - merge with clean_total_reads_factor & donor_group -> melt df
seqdat <- select(env, c("donor_id", "clean_total_reads_factor", "donor_group"))
ab <- abundances(dat) %>% t() %>% as.data.frame.matrix() %>% rownames_to_column() %>% cbind(seqdat)
abm <- melt(ab)
# histogram(abm$value)
# histogram(asin(sqrt(abm$value)))
# histogram(abm$value)
# histogram(log10(abm$value + 1))
# histogram(log10(abm$value + impt))


# Find smallest value greater than 0 in abundnace df
abun <- abundances(dat)
impt <- min(abun[abun > 0])/2

sqr <- expression(sqrt(Abundance)*'by Read Depth Quantile')

z1 <-ggplot(abm, aes(x=asin(sqrt(value)), color = donor_group)) +
  geom_density(alpha = .7) +
  facet_wrap(~clean_total_reads_factor, nrow = 1) +
  labs(x = paste0("ArcSin(Sqrt(Abundance))")) +
  # labs(x = expression(sqrt(Abundance)*'by Read Depth Quantile')) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
z1


z2 <-ggplot(abm, aes(x=variable, y=asin(sqrt(value)), color = donor_group)) +
  geom_beeswarm(aes(color = donor_group, fill = donor_group), shape = 21, size=.3, alpha = .7, cex = .5) +
  facet_wrap(~clean_total_reads_factor, nrow = 1) +
  labs(x = "Species by Read Depth Quantile", 
       y = paste0("ArcSin(Sqrt(Abundance))")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 60, size=0.5, vjust = 0.3))
z2


QC_speciesAbundanceVsreads <- cowplot::plot_grid(z1, z2, align = "h", nrow = 2)


ggsave(QC_speciesAbundanceVsreads, filename = "Group_Species_Abundance_by_Read_Depth.png",
       width = 10, height = 8, dpi = 1200)



