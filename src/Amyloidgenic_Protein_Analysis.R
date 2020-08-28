# Amyloidgenic Protein Analysis

############  load packages  ############

# rm(list = ls())

source("src/load_packages.R")
source("src/load_phyloseq_obj.R")
source("src/miscellaneous_funcs.R")
source("src/Metadata_prep_funcs.R")
source("src/Stats_funcs.R")



kos <- abundances(dat.KOs.slim) %>% as.data.frame() %>% rownames_to_column()

########### KOs ########### 
csga.df_a <- dat.KOs.slim %>% 
  microbiome::transform("compositional") %>% 
  abundances() %>% t() %>% 
  as.data.frame() %>% 
  dplyr::select(K04334) %>% 
  rownames_to_column("donor_id")

csga.df_a$K04334 <- asin(sqrt(csga.df_a$K04334))

csga.df <- group_col_from_ids(csga.df_a, csga.df_a$donor_id)
csga.df$group <- factor(csga.df$group, levels = c("HC", "PD","PC") )

########### Eggnogs - COGs ########### 
# eggnogs.df_a <- dat.Eggnogs.slim %>% 
#   abundances() %>% t() %>% 
#   as.data.frame() %>% 
#   dplyr::select(COG4625) %>% 
#   rownames_to_column("donor_id")
# eggnogs.df <- group_col_from_ids(eggnogs.df_a, eggnogs.df_a$donor_id)
# eggnogs.df$group <- factor(eggnogs.df$group, levels = c("HC", "PD","PC") )

# Basic Three Group plotting function
# cols <- c("PD"= "#bfbfbf", "PC" = "#ed7d31", "HC" = "#5b9bd5")
# title <- "CsgA - K04334"
# ylabel <- "Normalized Abundance"
# 
# 
# boxplot_all(csga.df, csga.df$group, csga.df$K04334, cols, title, ylabel)
# title <- "Pathogenicity - COG4625"
# boxplot_all(eggnogs.df, eggnogs.df$group, eggnogs.df$COG4625, cols, title, ylabel)



####################### csgA qPCR - SMS Analysis ####################### 

qpcr_dat <- 
  read.csv(file = "files/CsgA_qPCR.csv", 
           row.names = 1,  header= TRUE) %>% 
  rownames_to_column(var = "donor_id") 
qpcr_dat$donor_id <- gsub("_", ".", qpcr_dat$donor_id)

amyloids.ko <- left_join(csga.df, qpcr_dat, by = "donor_id") %>% 
  mutate(SMS_CsgA_Detected = if_else(K04334 > 0, "Positive", 
                                    "Negative" ))

qpcr_prev <- amyloids.ko %>% 
  dplyr::group_by(group, CsgA_Detected) %>% 
  tally()

HC.qpcr.pos <- filter(qpcr_prev, group == "HC" & CsgA_Detected == "Positive")$n
HC.qpcr.n <-sum(filter(qpcr_prev, group == "HC" & CsgA_Detected != "NA")$n)
PC.qpcr.pos <- filter(qpcr_prev, group == "PC" & CsgA_Detected == "Positive")$n
PC.qpcr.n <- sum(filter(qpcr_prev, group == "PC" & CsgA_Detected != "NA")$n)
PD.qpcr.pos <- filter(qpcr_prev, group == "PD" & CsgA_Detected == "Positive")$n
PD.qpcr.n <- sum(filter(qpcr_prev, group == "PD" & CsgA_Detected != "NA")$n)

qpcr_prev.HC <- HC.qpcr.pos/HC.qpcr.n
qpcr_prev.PC <- PC.qpcr.pos/PC.qpcr.n
qpcr_prev.PD <- PD.qpcr.pos/PD.qpcr.n

sms_prev <- amyloids.ko %>% 
  dplyr::group_by(group, SMS_CsgA_Detected) %>% 
  tally()

HC.sms.pos <- filter(sms_prev, group == "HC" & SMS_CsgA_Detected == "Positive")$n
HC.sms.n <-sum(filter(sms_prev, group == "HC" & SMS_CsgA_Detected != "NA")$n)
PC.sms.pos <- filter(sms_prev, group == "PC" & SMS_CsgA_Detected == "Positive")$n
PC.sms.n <- sum(filter(sms_prev, group == "PC" & SMS_CsgA_Detected != "NA")$n)
PD.sms.pos <- filter(sms_prev, group == "PD" & SMS_CsgA_Detected == "Positive")$n
PD.sms.n <- sum(filter(sms_prev, group == "PD" & SMS_CsgA_Detected != "NA")$n)

sms_prev.HC <- HC.sms.pos/HC.sms.n
sms_prev.PC <- PC.sms.pos/PC.sms.n
sms_prev.PD <- PD.sms.pos/PD.sms.n

csga_prev <- data.frame("group" = c("HC", "PC", "PD"), 
           "sms_prev" = c(sms_prev.HC, sms_prev.PC, sms_prev.PD),
           "qpcr_prev" = c(qpcr_prev.HC, qpcr_prev.PC, qpcr_prev.PD))

################ Stats ################ 
process_meta(dat)
amyloids.ko <- left_join(amyloids.ko, env, by = "donor_id")
shapiro.test(amyloids.ko$K04334) 
# distribution is non-normal 

###### PD vs HC tests
## Two-Proportions Z-Test : qPCR data
prop.test(x = c(HC.qpcr.pos, PD.qpcr.pos), 
          n = c(HC.qpcr.n, PD.qpcr.n),
          p = NULL, alternative = "two.sided",
          correct = TRUE)
## Two-Proportions Z-Test : SMS data
prop.test(x = c(HC.sms.pos, PD.sms.pos), 
          n = c(HC.sms.n, PD.sms.n),
          p = NULL, alternative = "two.sided",
          correct = TRUE)


######  PD vs PC tests
## Two-Proportions Z-Test : qPCR data
prop.test(x = c(PC.qpcr.pos, PD.qpcr.pos), 
          n = c(PC.qpcr.n, PD.qpcr.n),
          p = NULL, alternative = "two.sided",
          correct = TRUE)
## Two-Proportions Z-Test : SMS data
prop.test(x = c(PC.sms.pos, PD.sms.pos), 
          n = c(PC.sms.n, PD.sms.n),
          p = NULL, alternative = "two.sided",
          correct = TRUE)

## No significant Diferences in proportions of CsgA in either qPCR or SMS data


# Maas.pd.pc.csgA <- read_tsv(paste0("data/MaAsLin2_Analysis/KOs.slim_PDvPC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
#   filter(value == "Population Control") %>% filter(feature == "K04334")
# 
# Maas.pd.hc.csgA <- read_tsv(paste0("data/MaAsLin2_Analysis/KOs.slim_PDvHC_maaslin2_output/all_results.tsv"), col_names = T) %>% 
#   filter(value == "Household Control") %>% filter(feature == "K04334")



################ plot CsgA ################ 
#-------------------------------
# PLOT STATS
#-------------------------------
process_meta(dat)
csga.stats <- env
csga.stats$K04334 <- amyloids.ko$K04334
csga.stats$Paired.plot <- as.numeric(levels(csga.stats$Paired))[csga.stats$Paired]

observed.PdPC <- lm.PdPc(metadf=csga.stats, metric="K04334")
observed.PdHC <- lmm.PdHc(metadf=csga.stats, metric="K04334")
## Pull p-values
observed.PdPC.pval <- summary(observed.PdPC)$coefficients["descriptionPopulation Control","Pr(>|t|)"]
observed.PdHC.pval <- summary(observed.PdHC)$tTable["descriptionHousehold Control","p-value"]

cols <- c("PC"= "#bfbfbf", "PD" = "#ed7d31", "HC" = "#5b9bd5")
cols.qpcr =  c("Positive"= "#8b0000", "Negative" = "#838383", "NA" = "#ffffff")
amyloids.ko$group <- factor(amyloids.ko$group, levels = c("PC", "PD","HC") )

set.seed(42)
p1 <- ggplot(data=amyloids.ko, aes(x=group, y=K04334)) +
  geom_boxplot(aes(color = group), outlier.alpha = 0) +
  geom_point(aes(fill = CsgA_Detected), position = position_jitterdodge(jitter.width = 0.3),
             shape=21, size=3, alpha = 0.8) +
  geom_text(data=csga_prev, 
            aes(x=group, y=0.007, 
                label=paste0(round(sms_prev*100, 1), "%, ", round(qpcr_prev*100, 1), "%") ), 
            size = 4, check_overlap = F) +
  geom_signif(comparisons = list(c("PD", "HC")), tip_length = 0.02,
              annotations = sig_mapper(observed.PdHC.pval, porq = "p", symbols = F)) +
  geom_signif(comparisons = list(c("PC", "PD")), tip_length = 0.02,
              annotations = sig_mapper(observed.PdPC.pval, porq = "p", symbols = F)) +
  theme_classic() +
  ggtitle("csgA (K04334)") +
  ylab("WGS - Normalized Abundance") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
  scale_color_manual(values = cols, name ="Group") +
  scale_fill_manual(values = cols.qpcr, name ="qPCR csgA \n detected") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.line = element_line(colour="grey"),
        panel.grid = element_blank())

BTp <- filter(amyloids.ko, CsgA_Detected == "Positive")

p2 <- ggplot(BTp, aes(y=(K04334), x=2^(32.44-CT_Avg))) +
  geom_smooth(method="lm", se=F, colour = "grey", linetype = "dashed") +
  geom_point(aes(fill=group), shape=21, size=3, alpha =0.8) +
  ylab("WGS - Normalized Abundance") +
  xlab(expression(paste(2^"(32.44 - CT Average)"))) +
  theme_bw() +
  scale_fill_manual(values =cols) +
  scale_color_manual(values = cols) +
  ggtitle("csgA (K04334): SMS vs qPCR \n qPCR positive samples") +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 1250, label.y = 0.0025) +
  # stat_regline_equation(label.x = 1250, label.y = 0.0025) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_line(colour="grey"),
        panel.grid = element_blank())

p1.a <- p1 + theme(legend.position = "none")
p2.a <- p2 + theme(legend.position = "none")
c1 <- cowplot::plot_grid(p1.a, p2.a, ncol = 2, align="hv", labels = c("C", "D"))
c2 <- cowplot::plot_grid(c1, get_legend(p1), ncol = 2, 
                        rel_widths = c(6, 1), axis = "l", align="h")
c2

ggsave(c2, filename = "data/Amyloidgenic_Protein_Analysis/SMSvQPCR_data.svg", height = 4, width =12)



###################     Boxplots for all bacterial amyloids   #######################




