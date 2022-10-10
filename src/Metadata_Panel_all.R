# Study Design Panel

library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(ggrepel)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)
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
library("foreach")
library("doParallel")
library(ggbeeswarm)
library(FSA)
library(ggpubr)
library(ggsci)
library(microbiome)
library(ggridges)
library(future)
library(cowplot)

rm(list = ls())

# Load Phlyoseq Obj
load("files/Species_PhyloseqObj.RData")
source("src/Metadata_prep_funcs.R")
source("src/miscellaneous_funcs.R")

########################### Metadata Pre-processing & Selection ###########################


ps <- readRDS("files/Phyloseq_Merged/PhyloseqObj_slim_clean.rds")
dat <- ps$Species

# Run Metadata pre-processing function
process_meta_study_design_plot(dat)
trim_meta(env, 3/nsamples(dat))

env$paired <- as.character(env$paired)
env[which(env$paired != "No"), "paired"] <- "Yes"
env$paired <- factor(env$paired)


### Manual Colors:

# Yes/No
yn.na.Q <- c("Yes" = "#bfd3e6", "No" = "#810f7c", "NA" = "#ffffff")
# Frequency
freq.Q <-
  setNames(
    colorRampPalette(brewer.pal(9, "BuPu"))(6),
    c(
      "Never", "Rarely", "Less than half the time", "Not in the past 7 days",
      "Almost every day", "Every day"
    )
  )
# Numeric Ranks
num.Q <-
  setNames(
    colorRampPalette(brewer.pal(9, "BuPu"))(7),
    0:6
  )

# Donor group vars
donorgroup.Q <-
  setNames(
    viridis(3, alpha = 0.8),
    c("HC", "PC", "PD")
  )
# Sex vars
sex.Q <-
  setNames(
    colorRampPalette(brewer.pal(3, "BuPu"))(2),
    c("male", "female")
  )

race.Q <- c("White" = "#bfd3e6")

# Remove resveratrol and creatine - previously filtered
env.design <- env %>% dplyr::select(-c(description, host_age))

env.design$donor_group4melt <- env.design$donor_group
env.design <- melt(env.design, id = c("donor_id", "donor_group4melt"))
env.design$donor_group4melt <- factor(env.design$donor_group4melt, levels = c("HC", "PD", "PC"))

env.design <- env.design[order(env.design$donor_group4melt), ]
env.design$order_col <- 1:nrow(env.design)



## Order Variable Col
y.axis.order.genres <- c(
  other_metadata, grouping, Anthro, GI, nonstarters,
  general_drugs, pd_drugs, supplements, allergy, smoke, diet
)

y.axis.order <- c()
for (i in seq_along(y.axis.order.genres)) {
  x <- y.axis.order.genres[i]
  x.filtered <- x[x %in% env.design$variable]
  print(x.filtered)
  y.axis.order <- c(y.axis.order, x.filtered)
}
env.design$variable <- factor(env.design$variable, levels = y.axis.order)
unique(env.design$value)

p1 <- ggplot(data = env.design, aes(x = reorder(donor_id, order_col), y = variable, fill = value)) +
  geom_tile() +
  theme_minimal() +
  labs(x = "Donors") +
  scale_fill_manual(values = c(yn.na.Q, freq.Q, num.Q, donorgroup.Q, sex.Q, race.Q)) +
  theme(
    # axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=5),
    axis.text.x = element_blank(),
    # axis.title.x=element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "none"
  )
p1



yn.na.Q_legend <-
  env.design %>%
  filter(value %in% names(yn.na.Q)) %>%
  ggplot(aes(x = donor_id, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c(yn.na.Q))

num.Q_legend <-
  env.design %>%
  filter(value %in% names(num.Q)) %>%
  ggplot(aes(x = donor_id, y = variable, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c(num.Q))

donorgroup.Q_legend <-
  env.design %>%
  filter(value %in% names(donorgroup.Q)) %>%
  ggplot(aes(x = donor_id, y = variable, fill = value)) +
  geom_tile() +
  labs(fill = "Donor Group") +
  scale_fill_manual(values = c(donorgroup.Q))

sex.Q_legend <-
  env.design %>%
  filter(value %in% names(sex.Q)) %>%
  ggplot(aes(x = donor_id, y = variable, fill = value)) +
  geom_tile() +
  labs(fill = "Sex") +
  scale_fill_manual(values = c(sex.Q))

env.design.freq <- env.design %>%
  filter(value %in% names(freq.Q))
env.design.freq$value <- factor(env.design.freq$value, levels = c(
  "Never", "Rarely", "Less than half the time",
  "Not in the past 7 days", "Almost every day", "Every day"
))
freq.Q_legend <-
  env.design.freq %>%
  ggplot(aes(x = donor_id, y = variable, fill = value)) +
  geom_tile() +
  labs(fill = "Dietary Consumption \n Frequency") +
  scale_fill_manual(values = c(freq.Q))



Legends <- plot_grid(get_legend(freq.Q_legend), NULL,
  get_legend(num.Q_legend), get_legend(yn.na.Q_legend),
  get_legend(sex.Q_legend), get_legend(donorgroup.Q_legend),
  ncol = 2, align = "v"
)
Legends

heatmap <- plot_grid(NULL, Legends, p1, ncol = 3, rel_widths = c(0.25, 1, 5))
# heatmap


####################   Load Permanova dataframe from files

permdf <- read.csv(file = "files/permanova_data.csv") %>%
  filter(metacat != "AA_notmatched")
# rename description to donor_group
permdf[permdf$vars == "description", "vars"] <- "donor_group"

## 1) remove features in permanova not shared with meta-analysis

permdf.metaplot <- filter(permdf, vars %in% env.design$variable)

## 2) Re-order y-axis based on heatmap groups

y.axis.order <- factor(y.axis.order)
permdf.metaplot$vars <- factor(permdf.metaplot$vars, levels = y.axis.order)

## 3) Plot Bars
#################### All Metadata - METADATA CATEGORY FILL

permdf.metaplot$FDR.symbol <- sig.symbol.generator(permdf.metaplot$FDR)

p2 <- ggplot() +
  geom_col(data = permdf.metaplot, aes(x = log10(R2), y = vars, fill = metacat)) +
  theme_minimal() +
  labs(x = expression("Variance Explained - log"[10] * "(" ~ R^2 ~ "%)"), fill = "Category") +
  # ggtitle("PERMANOVA Covariate Independence") +
  annotate("text", y = unique(permdf.metaplot$vars), x = (log10(permdf.metaplot$R2) + 0.2), label = permdf.metaplot$FDR.symbol, size = 2) +
  # annotate("text", y = unique(permdf.metaplot$vars), x = (log10(permdf.metaplot$R2)+0.2), label = paste0("n =", permdf.metaplot$n_meta), size = 2) +
  scale_fill_simpsons() +
  theme(
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(.80, .50)
  )
p2

## 4) Cowplot merge

p2a <- p2 + theme(axis.text.y = element_blank())
p3 <- plot_grid(p1, p2a, ncol = 2, align = "h", rel_widths = c(2, 1))


# fig1a_supp <- plot_grid(NULL, Legends, p3,
#                         ncol = 3, rel_widths = c(0.25, 1, 5))
fig1a_supp <- plot_grid(p3, NULL, Legends,
  ncol = 3, rel_widths = c(5, 0.25, 0.75)
)
fig1a_supp



# ggsave(fig1a_supp,
#   filename = "figures/metadata_heatmap_permanova_variance.svg",
#   width = 14, height = 6
# )
