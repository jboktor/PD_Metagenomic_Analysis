##  Meta-Analysis Confounder Analysis ##  

# Packages
library("tidyverse"); library("coin"); library(ggplot2); library(ggrepel); library(microbiome)


### NOTES/UPDATES REQUIRED: 
# STATS: color by features significant in PDvPC &/or PDvHC comparison rather than using single maaslin2 model
# GRAPHICS: nest count table in each plot for n# of PD, PC, and HC available for metdata covariate 
  # Also add R2 for each covariate beside name


rm(list = ls())
######## Load Data & functions
source("src/load_phyloseq_obj.R")
dat_lg <- microbiome::transform(dat, "log10")


## Read in table of significant microbes - from MaAsLin2(Disease glm (no age/bmi/sex) covars)
Maas.test2 <- read_tsv("files/maaslin2-export_OLD/taxa_allsamples_PDeffect_maaslin2_output/all_results.tsv", col_names = T)
Maas.test2 <- Maas.test2[order(Maas.test2$feature),] # order df alphabetically by species
adj.p.val <- tibble(Maas.test2$qval)

# Set variables for Analysis
feat.all.dfm <- as.data.frame.matrix(abundances(dat_lg))
feat.all.df <- rownames_to_column(feat.all.dfm)
feat.all.df <- feat.all.df[order(feat.all.df$rowname),]
feat.all <- filter(feat.all.df, rowname %in% Maas.test2$feature) %>% column_to_rownames(var="rowname")

##############################################################################

# Readin PERMANOVA Analysis to pass significant metdata features;
#contains metadata pre-processing as well
permdf <- read.csv(file = 'files/permanova_data.csv') %>% 
  filter(metacat != "AA_notmatched")
meta <- env
alpha.meta <- 5e-02


###############################     Disease Status Variance   ###############################################
#  variance explained by disease status
ss.disease <- apply(feat.all, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=meta %>% pull(PD))

# calculate mean abundance - not trimmed mean
t.mean <- apply(feat.all, 1, mean) #, trim=0.1)

df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val = adj.p.val$`Maas.test2$qval`,
  meta.significance= adj.p.val < 0.05)



################################# Test possible Confounders - Essential/Foundational METADATA #############################################

permanova_essential <- permdfess$vars[-grep("PD", permdfsig$vars)]
df.list <- list()

for (meta.var in permanova_essential){
  cat('###############################', meta.var, '###############################\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$PD, meta.c %>% pull(meta.var)))
  feat.red <- feat.all[,meta.c$donor_id]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
    print(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    print(meta.c %>% pull(meta.var) %>% unique)
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      print(meta.c %>% pull(meta.var) %>% unique)
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}


##################### Facet Wrapped - Select Metadata Plot  ##############################################

g1 <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when(type=='sex' ~ 'Sex',
                         type=='host_age_factor' ~ 'Age',
                         type=='host_body_mass_index' ~ 'BMI',
                         type=='bristol_stool_scale' ~ 'BSS',
                         type=='description' ~ 'Donor Group',
                         type=='Paired' ~ 'Household Effect',
                         type=='antibiotic_use' ~ 'Antibiotics',
                         type=='bowl_movments_per_day' ~ 'Bowel Movement',
                         type=='laxatives' ~ 'Laxatives',
                         TRUE ~ type)) %>%
  ggplot(aes(x=disease, y=meta, size=t.mean+1e-08, col=meta.significance, alpha = meta.significance)) +
  geom_point(shape=19) +
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=3, scales = "free") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from=0, to=1, by=0.1)) +
  scale_colour_manual(values = c('black', '#CC071E'),
                      name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
  scale_alpha_manual(values=c("FALSE" = 0.1, "TRUE" = 0.75)) + 
  scale_size_area(name=expression('log'[10]*'(mean abundance)'),
                  breaks=c(0, 0.25, 0.5, 0.75))+
  guides( size = "legend", colour='legend')
g1 

# ggsave(g1,
#        filename = 'data/Confounder_Analysis/confounder_plot_KEY_metadata.svg',
#        width = 8,
#        height = 8)



################################# Test possible Confounders - PERMANOVA Significant METADATA #############################################

################################# Refresh df.plot.all
df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val = adj.p.val$`Maas.test2$qval`,
  meta.significance= adj.p.val < 0.05)
################################# 
df.list <- list()

permanova_significant <- permdfsig$vars[-grep("PD", permdfsig$vars)]

for (meta.var in permanova_significant ){
  cat('###############################', meta.var, '###############################\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$PD, meta.c %>% pull(meta.var)))
  feat.red <- feat.all[,meta.c$donor_id]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
    print(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    print(meta.c %>% pull(meta.var) %>% unique)
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      print(meta.c %>% pull(meta.var) %>% unique)
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}

##################### Facet Wrapped - PERMANOVA Significant Metadata Plot  ##############################################

g2 <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when(type=='PD' ~ 'Disease Status',
                         type=='Paired' ~ 'Household Effect',
                         type=='description' ~ 'Donor Group',
                         type=='bristol_stool_scale' ~ 'BSS',
                         type=='antibiotic_use' ~ 'Antibiotics',
                         type=='laxatives' ~ 'Laxatives',
                         type=='ssri_antidepressants' ~ 'SSRI',
                         TRUE ~ type)) %>%
  ggplot(aes(x=disease, y=meta, size=t.mean+1e-08, col=meta.significance, alpha = meta.significance)) +
  geom_point(shape=19) +
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=3, scales = "free") +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from=0, to=1, by=0.1)) +
  # scale_y_continuous(breaks=seq(from=0, to=1, by=0.2)) +
  scale_colour_manual(values = c('black', '#CC071E'),
                      name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
  scale_alpha_manual(values=c("FALSE" = 0.1, "TRUE" = 0.75)) + 
  scale_size_area(name=expression('log'[10]*'(mean abundance)'),
                  breaks=c(0, 0.25, 0.5, 0.75))+
  guides( size = "legend", colour='legend')
g2 

# ggsave(g2,
#        filename = 'data/Confounder_Analysis/confounder_plot_PERMANOVA_Signficant_metadata.svg',
#        width = 8,
#        height = 8)





################################# Test possible Confounders  - ALL METADATA #############################################
################################# Refresh df.plot.all
df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean,
  adj.p.val = adj.p.val$`Maas.test2$qval`,
  meta.significance= adj.p.val < 0.05)
################################# 

df.list <- list()

for (meta.var in colnames(meta)){
  cat('###############################', meta.var, '###############################\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$PD, meta.c %>% pull(meta.var)))
  feat.red <- feat.all[,meta.c$donor_id]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
    print(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    print(meta.c %>% pull(meta.var) %>% unique)
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      print(meta.c %>% pull(meta.var) %>% unique)
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}


##################### Facet Wrapped - All Metadata Plot  ##############################################

g3 <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  mutate(facet=case_when(type=='sex' ~ 'Sex',tag,' ',
                         type=='host_age_factor' ~ 'Age',
                         type=='host_body_mass_index' ~ 'BMI',
                         type=='bristol_stool_scale' ~ 'BSS',
                         type=='description' ~ 'Donor Group',
                         type=='Paired' ~ 'Household Effect',
                         type=='antibiotic_use' ~ 'Antibiotics',
                         type=='laxatives' ~ 'Laxatives',
                         type=='ssri_antidepressants' ~ 'SSRI',
                         TRUE ~ type)) %>%
  ggplot(aes(x=disease, y=meta, size=t.mean+1e-08, col=meta.significance, alpha = meta.significance)) +
  geom_point(shape=19) +
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  facet_wrap(~facet, ncol=8) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = seq(from=0, to=1, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0, to=1, by=0.2)) +
  scale_colour_manual(values = c('black', '#CC071E'),
                      name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
  scale_alpha_manual(values=c("FALSE" = 0.1, "TRUE" = 0.75)) + 
  scale_size_area(name=expression('log'[10]*'(mean abundance)'),
                  breaks=c(0, 0.25, 0.5, 0.75)) +
  guides( size = "legend", colour='legend')
g3

# ggsave(g3, 
#        filename = 'data/Confounder_Analysis/confounder_plot_all_metadata.svg', 
#        width = 15, 
#        height = 15)








##################################  Single Confounder vs Disease Variation Plots ############################################

################# INPUT Any confounder of interest below  #################

df.plot.study <- df.plot.all %>%
  gather(key=type, value=meta, -species, -disease,
         -t.mean, -adj.p.val, -meta.significance) %>%
  filter(!str_detect(type, '.significance')) %>%
  filter(complete.cases(.)) %>%
  filter(type=='donor_group')   ################# Input HERE:

alpha.meta <- 5e-02
g2 <- df.plot.study %>%
  ggplot(aes(x=disease, y=meta)) +
  geom_point(aes(size= (t.mean), fill=meta.significance, alpha = meta.significance), shape=21) +
  xlab(paste0('Variance explained by Disease Status\n average: ',
              formatC(mean(df.plot.study$disease)*100, digits=2), '%')) +
  ylab(paste0('Variance explained by Household Effect \n average: ',
              formatC(mean(df.plot.study$meta)*100, digits=2), '%')) +
  theme_bw() +
  theme(panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c('grey', '#CC071E'),
                      name=paste0('Significance\n(', alpha.meta, ' FDR)')) +
  scale_alpha_manual(values=c("FALSE" = 0.4, "TRUE" = 0.8)) +
  scale_size_area(name=expression('log'[10]*'(mean abundance)'),
                  breaks=c(0, 0.25, 0.5, 0.75)) +
  geom_text_repel(data=subset(df.plot.study, meta.significance), 
                  aes(label=species), direction = "x", angle = 90,
                  nudge_y = 0.01, vjust = 0,force = 3, size = 3)+
  guides( size = "legend", colour='legend')
g2


# ggsave(g2, 
#        filename = 'data/Confounder_Analysis/confounder_plot_DonorGroup_vs_PDstatus_metadata.svg', 
#        width = 7, 
#        height = 8)





#############  CODE JUNKYARD   #############
# 
# df.plot.study <- na.omit(select(df.plot.all, c("species", "disease", "t.mean", "adj.p.val", "meta.significance", "description")))
# df.plot.study.sig <- filter(df.plot.study, adj.p.val < 0.05)
# df.plot.study.sig$species <- gsub("s__", "", df.plot.study.sig$species)
# # df.plot.study.sig$species <- filter(df.plot.study, meta.significance)
# 
# g2 <- df.plot.study %>%
#   ggplot(aes(x=disease, y=description)) +
#   geom_point(aes(size= (t.mean), fill=meta.significance), shape=21,
#              col=alpha(c('black'), alpha=0.4)) +
#   xlab(paste0('Variance explained by Disease Status\n average: ',
#               formatC(mean(df.plot.study$disease)*100, digits=2), '%')) +
#   ylab(paste0('Variance explained by Donor Group \n average: ',
#               formatC(mean(df.plot.study$description)*100, digits=2), '%')) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank()) +
#   scale_x_continuous(breaks = seq(from=0, to=0.2, by=0.1)) +
#   scale_y_continuous(breaks=seq(from=0, to=0.6, by=0.1)) +
#   scale_fill_manual(values = alpha(c('grey', '#CC071E'),
#                                    alpha=c(0.4, .8)),
#                     name=paste0('FDR < ', alpha.meta)) +
#   # scale_size_area(name=expression('log'[10]*'(abundance)'),
#   #                 breaks=c(0, 0.25, 0.5, 0.75)) +
#   geom_text_repel(data=df.plot.study.sig, aes(label=species), direction = "x", angle = 90,
#                   nudge_y = 0.01, vjust = 0,force = 3, size = 3)+
#   guides( size = "legend", colour='legend')
# g2
# 
# g2a <- g2 + theme(legend.position = c(.80, .25))

# ggsave(g2a,
#        filename = paste0('confounder_plot_DiseaseStatus_by_DonorGroup.svg'),
#        width = (mean(df.plot.study$disease)*100*2), 
#        height = (mean(df.plot.study$description)*100*2))



# cat('Successfully computed confounder effects in',
#     proc.time()[1]-start.time, 'second...\n')

# #######################
# End of script
# #######################
