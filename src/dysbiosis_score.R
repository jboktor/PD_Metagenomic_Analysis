# Dysbiosis score

####### Load PhyloSeq objects  ####### 
source("src/beta_diversity.R")
load("files/low_quality_samples.RData")


#--------------------------------------------------------------------------------
#                                  TBC
#--------------------------------------------------------------------------------
remove_dats()
load_tbc()

dys.tbc <- dysbiosis_score(dat.EGGNOGs.slim, dist = "bray")
# dys.tbc <- dysbiosis_score(dat.species, dist = "Aitchisons")

met <- meta(dat.species) %>%
  dplyr::select(donor_id, cohort)
dys.df.tbc <- left_join(dys.tbc, met, by = "donor_id")

#--------------------------------------------------------------------------------
#                                  RUSH
#--------------------------------------------------------------------------------

remove_dats()
load_rush()

dys.rush <- dysbiosis_score(dat.EGGNOGs.slim, dist = "bray")
# dys.rush <- dysbiosis_score(dat.species, dist = "Aitchisons")

met <- meta(dat.species) %>%
  dplyr::select(donor_id, cohort)
dys.df.rush <- left_join(dys.rush, met, by = "donor_id")


#--------------------------------------------------------------------------------
#                                  Faceted Combined 
#--------------------------------------------------------------------------------

dys.df <- dplyr::bind_rows(dys.df.tbc, dys.df.rush)

dysbiosis_thres <- 
  dys.df %>% 
  filter(group == "PC") %>%
  dplyr::group_by(cohort) %>% 
  dplyr::summarise(quant = quantile(median, na.rm = T, probs = 0.9))

dys_faceted <- 
  dys.df %>% 
  ggplot() +
  geom_density(aes(median, color = group, fill = group), alpha = 0.6) +
  facet_wrap(~cohort, nrow = 2) +
  geom_vline(data = dysbiosis_thres, aes(xintercept = quant), size = 1) +
  geom_rect(data = dysbiosis_thres, aes(xmin = quant, xmax = Inf, ymin = -Inf, ymax = Inf), alpha = 0.2) +
  labs(x = "Dysbiosis score") +
  scale_fill_manual(values = cols.pdpchc) +
  scale_colour_manual(values = cols.pdpchc.rim) +
  theme_bw() +
  theme(panel.grid = element_blank())
dys_faceted

# ggsave(dys_faceted, filename = "data/Community_Composition/Dysbiosis/Dysbiosis_score_Species_Bray-Curtis.svg",
#        width = 6, height = 5)



dysbiosis_thres.tbc <- filter(dysbiosis_thres, cohort == "TBC")
dysbiosis_thres.rush <- filter(dysbiosis_thres, cohort == "Rush")

dys.df.classified <- 
  dys.df %>%
  mutate(
    classification = if_else(
      cohort == "TBC" &
        median > dysbiosis_thres.tbc$quant,
      "dysbiotic",
      if_else(
        cohort == "Rush" &
          median > dysbiosis_thres.rush$quant,
        "dysbiotic",
        "normal"
      )
    )
  )

dysbiosis.stats <- 
  dys.df.classified %>%
  dplyr::group_by(group, cohort, classification) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(percentage = n/sum(n)*100) 

dysbiosis.barplot <- 
  dysbiosis.stats %>% 
  ggplot(aes(x = group, y = percentage, fill = classification)) + 
  geom_bar(stat = "identity", width = 0.5) +
  facet_wrap(~cohort) +
  theme_bw() +
  labs(y = "Percentage %") +
  geom_text(aes(label=n), size=4, color = "white",
            position = position_stack(vjust = 0.5),
            vjust = 0) +
  scale_fill_manual(values = 
                      c("dysbiotic" = "#800000",
                        "normal" = "#959595")) +
  theme(axis.title.x = element_blank(),
        panel.grid = element_blank())
dysbiosis.barplot

# ggsave(dysbiosis.barplot, filename = "data/Community_Composition/Dysbiosis/Dysbiosis_barplot_Species_Bray-Curtis.svg",
#        width = 5, height = 5)


# #--------------------------------------------------------------------------------
# #                            Merged Dysbiosis 
# #--------------------------------------------------------------------------------
# 
# remove_dats()
# load_all_cohorts()
# 
# dys <- dysbiosis_score(dat.species, dist = "bray")
# met <- meta(dat.species) %>% 
#   dplyr::select(donor_id, cohort)
# dys.df <- left_join(dys, met, by = "donor_id")
# 
# dys_all <- 
#   dys.df %>% 
#   ggplot() +
#   geom_density(aes(median, color = group)) +
#   # facet_wrap(~cohort) +
#   labs(x = "Dysbiosis score") +
#   theme_classic()
# dys_all


# #--------------------------------------------------------------------------------
# #                                  TBC 
# #--------------------------------------------------------------------------------
# 
# dysbiosis_thres.tbc <-
#   dys.tbc %>%
#   filter(group == "PC") %>%
#   dplyr::select(median) %>%
#   quantile(na.rm = T, probs = 0.9)
# dys.plot.tbc <- 
#   dys.tbc %>% 
#   ggplot() +
#   geom_density(aes(median, color = group, fill = group), alpha = 0.7) +
#   geom_vline(xintercept = dysbiosis_thres.tbc, size = 1) +
#   labs(x = "Dysbiosis score", title =  "TBC") +
#   scale_fill_manual(values = cols.pdpchc) +
#   scale_colour_manual(values = cols.pdpchc.rim) +
#   annotate("rect", xmin = dysbiosis_thres.tbc, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# dys.plot.tbc
# 
# 
# #--------------------------------------------------------------------------------
# #                                  RUSH
# #--------------------------------------------------------------------------------
# 
# dysbiosis_thres.rush <-
#   dys.rush %>%
#   filter(group == "PC") %>%
#   dplyr::select(median) %>%
#   quantile(na.rm = T, probs = 0.9)
# dys.plot.rush <- 
#   dys.rush %>% 
#   ggplot() +
#   geom_density(aes(median, color = group, fill = group), alpha = 0.7) +
#   geom_vline(xintercept = dysbiosis_thres.rush, size = 1) +
#   labs(x = "Dysbiosis score", title =  "RUSH") +
#   scale_fill_manual(values = cols.pdpchc) +
#   scale_colour_manual(values = cols.pdpchc.rim) +
#   annotate("rect", xmin = dysbiosis_thres.rush, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
#   theme_classic() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# dys.plot.rush
# 
# cowplot::plot_grid(dys.plot.tbc, dys.plot.rush, nrow = 2, axis = "lr")



