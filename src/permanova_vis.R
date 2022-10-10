# PERMANOVA PLOTS

source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/metadata_prep_funcs.R")

permdf_all <-
  read.csv(file = "files/permanova_analysis_2022-09-29.csv", header = TRUE, row.names = NULL)

permdf_all %>%
  filter(vars == "paired") %>% 
  summarize(min = min(R2),
            max = max(R2))

#_______________________________________________________________________________
# Save Supplementary Table ----
perm_save <- list()
perm_save[["Merged"]] <- permdf_all %>% filter(study == "Merged")
perm_save[["TBC"]] <- permdf_all %>% filter(study == "TBC")
perm_save[["RUMC"]] <- permdf_all %>% filter(study == "Rush")

supp_loc <- "files/Supplementary Tables"
openxlsx::write.xlsx(perm_save,
                     file = glue("{supp_loc}/Table_S2 permanova_analysis_{Sys.Date()}.xlsx"))
#_______________________________________________________________________________

permdf_all %>% 
  filter(vars %in% c("state_of_residence", "quadrant_of_residence"), study == "TBC")

#' after visualizing both state and quadrant results, quadrant seems comparable 
#' in both variance and significance, and reduces single sample categorization
#' thus, we will move ahead using geographic quadrant for downstream modeling.

permdf <- permdf_all %>% 
  filter(study == "Merged") %>% 
  mutate(FDR.symbol = sig.symbol.generator(FDR, shh = T)) %>% 
  filter(vars != "state_of_residence") 


#-------------------------------------------------------------------------------
#####                           Plotting                                 #####
#-------------------------------------------------------------------------------

perm.aitchisons_0.1 <-
  permdf %>%
  dplyr::filter(FDR <= 0.1,
                data_type %ni% c("Enzymes.slim", "Pfams.slim")) %>% 
  dplyr::mutate(var_labs = paste0(vars, " (n=", n_meta, ")")) %>%
  arrange(FDR) %>%
  ggplot(aes(
    x = data_type,
    y = fct_reorder(var_labs, R2),
    fill = FDR
  )) +
  theme_bw() +
  geom_point(aes(size = R2), shape = 21, stroke = 0.2) +
  labs(
    x = NULL,
    y = NULL,
    fill = "FDR",
    size = expression(R ^ "2" ~ "(%)")
  ) +
  scale_y_discrete(position = "right") +
  scale_fill_viridis_c(
    option = "magma",
    begin = 0,
    end = 1,
    direction = -1,
    na.value = "transparent"
  ) +
  scale_x_discrete(labels = c("Eggnogs.slim" = "eggNOGs",
                              "KOs.slim" = "KOs")) +
  scale_size(breaks = c(1, 5, 10, 20), labels = c(1, 5, 10, 20)) +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "left",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
perm.aitchisons_0.1

ggsave(perm.aitchisons_0.1,
  filename = paste0("data/PERMANOVA_Analysis/PERMANOVA_Aitchisons_summary_q.1_", Sys.Date(), ".svg"),
  width = 5, height = 4
)

#  Here we additionally filter for assocations that also show significance 
# in at least one cohort independently (RUMC or TBC) to remove assocations that
# are driven by cohort to cohort differences

#' The reason for this approach and not modeling the PERMNOVA with a 
# blocking factor is because the majority of meta-data are specific to one cohort,
# which creates issues when modeling with a random effect


individual_cohort_assocations <- 
  permdf_all %>% 
  filter(FDR <= 0.1, study != "Merged") %>% 
  group_by(data_type, vars) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  mutate(data_type_var = paste0(data_type,"-", vars))


# Visualizing metadata with association in RUMC, TBC, and combined datasets
permdf_all_label <- permdf_all %>% 
  filter(FDR <= 0.1)

permdf_all %>% 
  filter(data_type %ni% c("Enzymes.slim", "Pfams.slim")) %>% 
  ggplot(aes(x=-log10(FDR), y=R2)) +
  geom_point(aes(color = data_type)) +
  facet_wrap(~study) +
  scale_y_continuous(limits = c(0, 40)) +
  geom_rect(aes(xmin = -Inf, xmax = -log10(0.1), ymin = -Inf, ymax = Inf), 
            alpha = .02, fill = "lightgrey") +
  geom_text_repel(data = permdf_all_label,
                  aes(x=-log10(FDR), y=R2, label = vars))


#_______________________________________________________________________________
#' Annotating associations as either Merged dataset specific, or shared in at least one
#' cohort specific analysis

perm.aitchisons_0.1_annot <- permdf %>%
  mutate(data_type_var = paste0(data_type,"-", vars)) %>% 
  mutate(stratified_association = case_when(
    data_type_var %in% individual_cohort_assocations$data_type_var ~ "Yes",
    TRUE ~ "No"
  )) %>% 
  dplyr::filter(FDR <= 0.1,
                data_type %ni% c("Enzymes.slim", "Pfams.slim")) %>% 
  dplyr::mutate(var_labs = paste0(vars, " (n=", n_meta, ")")) %>%
  arrange(FDR) %>%
  ggplot(aes(
    x = data_type,
    y = fct_reorder(var_labs, R2),
    fill = FDR
  )) +
  theme_bw() +
  geom_point(aes(size = R2, shape = stratified_association), stroke = 0.2) +
  labs(
    x = NULL,
    y = NULL,
    fill = "FDR",
    size = expression(R ^ "2" ~ "(%)")
  ) +
  scale_y_discrete(position = "right") +
  scale_fill_viridis_c(
    option = "magma",
    begin = 0,
    end = 1,
    direction = -1,
    na.value = "transparent"
  ) +
  scale_x_discrete(labels = c("eggNOGs.slim" = "eggNOGs",
                              "KOs.slim" = "KOs",
                              "Pathways.slim" = "Pathways")) +
  scale_size(breaks = c(1, 5, 10, 20), labels = c(1, 5, 10, 20)) +
  scale_shape_manual(values = c("Yes" = 24, "No" = 21)) +
  labs(shape = "Independently \nsignificant in \nRUMC or TBC") +
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "left",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

perm.aitchisons_0.1_annot

ggsave(perm.aitchisons_0.1_annot,
       filename = paste0("data/PERMANOVA_Analysis/PERMANOVA_Aitchisons_summary_q.1_stratified_association_", Sys.Date(), ".svg"),
       width = 6, height = 4.5
)
