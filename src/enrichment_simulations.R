
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")

# Enrichment calculation simulated dataset ----

total_lib_fake <- seq(from = 50, to = 50000, by = 5000)
total_sig_lib_fake <- seq(from = 20, to = 500, by = 20)
path_lib_fake <- seq(from = 10, to = 150, by = 10)
path_sig_lib_fake <- seq(from = 0, to = 50, by = 5)

simulation_df <- tibble()
start_time <- Sys.time()
for (total_lib in total_lib_fake) {
  for (total_sig_lib in total_sig_lib_fake) {
    for (path_lib in path_lib_fake) {
      for (path_sig_lib in path_sig_lib_fake) {
        if (total_lib > total_sig_lib &
          total_lib > path_lib &
          path_lib > path_sig_lib &
          total_lib_fake > path_lib_fake) {
          stats <-
            enrichment_formula(
              N = total_lib,
              n = total_sig_lib,
              m = path_lib,
              k = path_sig_lib
            )
          row2add <- data.frame(
            "N" = total_lib,
            "n" = total_sig_lib,
            "m" = path_lib,
            "k" = path_sig_lib,
            "enrichment" = stats
          )
          simulation_df <- rbind(simulation_df, row2add)
        }
      }
    }
  }
}
end_time <- Sys.time()
cat(
  "simulation complete : ",
  end_time - start_time,
  attr(end_time - start_time, "units"),
  "\n"
)

simulated_data <-
  simulation_df %>%
  filter(m < n)

simulated_data_plot <-
  simulated_data %>%
  ggpairs()

# ggsave(simulated_data_plot, filename = "figures/correlations/enrichment_simulations.png",
#        width = 12, height = 12)

simulated_data %>%
  ggplot(aes(x = k / m, y = enrichment)) +
  geom_density2d_filled() +
  my_clean_theme()

simulated_data %>%
  ggplot(aes(x = N / n, y = k)) +
  geom_density2d_filled() +
  my_clean_theme()

simulated_data %>%
  ggplot(aes(x = N / n, y = -log10(enrichment + 1))) +
  geom_density2d_filled() +
  my_clean_theme()

# write.xlsx(simulated_data, 'figures/correlations/enrichment_simulations.xlsx')
