
source("src/load_packages.R")
source("src/miscellaneous_funcs.R")
source("src/load_phyloseq_obj.R")
source("src/metadata_prep_funcs.R")

library(mlbench)
library(tidymodels)
library(kernlab)
library(furrr)
library(scales)


library(mlbench)
sim_data <- function(n) {
  tmp <- mlbench.friedman1(n, sd = 1)
  tmp <- cbind(tmp$x, tmp$y)
  tmp <- as.data.frame(tmp)
  names(tmp)[ncol(tmp)] <- "y"
  tmp
}

set.seed(9815)
train_dat <- sim_data(100)
large_dat <- sim_data(10^5)
# NESTED RESAMPLING
# To get started, the types of resampling methods need to be specified. This isn’t a large data set, so 5 repeats of 10-fold cross validation will be used as the outer resampling method for generating the estimate of overall performance. To tune the model, it would be good to have precise estimates for each of the values of the tuning parameter so let’s use 25 iterations of the bootstrap. This means that there will eventually be 5 * 10 * 25 = 1250 models that are fit to the data per tuning parameter. These models will be discarded once the performance of the model has been quantified.

# To create the tibble with the resampling specifications:


library(tidymodels)
results <- nested_cv(train_dat,
  outside = vfold_cv(repeats = 3),
  inside = bootstraps(times = 11)
)
results

# The splitting information for each resample is contained in the split objects. Focusing on the second fold of the first repeat:

results$splits[[2]]

# <90/10/100> indicates the number of observations in the analysis set, assessment set, and the original data.

# Each element of inner_resamples has its own tibble with the bootstrapping splits.

results$inner_resamples[[5]]

# These are self-contained, meaning that the bootstrap sample is aware that it is a sample of a specific 90% of the data:

results$inner_resamples[[5]]$splits[[1]]

# To start, we need to define how the model will be created and measured. Let’s use a radial basis support vector machine model via the function kernlab::ksvm. This model is generally considered to have two tuning parameters: the SVM cost value and the kernel parameter sigma. For illustration purposes here, only the cost value will be tuned and the function kernlab::sigest will be used to estimate sigma during each model fit. This is automatically done by ksvm.
#
# After the model is fit to the analysis set, the root-mean squared error (RMSE) is computed on the assessment set. One important note: for this model, it is critical to center and scale the predictors before computing dot products. We don’t do this operation here because mlbench.friedman1 simulates all of the predictors to be standardized uniform random variables.
#
# Our function to fit the model and compute the RMSE is:

library(kernlab)

# `object` will be an `rsplit` object from our `results` tibble
# `cost` is the tuning parameter
svm_rmse <- function(object, cost) {
  y_col <- ncol(object$data)
  mod <-
    svm_rbf(mode = "regression", cost = cost) %>%
    set_engine("kernlab") %>%
    fit(y ~ ., data = analysis(object))

  holdout_pred <-
    predict(mod, assessment(object) %>% dplyr::select(-y)) %>%
    bind_cols(assessment(object) %>% dplyr::select(y))
  rmse(holdout_pred, truth = y, estimate = .pred)$.estimate
}

# In some case, we want to parameterize the function over the tuning parameter:
rmse_wrapper <- function(cost, object) svm_rmse(object, cost)
# For the nested resampling, a model needs to be fit for each tuning parameter and each bootstrap split. To do this, create a wrapper:

# `object` will be an `rsplit` object for the bootstrap samples
tune_over_cost <- function(object) {
  tibble(cost = 2^seq(-2, 2, by = 1)) %>%
    mutate(RMSE = map_dbl(cost, rmse_wrapper, object = object))
}

# Since this will be called across the set of outer cross - validation splits, another wrapper is required:# `object` is an `rsplit` object in `results$inner_resamples`
summarize_tune_results <- function(object) {
  # Return row-bound tibble that has the 25 bootstrap results
  map_df(object$splits, tune_over_cost) %>%
    # For each value of the tuning parameter, compute the
    # average RMSE which is the inner bootstrap estimate.
    group_by(cost) %>%
    dplyr::summarize(
      mean_RMSE = mean(RMSE, na.rm = TRUE),
      n = length(RMSE),
      .groups = "drop"
    )
}
# Now that those functions are defined, we can execute all the inner resampling loops:

# tuning_results <- map(results$inner_resamples, summarize_tune_results)
# Alternatively, since these computations can be run in parallel, we can use the furrr package. Instead of using map(), the function future_map() parallelizes the iterations using the future package. The multisession plan uses the local cores to process the inner resampling loop. The end results are the same as the sequential computations.

library(furrr)
plan(multisession)
tuning_results <- future_map(results$inner_resamples,
  summarize_tune_results,
  .options = future_options(seed = TRUE),
  .progress = TRUE
)


# The object tuning_results is a list of data frames for each of the 50 outer resamples.

# Let’s make a plot of the averaged results to see what the relationship is between the RMSE and the tuning parameters for each of the inner bootstrapping operations:

library(scales)

pooled_inner <- tuning_results %>% bind_rows()

best_cost <- function(dat) dat[which.min(dat$mean_RMSE), ]

p <-
  ggplot(pooled_inner, aes(x = cost, y = mean_RMSE)) +
  scale_x_continuous(trans = "log2") +
  xlab("SVM Cost") +
  ylab("Inner RMSE")

for (i in 1:length(tuning_results)) {
  p <- p +
    geom_line(data = tuning_results[[i]], alpha = .2) +
    geom_point(data = best_cost(tuning_results[[i]]), pch = 16, alpha = 3 / 4)
}

p <- p + geom_smooth(data = pooled_inner, se = FALSE)
p


cost_vals <-
  tuning_results %>%
  map_df(best_cost) %>%
  select(cost)

results <-
  bind_cols(results, cost_vals) %>%
  mutate(cost = factor(cost, levels = paste(2^seq(-2, 8, by = 1))))

ggplot(results, aes(x = cost)) +
  geom_bar() +
  xlab("SVM Cost") +
  scale_x_discrete(drop = FALSE)


results <-
  results %>%
  mutate(RMSE = map2_dbl(splits, cost, svm_rmse))

summary(results$RMSE)


not_nested <-
  future_map(results$splits, tune_over_cost) %>%
  bind_rows()

outer_summary <- not_nested %>%
  group_by(cost) %>%
  dplyr::summarize(outer_RMSE = mean(RMSE), n = length(RMSE))
outer_summary


ggplot(outer_summary, aes(x = cost, y = outer_RMSE)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = "log2") +
  xlab("SVM Cost") +
  ylab("RMSE")
