# Load libraries 
library(tidyverse)
library(data.table)
library(randomForest)
library(parallel)
library(ggplot2)
library(viridis)
library(scales)

options(mc.cores = parallel::detectCores())
dir.create("../Figures/ensemble_cv/RF/", recursive = TRUE, showWarnings = FALSE)
dir.create("../Outputs/ensemble_cv/RF/", recursive = TRUE, showWarnings = FALSE)

# Load and prepare data 
data_raw <- readr::read_csv("../Data/final_dataset.csv")

covariates <- c("uwwLongitude","uwwLatitude","time",
                "pop_density","imd","bame","prop_ind","prop_64","prop_16")

# scaling params from unscaled rows that have outcome
train_stats <- data_raw %>%
  dplyr::filter(!is.na(logRNA)) %>%
  dplyr::summarise(
    m_pd = mean(pop_density), s_pd  = sd(pop_density),
    m_imd = mean(imd), s_imd = sd(imd),
    m_ba = mean(bame), s_ba  = sd(bame),
    m_pi = mean(prop_ind), s_pi  = sd(prop_ind),
    m_p64 = mean(prop_64), s_p64 = sd(prop_64),
    m_p16 = mean(prop_16), s_p16 = sd(prop_16)
  )
scale_with <- function(x, m, s) (x - m) / s

# build the scaled training frame using those params
data <- data_raw %>%
  dplyr::filter(!is.na(logRNA)) %>%
  dplyr::mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd,  train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  )

# quick guards
stopifnot(all(covariates %in% names(data)))

# CV folds 
set.seed(123)
data <- data %>%
  dplyr::mutate(row_id = dplyr::row_number(),
                fold = sample(rep(1:10, length.out = n())))

# CV
run_fold <- function(k) {
  cat("Running Fold", k, "\n")
  train <- data %>% dplyr::filter(fold != k)
  test  <- data %>% dplyr::filter(fold == k)
  
  rf_formula <- as.formula(
    paste("logRNA ~", paste(covariates, collapse = " + "))
  )
  
  # Base RF 1
  set.seed(100 + k*1)
  rf1_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750,
    mtry = 5
  )
  rf1_pred_train <- predict(rf1_fit, newdata = train)
  rf1_pred_test  <- predict(rf1_fit, newdata = test)
  
  # Base RF 2
  set.seed(200 + k*2)
  rf2_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750,
    mtry = 5
  )
  rf2_pred_train <- predict(rf2_fit, newdata = train)
  rf2_pred_test  <- predict(rf2_fit, newdata = test)
  
  # Base RF 3
  set.seed(300 + k*3)
  rf3_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750,
    mtry = 5
  )
  rf3_pred_train <- predict(rf3_fit, newdata = train)
  rf3_pred_test  <- predict(rf3_fit, newdata = test)
  
  # Meta-model
  meta_train <- train %>%
    dplyr::mutate(
      rf1_pred = rf1_pred_train,
      rf2_pred = rf2_pred_train,
      rf3_pred = rf3_pred_train
    )
  meta_fit <- lm(logRNA ~ rf1_pred + rf2_pred + rf3_pred, data = meta_train)
  
  meta_test <- test %>%
    dplyr::select(row_id, fold, logRNA, time, dplyr::all_of(covariates)) %>%
    dplyr::mutate(
      rf1_pred = rf1_pred_test,
      rf2_pred = rf2_pred_test,
      rf3_pred = rf3_pred_test
    )
  meta_test$stacked_pred <- predict(meta_fit, newdata = meta_test)
  
  return(meta_test)
}

# Parallelisation
cv_results <- parallel::mclapply(1:10, run_fold, mc.cores = 10)
cv_results <- dplyr::bind_rows(cv_results)

# Metrics 
cv_results <- cv_results %>%
  dplyr::mutate(Bias = logRNA - stacked_pred)

# Stacking
stopifnot(all(c("rf1_pred","rf2_pred","rf3_pred") %in% names(cv_results)))
meta_fit_oof <- lm(logRNA ~ rf1_pred + rf2_pred + rf3_pred, data = cv_results)

# Train three base learners on FULL training data
rf_formula <- as.formula(
  paste("logRNA ~", paste(c("uwwLongitude","uwwLatitude","time",
                            "pop_density","imd","bame","prop_ind","prop_64","prop_16"),
                          collapse = " + "))
)
set.seed(1101); rf1_full <- randomForest(rf_formula, data = data, ntree = 750, mtry = 5)
set.seed(2202); rf2_full <- randomForest(rf_formula, data = data, ntree = 750, mtry = 5)
set.seed(3303); rf3_full <- randomForest(rf_formula, data = data, ntree = 750, mtry = 5)

# Metrics
ensemble_metrics_by_fold <- cv_results %>%
  dplyr::group_by(fold) %>%
  dplyr::summarise(
    MAE = mean(abs(Bias)),
    RMSE = sqrt(mean(Bias^2)),
    Bias = mean(Bias),
    Correlation = cor(logRNA, stacked_pred),
    .groups = "drop"
  )
readr::write_csv(ensemble_metrics_by_fold, "../Outputs/ensemble_cv/RF/optimised_ensemble_rf_metrics_by_fold.csv")

overall_ensemble_metrics <- tibble::tibble(
  MAE = mean(abs(cv_results$Bias)),
  RMSE = sqrt(mean(cv_results$Bias^2)),
  Bias = mean(cv_results$Bias),
  Correlation = cor(cv_results$logRNA, cv_results$stacked_pred)
)
print("\nOverall ensemble metrics:")
print(overall_ensemble_metrics)
readr::write_csv(overall_ensemble_metrics, "../Outputs/ensemble_cv/RF/optimised_ensemble_rf_metrics_overall.csv")

# Empirical CI
ts_summary <- cv_results %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(
    actual_mean = mean(logRNA, na.rm=TRUE),
    pred_mean = mean(stacked_pred, na.rm=TRUE),
    ci_lower = quantile(stacked_pred, 0.025, na.rm=TRUE),
    ci_upper = quantile(stacked_pred, 0.975, na.rm=TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

# Final df
data_oos <- data %>%
  dplyr::mutate(time = as.numeric(time)) %>%
  dplyr::left_join(
    cv_results %>% dplyr::select(row_id, stacked_pred, fold, Bias, time),
    by = "row_id"
  )

# Clean up any duplicated or renamed 'time' columns
if ("time.y" %in% names(data_oos)) {
  data_oos$time <- data_oos$time.y
  data_oos$time.y <- NULL
}
if ("time.x" %in% names(data_oos)) {
  data_oos$time.x <- NULL
}

data_oos <- data_oos %>%
  dplyr::left_join(
    ts_summary %>% dplyr::select(time, ci_lower, ci_upper),
    by = "time"
  )

# Clean up any named/attributes for export
for (col in c("stacked_pred", "ci_lower", "ci_upper")) {
  if (col %in% names(data_oos)) {
    data_oos[[col]] <- as.numeric(data_oos[[col]])
    attributes(data_oos[[col]]) <- NULL
  }
}
readr::write_csv(data_oos, "../Outputs/ensemble_cv/RF/data_oos_optimised.csv")


# LSOA PREDICTIONS

# Build LSOA covariates aligned to training 
lsoa_raw <- readr::read_csv("../Data/lsoa_prediction_table.csv")

# Map weeks to the SAME time index used in training
time_levels <- sort(unique(data$time)) 
lsoa_pred <- lsoa_raw %>%
  dplyr::mutate(time = as.integer(factor(week, levels = time_levels))) %>%
  dplyr::filter(!is.na(time))  

# Apply the training scaling params captured earlier 
lsoa_pred <- lsoa_pred %>%
  dplyr::mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd,  train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  )

pred_cols <- c("uwwLongitude","uwwLatitude","time",
               "pop_density","imd","bame","prop_ind","prop_64","prop_16")
stopifnot(all(pred_cols %in% names(lsoa_pred)))
stopifnot(!anyNA(lsoa_pred[, pred_cols]))
lsoa_X <- lsoa_pred %>% dplyr::select(dplyr::all_of(pred_cols))

# Predict with full base learners and combine via OOF meta-model
lsoa_pred$rf1_pred <- predict(rf1_full, newdata = lsoa_X)
lsoa_pred$rf2_pred <- predict(rf2_full, newdata = lsoa_X)
lsoa_pred$rf3_pred <- predict(rf3_full, newdata = lsoa_X)
lsoa_pred$stacked_pred <- predict(meta_fit_oof, newdata = lsoa_pred)

# ===== 95% CIs and SD for the stacked predictions
rf_stats_batched <- function(model, newdata, probs = c(0.025, 0.5, 0.975),
                             batch_size = 20000) {
  n <- nrow(newdata)
  qmat <- matrix(NA_real_, nrow = n, ncol = length(probs))
  colnames(qmat) <- paste0("q", sprintf("%05.3f", probs)) # q0.025, q0.500, q0.975
  sds <- numeric(n)
  idx <- seq_len(n)
  batches <- split(idx, ceiling(seq_along(idx) / batch_size))
  for (b in batches) {
    pa <- predict(model, newdata = newdata[b, , drop = FALSE], predict.all = TRUE)
    trees <- if (!is.null(pa$individual)) pa$individual else pa$predictions
    sds[b] <- apply(trees, 1, sd)
    qmat[b, ] <- t(apply(trees, 1, function(x) stats::quantile(x, probs = probs, names = FALSE)))
  }
  list(quantiles = as.data.frame(qmat), sd = sds)
}

# Per-base quantiles & SDs 
stats_rf1 <- rf_stats_batched(rf1_full, lsoa_X)
stats_rf2 <- rf_stats_batched(rf2_full, lsoa_X)
stats_rf3 <- rf_stats_batched(rf3_full, lsoa_X)

# Combine base quantiles into stacked CI with sign-aware bounds
coefs <- coef(meta_fit_oof)  # (Intercept), rf1_pred, rf2_pred, rf3_pred
b0 <- unname(coefs[1]); b1 <- unname(coefs["rf1_pred"]); b2 <- unname(coefs["rf2_pred"]); b3 <- unname(coefs["rf3_pred"])

pick_lower <- function(b, q025, q975) if (b >= 0) q025 else q975
pick_upper <- function(b, q025, q975) if (b >= 0) q975 else q025

stack_lower <- b0 +
  b1 * mapply(pick_lower, b1, stats_rf1$quantiles$q0.025, stats_rf1$quantiles$q0.975) +
  b2 * mapply(pick_lower, b2, stats_rf2$quantiles$q0.025, stats_rf2$quantiles$q0.975) +
  b3 * mapply(pick_lower, b3, stats_rf3$quantiles$q0.025, stats_rf3$quantiles$q0.975)

stack_upper <- b0 +
  b1 * mapply(pick_upper, b1, stats_rf1$quantiles$q0.025, stats_rf1$quantiles$q0.975) +
  b2 * mapply(pick_upper, b2, stats_rf2$quantiles$q0.025, stats_rf2$quantiles$q0.975) +
  b3 * mapply(pick_upper, b3, stats_rf3$quantiles$q0.025, stats_rf3$quantiles$q0.975)

# meta model 
stack_sd <- sqrt(
  (b1 * stats_rf1$sd)^2 +
    (b2 * stats_rf2$sd)^2 +
    (b3 * stats_rf3$sd)^2
)

lsoa_pred$q0.025     <- as.numeric(stack_lower)
lsoa_pred$q0.975     <- as.numeric(stack_upper)
lsoa_pred$stacked_sd <- as.numeric(stack_sd)

# Save LSOA predictions with CIs and SD
readr::write_csv(lsoa_pred, "../Outputs/ensemble_cv/RF/lsoa_predictions_stacked_rf.csv")

