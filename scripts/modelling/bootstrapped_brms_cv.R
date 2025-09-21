# Load libraries
library(tidyverse)
library(data.table)
library(parallel)
library(scales)
library(brms)

options(mc.cores = parallel::detectCores())

# Load and prep data
data <- readr::read_csv("../Data/final_dataset.csv") %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::filter(!is.na(logRNA)) %>%
  dplyr::mutate(across(
    c(pop_density, imd, bame, prop_ind, prop_64, prop_16),
    ~ scale(.x)[,1]
  ))

# time indexing
if (!"time_index" %in% names(data)) {
  if ("time" %in% names(data)) {
    data <- data %>% mutate(time_index = as.integer(factor(time)))
    message("Added time_index column from 'time'.")
  } else {
    stop("Neither time_index nor time found in data! Please check.")
  }
}

covariates <- c("uwwLongitude", "uwwLatitude", "time_index",
                "pop_density", "imd", "bame", "prop_ind", "prop_64", "prop_16")

missing_vars <- setdiff(c("logRNA", covariates), names(data))
if (length(missing_vars) > 0) {
  stop("The following variables are missing from the data: ", paste(missing_vars, collapse = ", "))
}

# Spatially asign folds by site 
set.seed(123)
if (!"fold" %in% names(data)) {
  if ("Site code" %in% names(data)) {
    unique_sites <- unique(data$`Site code`)
    site_folds <- tibble(
      `Site code` = sample(unique_sites),
      fold = rep(1:10, length.out = length(unique_sites))
    )
    data <- dplyr::left_join(data, site_folds, by = "Site code")
    message("Assigned folds using Site code.")
  } else {
    data <- data %>%
      dplyr::mutate(fold = sample(rep(1:10, length.out = n())))
    message("Assigned folds randomly (no Site code found).")
  }
}

if (!"row_id" %in% names(data)) stop("row_id is missing from data!")

# CV and bootstrapping 
run_fold_brms <- function(k) {
  tryCatch({
    cat("Running Fold", k, "\n")
    train <- data %>% dplyr::filter(fold != k)
    test <- data %>% dplyr::filter(fold == k)
    if (nrow(test) == 0 || nrow(train) < 5) stop("Fold has insufficient train/test rows")
    
    # Fit BRMS
    brms_fit <- brm(
      logRNA ~ s(uwwLongitude, uwwLatitude) + s(time_index) +
        pop_density + imd + bame + prop_ind + prop_64 + prop_16,
      data = train,
      family = gaussian(),
      chains = 4, 
      iter = 2000,
      warmup = 1000,
      seed = 123,
      control = list(adapt_delta = 0.95),
      silent = TRUE, refresh = 0
    )
    pred_brms <- posterior_predict(brms_fit, newdata = test) 
    
    brms_mean <- colMeans(pred_brms)
    q0.025 <- apply(pred_brms, 2, quantile, probs = 0.025)
    q0.975 <- apply(pred_brms, 2, quantile, probs = 0.975)
    coverage <- as.numeric(test$logRNA >= q0.025 & test$logRNA <= q0.975)
    
    results <- tibble(
      row_id = test$row_id,
      fold = k,
      truth = test$logRNA,
      BRMS_pred  = brms_mean,
      q0.025 = q0.025,
      q0.975 = q0.975,
      coverage = coverage
    )
    
    eval_metrics <- function(obs, pred) {
      tibble(
        MAE = mean(abs(obs - pred)),
        RMSE = sqrt(mean((obs - pred)^2)),
        Bias = mean(obs - pred),
        Correlation = cor(obs, pred)
      )
    }
    
    fold_metrics <- eval_metrics(results$truth, results$BRMS_pred) %>%
      mutate(
        Model = "BRMS",
        Fold = k,
        Coverage = mean(coverage, na.rm = TRUE)
      )
    list(results = results, metrics = fold_metrics)
  }, error = function(e) {
    cat(sprintf("Fold %d failed: [%s] %s\n", k, class(e)[1], conditionMessage(e)))
    return(NULL)
  })
}

# Parrallelising folds
n_folds <- length(unique(data$fold))
all_folds <- parallel::mclapply(1:n_folds, run_fold_brms, mc.cores = min(10, n_folds))
good_folds <- Filter(Negate(is.null), all_folds)
if (length(good_folds) == 0) stop("All folds failed! Check variable names and data splits.")

all_results <- purrr::map(good_folds, "results")
all_metrics <- purrr::map(good_folds, "metrics")
final_results <- dplyr::bind_rows(all_results)
final_metrics <- dplyr::bind_rows(all_metrics)

overall_coverage <- mean(final_results$coverage, na.rm=TRUE)
cat(sprintf("\nOverall 95%% empirical coverage (BRMS): %.2f%%\n", 100 * overall_coverage))

# Merge predictions and fold back to data via row_id 
data_oos <- data %>%
  dplyr::left_join(
    final_results %>%
      dplyr::select(row_id, fold, BRMS_pred, q0.025, q0.975, coverage),
    by = "row_id"
  )

# Remoive duplicates (debugging)
if ("time.y" %in% names(data_oos)) {
  data_oos$time <- data_oos$time.y
  data_oos$time.y <- NULL
}
if ("time.x" %in% names(data_oos)) {
  data_oos$time.x <- NULL
}

# Debugging
badcols <- which(!sapply(data_oos, function(x) is.atomic(x) && is.vector(x)))
if (length(badcols) > 0) {
  print("Non-vector columns found in data_oos:")
  print(names(data_oos)[badcols])
  for (col in names(data_oos)[badcols]) {
    if (is.matrix(data_oos[[col]]) && ncol(data_oos[[col]]) == 1) {
      data_oos[[col]] <- as.numeric(data_oos[[col]][,1])
    } else if (is.matrix(data_oos[[col]])) {
      warning(paste("Column", col, "is a multi-column matrix! Only first column will be used."))
      data_oos[[col]] <- as.numeric(data_oos[[col]][,1])
    }
    if (is.list(data_oos[[col]])) {
      data_oos[[col]] <- unlist(data_oos[[col]])
    }
    attributes(data_oos[[col]]) <- NULL
  }
}

# Print output
cat("\n=== Final Metrics Across Folds ===\n")
print(final_metrics)

cat("\n=== Sample of Final Predictions ===\n")
print(head(data_oos))

cat("\nSummary Statistics:\n")
summary_stats <- final_metrics %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    MAE = mean(MAE),
    RMSE = mean(RMSE),
    Bias = mean(Bias),
    Correlation = mean(Correlation),
    Coverage = mean(Coverage, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_stats)

# Save Output to .csv
readr::write_csv(data_oos, "../Outputs/ML_cv/BRMS/final_predictions_with_folds.csv")
readr::write_csv(final_metrics, "../Outputs/ML_cv/BRMS/final_metrics_by_fold.csv")


################################################################################################
# LSOA Predictions using BRMS model ===

# Load Data 
data_raw <- readr::read_csv("../Data/final_dataset.csv")

train_stats <- data_raw %>%
  dplyr::filter(!is.na(logRNA)) %>%
  dplyr::summarise(
    m_pd = mean(pop_density), s_pd = sd(pop_density),
    m_imd = mean(imd), s_imd = sd(imd),
    m_ba = mean(bame), s_ba = sd(bame),
    m_pi = mean(prop_ind), s_pi = sd(prop_ind),
    m_p64 = mean(prop_64), s_p64 = sd(prop_64),
    m_p16 = mean(prop_16), s_p16 = sd(prop_16)
  )
scale_with <- function(x, m, s) (x - m) / s

# Load the LSOA covariates 
lsoa_raw <- readr::read_csv("../Data/lsoa_prediction_table.csv")

# Map LSOA weeks time index from model training 
time_levels <- sort(unique(data$time_index))
lsoa_pred <- lsoa_raw %>%
  mutate(time_index = as.integer(factor(week, levels = time_levels))) %>%
  filter(!is.na(time_index))

# Apply training scaling
lsoa_pred <- lsoa_pred %>%
  mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd, train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  )

# Keep predictors in same order as training
pred_cols <- c("uwwLongitude","uwwLatitude","time_index",
               "pop_density","imd","bame","prop_ind","prop_64","prop_16")
stopifnot(all(pred_cols %in% names(lsoa_pred)))

# Fit BRMS model on FULL training data
brms_full <- brm(
  logRNA ~ s(uwwLongitude, uwwLatitude) + s(time_index) +
    pop_density + imd + bame + prop_ind + prop_64 + prop_16,
  data = data, 
  family = gaussian(),
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

# Predict for all LSOA week rows
pred_draws <- posterior_predict(brms_full, newdata = lsoa_pred)
lsoa_pred$BRMS_pred <- colMeans(pred_draws)
lsoa_pred$q0.025 <- apply(pred_draws, 2, quantile, probs = 0.025)
lsoa_pred$q0.975 <- apply(pred_draws, 2, quantile, probs = 0.975)

# Save LSOA predictions as .csv
readr::write_csv(lsoa_pred, "../Outputs/ML_cv/BRMS/lsoa_predictions_brms.csv")



