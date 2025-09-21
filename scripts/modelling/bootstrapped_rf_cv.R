# Load libraries
library(tidyverse)
library(data.table)
library(parallel)
library(scales)
library(randomForest)
library(ranger)
library(stringr)


options(mc.cores = parallel::detectCores())

# Load and prepare data
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

# Assign folds spatially by site
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

# CV loop and bootstrapping
run_fold_rf <- function(k) {
  tryCatch({
    cat("Running Fold", k, "\n")
    train <- data %>% dplyr::filter(fold != k)
    test  <- data %>% dplyr::filter(fold == k)
    if (nrow(test) == 0 || nrow(train) < 5) stop("Fold has insufficient train/test rows")
    
    y_test <- test$logRNA
    rf_formula <- as.formula(
      paste("logRNA ~", paste(covariates, collapse = " + "))
    )
    
    # Fit RF 
    rf_fit <- ranger(
      formula = rf_formula,
      data = train,
      #ntree = 750,
      num.trees = 750,
      mtry = 5
    )
    
    rf_pred_all <- predict(rf_fit, data = test, predict.all = TRUE)
    RF_pred <- rowMeans(rf_pred_all$predictions)
    q0.025 <- apply(rf_pred_all$predictions, 1, quantile, probs = 0.025)
    q0.975 <- apply(rf_pred_all$predictions, 1, quantile, probs = 0.975)
    coverage <- as.numeric(y_test >= q0.025 & y_test <= q0.975)
    
    results <- tibble(
      row_id = test$row_id,
      fold = k,
      truth = y_test,
      RF_pred = RF_pred,
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
    
    fold_metrics <- eval_metrics(results$truth, results$RF_pred) %>%
      mutate(
        Model = "RF",
        Fold = k,
        Coverage = mean(coverage, na.rm = TRUE)
      )
    list(results = results, metrics = fold_metrics)
  }, error = function(e) {
    cat(sprintf("Fold %d failed: [%s] %s\n", k, class(e)[1], conditionMessage(e)))
    return(NULL)
  })
}


# Runs folds in parrallel 
n_folds <- length(unique(data$fold))
all_folds <- parallel::mclapply(1:n_folds, run_fold_rf, mc.cores = min(10, n_folds))
good_folds <- Filter(Negate(is.null), all_folds)
if (length(good_folds) == 0) stop("All folds failed! Check variable names and data splits.")

all_results <- purrr::map(good_folds, "results")
all_metrics <- purrr::map(good_folds, "metrics")
final_results <- dplyr::bind_rows(all_results)
final_metrics <- dplyr::bind_rows(all_metrics)

overall_coverage <- mean(final_results$coverage, na.rm=TRUE)
cat(sprintf("\nOverall 95%% empirical coverage (RF): %.2f%%\n", 100 * overall_coverage))

# Merge predictions and fold back to data via row_id
data_oos <- data %>%
  dplyr::left_join(
    final_results %>%
      dplyr::select(row_id, fold, RF_pred, q0.025, q0.975, coverage),
    by = "row_id"
  )

# Debugging
if ("time.y" %in% names(data_oos)) {
  data_oos$time <- data_oos$time.y
  data_oos$time.y <- NULL
}
if ("time.x" %in% names(data_oos)) {
  data_oos$time.x <- NULL
}

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

# Print outputs
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

# Save output to .csv
readr::write_csv(data_oos, "../Outputs/ML_cv/RF/final_predictions_with_folds.csv")
readr::write_csv(final_metrics, "../Outputs/ML_cv/RF/final_metrics_by_fold.csv")

#####################################################################################################################
# LSOA predictions 

# Load LSOA prediction table
lsoa.cov <- readr::read_csv("../Data/lsoa_prediction_table.csv") %>%
  dplyr::mutate(time_index = week)

# Check cols
req_cols <- c("uwwLongitude","uwwLatitude","time_index",
              "pop_density","imd","bame","prop_ind","prop_64","prop_16")
missing <- setdiff(req_cols, names(lsoa.cov))
if (length(missing) > 0) {
  stop("LSOA table is missing required columns: ", paste(missing, collapse = ", "))
}

# Align the LSOA time_index with training time index
time_levels <- sort(unique(data$time))  
lsoa.cov <- lsoa.cov %>%
  dplyr::mutate(time_index = as.integer(factor(time_index, levels = time_levels)))

# Drop weeks not seen in training 
bad_weeks <- lsoa.cov %>% dplyr::filter(is.na(time_index)) %>% dplyr::distinct(week) %>% dplyr::pull(week)
if (length(bad_weeks) > 0) {
  warning("Dropping ", length(bad_weeks), " LSOA rows with weeks not in training: ",
          paste(head(bad_weeks, 10), collapse = ", "),
          if (length(bad_weeks) > 10) " ..." else "")
  lsoa.cov <- lsoa.cov %>% dplyr::filter(!is.na(time_index))
}

# Apply training scaling to LSOA covariates
train_stats <- data %>% dplyr::summarise(
  m_pd = mean(pop_density), s_pd  = sd(pop_density),
  m_imd = mean(imd), s_imd = sd(imd),
  m_ba = mean(bame), s_ba = sd(bame),
  m_pi = mean(prop_ind), s_pi = sd(prop_ind),
  m_p64 = mean(prop_64), s_p64 = sd(prop_64),
  m_p16 = mean(prop_16), s_p16 = sd(prop_16)
)
scale_with <- function(x, m, s) (x - m) / s

lsoa.cov <- lsoa.cov %>%
  dplyr::mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd, train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  )

# Build the prediction covariate frame in the SAME order as training
lsoa_covariate_df <- lsoa.cov %>%
  dplyr::select(all_of(req_cols))

# Refit RF on all training data and calculate quantile prediction
rf_formula <- as.formula(paste("logRNA ~", paste(req_cols, collapse = " + ")))

rf_fit_full <- ranger(
  formula = rf_formula,
  data = data,
  num.trees = 750,
  mtry = 5,
  quantreg  = TRUE
)

# Predict median and 95% interval 
pred_q <- predict(
  rf_fit_full,
  data      = lsoa_covariate_df,
  type      = "quantiles",
  quantiles = c(0.025, 0.5, 0.975)
)

lsoa.cov$RF_pred <- pred_q$predictions[, 2]
lsoa.cov$q0.025  <- pred_q$predictions[, 1]
lsoa.cov$q0.975  <- pred_q$predictions[, 3]

# Save as .csv
readr::write_csv(lsoa.cov, "../Outputs/ML_cv/RF/lsoa_predictions.csv")

