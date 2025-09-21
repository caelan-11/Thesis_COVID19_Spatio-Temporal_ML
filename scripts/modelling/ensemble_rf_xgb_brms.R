# Load libraries
library(tidyverse)
library(data.table)
library(mgcv)
library(xgboost)
library(randomForest)
library(parallel)
library(ggplot2)
library(viridis)
library(scales)

options(mc.cores = parallel::detectCores())
dir.create("../Figures/ensemble_cv/", recursive = TRUE, showWarnings = FALSE)
dir.create("../Outputs/ensemble_cv/", recursive = TRUE, showWarnings = FALSE)

# Load and prepare data 
data <- readr::read_csv("../Data/final_dataset.csv")

covariates <- c("uwwLongitude", "uwwLatitude", "time",
                "pop_density", "imd", "bame", "prop_ind", "prop_64", "prop_16")

data <- data %>%
  dplyr::filter(!is.na(logRNA)) %>%
  dplyr::mutate(across(c(pop_density, imd, bame, prop_ind, prop_64, prop_16), ~ scale(.x)[,1]))

# Create folds
set.seed(123)
data <- data %>%
  dplyr::mutate(row_id = dplyr::row_number(),
                fold = sample(rep(1:10, length.out = n())))

# CV
run_fold <- function(k) {
  cat("Running Fold", k, "\n")
  train <- data %>% dplyr::filter(fold != k)
  test  <- data %>% dplyr::filter(fold == k)
  
  # GAM
  gam_fit <- mgcv::gam(
    logRNA ~ s(uwwLongitude, uwwLatitude) + s(time) +
      pop_density + imd + bame + prop_ind + prop_64 + prop_16,
    data = train
  )
  gam_pred <- predict(gam_fit, newdata = test)
  
  # RF
  rf_formula <- as.formula(
    paste("logRNA ~", paste(covariates, collapse = " + "))
  )
  rf_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750,
    mtry = 5
  )
  rf_pred <- predict(rf_fit, newdata = test)
  
  # XGBoost
  x_train <- model.matrix(~ . - logRNA, data = train[, c("logRNA", covariates)])
  y_train <- train$logRNA
  x_test <- model.matrix(~ . - logRNA, data = test[, c("logRNA", covariates)])
  
  
  # Optimised parameters via grid search
  xgb_fit <- xgboost(
    x = x_train,
    y = y_train,
    objective = "reg:squarederror",
    nrounds = 250,
    max_depth = 6,
    eta = 0.1,
    subsample = 0.8,
    colsample_bytree = 1,
    verbose = 0
  )
  
  xgb_pred <- predict(xgb_fit, newdata = x_test)
  
  # Meta-model (stack)
  meta_train <- train %>%
    dplyr::mutate(
      gam_pred = predict(gam_fit, newdata = train),
      rf_pred  = predict(rf_fit, newdata = train),
      xgb_pred = predict(xgb_fit, newdata = model.matrix(~ . - logRNA, data = train[, c("logRNA", covariates)]))
    )
  meta_fit <- lm(logRNA ~ gam_pred + rf_pred + xgb_pred, data = meta_train)
  
  # OOS predictions for held-out fold
  meta_test <- test %>%
    dplyr::select(row_id, fold, logRNA, time, dplyr::all_of(covariates)) %>%
    dplyr::mutate(
      gam_pred = gam_pred,
      rf_pred = rf_pred,
      xgb_pred = xgb_pred
    )
  meta_test$stacked_pred <- predict(meta_fit, newdata = meta_test)
  return(meta_test)
}

# Run folds in parallel ----
cv_results <- parallel::mclapply(1:10, run_fold, mc.cores = 10)
cv_results <- dplyr::bind_rows(cv_results)

# compute bias and other metrics 
cv_results <- cv_results %>%
  dplyr::mutate(Bias = logRNA - stacked_pred)

ensemble_metrics_by_fold <- cv_results %>%
  dplyr::group_by(fold) %>%
  dplyr::summarise(
    MAE = mean(abs(Bias)),
    RMSE = sqrt(mean(Bias^2)),
    Bias = mean(Bias),
    Correlation = cor(logRNA, stacked_pred)
  ) %>%
  dplyr::ungroup()

readr::write_csv(ensemble_metrics_by_fold, "../Outputs/ensemble_cv/RF_XGB_GAM/optimised_ensemble_metrics_by_fold.csv")

# Overall ensemble metrics 
overall_ensemble_metrics <- tibble::tibble(
  MAE = mean(abs(cv_results$Bias)),
  RMSE = sqrt(mean(cv_results$Bias^2)),
  Bias = mean(cv_results$Bias),
  Correlation = cor(cv_results$logRNA, cv_results$stacked_pred)
)

print("\nOverall ensemble metrics:")
print(overall_ensemble_metrics)
readr::write_csv(overall_ensemble_metrics, "../Outputs/ensemble_cv/RF_XGB_GAM/optimised_ensemble_metrics_overall.csv")

 # empirical CI
ts_summary <- cv_results %>%
  dplyr::group_by(time) %>%
  dplyr::summarise(
    actual_mean = mean(logRNA, na.rm=TRUE),
    pred_mean = mean(stacked_pred, na.rm=TRUE),
    ci_lower = quantile(stacked_pred, 0.025, na.rm=TRUE),
    ci_upper = quantile(stacked_pred, 0.975, na.rm=TRUE),
    n = dplyr::n()
  ) %>%
  dplyr::ungroup()

data_oos <- data %>%
  dplyr::mutate(time = as.numeric(time)) %>%
  dplyr::left_join(
    cv_results %>% dplyr::select(row_id, stacked_pred, fold, Bias, time),
    by = "row_id"
  )

# Clean up 
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

for (col in c("stacked_pred", "ci_lower", "ci_upper")) {
  if (col %in% names(data_oos)) {
    data_oos[[col]] <- as.numeric(data_oos[[col]])
    attributes(data_oos[[col]]) <- NULL
  }
}

# Save results to .csv
readr::write_csv(data_oos, "../Outputs/ensemble_cv/RF_XGB_GAM/data_oos_optimised.csv")

