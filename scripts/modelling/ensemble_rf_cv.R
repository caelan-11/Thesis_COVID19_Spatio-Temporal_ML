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

# CV with stakking 
run_fold <- function(k) {
  cat("Running Fold", k, "\n")
  train <- data %>% dplyr::filter(fold != k)
  test  <- data %>% dplyr::filter(fold == k)
  
  rf_formula <- as.formula(
    paste("logRNA ~", paste(covariates, collapse = " + "))
  )
  
  # Base RF 1 (different random seed)
  set.seed(100 + k*1)
  rf1_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750, # Default = 500
    mtry = 5 # Default = 3
  )
  rf1_pred_train <- predict(rf1_fit, newdata = train)
  rf1_pred_test  <- predict(rf1_fit, newdata = test)
  
  # Base RF 2 (different random seed)
  set.seed(200 + k*2)
  rf2_fit <- randomForest(
    formula = rf_formula,
    data = train,
    ntree = 750,
    mtry = 5
  )
  rf2_pred_train <- predict(rf2_fit, newdata = train)
  rf2_pred_test  <- predict(rf2_fit, newdata = test)
  
  # Base RF 3 (different random seed)
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

# Parallelisation of folds
cv_results <- parallel::mclapply(1:10, run_fold, mc.cores = 10)
cv_results <- dplyr::bind_rows(cv_results)

# Metrics 
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


# Empircal CI 
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

# Save to0 .csv
readr::write_csv(data_oos, "../Outputs/ensemble_cv/RF/data_oos_optimised.csv")


