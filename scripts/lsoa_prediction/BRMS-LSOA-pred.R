# Import Libraries
library(tidyverse)
library(brms)
library(parallel)
library(readr)
library(stringr)

options(mc.cores = parallel::detectCores())

# Load data
data_raw <- readr::read_csv("../Data/final_dataset.csv")

if (!"time_index" %in% names(data_raw)) {
  if ("time" %in% names(data_raw)) {
    data_raw <- data_raw %>% mutate(time_index = as.integer(factor(time)))
  } else {
    stop("Neither 'time_index' nor 'time' present in data.")
  }
}

# variables to scale
scale_vars <- c("pop_density","imd","bame","prop_ind","prop_64","prop_16")

# scaling params from rows with observed outcome
train_stats <- data_raw %>%
  filter(!is.na(logRNA)) %>%
  summarise(
    m_pd = mean(pop_density, na.rm=TRUE), s_pd  = sd(pop_density, na.rm=TRUE),
    m_imd = mean(imd, na.rm=TRUE), s_imd = sd(imd, na.rm=TRUE),
    m_ba = mean(bame, na.rm=TRUE), s_ba  = sd(bame, na.rm=TRUE),
    m_pi = mean(prop_ind, na.rm=TRUE), s_pi  = sd(prop_ind, na.rm=TRUE),
    m_p64 = mean(prop_64, na.rm=TRUE), s_p64 = sd(prop_64, na.rm=TRUE),
    m_p16 = mean(prop_16, na.rm=TRUE), s_p16 = sd(prop_16, na.rm=TRUE)
  )
scale_with <- function(x, m, s) (x - m) / s

# Build the scaled training frame
data <- data_raw %>%
  filter(!is.na(logRNA)) %>%
  mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd,  train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  ) %>%
  mutate(row_id = row_number())

# predictors used in the model
covariates <- c("uwwLongitude","uwwLatitude","time_index",
                "pop_density","imd","bame","prop_ind","prop_64","prop_16")
stopifnot(all(c("logRNA", covariates) %in% names(data)))

# Fit BRMS on all training data 
brms_full <- brm(
  formula = logRNA ~ s(uwwLongitude, uwwLatitude) + s(time_index) +
    pop_density + imd + bame + prop_ind + prop_64 + prop_16,
  data   = data,
  family  = gaussian(),
  chains = 4,
  cores = 1,           
  iter = 2000,
  warmup = 1000,
  seed = 123,
  control = list(adapt_delta = 0.95)
)

# In-sample fitted summaries 
epred_train <- posterior_epred(brms_full, newdata = data, ndraws = 1000)
fitted_train <- tibble(
  row_id   = data$row_id,
  BRMS_pred = colMeans(epred_train),
  BRMS_sd   = apply(epred_train, 2, sd)
)
# posterior_predict
ppc_train <- posterior_predict(brms_full, newdata = data, ndraws = 1000)
fitted_train$q0.025 <- apply(ppc_train, 2, quantile, probs = 0.025)
fitted_train$q0.975 <- apply(ppc_train, 2, quantile, probs = 0.975)

data_fitted <- data %>%
  left_join(fitted_train, by = "row_id")

# readr::write_csv(data_fitted, "../Outputs/ML_cv/BRMS/fitted_on_training.csv")

# LSOA predictions 
lsoa_raw <- readr::read_csv("../Data/lsoa_prediction_table.csv")

time_levels <- sort(unique(data$time_index))
lsoa_pred <- lsoa_raw %>%
  mutate(
    time_index = as.integer(factor(week, levels = time_levels))
  ) %>%
  filter(!is.na(time_index)) 

# Apply training scaling to covariates
lsoa_pred <- lsoa_pred %>%
  mutate(
    pop_density = scale_with(pop_density, train_stats$m_pd,  train_stats$s_pd),
    imd = scale_with(imd, train_stats$m_imd, train_stats$s_imd),
    bame = scale_with(bame, train_stats$m_ba, train_stats$s_ba),
    prop_ind = scale_with(prop_ind, train_stats$m_pi, train_stats$s_pi),
    prop_64 = scale_with(prop_64, train_stats$m_p64, train_stats$s_p64),
    prop_16 = scale_with(prop_16, train_stats$m_p16, train_stats$s_p16)
  )

# Sanity check
stopifnot(all(covariates %in% names(lsoa_pred)))

# Predict
mean_draws <- posterior_epred(brms_full, newdata = lsoa_pred, ndraws = 1000)
pred_draws <- posterior_predict(brms_full, newdata = lsoa_pred, ndraws = 1000)

lsoa_pred$BRMS_pred <- colMeans(mean_draws)
lsoa_pred$BRMS_sd <- apply(mean_draws, 2, sd)
lsoa_pred$q0.025 <- apply(pred_draws, 2, quantile, probs = 0.025)
lsoa_pred$q0.975 <- apply(pred_draws, 2, quantile, probs = 0.975)


readr::write_csv(lsoa_pred, "../Outputs/ML_cv/BRMS/lsoa_predictions_brms.csv")

cat("Done. Files written:\n",
    "- ../Outputs/ML_cv/BRMS/lsoa_predictions_brms.csv (LSOA predictions)\n")


