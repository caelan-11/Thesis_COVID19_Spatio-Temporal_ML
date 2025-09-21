# Libraries 
library(fmesher)
library(INLA)
library(Matrix)
library(Metrics)
library(sf)
library(tibble)
library(dplyr)
library(lubridate)
library(terra)
library(tidyr)
library(raster)
library(stringr)
library(inlabru)
library(zoo)
library(INLAspacetime)
library(blockCV)
library(ggplot2)
library(readr)

options(saveWorkspace = FALSE)
inla.setOption(num.threads = "6")

# Metric helpers
MSE <- function(z, zhat) mean((z - zhat)^2, na.rm = TRUE)
MAE <- function(z, zhat) mean(abs(z - zhat), na.rm = TRUE)
MAPE <- function(z, zhat) mean(abs((zhat - z)/z)[is.finite(z)], na.rm = TRUE) * 100
BIAS <- function(z, zhat) mean(zhat - z, na.rm = TRUE)
pBIAS<- function(z, zhat) mean((zhat - z)/z[is.finite(z)], na.rm = TRUE) * 100
CORR <- function(z, zhat) cor(z, zhat, use = "complete.obs", method = "spearman")
COV <- function(z, lower, upper) mean(z >= lower & z <= upper, na.rm = TRUE) * 100
PMCC <- function(z, zhat) sum((z - zhat)^2, na.rm = TRUE)

# Load & prep data
data <- read_csv("../Data/final_dataset.csv") %>%
  mutate(across(c(pop_density, imd, bame, prop_ind, prop_64, prop_16), ~ scale(.x)[,1])) %>%
  mutate(
    time_index = as.integer(factor(time)),
    logRNA = as.numeric(logRNA),
    site_code = as.factor(Site.code)
  )

train_time_levels <- sort(unique(data$time_index))


# Build sf objects
stations_sf <- st_as_sf(data, coords = c("uwwLongitude", "uwwLatitude"), crs = 4326) %>%
  st_transform(27700)

# Regions
regions_sf <- st_read(
  "../Data/Regions/Regions_December_2021_EN_BUC_2022_-8285331225851041707/RGN_DEC_2021_EN_BUC.shp",
  quiet = TRUE
) %>% st_transform(27700)

stations_sf <- st_join(stations_sf, regions_sf["RGN21NM"]) %>%
  mutate(region = as.factor(RGN21NM))


# Spatial CV folds 
sites <- stations_sf %>% dplyr::select(site_code, geometry) %>% distinct()

bb <- st_bbox(stations_sf)
r_t  <- rast(ext(bb$xmin, bb$xmax, bb$ymin, bb$ymax), crs = "EPSG:27700", resolution = 10000) # 10 km-ish blocks

folds <- blockCV::cv_spatial(x = sites, r = r_t, k = 10, seed = 12, plot = FALSE)

fold_blocks <- st_as_sf(folds$blocks) %>% st_transform(27700)
stations_sf <- st_join(stations_sf, fold_blocks, left = TRUE)

if ("folds" %in% names(stations_sf)) {
  stations_sf <- stations_sf %>% rename(cv_fold = folds)
} else if ("foldID" %in% names(stations_sf)) {
  stations_sf <- stations_sf %>% rename(cv_fold = foldID)
} else {
  stop("Could not find a fold id column after blockCV; check object names.")
}

# Mesh 
max.edge_km <- 30  
bound.outer_km <- diff(range(st_coordinates(stations_sf)[,1])) / 1000 / 5 
n_time <- length(train_time_levels)

# Results containers
results_all <- data.frame()
metrics_all <- data.frame()

# CV loop
for (k in 1:10) {
  cat("Starting fold", k, "at", Sys.time(), "\n")
  start_time <- Sys.time()
  
  train <- dplyr::filter(stations_sf, cv_fold != k)
  val <- dplyr::filter(stations_sf, cv_fold == k)
  
  train$site_code <- as.factor(train$site_code)
  val$site_code <- as.factor(val$site_code)
  train$region <- as.factor(train$region)
  val$region <- as.factor(val$region)
  
  coords_train_km <- st_coordinates(train) / 1000
  coords_val_km <- st_coordinates(val)   / 1000
  
  domain_km <- inla.nonconvex.hull(coords_train_km)
  mesh <- fm_mesh_2d_inla(
    boundary = domain_km,
    loc = coords_train_km,
    max.edge = c(1, 2) * max.edge_km,
    offset = c(max.edge_km, bound.outer_km),
    cutoff = max.edge_km / 5
  )
  
  # SPDE
  spde <- inla.spde2.pcmatern(
    mesh = mesh, alpha = 2,
    prior.range = c(5, 0.05),  
    prior.sigma = c(5, 0.05),  
    constr = TRUE
  )
  
  # A matrices
  A.train <- inla.spde.make.A(mesh, loc = coords_train_km, group = train$time_index, n.group = n_time)
  A.val <- inla.spde.make.A(mesh, loc = coords_val_km,   group = val$time_index,   n.group = n_time)
  
  s.index <- inla.spde.make.index("spatial.field", n.spde = mesh$n, n.group = n_time)
  
  # stacks
  stk.train <- inla.stack(
    data = list(y = train$logRNA),
    A = list(A.train, 1),
    effects = list(
      c(s.index, list(Intercept = 1)),
      data.frame(
        time = train$time_index,
        site_code = train$site_code,
        region = train$region,
        pop_density = train$pop_density,
        prop_64 = train$prop_64,
        prop_16 = train$prop_16,
        imd = train$imd,
        bame = train$bame,
        prop_ind = train$prop_ind
      )
    ),
    tag = "train"
  )
  
  stk.val <- inla.stack(
    data = list(y = NA),
    A = list(A.val, 1),
    effects = list(
      c(s.index, list(Intercept = 1)),
      data.frame(
        time = val$time_index,
        site_code = val$site_code,
        region = val$region,
        pop_density = val$pop_density,
        prop_64 = val$prop_64,
        prop_16 = val$prop_16,
        imd = val$imd,
        bame = val$bame,
        prop_ind = val$prop_ind
      )
    ),
    tag = "val"
  )
  
  stk.full <- inla.stack(stk.train, stk.val)
  
  # model formula
  formula <- y ~ -1 + Intercept +
    pop_density + prop_64 + prop_16 + imd + bame + prop_ind +
    f(region,   model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
    f(site_code, model = "iid", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
    f(time, model = "rw1", hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
    f(spatial.field, model = spde, group = spatial.field.group,
      control.group = list(model = "ar1", hyper = list(theta = list(prior = "pccor0", param = c(0.1, 0.9)))))
  
  # fit model
  fit <- inla(
    formula,
    data = inla.stack.data(stk.full),
    family = "gaussian",
    control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE, link = 1),
    control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    verbose = TRUE
  )

  idx.val <- inla.stack.index(stk.full, tag = "val")$data
  val$mean <- fit$summary.fitted.values[idx.val, "mean"]
  val$sd <- fit$summary.fitted.values[idx.val, "sd"]        
  val$q0.025 <- fit$summary.fitted.values[idx.val, "0.025quant"]
  val$q0.975 <- fit$summary.fitted.values[idx.val, "0.975quant"]
  val$Fold <- k
  
  results_all  <- rbind(results_all, val)
  
  # metrics
  z <- val$logRNA
  zhat <- val$mean
  
  metrics_all <- bind_rows(metrics_all, data.frame(
    Fold = k,
    MAE = MAE(z, zhat),
    MSE = MSE(z, zhat),
    MAPE = MAPE(z, zhat),
    BIAS = BIAS(z, zhat),
    pBIAS= pBIAS(z, zhat),
    CORR = CORR(z, zhat),
    COV = COV(z, lower = val$q0.025, upper = val$q0.975),
    PMCC = PMCC(z, zhat)
  ))
  
  cat("Fold", k, "completed in", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
}


# Save CV predictions & overall metrics 
results_all <- results_all %>% mutate(Bias = logRNA - mean)

overall_inla_metrics <- tibble(
  MAE = mean(abs(results_all$Bias), na.rm = TRUE),
  RMSE = sqrt(mean(results_all$Bias^2, na.rm = TRUE)),
  Bias = mean(results_all$Bias, na.rm = TRUE),
  Correlation = cor(results_all$logRNA, results_all$mean, use = "complete.obs"),
  Coverage = mean(results_all$logRNA >= results_all$q0.025 & results_all$logRNA <= results_all$q0.975, na.rm = TRUE) * 100,
  CI_Amplitude = mean(results_all$q0.975 - results_all$q0.025, na.rm = TRUE),
  PMCC = sum((results_all$logRNA - results_all$mean)^2, na.rm = TRUE)
)

write_csv(results_all, "../Outputs/INLA_cv/cv_predictions_all_folds.csv")
write_csv(metrics_all, "../Outputs/INLA_cv/cv_metrics_summary.csv")
write_csv(overall_inla_metrics, "../Outputs/INLA_cv/inla_metrics_overall.csv")

# Fit on ALL data and predict to LSOA table
coords_all_km <- st_coordinates(stations_sf) / 1000

domain_km <- inla.nonconvex.hull(coords_all_km)
mesh_all <- fm_mesh_2d_inla(
  boundary = domain_km,
  loc = coords_all_km,
  max.edge = c(1, 2) * max.edge_km,
  offset = c(max.edge_km, bound.outer_km),
  cutoff = max.edge_km / 5
)

spde_all <- inla.spde2.pcmatern(
  mesh = mesh_all, alpha = 2,
  prior.range = c(5, 0.05),
  prior.sigma = c(5, 0.05),
  constr = TRUE
)

A.all <- inla.spde.make.A(mesh_all, loc = coords_all_km,
                          group = stations_sf$time_index, n.group = n_time)
s.index.all <- inla.spde.make.index("spatial.field", n.spde = mesh_all$n, n.group = n_time)

stk.all <- inla.stack(
  data = list(y = stations_sf$logRNA),
  A = list(A.all, 1),
  effects = list(
    c(s.index.all, list(Intercept = 1)),
    data.frame(
      time = stations_sf$time_index,
      site_code = stations_sf$site_code,
      region = stations_sf$region,
      pop_density = stations_sf$pop_density,
      prop_64 = stations_sf$prop_64,
      prop_16 = stations_sf$prop_16,
      imd = stations_sf$imd,
      bame = stations_sf$bame,
      prop_ind = stations_sf$prop_ind
    )
  ),
  tag = "train"
)

lsoa.cov <- read_csv("../Data/lsoa_prediction_table.csv") %>%
  mutate(
    geometry = str_remove_all(geometry, "c\\(|\\)"),
    Easting = as.numeric(str_split_fixed(geometry, ",\\s*", 2)[,1]),
    Northing = as.numeric(str_split_fixed(geometry, ",\\s*", 2)[,2])
  )

coords.lsoa_km <- as.matrix(cbind(lsoa.cov$Easting, lsoa.cov$Northing)) / 1000

lsoa.cov$time_index <- as.integer(factor(lsoa.cov$week, levels = train_time_levels))
keep_pred <- !is.na(lsoa.cov$time_index)
lsoa.cov <- lsoa.cov[keep_pred, , drop = FALSE]
coords.lsoa_km <- coords.lsoa_km[keep_pred, , drop = FALSE]

A.pred <- inla.spde.make.A(
  mesh = mesh_all,
  loc = coords.lsoa_km,
  group = lsoa.cov$time_index,
  n.group = n_time
)

s.index.pred <- inla.spde.make.index("spatial.field", n.spde = mesh_all$n, n.group = n_time)

stk.pred <- inla.stack(
  data = list(y = NA),
  A = list(A.pred, 1),
  effects = list(
    c(s.index.pred, list(Intercept = 1)),
    data.frame(
      time = lsoa.cov$time_index,
      pop_density = lsoa.cov$pop_density,
      prop_64 = lsoa.cov$prop_64,
      prop_16 = lsoa.cov$prop_16,
      imd = lsoa.cov$imd,
      bame = lsoa.cov$bame,
      prop_ind  = lsoa.cov$prop_ind
    )
  ),
  tag = "pred"
)

stk.full <- inla.stack(stk.all, stk.pred)

fit_all <- inla(
  formula,
  data = inla.stack.data(stk.full),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(stk.full), compute = TRUE, link = 1),
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  verbose = TRUE
)

idx.pred <- inla.stack.index(stk.full, tag = "pred")$data
lsoa.cov$mean <- fit_all$summary.fitted.values[idx.pred, "mean"]
lsoa.cov$sd <- fit_all$summary.fitted.values[idx.pred, "sd"]
lsoa.cov$q0.025 <- fit_all$summary.fitted.values[idx.pred, "0.025quant"]
lsoa.cov$q0.975 <- fit_all$summary.fitted.values[idx.pred, "0.975quant"]

# save LSOA predictions
write_csv(lsoa.cov, "../Outputs/INLA_predictions/lsoa_predictions.csv")
