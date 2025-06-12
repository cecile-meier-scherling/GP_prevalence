library(tidyverse)
library(fields)
library(mvtnorm)
library(epitools)
library(lubridate)

set.seed(2)

# Read in data and filter data --------------------------------------------
dat <- read.csv("analysis/data_derived/validated_and_candidate_get_prevalence.csv")

# Filter for Uganda and Ethiopia
dat <- dat %>%
  filter(country_name %in% c("Uganda", "Ethiopia", "Rwanda", "Kenya", "Tanzania", "Democratic Republic of the Congo"),
         mutation == "k13:469:Y")

# transform raw numerator & denominator into "z" value using empirical logit
dat <- dat |>
  mutate(time = as.Date(collection_day, format = "%m/%d/%y"),
         z_est = qlogis((numerator + 0.5) / (denominator + 1.0)),
         logit_z = plogis(z_est))

# plot
dat |>
  arrange(-prevalence) |>
  ggplot() + theme_bw() +
  geom_point(aes(x = longitude, y = latitude, fill = plogis(z_est)),
             shape = 21, size = 2) +
  facet_wrap(~year(time)) +
  scale_fill_viridis_c(option = "magma")

# set space and time bounds
lat_min <- min(dat$latitude)
lat_max <- max(dat$latitude)
lon_min <- min(dat$longitude)
lon_max <- max(dat$longitude)

# create a data.frame spanning the full space-time grid
n_grid_space <- 41
df_grid <- expand_grid(lat = seq(lat_min, lat_max, l = n_grid_space),
                       lon = seq(lon_min, lon_max, l = n_grid_space),
                       year = min(dat$year):2024) |>
  mutate(time = as.Date(sprintf("%s-01-01", year)))

# get distance between observations in space and time
dist_space <- dat |>
  select(latitude, longitude) |>
  dist() |>
  as.matrix()

dist_time <- dist(dat$time) |>
  as.matrix()

# get distance between observations and prediction grid
crossdist_space <- fields::rdist(x1 = cbind(dat$latitude, dat$longitude),
                                 x2 = cbind(df_grid$lat, df_grid$lon))

crossdist_time <- fields::rdist(x1 = dat$time,
                                x2 = df_grid$time)

# Parameters of the Gaussian Process (GP) model
alpha <- 3           # marginal variance (controls overall scale of covariance)
beta  <- 0.2           # noise variance (nugget effect)
lambda_space <- 1  # spatial length‐scale (how quickly spatial correlation decays)
lambda_time  <- 200^2  # temporal length‐scale (how quickly temporal correlation decays)

# Define the GP covariance function
get_GP_model <- function(dist_space, dist_time,
                         alpha, beta,
                         lambda_space, lambda_time,
                         noise_on = TRUE) {
  # Base covariance: product of spatial & temporal kernels (squared‐exponential)
  ret <- alpha *
    exp(- dist_space^2 / (2 * lambda_space)) *
    exp(- dist_time^2  / (2 * lambda_time))

  # Add independent (diagonal) noise if requested
  if (noise_on) {
    n <- nrow(dist_space)
     ret <- ret + diag(beta, n)
  }

  ret
}

# Compute the full covariance matrix K (with noise) for observed sites
K <- get_GP_model(
  dist_space,        # pairwise spatial distances among data sites
  dist_time,         # pairwise temporal distances among data sites
  alpha, beta,
  lambda_space, lambda_time
)

# Cholesky factorization of K: K = L %*% t(L)
# This will be used for efficient solves and draws
L <- chol(K)

# Compute the cross‐covariance between prediction grid and observed sites,
# but without adding noise on the diagonal (we want pure covariance)
K_pred <- get_GP_model(
  crossdist_space,   # spatial distances between grid points & data sites
  crossdist_time,    # temporal distances between grid points & data sites
  alpha, beta,
  lambda_space, lambda_time,
  noise_on = FALSE   # no noise term for predictive covariances
)

# Perform GP prediction on empirical logit of observed prevalence. Predictions
# will be on real line, so transform these back to [0,1] range
GP_mean <- mean(dat$z_est)
df_pred <- df_grid |>
  mutate(z_pred = t(K_pred) %*% backsolve(L, forwardsolve(t(L), dat$z_est - GP_mean)),
         z_pred = as.vector(z_pred) + GP_mean,
         p_pred = plogis(z_pred))

# plot prediction
df_pred |>
  ggplot() + theme_bw() +
  geom_raster(aes(x = lon, y = lat, fill = plogis(z_pred))) +
  geom_point(aes(x = longitude, y = latitude, fill = plogis(z_est), size = denominator), shape = 21,
             colour = grey(0.4), data = mutate(dat, time = collection_day)) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma", limit = c(0, NA)) +
  scale_size(range = c(1, 3)) +
  ggtitle("Estimated prevalence surface")

