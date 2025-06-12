# title: "Gaussian process "
# author: "Shazia Ruybal-Pesántez"

# Load packages -----------------------------------------------------------
library(tidyverse)
library(fields)
library(mvtnorm)
library(epitools)
library(GenSA)
library(lubridate)
library(sf)

# Set the set for reproducibility
set.seed(2)
sf_use_s2(FALSE)  # Critical for planar operations

# Read in data and filter data --------------------------------------------
dat <- read.csv("analysis/data_derived/validated_and_candidate_get_prevalence.csv")
admin0_africa <- readRDS("analysis/data_derived/sf_admin0_africa.rds")  %>%
  st_make_valid()
africa_union <- st_union(admin0_africa)

# Filter for individual mutation
dat <- dat %>%
  filter(mutation == "k13:675:V")

# Transform raw numerator & denominator into "z" value using empirical logit
dat <- dat |>
  mutate(time = as.Date(collection_day, format = "%m/%d/%y"),
         z_est = qlogis((numerator + 0.5) / (denominator + 1.0)))

# Obtain the box from Africa
bbox_africa <- st_bbox(admin0_africa)

# Check how many data points have a denominator less than 20
dim(dat |> filter(denominator > 20))

## INTERMEDIATE SOLUTION: Filter out data points with denominator less than 20
dat <- dat |>
  filter(denominator > 20)

# Plot data --------------------------------------------
# Plot raw data
prevalence_admin0_plot <- ggplot() +
  geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = dat %>% arrange(prevalence),
             aes(x = longitude, y = latitude, fill = prevalence),
             shape = 21, size = 1) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma") +
  theme_bw()

ggsave("analysis/plots/prediction/data_prevalence.png", prevalence_admin0_plot)

# Plot logit transformed data
logit_admin0_plot <- ggplot() +
  geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
  geom_point(data = dat %>% arrange(prevalence),
             aes(x = longitude, y = latitude, fill = z_est),
             shape = 21, size = 1) +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "magma") +
  theme_bw()

ggsave("analysis/plots/prediction/logit_prevalence.png", logit_admin0_plot)

# Gaussian Process --------------------------------------------

# Create regular grid over Africa
grid <- st_make_grid(africa_union,
                     what = "centers",
                     cellsize = 0.5,   # ~50km at equator
                     square = TRUE)
grid_sf <- st_sf(geometry = grid)
# Filter grid to include only points within Africa
grid_africa <- st_filter(grid_sf, africa_union, .predicate = st_within)
# Extract coordinates (lon, lat)
coords <- st_coordinates(grid_africa) |>
  as_tibble() |>
  rename(lon = X, lat = Y)

# Create space-time grid by crossing with years
df_grid_africa <- expand_grid(coords, year = min(dat$year):2024) |>
  mutate(time = as.Date(sprintf("%s-01-01", year)))

## INTERMEDIATE SOLUTION: aggregate duplicates (or near duplicate) observations
duplicates <- dat %>%
  group_by(latitude, longitude, year) %>%
  group_split()
dat <- dat %>%
  group_by(latitude, longitude, year) %>%
  summarise(
    numerator = sum(numerator),
    denominator = sum(denominator),
    prevalence = numerator / denominator,
    .groups = "drop"
  ) %>%
  mutate(
    time = as.Date(sprintf("%s-01-01", year)),
    z_est = qlogis((numerator + 0.5) / (denominator + 1))
  )

# Get distance between observations in space and time
dist_space <- dat |>
  select(latitude, longitude) |>
  dist() |>
  as.matrix()
dist_time <- dist(dat$time) |>
  as.matrix()

# Get distance between observations and prediction grid
crossdist_space <- fields::rdist(x1 = cbind(dat$latitude, dat$longitude),
                                 x2 = cbind(df_grid_africa$lat, df_grid_africa$lon))
crossdist_time <- fields::rdist(x1 = dat$time,
                                x2 = df_grid_africa$time)

# Define the GP covariance function
get_GP_model <- function(dist_space, dist_time,
                         alpha, beta,
                         lambda_space, lambda_time,
                         noise_on = TRUE) {
  # Base covariance: product of spatial & temporal kernels (squared‐exponential)
  ret <- alpha *
    exp(- dist_space^2 / (2 * lambda_space^2)) *
    exp(- dist_time^2  / (2 * lambda_time^2))

  # Add independent (diagonal) noise if requested
  if (noise_on) {
    n <- nrow(dist_space)
    ret <- ret + diag(beta, n)
  }

  ret
}

# Calculate the Negative Log Likelihood
neg_log_likelihood <- function(par, dist_space, dist_time, z_est) {
  alpha <- exp(par[1])  # ensure positivity
  beta  <- exp(par[2])
  lambda_space <- exp(par[3])
  lambda_time  <- exp(par[4])

  # Covariance matrix
  K <- get_GP_model(dist_space, dist_time,
                    alpha = alpha,
                    beta = beta,
                    lambda_space = lambda_space,
                    lambda_time = lambda_time,
                    noise_on = TRUE)

  # Add numerical jitter if K is not positive definite
  jitter <- 1e-6
  diag(K) <- diag(K) + jitter

  # Cholesky decomposition
  L <- tryCatch(chol(K), error = function(e) return(NULL))
  if (is.null(L)) return(Inf)

  n <- length(z_est)
  GP_mean <- qlogis(0.01)
  z_centered <- z_est - GP_mean

  # Log determinant from Cholesky
  log_det_K <- 2 * sum(log(diag(L)))

  # Quadratic term
  alpha_term <- backsolve(L, forwardsolve(t(L), z_centered))
  quad_form <- sum(alpha_term^2)

  # Return negative log-likelihood
  0.5 * (log_det_K + quad_form + n * log(2 * pi))
}


# --- MULTI-RUN HYPERPARAMETER OPTIMIZATION ---------------------------------
n_runs <- 10  # You can increase this later
set.seed(42)  # Ensure reproducibility

results_list <- vector("list", n_runs)

# Define Fixed (Baseline) Values
alpha_fixed <- 0.5
beta_fixed  <- 0.98
lambda_time_fixed <- 880
lambda_space_fixed <- 2.50

# Define Sweep Grids
alpha_vals <- seq(0.01, 2, length.out = 100)
beta_vals <- seq(0.01, 2, length.out = 100)
lambda_time_vals <- seq(100, 2000, length.out = 100)
lambda_space_vals <- seq(0.01, 20, length.out = 100)

sweep_log_likelihood <- function(param_name, param_vals) {
  results <- vector("list", length(param_vals))

  for (i in seq_along(param_vals)) {
    current_val <- param_vals[i]

    # Assign parameters
    alpha_val        <- if (param_name == "alpha") current_val else alpha_fixed
    beta_val         <- if (param_name == "beta") current_val else beta_fixed
    lambda_space_val <- if (param_name == "lambda_space") current_val else lambda_space_fixed
    lambda_time_val  <- if (param_name == "lambda_time") current_val else lambda_time_fixed

    # Log-scale parameters for optimization
    par <- log(c(alpha_val, beta_val, lambda_space_val, lambda_time_val))

    nll <- neg_log_likelihood(par,
                              dist_space = dist_space,
                              dist_time = dist_time,
                              z_est = dat$z_est)

    results[[i]] <- tibble(
      param_name = param_name,
      alpha = alpha_val,
      beta = beta_val,
      lambda_space = lambda_space_val,
      lambda_time = lambda_time_val,
      neg_log_lik = nll
    )
  }

  bind_rows(results)
}

sweep_alpha        <- sweep_log_likelihood("alpha", alpha_vals)
sweep_beta         <- sweep_log_likelihood("beta", beta_vals)
sweep_lambda_space <- sweep_log_likelihood("lambda_space", lambda_space_vals)
sweep_lambda_time  <- sweep_log_likelihood("lambda_time", lambda_time_vals)

write.csv(sweep_alpha, "analysis/plots/prediction/sweep_alpha.csv")
write.csv(sweep_beta, "analysis/plots/prediction/sweep_beta.csv")
write.csv(sweep_lambda_space, "analysis/plots/prediction/sweep_lambda_space.csv")
write.csv(sweep_lambda_time, "analysis/plots/prediction/sweep_lambda_time.csv")

# Run and store best results from each sweep
best_alpha <- sweep_alpha %>% filter(neg_log_lik == min(neg_log_lik))
best_beta <- sweep_beta %>% filter(neg_log_lik == min(neg_log_lik))
best_lambda_space <- sweep_lambda_space %>% filter(neg_log_lik == min(neg_log_lik))
best_lambda_time <- sweep_lambda_time %>% filter(neg_log_lik == min(neg_log_lik))

plot_prediction <- function(best_run){
  # Extract best parameters
  alpha_opt <- best_run$alpha
  beta_opt <- best_run$beta
  lambda_space_opt <- best_run$lambda_space
  lambda_time_opt <- best_run$lambda_time

  cat("Optimized params:\n")
  cat("alpha =", alpha_opt, "\n")
  cat("beta =", beta_opt, "\n")
  cat("lambda_space =", lambda_space_opt, "\n")
  cat("lambda_time =", lambda_time_opt, "\n")

  # Compute the full covariance matrix K (with noise) for observed sites
  K <- get_GP_model(
    dist_space,        # pairwise spatial distances among data sites
    dist_time,         # pairwise temporal distances among data sites
    alpha_opt, beta_opt,
    lambda_space_opt, lambda_time_opt
  )

  # Cholesky factorization of K: K = L %*% t(L)
  # This will be used for efficient solves and draws
  L <- chol(K)

  # Compute the cross‐covariance between prediction grid and observed sites,
  # but without adding noise on the diagonal (we want pure covariance)
  K_pred <- get_GP_model(
    crossdist_space,   # spatial distances between grid points & data sites
    crossdist_time,    # temporal distances between grid points & data sites
    alpha_opt, beta_opt,
    lambda_space_opt, lambda_time_opt,
    noise_on = FALSE   # no noise term for predictive covariances
  )

  # Perform GP prediction on empirical logit of observed prevalence. Predictions
  # will be on real line, so transform these back to [0,1] range
  GP_mean <- mean(dat$z_est)
  df_pred <- df_grid_africa |>
    mutate(z_pred = t(K_pred) %*% backsolve(L, forwardsolve(t(L), dat$z_est - GP_mean)),
           z_pred = as.vector(z_pred) + GP_mean,
           p_pred = plogis(z_pred))

  df_pred <- df_pred |>
    mutate(
      prevalence = plogis(z_pred),
      prevalence_bin = cut(
        prevalence,
        breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1),
        labels = c("<1%", "1–5%", "5–10%", "10–25%", "25–50%", "50–75%", ">75%"),
        include.lowest = TRUE,
        right = FALSE
      )
    )

  # Predictive variance
  # Solve: v = L \ (K_pred)
  v <- forwardsolve(t(L), K_pred) # Lᵗ⁻¹ K_pred, dimensions: n x m
  pred_var <- alpha_opt - colSums(v^2)

  # Sanity check: variance must be >= 0
  pred_var[pred_var < 0] <- 0

  # Standard deviation on logit scale
  df_pred$z_sd <- sqrt(pred_var)

  df_pred <- df_pred %>%
    mutate(
      lower = plogis(z_pred - 1.96 * z_sd),
      upper = plogis(z_pred + 1.96 * z_sd)
    )

  # Convert lambda_time and lambda_space into years and km respectively
  best_run <- best_run |>
    mutate(lambda_time_year = lambda_time/365, # Convert days into years
           lambda_space_km = lambda_space*111) # Convert degrees into km using the approximation that 1 degree is about 111km

  # Plot prediction
  plot_predictions <- df_pred |>
    ggplot() + theme_bw() +
    geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
    geom_raster(aes(x = lon, y = lat, fill = plogis(z_pred))) +
    # geom_point(aes(x = longitude, y = latitude, fill = plogis(z_est), size = denominator), shape = 21,
    #            colour = grey(0.4), data = mutate(dat, time = time)) +
    facet_wrap(~year) +
    scale_fill_viridis_c(option = "magma", limit = c(0, NA)) +
    scale_size(range = c(1, 3)) +
    ggtitle(paste0("Estimated prevalence surface: alpha=", alpha_opt, " beta=", beta_opt, "\nlambdaSpace=", lambda_space_opt, " lambdaTime=", lambda_time_opt))

  prediction_plot_name = paste0("analysis/plots/prediction/test_prediction_plot_alpha", alpha_opt, "_beta", beta_opt, "_lambdaSpace", lambda_space_opt, "_lambdaTime", lambda_time_opt, ".png")
  ggsave(prediction_plot_name, plot_predictions)

  grouped_pred <- ggplot(df_pred) +
    theme_bw() +
    geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
    geom_raster(aes(x = lon, y = lat, fill = prevalence_bin)) +
    facet_wrap(~year) +
    scale_fill_viridis_d(option = "magma", direction = -1, name = "Prevalence") +
    ggtitle(paste0("Estimated prevalence surface: alpha=", alpha_opt, " beta=", beta_opt, "\nlambdaSpace=", lambda_space_opt, " lambdaTime=", lambda_time_opt))
  grouped_pred_name = paste0("analysis/plots/prediction/test_group_prediction_plot_alpha", alpha_opt, "_beta", beta_opt, "_lambdaSpace", lambda_space_opt, "_lambdaTime", lambda_time_opt, ".png")
  ggsave(grouped_pred_name, grouped_pred)

  confidence_plot <- ggplot(df_pred) +
    geom_sf(data = admin0_africa, fill = NA, color = "black", linewidth = 0.3) +
    geom_raster(aes(x = lon, y = lat, fill = z_sd)) +
    facet_wrap(~year) +
    scale_fill_viridis_c(option = "magma", name = "SD (logit)") +
  ggtitle(paste0("Posterior uncertainty (standard deviation on logit scale): alpha=", alpha_opt, " beta=", beta_opt, "\nlambdaSpace=", lambda_space_opt, " lambdaTime=", lambda_time_opt))
  confidence_name = paste0("analysis/plots/prediction/test_confidence_prediction_alpha", alpha_opt, "_beta", beta_opt, "_lambdaSpace", lambda_space_opt, "_lambdaTime", lambda_time_opt, ".png")
  ggsave(confidence_name, confidence_plot)
}

plot_prediction(best_alpha)
plot_prediction(best_beta)
plot_prediction(best_lambda_space)
plot_prediction(best_lambda_time)


