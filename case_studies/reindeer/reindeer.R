# =============================================================================
# Reindeer case study - Langevin BBIS
# =============================================================================
#   - import data + standardised covariates
#   - work in km for numerical stability
#   - single Langevin BBIS fit
#   - dt_max x M study (repeated refits over a grid of delta_max and M)
#   - summary plot of parameter estimates vs delta_max, by M
#
# Covariates come from reindeer_spatial_covariates.R (run that first).
#
# Individual subset: set `subset_ids` to restrict the fit to a subset of
# reindeer (see the SUBSET block below). NULL uses all individuals.
# =============================================================================

# prep workspace ####
source(here::here("functions/utility_functions.R"))
sourceDir("functions")
load_lib(here, dplyr, tidyr, mvnfast, parallel, terra, ggplot2, viridis,
         RColorBrewer, sf, oneimpact)

# ---- fitting / study settings -----------------------------------------------
ncores  <- 10
M       <- 25          # bridges for the single fit
dt_max  <- 3           # hours (reindeer fixes are ~3-hourly); single-fit value

# --------------------------------------------------------------------------- #
# SUBSET: choose which reindeer individuals to use                            #
# --------------------------------------------------------------------------- #
# subset_ids controls which individuals enter the fit. It can be either:
#   - NULL                -> use all individuals (default)
#   - a vector of ids     -> keep those `original_animal_id`s,
#                            e.g. c(3358, 3361, 3364)
#   - a single integer n  -> randomly sample n individuals
# Bridges are only generated within an `animal_year_id`, so tracks are grouped
# by animal-year regardless of the subset.
subset_ids <- NULL
set.seed(1)            # for reproducible random subsetting

# import data ####
data("reindeer")
reindeer <- sf::st_as_sf(reindeer, coords = c("x", "y"), crs = 25833)

# apply the individual subset
if (!is.null(subset_ids)) {
  ids_all <- sort(unique(reindeer$original_animal_id))
  if (length(subset_ids) == 1 && is.numeric(subset_ids) &&
      subset_ids == round(subset_ids) && !(subset_ids %in% ids_all)) {
    keep_ids <- sample(ids_all, min(subset_ids, length(ids_all)))
  } else {
    keep_ids <- subset_ids
  }
  reindeer <- reindeer[reindeer$original_animal_id %in% keep_ids, ]
  message("Using ", length(unique(reindeer$original_animal_id)),
          " individual(s): ", paste(sort(unique(reindeer$original_animal_id)),
                                     collapse = ", "))
}

# import standardised covariates (built by reindeer_spatial_covariates.R)
hbfull <- rast(here("case_studies/reindeer/data/reindeer_covariates.tif"))

# tracks: project to covariate CRS; ID = animal_year_id (bridge unit)
tracks <- reindeer |>
  sf::st_transform(crs(hbfull)) |>
  terra::vect()
tracks$time <- reindeer$t
tracks$ID   <- reindeer$animal_year_id

# Change the resolution and extent from m to km (numerical stability)
crs_km <- gsub("units=m", "units=km", crs(hbfull, proj = TRUE))
r <- rast(nrows = nrow(hbfull), ncols = ncol(hbfull),
          ext   = as.vector(ext(hbfull)) / 1000,
          crs   = crs_km)                 # define template raster

# standardise projections (to km grid)
hbfull <- project(hbfull, r)              # transform raster
tracks <- project(tracks, r) |>           # transform tracks
  as.data.frame(geom = "XY") |> 
  dplyr::arrange(ID, time) |>                 # order for within-ID bridging
  mutate(step = )


# ---- add polynomial (quadratic) terms ---------------------------------------
# Quadratic terms let selection for elevation, temperature and precipitation
# bend, so the fit can reveal an optimum rather than only monotone selection.
# Squares are taken on the standardised (z-scored, km-grid) layers, so each
# covariate and its square are exactly consistent with the fitted values.
poly_covs <- c("elevation", "summer_temperature", "summer_precipitation")
hbfull    <- add_poly_terms(hbfull, poly_covs, degree = 2)

#### fit Langevin BBIS - single fit ####
out <- fit_langevin_bbis(tracks, hbfull,
                         M = M,
                         dt_max = dt_max,
                         dt_units = "hours",
                         ncores = ncores,
                         fixed_sampling = FALSE)

# cache the single fit so reindeer_method_comparison.R can reuse it
make_path(here("case_studies/reindeer/fitted_estimates"))
saveRDS(out, here("case_studies/reindeer/fitted_estimates/reindeer_bbis_single_fit.rds"))

#### summary plot of estimated coefficients (single fit) ####
# bar plot of the habitat-selection coefficients from the single fit:
# covariates as categories on x, estimates as bar heights
coef_df <- data.frame(
  covariate = factor(names(hbfull), levels = names(hbfull)),
  estimate  = out$beta
)

coef_plot <- ggplot(coef_df, aes(x = covariate, y = estimate,
                                 fill = estimate > 0)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#d95f02")) +
  labs(x = NULL, y = "Estimated coefficient",
       title = sprintf("Langevin BBIS coefficients (M = %s, dt_max = %s h)",
                        M, dt_max)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(here("case_studies/reindeer", "reindeer_coef_estimates.png"), coef_plot,
       width = 7, height = 5)

#### predicted UD from estimated coefficients (single fit) ####
# The Langevin utilisation distribution is proportional to exp(covariates %*% beta).
# Build it from the (standardised, km-grid) starting rasters used in the fit.
lin_pred  <- sum(hbfull * out$beta)                          # linear predictor
ud        <- exp(lin_pred)
ud        <- ud / terra::global(ud, "sum", na.rm = TRUE)[[1]] # normalise to sum 1
names(ud) <- "UD"

# save the predicted UD raster
writeRaster(ud, here("case_studies/reindeer", "reindeer_predicted_UD.tif"),
            overwrite = TRUE)

# map of the predicted UD
ud_plot <- ggplot() +
  tidyterra::geom_spatraster(data = ud) +
  scale_fill_viridis_c(option = "inferno", na.value = "transparent", name = "UD") +
  coord_sf(expand = FALSE) +
  labs(title = "Predicted utilisation distribution (single fit)",
       x = NULL, y = NULL) +
  theme_minimal()

ggsave(here("case_studies/reindeer", "reindeer_predicted_UD.png"), ud_plot,
       width = 7, height = 6)

#### response curves for polynomial covariates (single fit) ####
# Marginal relative selection for each quadratic covariate: how exp(linear
# predictor) varies as that covariate moves across its observed (standardised)
# range while the others are held at their mean (0 on the z-score scale).
# A negative squared coefficient gives an interior optimum at -b1 / (2 b2).
beta_named <- setNames(out$beta, names(hbfull))

rng_of <- function(nm) c(terra::global(hbfull[[nm]], "min", na.rm = TRUE)[[1]],
                         terra::global(hbfull[[nm]], "max", na.rm = TRUE)[[1]])

resp_df <- dplyr::bind_rows(lapply(poly_covs, function(nm) {
  rng <- rng_of(nm)
  xs  <- seq(rng[1], rng[2], length.out = 200)
  b1  <- beta_named[[nm]]
  b2  <- beta_named[[paste0(nm, "_sq")]]
  lp  <- b1 * xs + b2 * xs^2
  data.frame(covariate = nm, x = xs, rel_sel = exp(lp - max(lp)))
}))

opt_df <- dplyr::bind_rows(lapply(poly_covs, function(nm) {
  rng   <- rng_of(nm)
  b1    <- beta_named[[nm]]; b2 <- beta_named[[paste0(nm, "_sq")]]
  xstar <- -b1 / (2 * b2)
  data.frame(covariate = nm, xstar = xstar,
             has_opt = is.finite(xstar) && b2 < 0 &&
                       xstar >= rng[1] && xstar <= rng[2])
}))

resp_plot <- ggplot(resp_df, aes(x = x, y = rel_sel)) +
  geom_line(linewidth = 0.8, colour = "#1b9e77") +
  geom_vline(data = subset(opt_df, has_opt), aes(xintercept = xstar),
             linetype = 2, colour = "grey40") +
  facet_wrap(~ covariate, scales = "free", ncol = 3) +
  labs(x = "Standardised covariate value (z-score)",
       y = "Relative selection  exp(linear predictor), peak = 1",
       title = "Langevin BBIS response curves (quadratic terms)") +
  theme_bw()

ggsave(here("case_studies/reindeer", "reindeer_response_curves.png"),
       resp_plot, width = 10, height = 4)

#### fit Langevin BBIS - dt_max & M refits ####
# define fitting criteria
ncores  <- 10
Ms      <- c(25, 50, 100)
deltas  <- exp(seq(log(0.5), log(24), length.out = 30))  # delta_max grid (hours)
nrefits <- 10

# number of pars
npar <- nlyr(hbfull) + 1
# add 1 column to track dt_max
params <- matrix(NA, ncol = npar + 1, nrow = nrefits * length(deltas))

for (M in Ms) {                       # for each number of bridges
  for (k in seq_along(deltas)) {      # for each delta_max
    for (i in 1:nrefits) {
      delta_max <- deltas[k]

      print(sprintf("Fitting M = %s, delta_max = %.4f, refit = %s",
                    M, delta_max, i))
      # fit
      out <- fit_langevin_bbis(tracks, hbfull,
                               M = M,
                               dt_max = delta_max,
                               dt_units = "hours",
                               ncores = ncores,
                               fixed_sampling = FALSE)

      # store outputs (par + delta_max)
      params[(k - 1L) * nrefits + i, ] <- c(out$par, delta_max)
    }
  }
  # convert to data.frame, name it dfM (df25 / df50 / df100), and save
  df <- as.data.frame(params) |>
    setNames(c(paste0("beta", seq_len(npar - 1L)), "sigma", "delta_max"))
  assign(sprintf("df%s", M), df)
  save(list = sprintf("df%s", M),
       file = sprintf("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=%s.Rda",
                      M))
}

#### generate summary plots ####
# import estimates
load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=25.Rda")
load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=50.Rda")
load("case_studies/reindeer/fitted_estimates/reindeer_deltamax_studyM=100.Rda")

# parameter columns (all betas + sigma), labelled for facet strips
beta_names <- paste0("beta", seq_len(nlyr(hbfull)))
cov_labels <- names(hbfull)
par_levels <- c(beta_names, "sigma")
par_labels <- c(setNames(sprintf("beta[%d]~(%s)", seq_along(cov_labels), cov_labels),
                         beta_names),
                sigma = "sigma")

# combine all data
df_all <- bind_rows(mutate(df25,  M = "M=25"),
                    mutate(df50,  M = "M=50"),
                    mutate(df100, M = "M=100")) |>
  pivot_longer(cols = all_of(par_levels),
               names_to = "par", values_to = "mu") |>
  mutate(par = factor(par, levels = par_levels))

# summarise estimates (median, sd, & confidence intervals)
z <- qnorm(0.975)
sum_all <- df_all |>
  dplyr::group_by(par, delta_max, M) |>
  dplyr::summarise(sd = sd(mu),
                   mu = median(mu),
                   .groups = "drop") |>
  dplyr::mutate(lo = mu - z * sd,
                hi = mu + z * sd)

# generate plot
plot <- ggplot(sum_all, aes(x = delta_max, y = mu,
                            color = factor(M, levels = c("M=25", "M=50", "M=100")))) +
  # BBIS estimates
  geom_point(data = df_all, alpha = 0.15, stroke = NA) +
  geom_line(aes(linetype = factor(M, levels = c("M=25", "M=50", "M=100"))),
            linewidth = 0.7) +
  # design
  facet_wrap(~ par, scales = "free",
             labeller = labeller(par = as_labeller(par_labels, label_parsed))) +
  scale_x_log10() +
  scale_color_brewer(palette = "Dark2") +
  labs(x = expression(Delta[max]), y = expression(Estimate),
       color = NULL, linetype = NULL) +
  theme_bw()

# save plot
ggsave(here("case_studies/reindeer", "reindeer_par_est_dmax_M.png"), plot,
       width = 9, height = 6)
