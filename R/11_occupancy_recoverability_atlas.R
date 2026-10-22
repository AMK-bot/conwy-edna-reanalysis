#!/usr/bin/env Rscript
# =============================================================================
# Phase-5 occupancy recoverability atlas
#
# Synthetic known-truth analysis only.  This script deliberately reads exactly
# two scientific inputs: the frozen preregistration and its machine-readable
# manifest.  It never reads a Conwy count table, observation mask, Layer-A fit,
# or other real-data object.  Every portable result row carries the scope tag
# `synthetic_known_truth_only`.
#
# Run from the repository root:
#   Rscript R/11_occupancy_recoverability_atlas.R
#   Rscript R/11_occupancy_recoverability_atlas.R --verify-rebuild
#
# Optional execution control (does not change the frozen design):
#   IMPACT_ATLAS_CORES=4 Rscript R/11_occupancy_recoverability_atlas.R
# =============================================================================

ATLAS_SCOPE <- "synthetic_known_truth_only"
ATLAS_ESTIMATORS <- c("naive_occurrence", "replication_aware", "depth_aware")
ATLAS_PARAMETER_BOUND <- 12
ATLAS_BOUNDARY_TOLERANCE <- 1e-5
ATLAS_INFORMATION_RELATIVE_TOLERANCE <- 1e-8
ATLAS_INFORMATION_ABSOLUTE_TOLERANCE <- 1e-8
ATLAS_INFORMATION_MAX_CONDITION <- 1e8
ATLAS_OPTIM_MAXIT <- 100L
ATLAS_OPTIM_PGTOL <- 1e-7
ATLAS_EXPECTED_PREREG_SHA256 <-
  "5d4ab9ce7d1f0dc695646f6e3c781b4931fad5fecc1513ddd79046c3b316e567"
ATLAS_EXPECTED_MANIFEST_SHA256 <-
  "0f6e6464ec54492cb23b8adbdc81b3b05a9ec98e73e1b04b413c764c03194843"

atlas_project_root <- function() {
  candidates <- c(normalizePath(".", mustWork = FALSE))
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    script <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
    candidates <- c(candidates, dirname(dirname(script)))
  }
  hit <- candidates[file.exists(file.path(candidates,
                                           "data/metadata/impact_sprint_manifest.json"))]
  if (!length(hit)) stop("Run the occupancy atlas from the repository root")
  normalizePath(hit[[1]])
}

ATLAS_ROOT <- atlas_project_root()
.atlas_input_audit <- new.env(parent = emptyenv())
.atlas_input_audit$scientific <- character()
.atlas_input_audit$generated <- character()

atlas_path <- function(relative_path) file.path(ATLAS_ROOT, relative_path)

atlas_record_input <- function(relative_path, type = c("scientific", "generated")) {
  type <- match.arg(type)
  normal <- gsub("\\\\", "/", relative_path)
  .atlas_input_audit[[type]] <- unique(c(.atlas_input_audit[[type]], normal))
  invisible(normal)
}

atlas_input_log <- function() {
  list(scientific_inputs_read = sort(.atlas_input_audit$scientific),
       generated_artifacts_read = sort(.atlas_input_audit$generated))
}

atlas_reset_input_log <- function() {
  .atlas_input_audit$scientific <- character()
  .atlas_input_audit$generated <- character()
  invisible(NULL)
}

atlas_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE))
    stop("The digest package is required for SHA-256 provenance")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

load_atlas_config <- function() {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("The jsonlite package is required")
  prereg_rel <- "docs/impact-sprint-preregistration.md"
  manifest_rel <- "data/metadata/impact_sprint_manifest.json"
  atlas_record_input(prereg_rel, "scientific")
  atlas_record_input(manifest_rel, "scientific")
  prereg_hash <- atlas_sha256(atlas_path(prereg_rel))
  manifest_hash <- atlas_sha256(atlas_path(manifest_rel))
  if (!identical(prereg_hash, ATLAS_EXPECTED_PREREG_SHA256))
    stop("Frozen preregistration hash mismatch; refusing to run")
  if (!identical(manifest_hash, ATLAS_EXPECTED_MANIFEST_SHA256))
    stop("Frozen impact manifest hash mismatch; refusing to run")
  manifest <- jsonlite::fromJSON(atlas_path(manifest_rel), simplifyVector = TRUE)
  occ <- manifest$occupancy_atlas
  stopifnot(
    identical(as.integer(occ$n_units), 120L),
    identical(as.integer(occ$simulations_per_scenario), 100L),
    identical(as.integer(occ$expected_scenarios), 648L),
    identical(as.character(occ$estimators), ATLAS_ESTIMATORS)
  )
  list(manifest = manifest, occ = occ,
       seed = as.integer(manifest$seeds$occupancy_atlas),
       prereg_sha256 = prereg_hash, manifest_sha256 = manifest_hash,
       prereg_path = prereg_rel, manifest_path = manifest_rel)
}

calibrate_occupancy_intercept <- function(target, river_position, slope) {
  if (!is.finite(target) || target <= 0 || target >= 1)
    stop("target occupancy must lie strictly between zero and one")
  stats::uniroot(function(intercept) {
    mean(stats::plogis(intercept + slope * river_position)) - target
  }, interval = c(-30, 30), tol = 1e-13)$root
}

build_atlas_scenarios <- function(config) {
  occ <- config$occ
  core <- expand.grid(
    target_mean_occupancy = as.numeric(occ$target_mean_occupancy),
    prevalence_structure = as.character(occ$prevalence_structure),
    baseline_detection = as.numeric(occ$baseline_detection),
    replicates = as.integer(occ$replicates),
    log_depth_sd = as.numeric(occ$log_depth_sd),
    depth_river_correlation = as.numeric(occ$depth_river_correlation),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  core$truth_scenario_id <- seq_len(nrow(core))
  river_position <- seq(-0.5, 0.5, length.out = as.integer(occ$n_units))
  core$occupancy_gradient_slope <- ifelse(
    core$prevalence_structure == "downstream_gradient",
    as.numeric(occ$occupancy_gradient_slope), 0)
  core$occupancy_intercept <- mapply(
    calibrate_occupancy_intercept,
    core$target_mean_occupancy,
    MoreArgs = list(river_position = river_position),
    core$occupancy_gradient_slope)
  core$calibrated_mean_occupancy <- mapply(function(a, b) {
    mean(stats::plogis(a + b * river_position))
  }, core$occupancy_intercept, core$occupancy_gradient_slope)

  mask_levels <- as.character(occ$mask_scenario)
  scenarios <- core[rep(seq_len(nrow(core)), each = length(mask_levels)), , drop = FALSE]
  scenarios$mask_scenario <- rep(mask_levels, times = nrow(core))
  scenarios$scenario_id <- seq_len(nrow(scenarios))
  scenarios$data_scope <- ATLAS_SCOPE
  scenarios <- scenarios[, c(
    "data_scope", "scenario_id", "truth_scenario_id",
    "target_mean_occupancy", "prevalence_structure",
    "occupancy_intercept", "occupancy_gradient_slope",
    "calibrated_mean_occupancy", "baseline_detection", "replicates",
    "log_depth_sd", "depth_river_correlation", "mask_scenario")]
  rownames(scenarios) <- NULL
  if (nrow(scenarios) != as.integer(occ$expected_scenarios))
    stop("Scenario construction does not match the frozen 648-scenario design")
  scenarios
}

atlas_simulation_seed <- function(base_seed, truth_scenario_id, simulation_id) {
  # Mask variants share a seed and therefore differ only in mask encoding.
  seed <- as.double(base_seed) + as.double(truth_scenario_id) * 1000 + simulation_id
  if (seed > .Machine$integer.max) stop("Derived seed exceeds R integer range")
  as.integer(seed)
}

set_atlas_rng <- function(seed) {
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion",
          sample.kind = "Rejection")
  suppressWarnings(RNGversion("4.0.0"))
  set.seed(seed)
}

simulate_atlas_truth <- function(core_row, simulation_id, config) {
  occ <- config$occ
  n_units <- as.integer(occ$n_units)
  k <- as.integer(core_row$replicates)
  seed <- atlas_simulation_seed(config$seed, core_row$truth_scenario_id,
                                simulation_id)
  set_atlas_rng(seed)

  river_position <- seq(-0.5, 0.5, length.out = n_units)
  psi <- stats::plogis(core_row$occupancy_intercept +
                         core_row$occupancy_gradient_slope * river_position)
  latent_occupancy <- stats::rbinom(n_units, 1L, psi)

  # Construct an exactly standardised, exactly orthogonal residual.  The sample
  # SD and correlation of planned log depth therefore equal the declared levels
  # (to floating-point tolerance), rather than merely doing so in expectation.
  river_cell <- rep(river_position, times = k)
  river_cell_z <- as.numeric(scale(river_cell))
  epsilon <- stats::rnorm(n_units * k)
  epsilon <- epsilon - mean(epsilon)
  epsilon <- epsilon - river_cell_z *
    sum(epsilon * river_cell_z) / sum(river_cell_z^2)
  epsilon <- epsilon / stats::sd(epsilon)
  rho <- as.numeric(core_row$depth_river_correlation)
  log_depth <- as.numeric(core_row$log_depth_sd) *
    (rho * river_cell_z + sqrt(1 - rho^2) * epsilon)
  log_depth_matrix <- matrix(log_depth, nrow = n_units, ncol = k)

  detection_probability <- stats::plogis(
    stats::qlogis(as.numeric(core_row$baseline_detection)) +
      as.numeric(occ$detection_depth_slope) * log_depth_matrix)
  generated_detection <- matrix(stats::rbinom(
    n_units * k, 1L,
    rep(latent_occupancy, times = k) * as.vector(detection_probability)),
    nrow = n_units, ncol = k)

  n_missing <- as.integer(round(as.numeric(occ$structured_missing_fraction) *
                                  n_units * k))
  unavailable <- matrix(FALSE, nrow = n_units, ncol = k)
  unavailable[order(as.vector(log_depth_matrix),
                    seq_len(n_units * k))[seq_len(n_missing)]] <- TRUE

  list(
    data_scope = ATLAS_SCOPE,
    truth_scenario_id = as.integer(core_row$truth_scenario_id),
    simulation_id = as.integer(simulation_id), simulation_seed = seed,
    river_position = river_position,
    true_occupancy_probability = psi,
    latent_occupancy = latent_occupancy,
    planned_log_depth = log_depth_matrix,
    detection_probability = detection_probability,
    generated_detection = generated_detection,
    depth_linked_unavailable = unavailable
  )
}

apply_atlas_mask <- function(simulation, mask_scenario) {
  unavailable <- simulation$depth_linked_unavailable
  surveyed_mask <- if (mask_scenario == "complete_correct") {
    matrix(TRUE, nrow = nrow(unavailable), ncol = ncol(unavailable))
  } else {
    !unavailable
  }
  analysis_mask <- if (mask_scenario == "structured_zero_fill") {
    matrix(TRUE, nrow = nrow(unavailable), ncol = ncol(unavailable))
  } else {
    surveyed_mask
  }
  observed_detection <- simulation$generated_detection
  observed_detection[!surveyed_mask] <- NA_integer_
  analysis_detection <- observed_detection
  if (mask_scenario == "structured_zero_fill")
    analysis_detection[!surveyed_mask] <- 0L
  list(
    surveyed_mask = surveyed_mask,
    analysis_mask = analysis_mask,
    observed_detection = observed_detection,
    analysis_detection = analysis_detection,
    erroneous_zero_fill = analysis_mask & !surveyed_mask
  )
}

wilson_interval <- function(successes, trials, level = 0.95) {
  if (trials <= 0) return(c(NA_real_, NA_real_))
  p <- successes / trials
  z <- stats::qnorm(1 - (1 - level) / 2)
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half <- z / denominator * sqrt(p * (1 - p) / trials +
                                   z^2 / (4 * trials^2))
  interval <- pmax(0, pmin(1, c(centre - half, centre + half)))
  # Remove floating residue at the exact binomial endpoints (for example,
  # 1.7e-18 instead of zero), which would otherwise mislabel a valid Wilson
  # interval as invalid when the point estimate is exactly zero or one.
  if (successes == 0) interval[[1]] <- 0
  if (successes == trials) interval[[2]] <- 1
  interval
}

fit_naive_occurrence <- function(masked) {
  analysable <- rowSums(masked$analysis_mask) > 0
  detected <- rowSums(masked$analysis_detection == 1L, na.rm = TRUE) > 0
  n <- sum(analysable)
  d <- sum(detected & analysable)
  estimate <- if (n) d / n else NA_real_
  interval <- wilson_interval(d, n)
  valid_interval <- n > 0 && all(is.finite(interval)) && interval[[1]] <= estimate &&
    estimate <= interval[[2]] && interval[[1]] < interval[[2]]
  list(estimate = estimate, lower95 = interval[[1]], upper95 = interval[[2]],
       fit_status = if (n) if (valid_interval) "ok" else "interval_invalid" else
         "no_analysable_units",
       fit_failure = n == 0, interval_failure = !valid_interval,
       boundary_failure = FALSE, singularity_failure = FALSE,
       structural_nonidentifiable = FALSE,
       optim_convergence = NA_integer_, successful_starts = NA_integer_,
       score_max_abs = NA_real_, information_min_eigen = NA_real_,
       information_condition_number = NA_real_, analysable_units = n)
}

prepare_likelihood_data <- function(simulation, masked, estimator) {
  observed_units <- which(rowSums(masked$analysis_mask) > 0)
  selected <- which(masked$analysis_mask[observed_units, , drop = FALSE],
                    arr.ind = TRUE)
  original_unit <- observed_units[selected[, "row"]]
  y <- masked$analysis_detection[cbind(original_unit, selected[, "col"])]
  unit <- match(original_unit, observed_units)
  x_all <- if (estimator == "replication_aware") {
    matrix(1, nrow = length(simulation$river_position), ncol = 1)
  } else {
    cbind(1, simulation$river_position)
  }
  w <- if (estimator == "replication_aware") {
    matrix(1, nrow = length(y), ncol = 1)
  } else {
    cbind(1, simulation$planned_log_depth[cbind(original_unit,
                                                 selected[, "col"])])
  }
  list(y = as.numeric(y), unit = as.integer(unit),
       x = x_all[observed_units, , drop = FALSE], w = w,
       x_prediction = x_all, observed_units = observed_units,
       replicates_by_unit = tabulate(unit, nbins = length(observed_units)))
}

.atlas_backend_ready <- FALSE

ensure_atlas_backend <- function() {
  if (isTRUE(.atlas_backend_ready) && exists("atlas_llgh_cpp", mode = "function"))
    return(invisible(TRUE))
  if (!requireNamespace("Rcpp", quietly = TRUE))
    stop("Rcpp is required for the exact full-grid likelihood analysis")
  # A directly invoked conda/micromamba R may not expose that environment's
  # compiler wrappers on PATH.  Resolve the wrapper directory from R.home so
  # the documented direct-Rscript route and `micromamba run` behave alike.
  r_env_bin <- normalizePath(file.path(R.home(), "..", "..", "bin"),
                             mustWork = TRUE)
  path_entries <- strsplit(Sys.getenv("PATH"), .Platform$path.sep,
                           fixed = TRUE)[[1]]
  if (!(r_env_bin %in% path_entries))
    Sys.setenv(PATH = paste(r_env_bin, Sys.getenv("PATH"),
                            sep = .Platform$path.sep))
  Rcpp::cppFunction(env = .GlobalEnv, code = '
  Rcpp::List atlas_llgh_cpp(Rcpp::NumericVector par,
                            Rcpp::NumericMatrix X,
                            Rcpp::NumericMatrix W,
                            Rcpp::NumericVector y,
                            Rcpp::IntegerVector unit,
                            bool need_hessian = false) {
    int n = X.nrow(), pb = X.ncol(), m = W.nrow(), pg = W.ncol();
    int pt = pb + pg;
    Rcpp::NumericVector psi(n), logpsi(n), log1mpsi(n);
    Rcpp::NumericVector p(m), logp(m), log1mp(m);
    Rcpp::NumericVector logq(n), detll(n), postz(n);
    Rcpp::IntegerVector anydet(n);
    for (int i = 0; i < n; ++i) {
      double eta = 0.0;
      for (int j = 0; j < pb; ++j) eta += X(i,j) * par[j];
      if (eta >= 0.0) {
        logpsi[i] = -std::log1p(std::exp(-eta));
        log1mpsi[i] = -eta - std::log1p(std::exp(-eta));
      } else {
        logpsi[i] = eta - std::log1p(std::exp(eta));
        log1mpsi[i] = -std::log1p(std::exp(eta));
      }
      psi[i] = std::exp(logpsi[i]);
    }
    for (int c = 0; c < m; ++c) {
      double eta = 0.0;
      for (int j = 0; j < pg; ++j) eta += W(c,j) * par[pb+j];
      if (eta >= 0.0) {
        logp[c] = -std::log1p(std::exp(-eta));
        log1mp[c] = -eta - std::log1p(std::exp(-eta));
      } else {
        logp[c] = eta - std::log1p(std::exp(eta));
        log1mp[c] = -std::log1p(std::exp(eta));
      }
      p[c] = std::exp(logp[c]);
      int u = unit[c] - 1;
      logq[u] += log1mp[c];
      detll[u] += y[c] * logp[c] + (1.0-y[c]) * log1mp[c];
      if (y[c] > 0.5) anydet[u] = 1;
    }
    double ll = 0.0;
    Rcpp::NumericVector grad(pt);
    for (int i = 0; i < n; ++i) {
      double lli;
      if (anydet[i]) {
        lli = logpsi[i] + detll[i];
        postz[i] = 1.0;
      } else {
        double a = log1mpsi[i], b = logpsi[i] + logq[i];
        double mx = a > b ? a : b;
        lli = mx + std::log(std::exp(a-mx) + std::exp(b-mx));
        postz[i] = std::exp(b-lli);
      }
      ll += lli;
      for (int j = 0; j < pb; ++j)
        grad[j] += X(i,j) * (postz[i] - psi[i]);
    }
    for (int c = 0; c < m; ++c) {
      int u = unit[c] - 1;
      for (int j = 0; j < pg; ++j)
        grad[pb+j] += W(c,j) * postz[u] * (y[c] - p[c]);
    }
    Rcpp::NumericMatrix hess(pt, pt);
    if (need_hessian) {
      Rcpp::NumericMatrix S(n, pg);
      for (int c = 0; c < m; ++c) {
        int u = unit[c] - 1;
        for (int j = 0; j < pg; ++j)
          S(u,j) += W(c,j) * (y[c] - p[c]);
      }
      for (int i = 0; i < n; ++i) {
        double vr = postz[i] * (1.0-postz[i]);
        double a = -psi[i] * (1.0-psi[i]) + vr;
        for (int j = 0; j < pb; ++j)
          for (int k = 0; k < pb; ++k)
            hess(j,k) += a * X(i,j) * X(i,k);
        for (int j = 0; j < pg; ++j)
          for (int k = 0; k < pg; ++k)
            hess(pb+j,pb+k) += vr * S(i,j) * S(i,k);
        for (int j = 0; j < pb; ++j)
          for (int k = 0; k < pg; ++k) {
            double cross = vr * X(i,j) * S(i,k);
            hess(j,pb+k) += cross;
            hess(pb+k,j) += cross;
          }
      }
      for (int c = 0; c < m; ++c) {
        int u = unit[c] - 1;
        double wt = postz[u] * p[c] * (1.0-p[c]);
        for (int j = 0; j < pg; ++j)
          for (int k = 0; k < pg; ++k)
            hess(pb+j,pb+k) -= wt * W(c,j) * W(c,k);
      }
    }
    return Rcpp::List::create(Rcpp::Named("loglik") = ll,
                              Rcpp::Named("gradient") = grad,
                              Rcpp::Named("hessian") = hess,
                              Rcpp::Named("posterior_z") = postz);
  }')
  .atlas_backend_ready <<- TRUE
  invisible(TRUE)
}

likelihood_starts <- function(dat, estimator) {
  detected_unit <- tabulate(dat$unit[dat$y > 0.5],
                            nbins = nrow(dat$x)) > 0
  occurrence <- (sum(detected_unit) + 0.5) / (length(detected_unit) + 1)
  p_starts <- c(0.2, 0.5, 0.8)
  slope_start <- 0
  if (estimator == "depth_aware" && stats::var(dat$x[, 2]) > 0) {
    denominator <- stats::var(dat$x[, 2]) * occurrence * (1 - occurrence)
    if (is.finite(denominator) && denominator > 1e-8)
      slope_start <- max(-3, min(3, stats::cov(dat$x[, 2],
                                               as.numeric(detected_unit)) / denominator))
  }
  lapply(seq_along(p_starts), function(i) {
    p0 <- p_starts[[i]]
    unit_detection <- mean(1 - (1 - p0)^dat$replicates_by_unit)
    psi0 <- max(0.01, min(0.99, occurrence / max(unit_detection, 1e-6)))
    if (estimator == "replication_aware") {
      c(stats::qlogis(psi0), stats::qlogis(p0))
    } else {
      c(stats::qlogis(psi0), if (i == 3L) 0 else slope_start,
        stats::qlogis(p0), if (i == 2L) 0.75 else 0)
    }
  })
}

fit_integrated_likelihood <- function(simulation, masked, estimator) {
  if (ncol(masked$analysis_mask) == 1L) {
    return(list(
      estimate = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      fit_status = "structural_nonidentifiable", fit_failure = TRUE,
      interval_failure = TRUE, boundary_failure = FALSE,
      singularity_failure = FALSE, structural_nonidentifiable = TRUE,
      optim_convergence = NA_integer_, successful_starts = 0L,
      score_max_abs = NA_real_, information_min_eigen = NA_real_,
      information_condition_number = NA_real_,
      analysable_units = sum(rowSums(masked$analysis_mask) > 0)))
  }
  ensure_atlas_backend()
  dat <- prepare_likelihood_data(simulation, masked, estimator)
  starts <- likelihood_starts(dat, estimator)
  lower <- rep(-ATLAS_PARAMETER_BOUND, length(starts[[1]]))
  upper <- rep(ATLAS_PARAMETER_BOUND, length(starts[[1]]))
  fits <- lapply(starts, function(start) {
    cache <- new.env(parent = emptyenv())
    evaluate <- function(par) {
      if (!is.null(cache$par) && identical(as.numeric(par), cache$par))
        return(cache$value)
      value <- atlas_llgh_cpp(par, dat$x, dat$w, dat$y, dat$unit, FALSE)
      cache$par <- as.numeric(par)
      cache$value <- value
      value
    }
    tryCatch(stats::optim(
      par = pmax(lower, pmin(upper, start)),
      fn = function(par) -evaluate(par)$loglik,
      gr = function(par) -evaluate(par)$gradient,
      method = "L-BFGS-B", lower = lower, upper = upper,
      control = list(maxit = ATLAS_OPTIM_MAXIT, pgtol = ATLAS_OPTIM_PGTOL,
                     factr = 1e7)),
      error = function(e) list(par = rep(NA_real_, length(start)), value = Inf,
                               convergence = 999L, message = conditionMessage(e)))
  })
  successful <- vapply(fits, function(x)
    identical(as.integer(x$convergence), 0L) && is.finite(x$value) &&
      all(is.finite(x$par)), logical(1))
  finite <- vapply(fits, function(x) is.finite(x$value) &&
                     all(is.finite(x$par)), logical(1))
  candidates <- if (any(successful)) which(successful) else which(finite)
  if (!length(candidates)) {
    return(list(
      estimate = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
      fit_status = "optimizer_failure", fit_failure = TRUE,
      interval_failure = TRUE, boundary_failure = FALSE,
      singularity_failure = FALSE, structural_nonidentifiable = FALSE,
      optim_convergence = 999L, successful_starts = 0L,
      score_max_abs = NA_real_, information_min_eigen = NA_real_,
      information_condition_number = NA_real_,
      analysable_units = length(dat$observed_units)))
  }
  best_index <- candidates[which.min(vapply(fits[candidates], `[[`, numeric(1),
                                            "value"))]
  best <- fits[[best_index]]
  final <- atlas_llgh_cpp(best$par, dat$x, dat$w, dat$y, dat$unit, TRUE)
  pb <- ncol(dat$x)
  beta <- best$par[seq_len(pb)]
  predicted <- stats::plogis(drop(dat$x_prediction %*% beta))
  estimate <- mean(predicted)
  information <- -as.matrix(final$hessian)
  information <- (information + t(information)) / 2
  eigenvalues <- tryCatch(eigen(information, symmetric = TRUE,
                                only.values = TRUE)$values,
                          error = function(e) rep(NA_real_, nrow(information)))
  min_eigen <- if (all(is.finite(eigenvalues))) min(eigenvalues) else NA_real_
  max_eigen <- if (all(is.finite(eigenvalues))) max(eigenvalues) else NA_real_
  condition_number <- if (is.finite(min_eigen) && min_eigen > 0)
    max_eigen / min_eigen else Inf
  boundary <- any(abs(best$par) >=
                    ATLAS_PARAMETER_BOUND - ATLAS_BOUNDARY_TOLERANCE)
  singular <- !all(is.finite(information)) || !is.finite(min_eigen) ||
    min_eigen <= max(ATLAS_INFORMATION_ABSOLUTE_TOLERANCE,
                     ATLAS_INFORMATION_RELATIVE_TOLERANCE * max_eigen) ||
    condition_number > ATLAS_INFORMATION_MAX_CONDITION
  optimizer_failure <- !any(successful)
  fit_failure <- optimizer_failure || boundary || singular
  lower95 <- upper95 <- NA_real_
  interval_failure <- TRUE
  if (!fit_failure && is.finite(estimate) && estimate > 0 && estimate < 1) {
    covariance <- tryCatch(solve(information), error = function(e) NULL)
    if (!is.null(covariance)) {
      gradient_mean <- colMeans(dat$x_prediction * predicted * (1 - predicted))
      gradient_logit <- gradient_mean / (estimate * (1 - estimate))
      variance_logit <- drop(t(gradient_logit) %*%
                               covariance[seq_len(pb), seq_len(pb), drop = FALSE] %*%
                               gradient_logit)
      if (is.finite(variance_logit) && variance_logit > 0) {
        se <- sqrt(variance_logit)
        lower95 <- stats::plogis(stats::qlogis(estimate) - 1.96 * se)
        upper95 <- stats::plogis(stats::qlogis(estimate) + 1.96 * se)
        interval_failure <- !(is.finite(lower95) && is.finite(upper95) &&
                                0 <= lower95 && lower95 < estimate &&
                                estimate < upper95 && upper95 <= 1)
      }
    }
  }
  status <- if (optimizer_failure) "optimizer_failure" else if (boundary)
    "boundary_failure" else if (singular) "singular_information" else if
      (interval_failure) "interval_invalid" else "ok"
  list(
    estimate = estimate, lower95 = lower95, upper95 = upper95,
    fit_status = status, fit_failure = fit_failure,
    interval_failure = interval_failure, boundary_failure = boundary,
    singularity_failure = singular, structural_nonidentifiable = FALSE,
    optim_convergence = as.integer(best$convergence),
    successful_starts = sum(successful),
    score_max_abs = max(abs(final$gradient)), information_min_eigen = min_eigen,
    information_condition_number = condition_number,
    analysable_units = length(dat$observed_units))
}

empty_fit_frame <- function(n) {
  data.frame(
    data_scope = rep(ATLAS_SCOPE, n), scenario_id = integer(n),
    truth_scenario_id = integer(n), simulation_id = integer(n),
    simulation_seed = integer(n), estimator = character(n),
    true_marginal_occupancy = numeric(n), estimate = rep(NA_real_, n),
    lower95 = rep(NA_real_, n), upper95 = rep(NA_real_, n),
    covered = rep(NA, n), fit_status = character(n),
    fit_failure = logical(n), interval_failure = logical(n),
    boundary_failure = logical(n), singularity_failure = logical(n),
    structural_nonidentifiable = logical(n),
    optim_convergence = rep(NA_integer_, n), successful_starts = integer(n),
    score_max_abs = rep(NA_real_, n),
    information_min_eigen = rep(NA_real_, n),
    information_condition_number = rep(NA_real_, n),
    stringsAsFactors = FALSE)
}

empty_draw_frame <- function(n) {
  data.frame(
    data_scope = rep(ATLAS_SCOPE, n), scenario_id = integer(n),
    truth_scenario_id = integer(n), simulation_id = integer(n),
    simulation_seed = integer(n), true_marginal_occupancy = numeric(n),
    realised_latent_prevalence = numeric(n), planned_cells = integer(n),
    truth_occupied_units = integer(n), structural_missing_cells = integer(n),
    surveyed_mask_cells = integer(n), analysis_mask_cells = integer(n),
    generated_detection_events = integer(n), observed_detection_events = integer(n),
    analysis_detection_events = integer(n), analysable_units = integer(n),
    erroneous_zero_filled_cells = integer(n), stringsAsFactors = FALSE)
}

make_generator_fixture <- function(simulation, masked, scenario_id, mask_scenario) {
  n <- length(simulation$generated_detection)
  data.frame(
    data_scope = rep(ATLAS_SCOPE, n), scenario_id = rep(as.integer(scenario_id), n),
    truth_scenario_id = rep(simulation$truth_scenario_id, n),
    simulation_id = rep(simulation$simulation_id, n),
    simulation_seed = rep(simulation$simulation_seed, n),
    mask_scenario = rep(mask_scenario, n),
    unit_id = rep(seq_along(simulation$river_position),
                  times = ncol(simulation$generated_detection)),
    replicate_id = rep(seq_len(ncol(simulation$generated_detection)),
                       each = nrow(simulation$generated_detection)),
    river_position = rep(simulation$river_position,
                         times = ncol(simulation$generated_detection)),
    true_occupancy_probability = rep(simulation$true_occupancy_probability,
                                     times = ncol(simulation$generated_detection)),
    latent_occupancy = rep(simulation$latent_occupancy,
                           times = ncol(simulation$generated_detection)),
    planned_log_depth = as.vector(simulation$planned_log_depth),
    detection_probability = as.vector(simulation$detection_probability),
    surveyed_mask = as.vector(masked$surveyed_mask),
    analysis_mask = as.vector(masked$analysis_mask),
    generated_detection = as.vector(simulation$generated_detection),
    observed_detection = as.vector(masked$observed_detection),
    analysis_detection = as.vector(masked$analysis_detection),
    erroneous_zero_fill = as.vector(masked$erroneous_zero_fill),
    stringsAsFactors = FALSE)
}

run_truth_scenario <- function(core_row, scenarios, config, n_simulations,
                               fixture_truth_id) {
  mask_levels <- as.character(config$occ$mask_scenario)
  n_draw <- n_simulations * length(mask_levels)
  n_fit <- n_draw * length(ATLAS_ESTIMATORS)
  draws <- empty_draw_frame(n_draw)
  fits <- empty_fit_frame(n_fit)
  fixtures <- list()
  for (simulation_id in seq_len(n_simulations)) {
    simulation <- simulate_atlas_truth(core_row, simulation_id, config)
    true_marginal <- mean(simulation$true_occupancy_probability)
    for (mask_index in seq_along(mask_levels)) {
      mask_scenario <- mask_levels[[mask_index]]
      scenario_id <- (as.integer(core_row$truth_scenario_id) - 1L) *
        length(mask_levels) + mask_index
      masked <- apply_atlas_mask(simulation, mask_scenario)
      draw_index <- (simulation_id - 1L) * length(mask_levels) + mask_index
      draws$scenario_id[[draw_index]] <- scenario_id
      draws$truth_scenario_id[[draw_index]] <- core_row$truth_scenario_id
      draws$simulation_id[[draw_index]] <- simulation_id
      draws$simulation_seed[[draw_index]] <- simulation$simulation_seed
      draws$true_marginal_occupancy[[draw_index]] <- true_marginal
      draws$realised_latent_prevalence[[draw_index]] <-
        mean(simulation$latent_occupancy)
      draws$planned_cells[[draw_index]] <- length(simulation$generated_detection)
      draws$truth_occupied_units[[draw_index]] <- sum(simulation$latent_occupancy)
      draws$structural_missing_cells[[draw_index]] <- sum(!masked$surveyed_mask)
      draws$surveyed_mask_cells[[draw_index]] <- sum(masked$surveyed_mask)
      draws$analysis_mask_cells[[draw_index]] <- sum(masked$analysis_mask)
      draws$generated_detection_events[[draw_index]] <-
        sum(simulation$generated_detection)
      draws$observed_detection_events[[draw_index]] <-
        sum(masked$observed_detection, na.rm = TRUE)
      draws$analysis_detection_events[[draw_index]] <-
        sum(masked$analysis_detection, na.rm = TRUE)
      draws$analysable_units[[draw_index]] <-
        sum(rowSums(masked$analysis_mask) > 0)
      draws$erroneous_zero_filled_cells[[draw_index]] <-
        sum(masked$erroneous_zero_fill)

      if (as.integer(core_row$truth_scenario_id) == fixture_truth_id &&
          simulation_id == 1L) {
        fixtures[[length(fixtures) + 1L]] <- make_generator_fixture(
          simulation, masked, scenario_id, mask_scenario)
      }

      for (estimator_index in seq_along(ATLAS_ESTIMATORS)) {
        estimator <- ATLAS_ESTIMATORS[[estimator_index]]
        fit <- if (estimator == "naive_occurrence") {
          fit_naive_occurrence(masked)
        } else {
          fit_integrated_likelihood(simulation, masked, estimator)
        }
        fit_index <- (draw_index - 1L) * length(ATLAS_ESTIMATORS) + estimator_index
        fits$scenario_id[[fit_index]] <- scenario_id
        fits$truth_scenario_id[[fit_index]] <- core_row$truth_scenario_id
        fits$simulation_id[[fit_index]] <- simulation_id
        fits$simulation_seed[[fit_index]] <- simulation$simulation_seed
        fits$estimator[[fit_index]] <- estimator
        fits$true_marginal_occupancy[[fit_index]] <- true_marginal
        fits$estimate[[fit_index]] <- fit$estimate
        fits$lower95[[fit_index]] <- fit$lower95
        fits$upper95[[fit_index]] <- fit$upper95
        fits$covered[[fit_index]] <- if (!fit$interval_failure)
          true_marginal >= fit$lower95 && true_marginal <= fit$upper95 else NA
        fits$fit_status[[fit_index]] <- fit$fit_status
        fits$fit_failure[[fit_index]] <- fit$fit_failure
        fits$interval_failure[[fit_index]] <- fit$interval_failure
        fits$boundary_failure[[fit_index]] <- fit$boundary_failure
        fits$singularity_failure[[fit_index]] <- fit$singularity_failure
        fits$structural_nonidentifiable[[fit_index]] <-
          fit$structural_nonidentifiable
        fits$optim_convergence[[fit_index]] <- fit$optim_convergence
        fits$successful_starts[[fit_index]] <- fit$successful_starts
        fits$score_max_abs[[fit_index]] <- fit$score_max_abs
        fits$information_min_eigen[[fit_index]] <- fit$information_min_eigen
        fits$information_condition_number[[fit_index]] <-
          fit$information_condition_number
      }
    }
  }
  list(draws = draws, fits = fits,
       fixture = if (length(fixtures)) do.call(rbind, fixtures) else NULL)
}

classify_atlas_gate <- function(absolute_bias, rmse, coverage,
                                fit_failure_rate, interval_failure_rate,
                                combined_identification_failure_rate,
                                structural_nonidentifiable_rate, gate) {
  if (is.finite(structural_nonidentifiable_rate) &&
      structural_nonidentifiable_rate > 0)
    return("non-identifiable")
  if (is.finite(combined_identification_failure_rate) &&
      combined_identification_failure_rate >=
      as.numeric(gate$nonidentifiable_failure_rate))
    return("non-identifiable")
  pass <- all(is.finite(c(absolute_bias, rmse, coverage, fit_failure_rate,
                          interval_failure_rate))) &&
    absolute_bias <= as.numeric(gate$max_absolute_bias) + 1e-12 &&
    rmse <= as.numeric(gate$max_rmse) + 1e-12 &&
    coverage >= as.numeric(gate$coverage_min) - 1e-12 &&
    coverage <= as.numeric(gate$coverage_max) + 1e-12 &&
    fit_failure_rate <= as.numeric(gate$max_fit_failure_rate) + 1e-12 &&
    interval_failure_rate <= as.numeric(gate$max_interval_failure_rate) + 1e-12
  if (pass) "recoverable" else "unreliable"
}

summarise_atlas_fits <- function(fits, scenarios, config) {
  n_groups <- nrow(scenarios) * length(ATLAS_ESTIMATORS)
  out <- data.frame(
    data_scope = rep(ATLAS_SCOPE, n_groups), scenario_id = integer(n_groups),
    estimator = character(n_groups), simulations = integer(n_groups),
    estimate_denominator = integer(n_groups), bias = numeric(n_groups),
    absolute_bias = numeric(n_groups), rmse = numeric(n_groups),
    coverage_numerator = integer(n_groups), coverage_denominator = integer(n_groups),
    coverage = numeric(n_groups), optimizer_failure_n = integer(n_groups),
    optimizer_failure_rate = numeric(n_groups), boundary_failure_n = integer(n_groups),
    boundary_failure_rate = numeric(n_groups), singularity_failure_n = integer(n_groups),
    singularity_failure_rate = numeric(n_groups), fit_failure_n = integer(n_groups),
    fit_failure_rate = numeric(n_groups), interval_failure_n = integer(n_groups),
    interval_failure_rate = numeric(n_groups),
    structural_nonidentifiable_n = integer(n_groups),
    structural_nonidentifiable_rate = numeric(n_groups),
    combined_identification_failure_n = integer(n_groups),
    combined_identification_failure_rate = numeric(n_groups),
    gate_class = character(n_groups), stringsAsFactors = FALSE)
  out[c("bias", "absolute_bias", "rmse", "coverage")] <- NA_real_
  key <- paste(fits$scenario_id, fits$estimator, sep = "|")
  indices <- split(seq_len(nrow(fits)), key)
  row_index <- 0L
  for (scenario_id in scenarios$scenario_id) {
    for (estimator in ATLAS_ESTIMATORS) {
      row_index <- row_index + 1L
      idx <- indices[[paste(scenario_id, estimator, sep = "|")]]
      sub <- fits[idx, , drop = FALSE]
      n <- nrow(sub)
      optimizer_failure <- sub$fit_status == "optimizer_failure"
      combined <- optimizer_failure | sub$boundary_failure |
        sub$singularity_failure | sub$structural_nonidentifiable
      valid_estimate <- !sub$fit_failure & is.finite(sub$estimate)
      valid_interval <- valid_estimate & !sub$interval_failure &
        is.finite(sub$lower95) & is.finite(sub$upper95) & !is.na(sub$covered)
      bias <- if (any(valid_estimate))
        mean(sub$estimate[valid_estimate] -
               sub$true_marginal_occupancy[valid_estimate]) else NA_real_
      rmse <- if (any(valid_estimate)) sqrt(mean(
        (sub$estimate[valid_estimate] -
           sub$true_marginal_occupancy[valid_estimate])^2)) else NA_real_
      coverage_numerator <- sum(sub$covered[valid_interval])
      coverage_denominator <- sum(valid_interval)
      coverage <- if (coverage_denominator)
        coverage_numerator / coverage_denominator else NA_real_
      values <- list(
        simulations = n, estimate_denominator = sum(valid_estimate), bias = bias,
        absolute_bias = abs(bias), rmse = rmse,
        coverage_numerator = coverage_numerator,
        coverage_denominator = coverage_denominator, coverage = coverage,
        optimizer_failure_n = sum(optimizer_failure),
        optimizer_failure_rate = mean(optimizer_failure),
        boundary_failure_n = sum(sub$boundary_failure),
        boundary_failure_rate = mean(sub$boundary_failure),
        singularity_failure_n = sum(sub$singularity_failure),
        singularity_failure_rate = mean(sub$singularity_failure),
        fit_failure_n = sum(sub$fit_failure),
        fit_failure_rate = mean(sub$fit_failure),
        interval_failure_n = sum(sub$interval_failure),
        interval_failure_rate = mean(sub$interval_failure),
        structural_nonidentifiable_n = sum(sub$structural_nonidentifiable),
        structural_nonidentifiable_rate = mean(sub$structural_nonidentifiable),
        combined_identification_failure_n = sum(combined),
        combined_identification_failure_rate = mean(combined))
      out$scenario_id[[row_index]] <- scenario_id
      out$estimator[[row_index]] <- estimator
      for (name in names(values)) out[[name]][[row_index]] <- values[[name]]
      out$gate_class[[row_index]] <- classify_atlas_gate(
        values$absolute_bias, values$rmse, values$coverage,
        values$fit_failure_rate, values$interval_failure_rate,
        values$combined_identification_failure_rate,
        values$structural_nonidentifiable_rate, config$occ$gate)
    }
  }
  merged <- merge(scenarios, out, by = c("data_scope", "scenario_id"),
                  sort = FALSE)
  merged$estimator <- factor(merged$estimator, levels = ATLAS_ESTIMATORS)
  merged <- merged[order(merged$scenario_id, merged$estimator), , drop = FALSE]
  merged$estimator <- as.character(merged$estimator)
  rownames(merged) <- NULL
  merged
}

run_occupancy_atlas <- function(config, scenarios = build_atlas_scenarios(config),
                                n_simulations = as.integer(config$occ$simulations_per_scenario),
                                cores = 1L) {
  ensure_atlas_backend()
  core <- scenarios[!duplicated(scenarios$truth_scenario_id), , drop = FALSE]
  fixture_truth_id <- with(core, truth_scenario_id[
    target_mean_occupancy == 0.35 &
      prevalence_structure == "downstream_gradient" &
      baseline_detection == 0.5 & replicates == 3L &
      log_depth_sd == 0.75 & depth_river_correlation == 0.7][1])
  # Small unit-test subsets may not contain the full-run fixture; their first
  # paired truth scenario becomes the fixture without affecting canonical runs.
  if (!length(fixture_truth_id) || !is.finite(fixture_truth_id))
    fixture_truth_id <- core$truth_scenario_id[[1]]
  cores <- max(1L, min(as.integer(cores), nrow(core)))
  batches <- split(seq_len(nrow(core)),
                   ceiling(seq_len(nrow(core)) / max(1L, 4L * cores)))
  chunks <- vector("list", nrow(core))
  completed <- 0L
  for (batch in batches) {
    worker <- function(i) run_truth_scenario(core[i, , drop = FALSE], scenarios,
                                             config, n_simulations,
                                             fixture_truth_id)
    result <- if (cores > 1L && .Platform$OS.type != "windows") {
      parallel::mclapply(batch, worker, mc.cores = min(cores, length(batch)),
                         mc.preschedule = TRUE, mc.set.seed = FALSE)
    } else lapply(batch, worker)
    chunks[batch] <- result
    completed <- completed + length(batch)
    message(sprintf("[impact atlas] completed %d/%d paired truth scenarios",
                    completed, nrow(core)))
  }
  draws <- do.call(rbind, lapply(chunks, `[[`, "draws"))
  fits <- do.call(rbind, lapply(chunks, `[[`, "fits"))
  fixture_chunks <- Filter(Negate(is.null), lapply(chunks, `[[`, "fixture"))
  fixture <- do.call(rbind, fixture_chunks)
  fits$estimator <- factor(fits$estimator, levels = ATLAS_ESTIMATORS)
  fits <- fits[order(fits$scenario_id, fits$simulation_id, fits$estimator), , drop = FALSE]
  fits$estimator <- as.character(fits$estimator)
  draws <- draws[order(draws$scenario_id, draws$simulation_id), , drop = FALSE]
  fixture <- fixture[order(fixture$scenario_id, fixture$unit_id,
                           fixture$replicate_id), , drop = FALSE]
  rownames(fits) <- rownames(draws) <- rownames(fixture) <- NULL
  atlas <- summarise_atlas_fits(fits, scenarios, config)
  list(scenarios = scenarios, draws = draws, fits = fits,
       fixture = fixture, atlas = atlas)
}

canonicalise_frame <- function(x, digits = 10L) {
  out <- x
  for (name in names(out)) {
    if (is.numeric(out[[name]]) && !is.integer(out[[name]])) {
      value <- out[[name]]
      out[[name]] <- ifelse(is.na(value), NA_character_,
                            formatC(value, format = "f", digits = digits))
    } else if (is.logical(out[[name]])) {
      out[[name]] <- ifelse(is.na(out[[name]]), NA_character_,
                            ifelse(out[[name]], "true", "false"))
    }
  }
  out
}

write_canonical_tsv <- function(x, relative_path) {
  path <- atlas_path(relative_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- tempfile(pattern = "atlas-", tmpdir = dirname(path))
  utils::write.table(canonicalise_frame(x), temporary, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE, na = "NA",
                     fileEncoding = "UTF-8")
  if (!file.rename(temporary, path)) stop("Could not atomically write ", path)
  invisible(path)
}

plot_occupancy_atlas <- function(atlas, png_path, svg_path) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required")
  d <- atlas
  target_levels <- sort(unique(d$target_mean_occupancy))
  detection_levels <- sort(unique(d$baseline_detection))
  sd_levels <- sort(unique(d$log_depth_sd))
  rho_levels <- sort(unique(d$depth_river_correlation))
  prevalence_levels <- c("uniform", "downstream_gradient")
  d$x_index <- (match(d$replicates, sort(unique(d$replicates))) - 1L) *
    length(detection_levels) * length(rho_levels) +
    (match(d$baseline_detection, detection_levels) - 1L) * length(rho_levels) +
    match(d$depth_river_correlation, rho_levels)
  d$y_index <- (match(d$target_mean_occupancy, target_levels) - 1L) *
    length(prevalence_levels) * length(sd_levels) +
    (match(d$prevalence_structure, prevalence_levels) - 1L) * length(sd_levels) +
    match(d$log_depth_sd, sd_levels)
  x_lookup <- unique(d[, c("x_index", "replicates", "baseline_detection",
                           "depth_river_correlation")])
  x_lookup <- x_lookup[order(x_lookup$x_index), ]
  x_labels <- sprintf("K%d\np%.1f\nr%s", x_lookup$replicates,
                      x_lookup$baseline_detection,
                      sub("^0\\.", ".", x_lookup$depth_river_correlation))
  y_lookup <- unique(d[, c("y_index", "target_mean_occupancy",
                           "prevalence_structure", "log_depth_sd")])
  y_lookup <- y_lookup[order(y_lookup$y_index), ]
  y_labels <- sprintf("psi %.2f | %s | sd %.2f",
                      y_lookup$target_mean_occupancy,
                      ifelse(y_lookup$prevalence_structure == "uniform", "U", "G"),
                      y_lookup$log_depth_sd)
  d$gate_class <- factor(d$gate_class,
                         levels = c("recoverable", "unreliable", "non-identifiable"))
  d$estimator_panel <- factor(d$estimator, levels = ATLAS_ESTIMATORS,
                              labels = c("Naive occurrence",
                                         "Replication-aware likelihood",
                                         "Depth-aware likelihood"))
  d$mask_panel <- factor(d$mask_scenario,
                         levels = c("complete_correct", "structured_correct",
                                    "structured_zero_fill"),
                         labels = c("Complete / correct",
                                    "20% unavailable / correct",
                                    "20% unavailable / zero-filled"))
  failure_points <- d[d$combined_identification_failure_rate > 0, , drop = FALSE]
  counts <- table(d$gate_class)
  subtitle <- sprintf(
    paste0("All 1,944 scenario-estimator cells | recoverable %d | unreliable %d | ",
           "non-identifiable %d | U = uniform, G = downstream gradient"),
    counts[["recoverable"]], counts[["unreliable"]], counts[["non-identifiable"]])
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x_index, y = y_index)) +
    ggplot2::geom_tile(ggplot2::aes(fill = gate_class), colour = "white",
                       linewidth = 0.12) +
    ggplot2::geom_vline(xintercept = c(6.5, 12.5), colour = "#334155",
                        linewidth = 0.18) +
    ggplot2::geom_hline(yintercept = c(4.5, 8.5), colour = "#334155",
                        linewidth = 0.18) +
    ggplot2::geom_point(data = failure_points,
                        ggplot2::aes(size = combined_identification_failure_rate),
                        shape = 21, fill = "white", colour = "#111827",
                        stroke = 0.25, alpha = 0.85) +
    ggplot2::facet_grid(estimator_panel ~ mask_panel) +
    ggplot2::scale_fill_manual(
      values = c("recoverable" = "#159D82", "unreliable" = "#E5A832",
                 "non-identifiable" = "#7C4D9E"), drop = FALSE) +
    ggplot2::scale_size_area(name = "Identification\nfailure rate",
                             limits = c(0, 1), max_size = 3.2,
                             breaks = c(0.2, 0.5, 1)) +
    ggplot2::scale_x_continuous(breaks = x_lookup$x_index, labels = x_labels,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = y_lookup$y_index, labels = y_labels,
                                expand = c(0, 0)) +
    ggplot2::labs(
      title = "Occupancy recoverability is a property of the observation design",
      subtitle = subtitle,
      x = "planned replication (K) x baseline detection (p) x depth-river correlation (r)",
      y = "target occupancy x prevalence structure x log-depth SD",
      fill = "Gate class",
      caption = paste0(
        "Each tile is one preregistered scenario; colours apply the joint bias, RMSE, ",
        "coverage and failure gate across 100 simulations. White circles scale with ",
        "combined optimisation/boundary/singularity/non-identifiability failure. ",
        "Synthetic known truth only; this does not estimate or repair Conwy occupancy.")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18,
                                         colour = "#132238"),
      plot.subtitle = ggplot2::element_text(size = 11, colour = "#475569"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      strip.background = ggplot2::element_rect(fill = "#E8EEF4", colour = NA),
      panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 5.5),
      axis.text.y = ggplot2::element_text(size = 6.3),
      legend.position = "bottom", legend.box = "horizontal",
      plot.caption = ggplot2::element_text(size = 8, colour = "#475569", hjust = 0),
      plot.margin = ggplot2::margin(12, 15, 12, 12))
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(png_path, p, width = 24, height = 16, units = "in", dpi = 160,
                  bg = "white")
  # Use base Cairo SVG rather than svglite so the frozen analysis needs no
  # package beyond the declared environment.
  grDevices::svg(svg_path, width = 24, height = 16, bg = "white",
                 onefile = TRUE)
  print(p)
  grDevices::dev.off()
  invisible(p)
}

atlas_output_paths <- function() c(
  scenarios = "data/derived/impact_occupancy_scenarios.tsv",
  simulation_audit = "data/derived/impact_occupancy_simulation_audit.tsv",
  fit_results = "data/derived/impact_occupancy_fit_results.tsv",
  atlas = "data/derived/impact_occupancy_atlas.tsv",
  generator_fixture = "data/derived/impact_occupancy_generator_fixture.tsv",
  figure_png = "figures/impact_occupancy_atlas.png",
  figure_svg = "figures/impact_occupancy_atlas.svg")

write_atlas_outputs <- function(result, config) {
  paths <- atlas_output_paths()
  write_canonical_tsv(result$scenarios, paths[["scenarios"]])
  write_canonical_tsv(result$draws, paths[["simulation_audit"]])
  write_canonical_tsv(result$fits, paths[["fit_results"]])
  write_canonical_tsv(result$atlas, paths[["atlas"]])
  write_canonical_tsv(result$fixture, paths[["generator_fixture"]])
  plot_occupancy_atlas(result$atlas, atlas_path(paths[["figure_png"]]),
                       atlas_path(paths[["figure_svg"]]))
  hashes <- vapply(paths, function(path) atlas_sha256(atlas_path(path)), character(1))
  class_counts <- table(factor(result$atlas$gate_class,
                               levels = c("recoverable", "unreliable",
                                          "non-identifiable")))
  input_log <- atlas_input_log()
  scientific_inputs <- input_log$scientific_inputs_read
  forbidden <- scientific_inputs[grepl(
    "layerA|phase1|occumb|msocc|occasions|count|seqtab|derived", scientific_inputs,
    ignore.case = TRUE)]
  if (length(forbidden))
    stop("Forbidden real Layer-A input was read: ", paste(forbidden, collapse = ", "))
  manifest <- list(
    schema_version = 1L, data_scope = ATLAS_SCOPE,
    evidence_scope = "synthetic known-truth design conclusions only",
    claim_boundary = paste0(
      "Cannot estimate Conwy occupancy, repair its observation mask, or revive the ",
      "withheld unconverged Layer-A correction."),
    frozen_sources = list(
      list(path = config$prereg_path, sha256 = config$prereg_sha256),
      list(path = config$manifest_path, sha256 = config$manifest_sha256)),
    scientific_inputs_read = as.list(scientific_inputs),
    generated_artifacts_read = as.list(input_log$generated_artifacts_read),
    layerA_real_inputs_read = list(),
    design = list(
      scenarios = nrow(result$scenarios), simulations_per_scenario =
        as.integer(config$occ$simulations_per_scenario),
      estimators = as.list(ATLAS_ESTIMATORS),
      expected_fit_rows = nrow(result$fits), paired_truth_scenarios =
        length(unique(result$scenarios$truth_scenario_id)),
      fixed_seed_base = config$seed,
      seed_formula = "base + 1000 * truth_scenario_id + simulation_id",
      rng = "R Mersenne-Twister / Inversion / Rejection; RNGversion 4.0.0"),
    likelihood = list(
      parameter_bounds = c(-ATLAS_PARAMETER_BOUND, ATLAS_PARAMETER_BOUND),
      deterministic_starts = 3L, hessian = "analytic observed information (Louis identity)",
      interval = "95% logit-scale delta/Wald interval for marginal occupancy",
      one_replicate_corrections = "structural_nonidentifiable"),
    output_rows = list(scenarios = nrow(result$scenarios),
                       simulation_audit = nrow(result$draws),
                       fit_results = nrow(result$fits), atlas = nrow(result$atlas),
                       generator_fixture = nrow(result$fixture)),
    gate_class_counts = as.list(as.integer(class_counts)),
    gate_class_count_names = as.list(names(class_counts)),
    canonical_output_sha256 = as.list(hashes))
  manifest_path <- atlas_path("data/derived/impact_occupancy_manifest.json")
  jsonlite::write_json(manifest, manifest_path, pretty = TRUE, auto_unbox = TRUE,
                       digits = 16, null = "null")
  c(hashes, impact_manifest = atlas_sha256(manifest_path))
}

read_prior_generated_receipts <- function(verify_rebuild) {
  if (!verify_rebuild) return(list(manifest = NULL, receipt = NULL))
  manifest_rel <- "data/derived/impact_occupancy_manifest.json"
  receipt_rel <- "data/derived/impact_occupancy_receipt.json"
  prior_manifest <- prior_receipt <- NULL
  if (file.exists(atlas_path(manifest_rel))) {
    atlas_record_input(manifest_rel, "generated")
    prior_manifest <- jsonlite::fromJSON(atlas_path(manifest_rel), simplifyVector = TRUE)
  }
  if (file.exists(atlas_path(receipt_rel))) {
    atlas_record_input(receipt_rel, "generated")
    prior_receipt <- jsonlite::fromJSON(atlas_path(receipt_rel), simplifyVector = TRUE)
  }
  list(manifest = prior_manifest, receipt = prior_receipt)
}

write_execution_receipt <- function(config, result, hashes, elapsed_seconds, cores,
                                    prior, verify_rebuild) {
  current_core_hashes <- hashes[names(hashes) %in%
                                  c(names(atlas_output_paths()), "impact_manifest")]
  prior_hashes <- if (!is.null(prior$manifest)) {
    x <- unlist(prior$manifest$canonical_output_sha256, use.names = TRUE)
    c(x, impact_manifest = if (!is.null(prior$receipt$impact_manifest_sha256))
      as.character(prior$receipt$impact_manifest_sha256) else NA_character_)
  } else character()
  comparable <- intersect(names(current_core_hashes), names(prior_hashes))
  # The manifest embeds generated-input audit state, so its own digest can differ
  # between baseline and verification runs. Canonical tables and figures must not.
  comparable <- setdiff(comparable, "impact_manifest")
  rebuild_pass <- verify_rebuild && length(comparable) == length(atlas_output_paths()) &&
    all(current_core_hashes[comparable] == prior_hashes[comparable])
  baseline_runtime <- if (!is.null(prior$receipt$analysis_runtime_seconds))
    as.numeric(prior$receipt$analysis_runtime_seconds) else NA_real_
  receipt <- list(
    schema_version = 1L, data_scope = ATLAS_SCOPE,
    analysis_runtime_seconds = unname(elapsed_seconds),
    baseline_analysis_runtime_seconds = baseline_runtime,
    execution_cores = as.integer(cores),
    full_design_executed = nrow(result$scenarios) == 648L &&
      nrow(result$fits) == 648L * 3L * 100L,
    deterministic_rebuild_requested = verify_rebuild,
    deterministic_rebuild_pass = rebuild_pass,
    deterministic_files_compared = as.list(comparable),
    impact_manifest_sha256 = unname(hashes[["impact_manifest"]]),
    canonical_output_sha256 = as.list(current_core_hashes),
    layerA_real_inputs_read = list(),
    test_command = paste("Rscript R/test_impact_atlas.R",
                         "(run separately; see impact_occupancy_test_receipt.tsv)"))
  jsonlite::write_json(receipt,
                       atlas_path("data/derived/impact_occupancy_receipt.json"),
                       pretty = TRUE, auto_unbox = TRUE, digits = 16, null = "null")
  receipt
}

atlas_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  verify_rebuild <- "--verify-rebuild" %in% args
  unknown <- setdiff(args, "--verify-rebuild")
  if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
  atlas_reset_input_log()
  config <- load_atlas_config()
  prior <- read_prior_generated_receipts(verify_rebuild)
  scenarios <- build_atlas_scenarios(config)
  cores <- suppressWarnings(as.integer(Sys.getenv("IMPACT_ATLAS_CORES", "4")))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  message(sprintf(
    "[impact atlas] full frozen design: %d scenarios x %d simulations x 3 estimators",
    nrow(scenarios), as.integer(config$occ$simulations_per_scenario)))
  message(sprintf("[impact atlas] %d paired truth scenarios; %d worker(s)",
                  length(unique(scenarios$truth_scenario_id)), cores))
  start <- proc.time()[["elapsed"]]
  result <- run_occupancy_atlas(config, scenarios, cores = cores)
  hashes <- write_atlas_outputs(result, config)
  elapsed <- proc.time()[["elapsed"]] - start
  receipt <- write_execution_receipt(config, result, hashes, elapsed, cores,
                                     prior, verify_rebuild)
  counts <- table(factor(result$atlas$gate_class,
                         levels = c("recoverable", "unreliable", "non-identifiable")))
  cat(sprintf(paste0(
    "\n[impact atlas] COMPLETE in %.3f seconds\n",
    "[impact atlas] gate cells: recoverable=%d, unreliable=%d, non-identifiable=%d\n",
    "[impact atlas] deterministic rebuild: %s\n"),
    elapsed, counts[[1]], counts[[2]], counts[[3]],
    if (isTRUE(receipt$deterministic_rebuild_pass)) "PASS" else
      if (verify_rebuild) "FAIL" else "not requested"))
  if (verify_rebuild && !isTRUE(receipt$deterministic_rebuild_pass))
    stop("Deterministic rebuild comparison failed")
  invisible(result)
}

if (!identical(Sys.getenv("IMPACT_ATLAS_SOURCE_ONLY"), "1")) atlas_main()
