#!/usr/bin/env Rscript

# Phase-5 functional-trait broadening null model.
# This script is governed by docs/impact-sprint-preregistration.md and
# data/metadata/impact_sprint_manifest.json. It must not change those files.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript R/14_trait_null_model.R [--output-root=PATH] [--validate-only]\n",
    "  default          Run the frozen 9,999-draw primary and strict null models.\n",
    "  --output-root    Write generated data/figures below PATH (inputs still come from repo).\n",
    "  --validate-only  Validate existing default outputs without rewriting them.\n"))
  quit(status = 0)
}

output_arg <- grep("^--output-root=", args, value = TRUE)
if (length(output_arg) > 1L) stop("Only one --output-root may be supplied", call. = FALSE)
validate_only <- "--validate-only" %in% args
known <- c("--validate-only", output_arg)
unknown <- setdiff(args, known)
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_root <- if (length(output_arg)) {
  path.expand(sub("^--output-root=", "", output_arg[[1]]))
} else repo_root
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE, OutDec = ".")
assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
in_path <- function(...) file.path(repo_root, ...)
out_path <- function(...) file.path(output_root, ...)
read_tsv <- function(path, ...) read.delim(path, sep = "\t", quote = "", comment.char = "",
                                           check.names = FALSE, ...)
write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
              na = "NA", fileEncoding = "UTF-8")
}

manifest_path <- in_path("data", "metadata", "impact_sprint_manifest.json")
assert(file.exists(manifest_path), "Missing frozen impact-sprint manifest")
manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = TRUE)
cfg <- manifest$trait_null
assert(identical(as.integer(cfg$balanced_units), 154L), "Frozen trait unit count changed")
assert(identical(as.integer(cfg$null_draws), 9999L), "Frozen null-draw count changed")
assert(identical(as.integer(manifest$seeds$trait_null), 2026071903L), "Frozen trait seed changed")

paths <- list(
  counts = in_path("data", "derived", "species_x_sample_counts_clean.rds"),
  covariates = in_path("data", "metadata", "conwy_covariates.tsv"),
  groups = in_path("data", "metadata", "fish_trait_analysis_groups.tsv"),
  gower = out_path("data", "derived", "impact_trait_gower.tsv"),
  units = out_path("data", "derived", "impact_trait_unit_rao.tsv"),
  draws = out_path("data", "derived", "impact_trait_null_draws.tsv"),
  summary = out_path("data", "derived", "impact_trait_null_summary.tsv"),
  receipt = out_path("data", "metadata", "impact_trait_null_receipt.tsv"),
  png = out_path("figures", "impact_trait_null_model.png"),
  svg = out_path("figures", "impact_trait_null_model.svg")
)

required_outputs <- unlist(paths[c("gower", "units", "draws", "summary", "receipt", "png", "svg")])
if (validate_only) {
  assert(all(file.exists(required_outputs)),
         paste("Missing trait-null outputs:", paste(required_outputs[!file.exists(required_outputs)], collapse = ", ")))
  got_summary <- read_tsv(paths$summary)
  got_draws <- read_tsv(paths$draws)
  got_gower <- read_tsv(paths$gower)
  assert(nrow(got_summary) == 2L, "Trait summary must have primary and strict rows")
  assert(nrow(got_draws) == 2L * as.integer(cfg$null_draws), "Trait null draw count changed")
  assert(all(got_draws$margin_pass %in% TRUE), "A persisted null draw failed its margin invariant")
  assert(all(got_gower$distance >= 0 & got_gower$distance <= 1), "Gower distances outside [0,1]")
  assert(all(got_summary$monte_carlo_p > 0 & got_summary$monte_carlo_p <= 1),
         "Invalid Monte Carlo p-value")
  cat(jsonlite::toJSON(list(status = "PASS", summary_rows = nrow(got_summary),
                            null_draw_rows = nrow(got_draws)), auto_unbox = TRUE, pretty = TRUE), "\n")
  quit(status = 0)
}

assert(file.exists(paths$counts), "Missing species count RDS")
assert(file.exists(paths$covariates), "Missing covariate table")
assert(file.exists(paths$groups), "Missing trait mapping")

counts <- round(readRDS(paths$counts))
groups <- read_tsv(paths$groups, na.strings = character())
cov <- read_tsv(paths$covariates)
assert(is.matrix(counts) && identical(dim(counts), c(647L, 14L)),
       "Species count matrix is not 647 x 14")
assert(identical(colnames(counts), groups$analysis_label),
       "Trait mapping no longer aligns with the species matrix")
assert(all(rownames(counts) %in% cov$ID), "Covariates do not cover every count-matrix library")

md <- cov[match(rownames(counts), cov$ID), , drop = FALSE]
md$site <- sub("_.*", "", md$ID)
md$timepoint <- sub(".*_", "", md$ID)
md$unit <- paste(md$site, md$timepoint, sep = "_")
md$section <- ifelse(md$site %in% sprintf("E%02d", 1:5), "upper",
                     ifelse(md$site %in% sprintf("E%02d", 6:10), "middle", "lower"))

unit_counts <- rowsum(counts, md$unit, reorder = TRUE)
first <- match(rownames(unit_counts), md$unit)
unit_md <- data.frame(
  unit = rownames(unit_counts), site = md$site[first], timepoint = md$timepoint[first],
  section = md$section[first], stringsAsFactors = FALSE)
grid <- xtabs(~ site + timepoint, data = unit_md) > 0
common_timepoints <- colnames(grid)[colSums(grid) == nrow(grid)]
keep <- unit_md$timepoint %in% common_timepoints
unit_counts <- unit_counts[keep, , drop = FALSE]
unit_md <- unit_md[keep, , drop = FALSE]
ord <- order(as.integer(sub("E", "", unit_md$site)),
             as.integer(sub("T", "", unit_md$timepoint)))
unit_counts <- unit_counts[ord, , drop = FALSE]
unit_md <- unit_md[ord, , drop = FALSE]
assert(nrow(unit_counts) == as.integer(cfg$balanced_units), "Balanced panel is not 154 units")
assert(length(unique(unit_md$site)) == as.integer(cfg$balanced_sites), "Balanced panel is not 14 sites")
assert(length(unique(unit_md$timepoint)) == as.integer(cfg$balanced_timepoints),
       "Balanced panel is not 11 timepoints")

primary_labels <- groups$analysis_label[groups$include_in_descriptive_summary == "yes"]
strict_labels <- as.character(cfg$strict_labels)
assert(length(primary_labels) == as.integer(cfg$primary_label_count), "Primary trait label count changed")
assert(identical(strict_labels, c("Anguilla_anguilla", "Barbatula_barbatula",
                                  "Gasterosteus_aculeatus", "Phoxinus_ujmonensis",
                                  "Salmo_trutta")), "Strict trait list changed")
assert(all(strict_labels %in% primary_labels), "Strict trait labels are not a primary subset")

trait_axes <- as.character(cfg$gower_axes)
assert(all(trait_axes %in% names(groups)), "Frozen Gower axis missing from trait mapping")
assert(!any(groups[match(primary_labels, groups$analysis_label), trait_axes] %in% c("", "NA")),
       "Primary trait panel contains missing machine-coded traits")

gower_for <- function(labels) {
  x <- groups[match(labels, groups$analysis_label), trait_axes, drop = FALSE]
  x[] <- lapply(x, factor)
  d <- as.matrix(cluster::daisy(x, metric = "gower", weights = rep(1, length(trait_axes))))
  dimnames(d) <- list(labels, labels)
  assert(max(abs(d - t(d))) < 1e-12, "Gower matrix is not symmetric")
  assert(max(abs(diag(d))) < 1e-12, "Gower diagonal is not zero")
  assert(all(d >= -1e-12 & d <= 1 + 1e-12), "Gower matrix outside [0,1]")
  d
}

rao_q <- function(x, d) {
  richness <- rowSums(x)
  numer <- rowSums((x %*% d) * x)
  out <- numer / pmax(richness^2, 1)
  out[richness == 0] <- NA_real_
  out[richness == 1] <- 0
  out
}

site_contrast <- function(q) {
  site_mean <- tapply(q, unit_md$site, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
  sec <- setNames(unit_md$section[match(names(site_mean), unit_md$site)], names(site_mean))
  if (any(vapply(c("upper", "middle", "lower"), function(s) !any(sec == s & is.finite(site_mean)), logical(1)))) {
    return(list(contrast = NA_real_, site_mean = site_mean))
  }
  list(contrast = mean(site_mean[sec == "lower"], na.rm = TRUE) -
                  mean(site_mean[sec == "upper"], na.rm = TRUE),
       site_mean = site_mean)
}

has_checkerboard <- function(m) {
  if (nrow(m) < 2L || ncol(m) < 2L) return(FALSE)
  rr <- combn(seq_len(nrow(m)), 2L)
  for (i in seq_len(ncol(rr))) {
    a <- m[rr[1, i], ]; b <- m[rr[2, i], ]
    if (any(a == 1L & b == 0L) && any(a == 0L & b == 1L)) return(TRUE)
  }
  FALSE
}

assert(requireNamespace("Rcpp", quietly = TRUE),
       "Rcpp is required for the preregistered 22-million-swap ordinary-laptop null run")
# Direct invocation of an isolated conda R does not always prepend that environment's
# compiler wrappers to PATH (micromamba run does). Resolve it from R.home for both routes.
r_env_bin <- normalizePath(file.path(R.home(), "..", "..", "bin"), mustWork = TRUE)
Sys.setenv(PATH = paste(r_env_bin, Sys.getenv("PATH"), sep = .Platform$path.sep))
Rcpp::cppFunction(code = '
Rcpp::List checker_swap_cpp(Rcpp::IntegerMatrix input, int accepted_target, int max_attempts) {
  Rcpp::RNGScope scope;
  Rcpp::IntegerMatrix m = Rcpp::clone(input);
  const int nr = m.nrow();
  const int nc = m.ncol();
  int accepted = 0;
  int attempts = 0;
  std::vector<int> plus;
  std::vector<int> minus;
  plus.reserve(nc);
  minus.reserve(nc);
  while (accepted < accepted_target && attempts < max_attempts) {
    int r1 = static_cast<int>(R::runif(0.0, static_cast<double>(nr)));
    int r2 = static_cast<int>(R::runif(0.0, static_cast<double>(nr - 1)));
    if (r2 >= r1) r2 += 1;
    plus.clear();
    minus.clear();
    for (int c = 0; c < nc; ++c) {
      if (m(r1, c) == 1 && m(r2, c) == 0) plus.push_back(c);
      if (m(r1, c) == 0 && m(r2, c) == 1) minus.push_back(c);
    }
    if (!plus.empty() && !minus.empty()) {
      int c1 = plus[static_cast<int>(R::runif(0.0, static_cast<double>(plus.size())))];
      int c2 = minus[static_cast<int>(R::runif(0.0, static_cast<double>(minus.size())))];
      m(r1, c1) = 0; m(r1, c2) = 1;
      m(r2, c1) = 1; m(r2, c2) = 0;
      accepted += 1;
    }
    attempts += 1;
  }
  return Rcpp::List::create(Rcpp::Named("matrix") = m,
                            Rcpp::Named("accepted") = accepted,
                            Rcpp::Named("attempts") = attempts);
}', depends = "Rcpp")

checker_swap <- function(m, accepted_target) {
  if (accepted_target <= 0L) {
    return(list(matrix = m, accepted = 0L, attempts = 0L, structural_fixed = !has_checkerboard(m)))
  }
  if (!has_checkerboard(m)) {
    return(list(matrix = m, accepted = 0L, attempts = 0L, structural_fixed = TRUE))
  }
  max_attempts <- max(100000L, as.integer(accepted_target * 500L))
  z <- checker_swap_cpp(matrix(as.integer(m), nrow = nrow(m), ncol = ncol(m)),
                        as.integer(accepted_target), as.integer(max_attempts))
  assert(z$accepted == accepted_target,
         sprintf("Checkerboard chain stalled at %d/%d accepted swaps", z$accepted, accepted_target))
  list(matrix = z$matrix, accepted = z$accepted, attempts = z$attempts,
       structural_fixed = FALSE)
}

run_null <- function(analysis, labels, seed_offset) {
  x_obs <- 1L * (unit_counts[, labels, drop = FALSE] > 0)
  d <- gower_for(labels)
  q_obs <- rao_q(x_obs, d)
  observed <- site_contrast(q_obs)
  assert(is.finite(observed$contrast), paste("Observed", analysis, "contrast is undefined"))

  tp_levels <- sort(unique(unit_md$timepoint))
  block_rows <- lapply(tp_levels, function(tp) which(unit_md$timepoint == tp))
  blocks <- lapply(block_rows, function(ii) x_obs[ii, , drop = FALSE])
  row_targets <- lapply(blocks, rowSums)
  col_targets <- lapply(blocks, colSums)

  set.seed(as.integer(manifest$seeds$trait_null) + seed_offset)
  burnin <- as.integer(cfg$burnin_accepted_swaps_per_timepoint)
  between <- as.integer(cfg$accepted_swaps_between_draws_per_timepoint)
  burn_attempts <- integer(length(blocks))
  structural_fixed <- logical(length(blocks))
  for (i in seq_along(blocks)) {
    z <- checker_swap(blocks[[i]], burnin)
    blocks[[i]] <- z$matrix
    burn_attempts[[i]] <- z$attempts
    structural_fixed[[i]] <- z$structural_fixed
  }

  n_draw <- as.integer(cfg$null_draws)
  null_contrast <- numeric(n_draw)
  margin_pass <- logical(n_draw)
  accepted_swaps <- integer(n_draw)
  attempts <- integer(n_draw)
  swap_target_pass <- logical(n_draw)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  x_null <- x_obs

  for (draw in seq_len(n_draw)) {
    for (i in seq_along(blocks)) {
      z <- checker_swap(blocks[[i]], between)
      blocks[[i]] <- z$matrix
      accepted_swaps[[draw]] <- accepted_swaps[[draw]] + z$accepted
      attempts[[draw]] <- attempts[[draw]] + z$attempts
      x_null[block_rows[[i]], ] <- blocks[[i]]
    }
    swap_target_pass[[draw]] <- accepted_swaps[[draw]] == between * sum(!structural_fixed)
    ok <- all(vapply(seq_along(blocks), function(i) {
      identical(as.integer(rowSums(blocks[[i]])), as.integer(row_targets[[i]])) &&
        identical(as.integer(colSums(blocks[[i]])), as.integer(col_targets[[i]]))
    }, logical(1)))
    margin_pass[[draw]] <- ok
    assert(ok, paste("Null margins changed in", analysis, "draw", draw))
    null_contrast[[draw]] <- site_contrast(rao_q(x_null, d))$contrast
    assert(is.finite(null_contrast[[draw]]), paste("Null contrast undefined in", analysis, "draw", draw))
    key <- paste0(as.integer(x_null), collapse = "")
    assign(key, TRUE, envir = seen)
    if (draw %% 1000L == 0L) message(sprintf("[trait-null] %s draw %d/%d", analysis, draw, n_draw))
  }

  unique_null <- length(ls(seen, all.names = TRUE))
  p <- (1 + sum(null_contrast >= observed$contrast)) / (n_draw + 1)
  null_sd <- sd(null_contrast)
  defined_sections <- all(vapply(c("upper", "middle", "lower"), function(s) {
    site_names <- unique(unit_md$site[unit_md$section == s])
    any(is.finite(observed$site_mean[site_names]))
  }, logical(1)))
  gate <- is.finite(observed$contrast) && observed$contrast > 0 &&
    p <= as.numeric(cfg$one_sided_alpha) &&
    unique_null >= as.integer(cfg$min_unique_null_matrices) &&
    is.finite(null_sd) && null_sd > 0 && defined_sections && all(swap_target_pass)

  gower_long <- do.call(rbind, lapply(seq_along(labels), function(i) {
    data.frame(analysis = analysis, label_a = labels[[i]], label_b = labels,
               distance = as.numeric(d[i, ]), stringsAsFactors = FALSE)
  }))
  unit_rows <- data.frame(
    analysis = analysis, unit = unit_md$unit, site = unit_md$site,
    timepoint = unit_md$timepoint, section = unit_md$section,
    richness = rowSums(x_obs), rao_q = q_obs, stringsAsFactors = FALSE)
  draw_rows <- data.frame(
    analysis = analysis, draw = seq_len(n_draw), contrast_lower_minus_upper = null_contrast,
    margin_pass = margin_pass, accepted_swaps = accepted_swaps, attempts = attempts,
    structural_fixed_timepoints = sum(structural_fixed),
    swap_target_pass = swap_target_pass,
    stringsAsFactors = FALSE)
  summary_row <- data.frame(
    analysis = analysis, n_labels = length(labels), n_units = nrow(x_obs),
    n_sites = length(unique(unit_md$site)), n_timepoints = length(unique(unit_md$timepoint)),
    observed_contrast = observed$contrast, null_mean = mean(null_contrast), null_sd = null_sd,
    null_q025 = unname(quantile(null_contrast, 0.025)),
    null_q975 = unname(quantile(null_contrast, 0.975)), monte_carlo_p = p,
    null_draws = n_draw, unique_null_matrices = unique_null,
    every_margin_pass = all(margin_pass), defined_sections = defined_sections,
    structural_fixed_timepoints = sum(structural_fixed),
    every_movable_block_met_swap_target = all(swap_target_pass),
    gate_pass = gate, gate_status = if (gate) "PASSED" else "FAILED_WITHHOLD_BEYOND_RICHNESS",
    stringsAsFactors = FALSE)
  list(gower = gower_long, units = unit_rows, draws = draw_rows,
       summary = summary_row, site_mean = observed$site_mean,
       burn_attempts = sum(burn_attempts), structural_fixed_timepoints = sum(structural_fixed))
}

t0 <- proc.time()[["elapsed"]]
primary <- run_null("primary_12_labels", primary_labels, 0L)
strict <- run_null("strict_5_labels", strict_labels, 100000L)
elapsed <- proc.time()[["elapsed"]] - t0

gower <- rbind(primary$gower, strict$gower)
units <- rbind(primary$units, strict$units)
draws <- rbind(primary$draws, strict$draws)
summary <- rbind(primary$summary, strict$summary)
both_pass <- all(summary$gate_pass)
summary$claim_gate <- if (both_pass) "PASSED_BOTH" else "FAILED_WITHHOLD_BEYOND_RICHNESS"

write_tsv(gower, paths$gower)
write_tsv(units, paths$units)
write_tsv(draws, paths$draws)
write_tsv(summary, paths$summary)
receipt <- data.frame(
  analysis = "Phase 5 functional-trait null model",
  evidence_scope = "analysis-label occurrence and cited trait proxies",
  seed = as.integer(manifest$seeds$trait_null), null_draws_per_analysis = as.integer(cfg$null_draws),
  balanced_units = nrow(unit_counts), common_timepoints = paste(common_timepoints, collapse = ","),
  primary_labels = length(primary_labels), strict_labels = length(strict_labels),
  primary_burnin_attempts = primary$burn_attempts,
  strict_burnin_attempts = strict$burn_attempts,
  primary_structural_fixed_timepoints = primary$structural_fixed_timepoints,
  strict_structural_fixed_timepoints = strict$structural_fixed_timepoints,
  canonical_runtime_excluded = TRUE,
  claim_gate = if (both_pass) "PASSED_BOTH" else "FAILED_WITHHOLD_BEYOND_RICHNESS",
  allowed_claim = if (both_pass)
    "Conditional functional broadening exceeds the fixed-richness/fixed-prevalence null in both panels"
  else
    "Report conditional trait-space contrast and null interval; do not claim broadening beyond richness/prevalence",
  stringsAsFactors = FALSE)
write_tsv(receipt, paths$receipt)

draw_figure <- function(device) {
  if (device == "png") {
    dir.create(dirname(paths$png), recursive = TRUE, showWarnings = FALSE)
    png(paths$png, width = 1800, height = 1050, res = 150)
  } else {
    dir.create(dirname(paths$svg), recursive = TRUE, showWarnings = FALSE)
    svg(paths$svg, width = 12, height = 7)
  }
  op <- par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.0, 1.0), oma = c(0, 0, 2.0, 0),
            family = "sans")
  on.exit({par(op); dev.off()}, add = TRUE)
  cols <- c(upper = "#36688D", middle = "#E6A33D", lower = "#A23E48")
  for (obj in list(primary, strict)) {
    analysis <- obj$summary$analysis
    labels <- if (analysis == "primary_12_labels") "Primary: 12 analysis-label proxies" else
      "Strict: 5 better-supported proxies"
    sm <- obj$site_mean
    sites <- sprintf("E%02d", 1:14)
    sec <- ifelse(as.integer(sub("E", "", sites)) <= 5, "upper",
                  ifelse(as.integer(sub("E", "", sites)) <= 10, "middle", "lower"))
    plot(seq_along(sites), sm[sites], type = "b", pch = 21, bg = cols[sec], col = "#333333",
         xaxt = "n", xlab = "Ordered River Conwy site", ylab = "Mean Rao Q across 11 timepoints",
         main = labels, ylim = range(c(0, sm), finite = TRUE))
    axis(1, at = seq_along(sites), labels = sites, las = 2, cex.axis = 0.8)
    abline(v = c(5.5, 10.5), col = "#777777", lty = 3)
    legend("topleft", legend = c("upper", "middle", "lower"), pt.bg = cols, pch = 21,
           bty = "n", horiz = TRUE, cex = 0.8)
  }
  for (obj in list(primary, strict)) {
    z <- obj$draws$contrast_lower_minus_upper
    s <- obj$summary
    hist(z, breaks = 45, freq = FALSE, col = "#D9E4EA", border = "#FFFFFF",
         xlab = "Null lower - upper Rao Q", ylab = "Density",
         main = paste0(if (s$analysis == "primary_12_labels") "Primary" else "Strict",
                       " fixed-margin null"))
    if (sd(z) > 0) lines(density(z), col = "#36688D", lwd = 2)
    abline(v = s$observed_contrast, col = "#A23E48", lwd = 3)
    legend("topright",
           legend = c(sprintf("Observed = %.3f", s$observed_contrast),
                      sprintf("one-sided p = %.4f", s$monte_carlo_p),
                      paste("Gate:", if (s$gate_pass) "PASS" else "WITHHOLD")),
           col = c("#A23E48", NA, NA), lwd = c(3, NA, NA), bty = "n", cex = 0.85)
  }
  mtext("Does downstream functional broadening exceed fixed richness and label prevalence?",
        outer = TRUE, cex = 1.25, font = 2)
}
draw_figure("png")
draw_figure("svg")

assert(nrow(summary) == 2L, "Expected two trait-null summary rows")
assert(nrow(draws) == 2L * as.integer(cfg$null_draws), "Null draw denominator changed")
assert(all(draws$margin_pass), "A null margin check failed")
assert(all(gower$distance >= 0 & gower$distance <= 1), "Persisted Gower value outside [0,1]")
assert(all(file.exists(required_outputs)), "A trait-null output was not written")

cat(jsonlite::toJSON(list(
  status = "PASS", runtime_seconds = round(elapsed, 3), null_draw_rows = nrow(draws),
  primary_gate = primary$summary$gate_status, strict_gate = strict$summary$gate_status,
  claim_gate = receipt$claim_gate), auto_unbox = TRUE, pretty = TRUE), "\n")
