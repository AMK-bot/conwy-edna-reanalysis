#!/usr/bin/env Rscript

# Phase-5 cluster-aware source-positive/reprocessed-negative disagreement analysis.
# The source table is a comparator, not biological truth.

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript R/15_pipeline_disagreement.R [--rebuild-input] [--output-root=PATH] [--validate-only]\n",
    "  default          Analyse the committed portable source-positive cell table.\n",
    "  --rebuild-input  Reconstruct and overwrite that table from the local source/RDS inputs.\n",
    "  --output-root    Write generated outputs under PATH; portable input still comes from repo.\n",
    "  --validate-only  Validate existing default outputs without rewriting them.\n"))
  quit(status = 0)
}
output_arg <- grep("^--output-root=", args, value = TRUE)
if (length(output_arg) > 1L) stop("Only one --output-root may be supplied", call. = FALSE)
rebuild_input <- "--rebuild-input" %in% args
validate_only <- "--validate-only" %in% args
known <- c("--rebuild-input", "--validate-only", output_arg)
unknown <- setdiff(args, known)
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_root <- if (length(output_arg)) path.expand(sub("^--output-root=", "", output_arg[[1]])) else repo_root
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
options(stringsAsFactors = FALSE, OutDec = ".", warn = 1)

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
zscore <- function(x) {
  s <- sd(x)
  assert(is.finite(s) && s > 0, "Cannot z-score a constant/non-finite predictor")
  as.numeric((x - mean(x)) / s)
}

manifest <- jsonlite::fromJSON(in_path("data", "metadata", "impact_sprint_manifest.json"),
                               simplifyVector = TRUE)
cfg <- manifest$pipeline_disagreement
assert(as.integer(cfg$source_positive_cells) == 2344L, "Frozen disagreement denominator changed")
assert(as.integer(cfg$bootstrap_draws) == 999L, "Frozen bootstrap count changed")

paths <- list(
  portable = in_path("data", "metadata", "impact_pipeline_disagreement_cells.tsv"),
  model = out_path("data", "derived", "impact_pipeline_model.tsv"),
  bootstrap = out_path("data", "derived", "impact_pipeline_bootstrap.tsv"),
  leaveout = out_path("data", "derived", "impact_pipeline_leaveout.tsv"),
  rates = out_path("data", "derived", "impact_pipeline_rates.tsv"),
  summary = out_path("data", "derived", "impact_pipeline_summary.tsv"),
  receipt = out_path("data", "metadata", "impact_pipeline_receipt.tsv"),
  png = out_path("figures", "impact_pipeline_disagreement.png"),
  svg = out_path("figures", "impact_pipeline_disagreement.svg")
)

required_outputs <- unlist(paths[c("model", "bootstrap", "leaveout", "rates", "summary", "receipt", "png", "svg")])
if (validate_only) {
  assert(file.exists(paths$portable), "Missing portable pipeline-disagreement input")
  assert(all(file.exists(required_outputs)), "Missing pipeline-disagreement output")
  cells <- read_tsv(paths$portable)
  boot <- read_tsv(paths$bootstrap)
  loo <- read_tsv(paths$leaveout)
  assert(nrow(cells) == as.integer(cfg$source_positive_cells), "Portable cell denominator changed")
  assert(sum(cells$miss) == as.integer(cfg$nonrecovered_cells), "Portable miss denominator changed")
  assert(all(c("reprocessed_assigned_fish_read_total",
               "reprocessed_assigned_fish_reads_excluding_focal_label") %in% names(cells)),
         "Portable input lacks focal-excluded assigned-fish totals")
  assert(length(unique(boot$draw)) == as.integer(cfg$bootstrap_draws), "Bootstrap draw count changed")
  assert(all(c("focal_excluded_model_term", "focal_excluded_estimate",
               "focal_excluded_fit_ok", "focal_excluded_failure_reason") %in% names(boot)),
         "Bootstrap output lacks the paired focal-excluded sensitivity")
  assert(length(unique(loo$group_id[loo$group_type == "taxon"])) == as.integer(cfg$leaveout_taxa),
         "Taxon leave-out denominator changed")
  cat(jsonlite::toJSON(list(status = "PASS", cells = nrow(cells), bootstrap_draws = length(unique(boot$draw)),
                            leaveout_rows = nrow(loo)), auto_unbox = TRUE, pretty = TRUE), "\n")
  quit(status = 0)
}

build_portable <- function() {
  source_path <- in_path("data", "metadata", "source_data_12S_normalised.tsv")
  raw_path <- in_path("data", "derived", "species_x_sample_counts_rawnames.rds")
  cov_path <- in_path("data", "metadata", "conwy_covariates.tsv")
  assert(file.exists(source_path) && file.exists(raw_path) && file.exists(cov_path),
         "Cannot rebuild portable disagreement input: a local source/RDS input is missing")
  A <- read_tsv(source_path)
  DROP_A <- c("Gadus_macrocephalus", "Hippoglossus_stenolepis", "Ammodytes_personatus",
              "Clupea_harengus", "Clupea_pallasii", "Copadichromis_virginalis",
              "Favonigobius_gymnauchen", "Pholis_crassispina", "Melanogrammus_aeglefinus",
              "Platichthys_stellatus", "Scomber_scombrus")
  RELAB_A <- c(Salmo_obtusirostris = "Salmo_trutta",
               Anguilla_rostrata = "Anguilla_anguilla",
               Salvelinus_fontinalis_x_Salvelinus_malma = "Salvelinus_malma")
  assert(sum(DROP_A %in% names(A)) == 10L && all(names(RELAB_A) %in% names(A)),
         "Source-study cleaning contract changed")
  Am <- as.matrix(A[, -1, drop = FALSE]); storage.mode(Am) <- "double"; rownames(Am) <- A$ID
  Am <- Am[, !colnames(Am) %in% DROP_A, drop = FALSE]
  for (src in names(RELAB_A)) {
    tgt <- RELAB_A[[src]]
    Am[, tgt] <- Am[, tgt] + Am[, src]
    Am <- Am[, colnames(Am) != src, drop = FALSE]
  }
  assert(ncol(Am) == 55L, "Source cleaned universe is not 55 labels")

  M <- readRDS(raw_path)
  MAP <- c(
    Salmo_macrostigma = "Salmo_trutta", Phoxinus_phoxinus = "Phoxinus_ujmonensis",
    Salmo_salar = "Salmo_salar", Salmo_trutta = "Salmo_trutta",
    `Barbatula_aff._barbatula_lineage_I` = "Barbatula_barbatula",
    Anguilla_anguilla = "Anguilla_anguilla", Gasterosteus_aculeatus = "Gasterosteus_aculeatus",
    Pungitius_pungitius = "Pungitius_pungitius",
    Oncorhynchus_mykiss = "Oncorhynchus_mykiss_x_Salmo_salar",
    Pungitius_kaibarae = "Pungitius_pungitius", Perca_fluviatilis = "Perca_fluviatilis",
    Salvelinus_umbla = "Salvelinus_malma", Cottus_gobio = "Cottus_perifretum",
    `Phoxinus_sp._Morella` = "Phoxinus_ujmonensis", Barbatula_barbatula = "Barbatula_barbatula",
    Anguilla_rostrata = "Anguilla_anguilla", Phoxinus_bigerri = "Phoxinus_ujmonensis",
    Abramis_brama = "Abramis_brama", Salvelinus_fontinalis = "Salvelinus_malma",
    Perca_flavescens = "Perca_flavescens", Barbus_barbus = "Barbus_barbus")
  DROP_M <- c("Clupea_harengus", "Gadus_morhua", "Myoxocephalus_thompsonii", "Sprattus_sprattus",
              "Microgobius_gulosus", "Pholis_gunnellus", "Taurulus_bubalis", "Ciliata_mustela",
              "Hippoglossus_stenolepis", "Rhamdia_cf._jequitinhonha_NS_3056B", "Merluccius_merluccius",
              "Melanogrammus_aeglefinus", "Hyperoplus_lanceolatus", "Unassigned")
  assert(setequal(colnames(M), c(names(MAP), DROP_M)), "Reprocessed raw-name universe changed")
  Mk <- M[, colnames(M) %in% names(MAP), drop = FALSE]
  tgt <- MAP[colnames(Mk)]
  harmonised <- unique(unname(tgt))
  Mh <- vapply(harmonised, function(x) rowSums(Mk[, tgt == x, drop = FALSE]), numeric(nrow(Mk)))
  rownames(Mh) <- rownames(M)
  shared <- intersect(rownames(Mh), rownames(Am))
  common <- intersect(colnames(Mh), colnames(Am))
  assert(length(shared) == as.integer(cfg$shared_libraries), "Shared-library denominator changed")
  assert(length(common) == as.integer(cfg$harmonised_labels), "Common-label denominator changed")
  Mh <- Mh[shared, common, drop = FALSE]
  Am <- Am[shared, common, drop = FALSE]
  raw_depth <- rowSums(M[shared, , drop = FALSE])
  assigned_fish_columns <- setdiff(colnames(M), "Unassigned")
  assert(length(assigned_fish_columns) == 34L,
         "Assigned raw-name fish-column denominator changed")
  assigned_fish_total <- rowSums(M[shared, assigned_fish_columns, drop = FALSE])
  assert(all(assigned_fish_total >= rowSums(Mh)),
         "Harmonised focal counts exceed their assigned-fish library totals")
  source_prev <- colMeans(Am > 0)
  cov <- read_tsv(cov_path)
  assert(all(shared %in% cov$ID), "Covariates do not cover all shared libraries")
  md <- cov[match(shared, cov$ID), , drop = FALSE]

  rows <- vector("list", sum(Am > 0)); k <- 1L
  for (i in seq_along(shared)) for (j in seq_along(common)) if (Am[i, j] > 0) {
    id <- shared[[i]]; taxon <- common[[j]]
    rows[[k]] <- data.frame(
      ID = id, site = sub("_.*", "", id), timepoint = sub(".*_", "", id),
      replicate = sub(".*_(R[0-9]+)_.*", "\\1", id),
      unit = paste(sub("_.*", "", id), sub(".*_", "", id), sep = "_"),
      section = as.character(md$river_section[[i]]), season = as.character(md$season[[i]]),
      taxon = taxon, source_relative_abundance = Am[i, j],
      reprocessed_count = Mh[i, j], miss = as.integer(Mh[i, j] == 0),
      reprocessed_raw_name_depth = raw_depth[[i]],
      reprocessed_assigned_fish_read_total = assigned_fish_total[[i]],
      reprocessed_assigned_fish_reads_excluding_focal_label =
        assigned_fish_total[[i]] - Mh[i, j],
      source_label_prevalence = source_prev[[j]],
      evidence_scope = "pipeline disagreement relative to source table; source is not truth",
      stringsAsFactors = FALSE)
    k <- k + 1L
  }
  cells <- do.call(rbind, rows)
  cells <- cells[order(cells$ID, cells$taxon), , drop = FALSE]
  rownames(cells) <- NULL
  assert(nrow(cells) == as.integer(cfg$source_positive_cells), "Source-positive cell denominator changed")
  assert(sum(cells$miss) == as.integer(cfg$nonrecovered_cells), "Non-recovered cell denominator changed")
  assert(sum(cells$miss == 0) == as.integer(cfg$recovered_cells), "Recovered cell denominator changed")
  assert(all(cells$reprocessed_assigned_fish_reads_excluding_focal_label >= 0) &&
           all(cells$reprocessed_assigned_fish_reads_excluding_focal_label +
                 cells$reprocessed_count == cells$reprocessed_assigned_fish_read_total),
         "Focal-excluded assigned-fish total is not an exact part-whole removal")
  assert(length(unique(cells$taxon)) == as.integer(cfg$common_source_positive_labels),
         "Source-positive common-label denominator changed")
  write_tsv(cells, paths$portable)
  # Analyse the serialized portable representation even on rebuild so
  # --rebuild-input and default runs share the same numeric round-trip.
  read_tsv(paths$portable)
}

if (rebuild_input || !file.exists(paths$portable)) cells <- build_portable() else cells <- read_tsv(paths$portable)
assert(nrow(cells) == as.integer(cfg$source_positive_cells), "Portable cell table is not 2,344 rows")
assert(sum(cells$miss) == as.integer(cfg$nonrecovered_cells), "Portable table does not contain 382 misses")
assert(all(cells$source_relative_abundance > 0), "Portable table includes a source-negative cell")
assert(all(cells$reprocessed_raw_name_depth > 0), "Raw-name depth must be positive in shared libraries")
assert(all(c("reprocessed_assigned_fish_read_total",
             "reprocessed_assigned_fish_reads_excluding_focal_label") %in% names(cells)),
       "Portable table lacks focal-excluded assigned-fish totals; rerun with --rebuild-input")
assert(all(cells$reprocessed_assigned_fish_reads_excluding_focal_label >= 0) &&
         all(cells$reprocessed_assigned_fish_reads_excluding_focal_label +
               cells$reprocessed_count == cells$reprocessed_assigned_fish_read_total),
       "Focal-excluded assigned-fish total does not exactly remove the focal count")
assert(!any(grepl("false-negative|ground truth", cells$evidence_scope, ignore.case = TRUE)),
       "Portable evidence scope crosses the claim boundary")

cells$section <- factor(cells$section, levels = c("upper", "middle", "lower"))
cells$season <- factor(cells$season, levels = c("Spring", "Summer", "Autumn", "Winter"))
cells$z_log_source_ra <- zscore(log10(cells$source_relative_abundance))
cells$z_log_depth <- zscore(log10(1 + cells$reprocessed_raw_name_depth))
cells$z_log_focal_excluded_assigned_fish_reads <- zscore(log10(
  1 + cells$reprocessed_assigned_fish_reads_excluding_focal_label))
cells$z_label_prevalence <- zscore(cells$source_label_prevalence)

full_formula <- miss ~ z_log_source_ra + z_log_depth + z_label_prevalence + section + season
primary_formula <- full_formula
mm <- model.matrix(full_formula, data = cells)
continuous <- cells[, c("z_log_source_ra", "z_log_depth", "z_label_prevalence")]
max_abs_cor <- max(abs(cor(continuous)[upper.tri(cor(continuous))]))
condition_number <- kappa(mm[, -1, drop = FALSE], exact = TRUE)
collinearity_pass <- max_abs_cor <= as.numeric(cfg$max_absolute_predictor_correlation) &&
  condition_number <= as.numeric(cfg$max_design_condition_number)
if (!collinearity_pass) primary_formula <- miss ~ z_log_source_ra + z_log_depth + section + season
sensitivity_formula <- miss ~ z_log_source_ra +
  z_log_focal_excluded_assigned_fish_reads + section + season

safe_glm <- function(formula, data, weights = NULL, required_terms = NULL) {
  warnings <- character()
  fit <- withCallingHandlers(
    tryCatch({
      if (is.null(weights)) {
        glm(formula, family = binomial(), data = data, control = glm.control(maxit = 100))
      } else {
        data$.analysis_weight <- as.numeric(weights)
        glm(formula, family = binomial(), data = data, weights = .analysis_weight,
            control = glm.control(maxit = 100))
      }
    }, error = function(e) e),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  if (inherits(fit, "error")) return(list(ok = FALSE, reason = paste("error", conditionMessage(fit)), fit = NULL))
  co <- coef(fit)
  if (!isTRUE(fit$converged)) return(list(ok = FALSE, reason = "not_converged", fit = fit))
  if (is.null(required_terms)) required_terms <- names(co)
  if (!all(required_terms %in% names(co)) || any(!is.finite(co[required_terms]))) {
    return(list(ok = FALSE, reason = "nonfinite_required_coefficient", fit = fit))
  }
  if (any(abs(co[required_terms]) > 20)) return(list(ok = FALSE, reason = "separation_or_boundary", fit = fit))
  aliased <- names(co)[!is.finite(co)]
  reason <- if (length(aliased)) paste0("nuisance_alias:", paste(aliased, collapse = ",")) else
    if (length(warnings)) paste(unique(warnings), collapse = " | ") else "ok"
  list(ok = TRUE, reason = reason, fit = fit)
}

primary <- safe_glm(primary_formula, cells)
assert(primary$ok, paste("Primary disagreement model failed:", primary$reason))
sensitivity <- safe_glm(sensitivity_formula, cells)
assert(sensitivity$ok, paste("Focal-excluded depth sensitivity failed:", sensitivity$reason))
secondary <- if (!collinearity_pass) safe_glm(miss ~ z_label_prevalence + section + season, cells) else NULL
assert(is.null(secondary) || secondary$ok, "Predeclared prevalence-only nested model failed")

model_rows <- function(fit, model_name) {
  sm <- summary(fit)$coefficients
  data.frame(model = model_name, term = rownames(sm), estimate = sm[, "Estimate"],
             model_se = sm[, "Std. Error"], z_value = sm[, "z value"], model_p = sm[, "Pr(>|z|)"],
             stringsAsFactors = FALSE, row.names = NULL)
}
primary_model_name <- if (collinearity_pass) "full_primary" else "abundance_depth_primary"
sensitivity_model_name <- "focal_excluded_assigned_fish_sensitivity"
model_table <- model_rows(primary$fit, primary_model_name)
if (!is.null(secondary)) model_table <- rbind(model_table, model_rows(secondary$fit, "prevalence_nested"))

sites_by_section <- split(unique(cells[, c("site", "section")])$site,
                          unique(cells[, c("site", "section")])$section)
time_levels <- sort(unique(cells$timepoint))
taxa <- sort(unique(cells$taxon))
block_len <- as.integer(cfg$time_block_length)

bootstrap_weights <- function() {
  site_mult <- setNames(integer(length(unique(cells$site))), sort(unique(cells$site)))
  for (s in names(sites_by_section)) {
    pool <- sites_by_section[[s]]
    tab <- table(sample(pool, length(pool), replace = TRUE))
    site_mult[names(tab)] <- as.integer(tab)
  }
  n_blocks <- ceiling(length(time_levels) / block_len)
  starts <- sample.int(length(time_levels), n_blocks, replace = TRUE)
  sampled_time <- unlist(lapply(starts, function(st) {
    time_levels[((st - 1L + seq_len(block_len) - 1L) %% length(time_levels)) + 1L]
  }), use.names = FALSE)[seq_along(time_levels)]
  time_mult <- table(factor(sampled_time, levels = time_levels))
  tax_mult <- table(factor(sample(taxa, length(taxa), replace = TRUE), levels = taxa))
  as.numeric(site_mult[cells$site]) * as.numeric(time_mult[cells$timepoint]) *
    as.numeric(tax_mult[cells$taxon])
}

set.seed(as.integer(manifest$seeds$pipeline_bootstrap))
B <- as.integer(cfg$bootstrap_draws)
coef_names <- names(coef(primary$fit))
sensitivity_coef_names <- names(coef(sensitivity$fit))
sensitivity_term_map <- setNames(
  ifelse(coef_names == "z_log_depth",
         "z_log_focal_excluded_assigned_fish_reads", coef_names),
  coef_names)
assert(identical(unname(sensitivity_term_map), sensitivity_coef_names),
       "Primary and focal-excluded coefficient grids cannot be paired")
required_boot_terms <- intersect(c("z_log_source_ra", "z_log_depth", "z_label_prevalence"), coef_names)
required_sensitivity_boot_terms <- intersect(
  c("z_log_source_ra", "z_log_focal_excluded_assigned_fish_reads"),
  sensitivity_coef_names)
boot_rows <- vector("list", B * length(coef_names)); kk <- 1L
fit_ok <- logical(B); fit_reason <- character(B)
sensitivity_fit_ok <- logical(B); sensitivity_fit_reason <- character(B)
t0 <- proc.time()[["elapsed"]]
for (b in seq_len(B)) {
  w <- bootstrap_weights()
  fit <- safe_glm(primary_formula, cells, weights = w, required_terms = required_boot_terms)
  sensitivity_fit <- safe_glm(sensitivity_formula, cells, weights = w,
                              required_terms = required_sensitivity_boot_terms)
  fit_ok[[b]] <- fit$ok; fit_reason[[b]] <- fit$reason
  sensitivity_fit_ok[[b]] <- sensitivity_fit$ok
  sensitivity_fit_reason[[b]] <- sensitivity_fit$reason
  co <- if (is.null(fit$fit)) setNames(rep(NA_real_, length(coef_names)), coef_names) else {
    z <- coef(fit$fit); z[coef_names]
  }
  sensitivity_co <- if (is.null(sensitivity_fit$fit)) {
    setNames(rep(NA_real_, length(sensitivity_coef_names)), sensitivity_coef_names)
  } else {
    z <- coef(sensitivity_fit$fit); z[sensitivity_coef_names]
  }
  for (term in coef_names) {
    sensitivity_term <- sensitivity_term_map[[term]]
    boot_rows[[kk]] <- data.frame(draw = b, term = term, estimate = unname(co[[term]]),
                                  fit_ok = fit$ok, failure_reason = fit$reason,
                                  focal_excluded_model_term = sensitivity_term,
                                  focal_excluded_estimate = unname(sensitivity_co[[sensitivity_term]]),
                                  focal_excluded_fit_ok = sensitivity_fit$ok,
                                  focal_excluded_failure_reason = sensitivity_fit$reason,
                                  stringsAsFactors = FALSE)
    kk <- kk + 1L
  }
  if (b %% 100L == 0L) message(sprintf("[pipeline] bootstrap %d/%d", b, B))
}
bootstrap <- do.call(rbind, boot_rows)
finite_fraction <- mean(fit_ok)
sensitivity_finite_fraction <- mean(sensitivity_fit_ok)

leaveout_specs <- rbind(
  data.frame(group_type = "taxon", group_id = taxa),
  data.frame(group_type = "site", group_id = sort(unique(cells$site))),
  data.frame(group_type = "timepoint", group_id = time_levels))
leave_rows <- list(); kk <- 1L
for (i in seq_len(nrow(leaveout_specs))) {
  typ <- leaveout_specs$group_type[[i]]; id <- leaveout_specs$group_id[[i]]
  keep <- if (typ == "taxon") cells$taxon != id else if (typ == "site") cells$site != id else cells$timepoint != id
  fit <- safe_glm(primary_formula, cells[keep, , drop = FALSE])
  sensitivity_fit <- safe_glm(sensitivity_formula, cells[keep, , drop = FALSE])
  co <- if (is.null(fit$fit)) setNames(rep(NA_real_, length(coef_names)), coef_names) else coef(fit$fit)[coef_names]
  sensitivity_co <- if (is.null(sensitivity_fit$fit)) {
    setNames(rep(NA_real_, length(sensitivity_coef_names)), sensitivity_coef_names)
  } else {
    coef(sensitivity_fit$fit)[sensitivity_coef_names]
  }
  for (term in coef_names) {
    sensitivity_term <- sensitivity_term_map[[term]]
    leave_rows[[kk]] <- data.frame(group_type = typ, group_id = id, term = term,
                                   estimate = unname(co[[term]]), fit_ok = fit$ok,
                                   failure_reason = fit$reason,
                                   focal_excluded_model_term = sensitivity_term,
                                   focal_excluded_estimate = unname(sensitivity_co[[sensitivity_term]]),
                                   focal_excluded_fit_ok = sensitivity_fit$ok,
                                   focal_excluded_failure_reason = sensitivity_fit$reason,
                                   stringsAsFactors = FALSE)
    kk <- kk + 1L
  }
}
leaveout <- do.call(rbind, leave_rows)

terms_of_interest <- intersect(c("z_log_source_ra", "z_log_depth", "z_label_prevalence"), coef_names)
effect_rows <- lapply(terms_of_interest, function(term) {
  z <- bootstrap$estimate[bootstrap$term == term & bootstrap$fit_ok & is.finite(bootstrap$estimate)]
  loo <- leaveout$estimate[leaveout$term == term & leaveout$fit_ok & is.finite(leaveout$estimate)]
  estimate <- coef(primary$fit)[[term]]
  lo <- if (length(z)) unname(quantile(z, 0.025)) else NA_real_
  hi <- if (length(z)) unname(quantile(z, 0.975)) else NA_real_
  sign_stability <- if (length(loo)) mean(sign(loo) == sign(estimate)) else 0
  no_reversal <- length(loo) == nrow(leaveout_specs) && all(sign(loo) == sign(estimate))
  gate <- finite_fraction >= as.numeric(cfg$min_finite_bootstrap_fraction) &&
    is.finite(lo) && is.finite(hi) && (lo > 0 || hi < 0) &&
    sign_stability >= as.numeric(cfg$min_leaveout_sign_stability) && no_reversal
  data.frame(term = term, estimate = estimate, bootstrap_q025 = lo, bootstrap_q975 = hi,
             bootstrap_finite_fraction = finite_fraction,
             leaveout_min = if (length(loo)) min(loo) else NA_real_,
             leaveout_max = if (length(loo)) max(loo) else NA_real_,
             leaveout_sign_stability = sign_stability, no_leaveout_reversal = no_reversal,
             gate_pass = gate, gate_status = if (gate) "ASSOCIATION_PASSED" else "DESCRIPTIVE_ONLY",
             stringsAsFactors = FALSE)
})
effects <- do.call(rbind, effect_rows)
sensitivity_terms_of_interest <- intersect(
  c("z_log_source_ra", "z_log_focal_excluded_assigned_fish_reads"),
  sensitivity_coef_names)
sensitivity_effect_rows <- lapply(sensitivity_terms_of_interest, function(term) {
  z <- bootstrap$focal_excluded_estimate[
    bootstrap$focal_excluded_model_term == term & bootstrap$focal_excluded_fit_ok &
      is.finite(bootstrap$focal_excluded_estimate)]
  loo <- leaveout$focal_excluded_estimate[
    leaveout$focal_excluded_model_term == term & leaveout$focal_excluded_fit_ok &
      is.finite(leaveout$focal_excluded_estimate)]
  estimate <- coef(sensitivity$fit)[[term]]
  lo <- if (length(z)) unname(quantile(z, 0.025)) else NA_real_
  hi <- if (length(z)) unname(quantile(z, 0.975)) else NA_real_
  sign_stability <- if (length(loo)) mean(sign(loo) == sign(estimate)) else 0
  no_reversal <- length(loo) == nrow(leaveout_specs) && all(sign(loo) == sign(estimate))
  gate <- sensitivity_finite_fraction >= as.numeric(cfg$min_finite_bootstrap_fraction) &&
    is.finite(lo) && is.finite(hi) && (lo > 0 || hi < 0) &&
    sign_stability >= as.numeric(cfg$min_leaveout_sign_stability) && no_reversal
  data.frame(term = term, estimate = estimate, bootstrap_q025 = lo, bootstrap_q975 = hi,
             bootstrap_finite_fraction = sensitivity_finite_fraction,
             leaveout_min = if (length(loo)) min(loo) else NA_real_,
             leaveout_max = if (length(loo)) max(loo) else NA_real_,
             leaveout_sign_stability = sign_stability, no_leaveout_reversal = no_reversal,
             gate_pass = gate, gate_status = if (gate) "ASSOCIATION_PASSED" else "DESCRIPTIVE_ONLY",
             stringsAsFactors = FALSE)
})
sensitivity_effects <- do.call(rbind, sensitivity_effect_rows)
model_table <- merge(model_table, effects, by = "term", all.x = TRUE, sort = FALSE)
sensitivity_model_table <- merge(
  model_rows(sensitivity$fit, sensitivity_model_name), sensitivity_effects,
  by = "term", all.x = TRUE, sort = FALSE)
model_table <- rbind(model_table, sensitivity_model_table)

rate_group <- function(var) {
  f <- factor(cells[[var]], levels = sort(unique(cells[[var]])))
  n <- as.vector(table(f)); miss <- as.vector(tapply(cells$miss, f, sum))
  data.frame(group_type = var, group_id = levels(f), source_positive_cells = n,
             nonrecovered_cells = miss, nonrecovery_fraction = miss / n, stringsAsFactors = FALSE)
}
rates <- rbind(rate_group("section"), rate_group("season"), rate_group("taxon"))

primary_depth_effect <- effects[effects$term == "z_log_depth", , drop = FALSE]
sensitivity_depth_effect <- sensitivity_effects[
  sensitivity_effects$term == "z_log_focal_excluded_assigned_fish_reads", , drop = FALSE]
assert(nrow(primary_depth_effect) == 1L && nrow(sensitivity_depth_effect) == 1L,
       "Depth effect rows are not unique")
depth_survives_focal_exclusion <- sensitivity_depth_effect$gate_pass[[1]]
depth_sensitivity_status <- if (depth_survives_focal_exclusion) {
  "SURVIVES_FOCAL_EXCLUSION"
} else {
  "DOES_NOT_SURVIVE_FOCAL_EXCLUSION"
}

summary <- data.frame(
  analysis = "cluster-aware pipeline disagreement",
  shared_libraries = length(unique(cells$ID)), harmonised_labels = as.integer(cfg$harmonised_labels),
  source_positive_labels = length(unique(cells$taxon)), source_positive_cells = nrow(cells),
  recovered_cells = sum(cells$miss == 0), nonrecovered_cells = sum(cells$miss),
  max_absolute_continuous_predictor_correlation = max_abs_cor,
  design_condition_number = condition_number, collinearity_gate_pass = collinearity_pass,
  primary_model = paste(deparse(primary_formula), collapse = " "),
  primary_depth_definition = paste0(
    "sum of all 35 raw-name columns, including the focal harmonised-label count and Unassigned; ",
    "retained as the preregistered primary predictor"),
  bootstrap_draws = B, finite_bootstrap_fits = sum(fit_ok), finite_bootstrap_fraction = finite_fraction,
  leaveout_groups = nrow(leaveout_specs), promoted_associations = sum(effects$gate_pass),
  sensitivity_model = paste(deparse(sensitivity_formula), collapse = " "),
  sensitivity_depth_definition = paste0(
    "sum of 34 assigned raw-name fish columns excluding Unassigned, minus every raw-name count ",
    "mapped to the focal harmonised label"),
  sensitivity_bootstrap_draws = B,
  sensitivity_finite_bootstrap_fits = sum(sensitivity_fit_ok),
  sensitivity_finite_bootstrap_fraction = sensitivity_finite_fraction,
  sensitivity_leaveout_groups = nrow(leaveout_specs),
  sensitivity_supported_associations = sum(sensitivity_effects$gate_pass),
  depth_association_survives_focal_exclusion = depth_survives_focal_exclusion,
  depth_association_sensitivity_status = depth_sensitivity_status,
  focal_excluded_depth_estimate = sensitivity_depth_effect$estimate[[1]],
  focal_excluded_depth_bootstrap_q025 = sensitivity_depth_effect$bootstrap_q025[[1]],
  focal_excluded_depth_bootstrap_q975 = sensitivity_depth_effect$bootstrap_q975[[1]],
  focal_excluded_depth_leaveout_sign_stability =
    sensitivity_depth_effect$leaveout_sign_stability[[1]],
  focal_excluded_depth_no_leaveout_reversal =
    sensitivity_depth_effect$no_leaveout_reversal[[1]],
  overall_gate = if (any(sensitivity_effects$gate_pass)) "ASSOCIATION_EVIDENCE" else
    "DESCRIPTIVE_DISAGREEMENT_ONLY",
  evidence_scope = "pipeline non-recovery relative to source table; neither pipeline is biological truth",
  stringsAsFactors = FALSE)
receipt <- data.frame(
  seed = as.integer(manifest$seeds$pipeline_bootstrap), bootstrap_draws = B,
  site_resampling = "within section", time_resampling = "circular moving blocks length 3",
  taxon_resampling = "whole taxon clusters", source_positive_cells = nrow(cells),
  nonrecovered_cells = sum(cells$miss), canonical_runtime_excluded = TRUE,
  primary_depth_part_whole_coupling = "includes focal reprocessed count used to define miss",
  sensitivity_bootstrap_draws = B,
  sensitivity_finite_bootstrap_fits = sum(sensitivity_fit_ok),
  sensitivity_leaveout_groups = nrow(leaveout_specs),
  depth_association_sensitivity_status = depth_sensitivity_status,
  gate_status = summary$overall_gate, stringsAsFactors = FALSE)

write_tsv(model_table, paths$model)
write_tsv(bootstrap, paths$bootstrap)
write_tsv(leaveout, paths$leaveout)
write_tsv(rates, paths$rates)
write_tsv(summary, paths$summary)
write_tsv(receipt, paths$receipt)

draw_figure <- function(device) {
  if (device == "png") {
    dir.create(dirname(paths$png), recursive = TRUE, showWarnings = FALSE)
    png(paths$png, width = 1800, height = 1400, res = 150)
  } else {
    dir.create(dirname(paths$svg), recursive = TRUE, showWarnings = FALSE)
    svg(paths$svg, width = 12, height = 9.3)
  }
  layout(matrix(1:2, ncol = 1), heights = c(1.15, 1.25))
  op <- par(oma = c(0, 0, 2.0, 0), family = "sans")
  on.exit({par(op); dev.off()}, add = TRUE)
  primary_figure <- effects[effects$term %in% c("z_log_source_ra", "z_log_depth"), , drop = FALSE]
  primary_figure$display_label <- c(
    z_log_source_ra = "Primary: source relative abundance",
    z_log_depth = "Primary: 35-column raw-name total"
  )[primary_figure$term]
  sensitivity_figure <- sensitivity_effects[
    sensitivity_effects$term %in%
      c("z_log_source_ra", "z_log_focal_excluded_assigned_fish_reads"), , drop = FALSE]
  sensitivity_figure$display_label <- c(
    z_log_source_ra = "Sensitivity: source relative abundance",
    z_log_focal_excluded_assigned_fish_reads =
      "Sensitivity: focal-excluded fish total"
  )[sensitivity_figure$term]
  ef <- rbind(primary_figure, sensitivity_figure)
  ef <- ef[nrow(ef):1, , drop = FALSE]
  y <- seq_len(nrow(ef))
  xlim <- range(c(ef$bootstrap_q025, ef$bootstrap_q975, ef$leaveout_min, ef$leaveout_max, 0), na.rm = TRUE)
  par(mar = c(4.5, 12.5, 4.2, 2.0))
  plot(ef$estimate, y, xlim = xlim, ylim = c(0.5, nrow(ef) + 0.5), yaxt = "n",
       pch = 21, bg = ifelse(ef$gate_pass, "#2C7A5B", "#D58B32"), cex = 1.5,
       xlab = "Log-odds coefficient (z-scored predictor)", ylab = "",
       main = "Primary and focal-excluded cluster-aware associations")
  axis(2, at = y, labels = ef$display_label, las = 1, cex.axis = 0.82)
  segments(ef$bootstrap_q025, y, ef$bootstrap_q975, y, lwd = 4, col = "#3A506B")
  segments(ef$leaveout_min, y - 0.13, ef$leaveout_max, y - 0.13, lwd = 1.5, col = "#A23E48")
  abline(v = 0, lty = 2, col = "#777777")
  mtext("Primary raw-name total contains focal count + Unassigned; sensitivity excludes both.",
        side = 3, line = 0.45, cex = 0.72)
  legend("topright",
         c("95% block-bootstrap interval", "leave-one-group range", "gate passed", "gate not passed"),
         col = c("#3A506B", "#A23E48", "#2C7A5B", "#D58B32"),
         lwd = c(4, 1.5, NA, NA), pch = c(NA, NA, 21, 21),
         pt.bg = c(NA, NA, "#2C7A5B", "#D58B32"), bty = "n", cex = 0.75)

  tr <- rates[rates$group_type == "taxon", ]
  tr <- tr[order(tr$nonrecovery_fraction), ]
  taxon_labels <- gsub("_", " ", tr$group_id)
  taxon_labels[taxon_labels == "Oncorhynchus mykiss x Salmo salar"] <- "Oncorhynchus mykiss proxy"
  par(mar = c(4.5, 12.0, 3.2, 2.0))
  bp <- barplot(tr$nonrecovery_fraction, names.arg = taxon_labels, horiz = TRUE,
                las = 1, col = "#6F8FA6", border = NA, xlim = c(0, max(tr$nonrecovery_fraction) * 1.2),
                xlab = "Non-recovery fraction among source-positive cells", main = "Observed comparator disagreement by label")
  text(tr$nonrecovery_fraction, bp, labels = sprintf("%d/%d", tr$nonrecovered_cells, tr$source_positive_cells),
       pos = 4, cex = 0.75)
  mtext("Where do the source and reprocessed tables disagree?", outer = TRUE, cex = 1.25, font = 2)
}
draw_figure("png")
draw_figure("svg")

elapsed <- proc.time()[["elapsed"]] - t0
assert(nrow(cells) == 2344L && sum(cells$miss) == 382L, "Final disagreement denominator drift")
assert(length(unique(bootstrap$draw)) == B, "Bootstrap draw denominator drift")
assert(all(table(bootstrap$focal_excluded_model_term) == B),
       "Focal-excluded bootstrap coefficient grid drift")
assert(nrow(leaveout_specs) == 43L, "Leave-out group denominator drift")
assert(all(leaveout$focal_excluded_fit_ok), "A focal-excluded leave-out fit failed")
assert(all(file.exists(required_outputs)), "A pipeline-disagreement output was not written")
cat(jsonlite::toJSON(list(status = "PASS", runtime_seconds = round(elapsed, 3),
                          finite_bootstrap_fraction = finite_fraction,
                          sensitivity_finite_bootstrap_fraction = sensitivity_finite_fraction,
                          depth_association_sensitivity_status = depth_sensitivity_status,
                          promoted_associations = sum(effects$gate_pass),
                          overall_gate = summary$overall_gate), auto_unbox = TRUE, pretty = TRUE), "\n")
