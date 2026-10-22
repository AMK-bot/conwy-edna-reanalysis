#!/usr/bin/env Rscript
# Unit, invariant, provenance and regression tests for the synthetic occupancy
# recoverability atlas. Run --unit-only before the full analysis if desired.

test_start <- proc.time()[["elapsed"]]
tests_passed <- 0L

pass <- function(label) {
  tests_passed <<- tests_passed + 1L
  cat(sprintf("PASS %02d: %s\n", tests_passed, label))
  invisible(TRUE)
}

check <- function(condition, label) {
  if (length(condition) != 1L || is.na(condition) || !isTRUE(condition))
    stop("FAIL: ", label, call. = FALSE)
  pass(label)
}

check_equal <- function(actual, expected, tolerance = 0, label) {
  okay <- if (is.numeric(actual) || is.numeric(expected)) {
    length(actual) == length(expected) &&
      all(is.na(actual) == is.na(expected)) &&
      all(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) <= tolerance)
  } else identical(actual, expected)
  check(okay, label)
}

Sys.setenv(IMPACT_ATLAS_SOURCE_ONLY = "1")
source("R/11_occupancy_recoverability_atlas.R")
unit_only <- "--unit-only" %in% commandArgs(trailingOnly = TRUE)

# ---- frozen inputs and scenario construction -------------------------------
atlas_reset_input_log()
config <- load_atlas_config()
log <- atlas_input_log()
check_equal(log$scientific_inputs_read,
            sort(c("data/metadata/impact_sprint_manifest.json",
                   "docs/impact-sprint-preregistration.md")),
            label = "only the two frozen scientific inputs are read")
check(length(log$generated_artifacts_read) == 0L,
      "config load reads no generated or real Layer-A artifact")

scenarios <- build_atlas_scenarios(config)
check(nrow(scenarios) == 648L &&
        length(unique(scenarios$truth_scenario_id)) == 216L,
      "frozen grid is 648 mask scenarios over 216 paired truth scenarios")
check(!anyDuplicated(scenarios$scenario_id), "scenario IDs are unique")
check(all(abs(scenarios$calibrated_mean_occupancy -
                scenarios$target_mean_occupancy) < 1e-10),
      "gradient intercept calibration reproduces every target mean")
check(all(table(scenarios$truth_scenario_id) == 3L),
      "each truth scenario has exactly three paired mask encodings")
check(all(vapply(split(scenarios$mask_scenario, scenarios$truth_scenario_id),
                 function(x) identical(as.character(x),
                                       as.character(config$occ$mask_scenario)),
                 logical(1))), "mask order is frozen and complete within each pair")

# ---- deterministic generator and explicit mask semantics -------------------
core <- scenarios[!duplicated(scenarios$truth_scenario_id), , drop = FALSE]
fixture_core <- subset(core,
                       target_mean_occupancy == 0.35 &
                         prevalence_structure == "downstream_gradient" &
                         baseline_detection == 0.5 & replicates == 3L &
                         log_depth_sd == 0.75 & depth_river_correlation == 0.7)[1, ]
sim_a <- simulate_atlas_truth(fixture_core, 17L, config)
sim_b <- simulate_atlas_truth(fixture_core, 17L, config)
check(identical(sim_a, sim_b), "fixed seed reproduces the complete latent generator state")
check_equal(sim_a$simulation_seed, 2026280918L,
            label = "derived seed regression value is stable")
check_equal(sum(sim_a$latent_occupancy), 47L,
            label = "latent-occupancy regression count is stable")
check_equal(sum(sim_a$generated_detection), 87L,
            label = "generated-detection regression count is stable")
check_equal(as.vector(sim_a$planned_log_depth)[1:8],
            c(-1.82050246, -1.56459320, -0.11239406, -0.18837457,
              -0.49976540, -0.69782152, -1.03818511, -0.74585942),
            tolerance = 5e-9, label = "planned-depth numeric regression fixture is stable")
check_equal(stats::sd(sim_a$planned_log_depth), 0.75, tolerance = 1e-12,
            label = "realised log-depth SD equals its declared level exactly")
check_equal(stats::cor(as.vector(sim_a$planned_log_depth),
                       rep(sim_a$river_position, times = 3L)),
            0.7, tolerance = 1e-12,
            label = "realised depth-river correlation equals its declared level exactly")

complete <- apply_atlas_mask(sim_a, "complete_correct")
correct <- apply_atlas_mask(sim_a, "structured_correct")
zero_fill <- apply_atlas_mask(sim_a, "structured_zero_fill")
check(all(complete$surveyed_mask) && all(complete$analysis_mask),
      "complete design surveys and analyses every planned cell")
check(sum(!correct$surveyed_mask) == 72L &&
        sum(correct$surveyed_mask) == 288L,
      "structured missingness removes exactly 20 percent of K=3 cells")
check(identical(correct$analysis_mask, correct$surveyed_mask),
      "correct structured mask keeps unavailable cells out of analysis")
check(all(zero_fill$analysis_mask) &&
        identical(zero_fill$surveyed_mask, correct$surveyed_mask),
      "zero-fill uses the exact same unavailability pattern but marks it analysable")
check(all(is.na(correct$analysis_detection[!correct$surveyed_mask])) &&
        all(zero_fill$analysis_detection[!zero_fill$surveyed_mask] == 0L),
      "structural absence remains NA unless the declared mask error zero-fills it")
check(sum(correct$observed_detection, na.rm = TRUE) == 79L &&
        sum(zero_fill$analysis_detection) == 79L,
      "paired mask regression preserves observed positives and changes denominators only")

fixture <- make_generator_fixture(sim_a, zero_fill, 999L,
                                  "structured_zero_fill")
required_fields <- c("true_occupancy_probability", "latent_occupancy",
                     "surveyed_mask", "analysis_mask", "generated_detection",
                     "observed_detection", "analysis_detection",
                     "erroneous_zero_fill")
check(all(required_fields %in% names(fixture)) && nrow(fixture) == 360L,
      "portable fixture persists separate truth, masks and detection fields")
check(all(fixture$data_scope == ATLAS_SCOPE), "every generator fixture row is scope-tagged")

# ---- estimators, derivatives, identifiability and gate boundaries -----------
check_equal(wilson_interval(0, 10), c(0, 0.2775327998628892), tolerance = 1e-12,
            label = "Wilson interval regression value is stable")
naive <- fit_naive_occurrence(correct)
check_equal(naive$estimate, 0.3652173913043478, tolerance = 1e-12,
            label = "naive occurrence regression estimate is stable")
check_equal(c(naive$lower95, naive$upper95),
            c(0.282897577416, 0.456250693670), tolerance = 1e-11,
            label = "naive Wilson interval regression values are stable")

ensure_atlas_backend()
likelihood_data <- prepare_likelihood_data(sim_a, correct, "depth_aware")
probe <- c(-0.5, 1.0, 0.2, 0.75)
analytic <- atlas_llgh_cpp(probe, likelihood_data$x, likelihood_data$w,
                           likelihood_data$y, likelihood_data$unit, TRUE)
fn <- function(x) atlas_llgh_cpp(x, likelihood_data$x, likelihood_data$w,
                                 likelihood_data$y, likelihood_data$unit,
                                 FALSE)$loglik
numeric_gradient <- numDeriv::grad(fn, probe)
numeric_hessian <- numDeriv::hessian(fn, probe)
check(max(abs(analytic$gradient - numeric_gradient)) < 1e-6,
      "analytic integrated-likelihood score matches numerical differentiation")
check(max(abs(analytic$hessian - numeric_hessian)) < 1e-5,
      "analytic observed Hessian matches numerical differentiation")

rep_fit_a <- fit_integrated_likelihood(sim_a, correct, "replication_aware")
rep_fit_b <- fit_integrated_likelihood(sim_a, correct, "replication_aware")
depth_fit <- fit_integrated_likelihood(sim_a, correct, "depth_aware")
check_equal(rep_fit_a$estimate, rep_fit_b$estimate, tolerance = 0,
            label = "bounded multistart likelihood fit is bitwise deterministic")
check(rep_fit_a$successful_starts == 3L && depth_fit$successful_starts == 3L,
      "both correction estimators execute all three deterministic starts")
check_equal(rep_fit_a$estimate, 0.407190545851, tolerance = 2e-8,
            label = "replication-aware fit regression estimate is stable")
check_equal(depth_fit$estimate, 0.394084081300, tolerance = 2e-8,
            label = "depth-aware fit regression estimate is stable")

k1_core <- core[core$replicates == 1L, , drop = FALSE][1, ]
k1_sim <- simulate_atlas_truth(k1_core, 1L, config)
k1_mask <- apply_atlas_mask(k1_sim, "complete_correct")
k1_fit <- fit_integrated_likelihood(k1_sim, k1_mask, "depth_aware")
check(k1_fit$structural_nonidentifiable && k1_fit$fit_failure &&
        k1_fit$interval_failure && is.na(k1_fit$estimate),
      "one-replicate correction is labelled structurally non-identifiable without fitting")

gate <- config$occ$gate
check(classify_atlas_gate(0.05, 0.10, 0.88, 0.05, 0.05, 0.05, 0, gate) ==
        "recoverable", "inclusive gate boundaries classify as recoverable")
check(classify_atlas_gate(0.0501, 0.10, 0.95, 0, 0, 0, 0, gate) ==
        "unreliable", "bias just beyond the gate classifies as unreliable")
check(classify_atlas_gate(0, 0, 0.95, 0.2, 0.2, 0.2, 0, gate) ==
        "non-identifiable", "20-percent identification-failure boundary is inclusive")
check(classify_atlas_gate(NA, NA, NA, 1, 1, 1, 1, gate) ==
        "non-identifiable", "structural non-identifiability takes gate precedence")

# Miniature end-to-end invariant exercise. This does not replace or shrink the
# canonical run; it tests assembly on two K=1 and two K=2 paired truth scenarios.
mini_ids <- c(head(core$truth_scenario_id[core$replicates == 1L], 2L),
              head(core$truth_scenario_id[core$replicates == 2L], 2L))
mini_scenarios <- scenarios[scenarios$truth_scenario_id %in% mini_ids, ]
mini <- run_occupancy_atlas(config, mini_scenarios, n_simulations = 2L, cores = 1L)
check(nrow(mini$fits) == 4L * 3L * 2L * 3L &&
        nrow(mini$atlas) == 4L * 3L * 3L,
      "miniature assembly preserves scenario x simulation x estimator dimensions")
check(all(mini$fits$data_scope == ATLAS_SCOPE) &&
        all(mini$atlas$data_scope == ATLAS_SCOPE),
      "miniature result rows remain explicitly synthetic")
check(all(mini$atlas$simulations == 2L) &&
        all(mini$atlas$coverage_denominator <= mini$atlas$simulations),
      "coverage denominators are conditional and never exceed simulations")

# Static I/O guard complements the runtime audit: the analysis source contains
# no generic reader capable of loading a real count/mask object.
script_text <- readLines("R/11_occupancy_recoverability_atlas.R", warn = FALSE)
forbidden_reader <- "(readRDS|read\\.csv|read\\.delim|fread)[[:space:]]*\\("
check(!any(grepl(forbidden_reader, script_text, ignore.case = TRUE)),
      "static I/O guard finds no RDS/CSV/TSV data reader in the analysis")

if (unit_only) {
  elapsed <- proc.time()[["elapsed"]] - test_start
  cat(sprintf("\nALL %d UNIT/INVARIANT TESTS PASSED in %.3f seconds\n",
              tests_passed, elapsed))
  quit(save = "no", status = 0)
}

# ---- canonical full-run output and provenance checks -----------------------
paths <- atlas_output_paths()
required_outputs <- c(paths,
                      impact_manifest = "data/derived/impact_occupancy_manifest.json",
                      receipt = "data/derived/impact_occupancy_receipt.json")
check(all(file.exists(vapply(required_outputs, atlas_path, character(1)))),
      "all canonical tables, figures, manifest and execution receipt exist")
if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table required for full checks")
read_tsv <- function(relative_path) data.table::fread(
  atlas_path(relative_path), data.table = FALSE, na.strings = "NA",
  showProgress = FALSE)
scenario_out <- read_tsv(paths[["scenarios"]])
draw_out <- read_tsv(paths[["simulation_audit"]])
fit_out <- read_tsv(paths[["fit_results"]])
atlas_out <- read_tsv(paths[["atlas"]])
fixture_out <- read_tsv(paths[["generator_fixture"]])
check(nrow(scenario_out) == 648L && nrow(draw_out) == 64800L &&
        nrow(fit_out) == 194400L && nrow(atlas_out) == 1944L,
      "full declared 648 x 100 x 3 analysis dimensions are exact")
check(!anyDuplicated(paste(fit_out$scenario_id, fit_out$simulation_id,
                           fit_out$estimator)) &&
        !anyDuplicated(paste(atlas_out$scenario_id, atlas_out$estimator)),
      "fit-level and atlas-level keys are unique")
check(all(vapply(list(scenario_out, draw_out, fit_out, atlas_out, fixture_out),
                 function(x) all(x$data_scope == ATLAS_SCOPE), logical(1))),
      "every canonical result and audit row is tagged synthetic known truth")
check(all(table(fit_out$scenario_id, fit_out$estimator) == 100L),
      "every scenario-estimator cell retains all 100 simulations")
check(all(atlas_out$simulations == 100L) &&
        all(atlas_out$estimate_denominator <= 100L) &&
        all(atlas_out$coverage_denominator <= atlas_out$estimate_denominator) &&
        all(atlas_out$coverage_numerator <= atlas_out$coverage_denominator),
      "estimate, interval and coverage denominators are nested and explicit")
check(all(atlas_out$fit_failure_n + (100L - atlas_out$fit_failure_n) == 100L) &&
        all(atlas_out$interval_failure_n +
              (100L - atlas_out$interval_failure_n) == 100L),
      "fit and interval failures remain in unconditional denominators")

k1_summary <- atlas_out[atlas_out$replicates == 1L, ]
check(all(k1_summary$gate_class[k1_summary$estimator != "naive_occurrence"] ==
            "non-identifiable") &&
        all(k1_summary$structural_nonidentifiable_rate[
          k1_summary$estimator != "naive_occurrence"] == 1),
      "all one-replicate correction cells remain structurally non-identifiable")
check(all(k1_summary$structural_nonidentifiable_rate[
  k1_summary$estimator == "naive_occurrence"] == 0),
      "one-replicate rule does not incorrectly suppress naive occurrence")

# Recompute raw failure/coverage denominators independently from fit rows.
fit_dt <- data.table::as.data.table(fit_out)
audit <- fit_dt[, .(
  fit_failure_n_check = sum(fit_failure),
  interval_failure_n_check = sum(interval_failure),
  boundary_failure_n_check = sum(boundary_failure),
  singularity_failure_n_check = sum(singularity_failure),
  structural_n_check = sum(structural_nonidentifiable),
  coverage_denominator_check = sum(!interval_failure & !fit_failure & !is.na(covered)),
  coverage_numerator_check = sum(covered[!interval_failure & !fit_failure], na.rm = TRUE)
), by = .(scenario_id, estimator)]
atlas_dt <- data.table::as.data.table(atlas_out)
joined <- merge(atlas_dt, audit, by = c("scenario_id", "estimator"), sort = FALSE)
check(all(joined$fit_failure_n == joined$fit_failure_n_check) &&
        all(joined$interval_failure_n == joined$interval_failure_n_check) &&
        all(joined$boundary_failure_n == joined$boundary_failure_n_check) &&
        all(joined$singularity_failure_n == joined$singularity_failure_n_check) &&
        all(joined$structural_nonidentifiable_n == joined$structural_n_check) &&
        all(joined$coverage_denominator == joined$coverage_denominator_check) &&
        all(joined$coverage_numerator == joined$coverage_numerator_check),
      "summary failure and conditional-coverage denominators reproduce raw fit rows")

reclassified <- mapply(
  classify_atlas_gate, atlas_out$absolute_bias, atlas_out$rmse, atlas_out$coverage,
  atlas_out$fit_failure_rate, atlas_out$interval_failure_rate,
  atlas_out$combined_identification_failure_rate,
  atlas_out$structural_nonidentifiable_rate,
  MoreArgs = list(gate = config$occ$gate), USE.NAMES = FALSE)
check(identical(reclassified, atlas_out$gate_class),
      "every atlas cell reproduces the frozen gate classification")

# Paired-mask audit: generator outcomes match within every truth draw; only the
# surveyed/analysis masks and zero-fill denominator change.
draw_joined <- merge(draw_out, scenario_out[, c("scenario_id", "mask_scenario")],
                     by = "scenario_id", sort = FALSE)
wide <- reshape(draw_joined[, c(
  "truth_scenario_id", "simulation_id", "mask_scenario", "simulation_seed",
  "true_marginal_occupancy", "realised_latent_prevalence", "truth_occupied_units",
  "generated_detection_events", "surveyed_mask_cells", "analysis_mask_cells",
  "erroneous_zero_filled_cells")],
  idvar = c("truth_scenario_id", "simulation_id"), timevar = "mask_scenario",
  direction = "wide")
paired_prefixes <- c("simulation_seed", "true_marginal_occupancy",
                     "realised_latent_prevalence", "truth_occupied_units",
                     "generated_detection_events")
paired_ok <- vapply(paired_prefixes, function(prefix) {
  a <- wide[[paste0(prefix, ".complete_correct")]]
  b <- wide[[paste0(prefix, ".structured_correct")]]
  c <- wide[[paste0(prefix, ".structured_zero_fill")]]
  isTRUE(all.equal(a, b, tolerance = 0)) && isTRUE(all.equal(b, c, tolerance = 0))
}, logical(1))
check(all(paired_ok), "all three masks share identical truth, depth seed and detections")
check(all(wide$surveyed_mask_cells.structured_correct ==
            wide$surveyed_mask_cells.structured_zero_fill) &&
        all(wide$analysis_mask_cells.structured_zero_fill ==
              wide$analysis_mask_cells.complete_correct) &&
        all(wide$erroneous_zero_filled_cells.structured_zero_fill > 0),
      "paired structured masks differ only through the declared zero-fill analysis mask")

manifest <- jsonlite::fromJSON(atlas_path(
  "data/derived/impact_occupancy_manifest.json"), simplifyVector = TRUE)
receipt <- jsonlite::fromJSON(atlas_path(
  "data/derived/impact_occupancy_receipt.json"), simplifyVector = TRUE)
check(length(manifest$layerA_real_inputs_read) == 0L &&
        identical(manifest$data_scope, ATLAS_SCOPE) &&
        length(manifest$scientific_inputs_read) == 2L,
      "manifest proves zero real Layer-A inputs and exactly two frozen inputs")
manifest_hashes <- unlist(manifest$canonical_output_sha256, use.names = TRUE)
actual_hashes <- vapply(names(paths), function(name)
  atlas_sha256(atlas_path(paths[[name]])), character(1))
check(identical(unname(manifest_hashes[names(actual_hashes)]),
                unname(actual_hashes)),
      "manifest SHA-256 values reproduce every canonical table and figure")
check(isTRUE(receipt$full_design_executed) &&
        isTRUE(receipt$deterministic_rebuild_requested) &&
        isTRUE(receipt$deterministic_rebuild_pass),
      "execution receipt proves the full design and second deterministic build")
check(length(receipt$layerA_real_inputs_read) == 0L,
      "execution receipt independently records no real Layer-A input")

if (!requireNamespace("png", quietly = TRUE)) stop("png package required for figure check")
image <- png::readPNG(atlas_path(paths[["figure_png"]]), native = FALSE)
check(dim(image)[[1]] >= 2500L && dim(image)[[2]] >= 3800L,
      "PNG atlas has high-information publication-scale dimensions")
svg_text <- readLines(atlas_path(paths[["figure_svg"]]), warn = FALSE)
check(length(svg_text) > 100L && any(grepl("<svg", svg_text, fixed = TRUE)),
      "SVG atlas exists as a substantive vector rendering")

elapsed <- proc.time()[["elapsed"]] - test_start
test_receipt <- data.frame(
  data_scope = ATLAS_SCOPE, tests_passed = tests_passed,
  unit_invariant_regression_tests = "PASS",
  full_dimension_denominator_provenance_tests = "PASS",
  deterministic_rebuild = "PASS",
  test_runtime_seconds = sprintf("%.3f", elapsed), stringsAsFactors = FALSE)
write.table(test_receipt,
            "data/derived/impact_occupancy_test_receipt.tsv", sep = "\t",
            quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("\nALL %d ATLAS TESTS PASSED in %.3f seconds\n", tests_passed, elapsed))
