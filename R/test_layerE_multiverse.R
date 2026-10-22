#!/usr/bin/env Rscript

# Unit, invariant, and regression tests for R/13_layerE_multiverse.R.
# Run after two complete analysis builds so the determinism receipt can be tested.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)
source("R/13_layerE_multiverse.R")

test_count <- 0L
check <- function(ok, message) {
  test_count <<- test_count + 1L
  if (length(ok) != 1L || is.na(ok) || !ok) {
    stop(sprintf("Layer-E test %d failed: %s", test_count, message), call. = FALSE)
  }
}

paths <- c(
  manifest = "data/metadata/impact_sprint_manifest.json",
  input = "data/derived/impact_library_counts.tsv",
  specifications = "data/derived/impact_layerE_specifications.tsv",
  library_draws = "data/derived/impact_layerE_library_draws.tsv",
  draw_results = "data/derived/impact_layerE_draw_results.tsv",
  results = "data/derived/impact_layerE_results.tsv",
  influence = "data/derived/impact_layerE_influence.tsv",
  gate = "data/derived/impact_layerE_gate.tsv",
  determinism = "data/metadata/impact_layerE_determinism_receipt.tsv",
  run_receipt = "data/metadata/impact_layerE_run_receipt.tsv",
  figure_png = "figures/impact_layerE_multiverse.png",
  figure_svg = "figures/impact_layerE_multiverse.svg",
  legacy_permanova = "data/derived/layerE_permanova.tsv"
)
check(all(file.exists(paths)), paste("missing output(s):",
                                    paste(paths[!file.exists(paths)], collapse = ", ")))
check(requireNamespace("jsonlite", quietly = TRUE), "jsonlite is unavailable")
manifest <- jsonlite::fromJSON(paths[["manifest"]], simplifyVector = TRUE)
cfg <- manifest$layerE_multiverse

## Portable input and observation-mask invariants.
lib <- load_impact_library(paths[["input"]])
label_cols <- grep("^analysis_label__", names(lib), value = TRUE)
genus_cols <- grep("^genus__", names(lib), value = TRUE)
check(nrow(lib) == 647L && length(label_cols) == 14L && length(genus_cols) == 12L,
      "portable dimensions are not 647 x (14 labels + 12 genera)")
check(sum(lib$observed_zero) == 18L &&
        !any(lib$structural_absence) && !any(lib$analytical_exclusion),
      "mask fields conflate observed zeros with structural/excluded cells")
check(all(rowSums(lib[, label_cols, drop = FALSE]) == lib$fish_read_depth) &&
        all(rowSums(lib[, genus_cols, drop = FALSE]) == lib$fish_read_depth),
      "portable label/genus read depths differ")
design <- design_from_library(lib)
check(nrow(design$unit_md) == 234L && length(design$common_timepoints) == 11L,
      "portable design is not 234 realized units / 11 common timepoints")

## Exact finite specification grid.
specifications <- read_tsv(paths[["specifications"]])
check(nrow(specifications) == 72L &&
        identical(specifications$spec_id, sprintf("E%03d", 1:72)),
      "specification table is not the complete stable 72-family grid")
check(setequal(specifications$panel, cfg$panels) &&
        setequal(specifications$taxonomy, cfg$taxonomy) &&
        setequal(specifications$library_detection_threshold_reads,
                 cfg$library_detection_threshold_reads) &&
        setequal(specifications$effort, cfg$effort) &&
        setequal(specifications$distance, cfg$distances),
      "a manifest factor level is absent from the specification grid")
check(!anyDuplicated(specifications[, c("panel", "taxonomy",
                                        "library_detection_threshold_reads",
                                        "effort", "distance")]),
      "specification table contains a duplicate family")
check(all(specifications$spatial_nonidentity_permutations == 27L) &&
        all(specifications$exact_p_denominator_with_observed == 28L),
      "specifications do not preserve exact 27+observed resolution")

## Complete dihedral permutations and whole-site blocks.
perms <- series_permutation_indices(14L)
check(identical(dim(perms), c(27L, 14L)) &&
        !anyDuplicated(apply(perms, 1L, paste, collapse = ",")),
      "restricted river permutation set is not 27 unique arrangements")
site_time_panel <- make_panel(
  aggregate_selected_units(lib, design, taxonomy_matrix(lib, "analysis_label"),
                           lib$library_id, 1L), design, "balanced_site_time")
site_time_groups <- permutation_groups(site_time_panel$row_sites, design, perms)
block_ok <- all(vapply(seq_len(nrow(site_time_groups)), function(i) {
  all(vapply(design$site_order, function(site) {
    length(unique(site_time_groups[i, site_time_panel$row_sites == site])) == 1L
  }, logical(1)))
}, logical(1)))
check(block_ok, "site x timepoint rows do not move as whole-site permutation blocks")

## Draw selections and full retained draw results.
library_draws <- read_tsv(paths[["library_draws"]])
check(nrow(library_draws) == 23400L &&
        all(table(library_draws$draw_id) == 234L) &&
        !anyDuplicated(library_draws[, c("draw_id", "unit_id")]),
      "100 one-library selection sets are incomplete")
check(all(library_draws$seed == as.integer(manifest$seeds$layerE_multiverse) +
            library_draws$draw_id),
      "library draw seed derivation differs from base seed + draw ID")
selected_join <- lib[match(library_draws$selected_library_id, lib$library_id),
                     c("library_id", "unit_id"), drop = FALSE]
check(!anyNA(selected_join$library_id) &&
        identical(selected_join$unit_id, library_draws$unit_id),
      "a selected library is absent from, or outside, its biological unit")

draw_results <- read_tsv(paths[["draw_results"]])
one_specs <- specifications$spec_id[
  specifications$effort == "one_library_100_draw_ensemble"]
check(nrow(draw_results) == 3600L &&
        setequal(draw_results$spec_id, one_specs) &&
        all(table(draw_results$spec_id) == 100L),
      "draw-result table does not retain all 36 x 100 declared rows")
check(all(draw_results$selected_library_count == 234L) &&
        all(draw_results$n_spatial_nonidentity_permutations == 27L),
      "draw rows lost the equalised library count or exact permutation count")
finite_draw_p <- draw_results$exact_p[is.finite(draw_results$exact_p)]
check(all(abs(finite_draw_p * 28 - round(finite_draw_p * 28)) < 1e-8),
      "a draw exact p-value is not on the 1/28 support")
check(all(nzchar(draw_results$failure_reason[!draw_results$valid])) &&
        sum(draw_results$valid) + sum(!draw_results$valid) == 3600L,
      "failed draw rows were dropped or lack a reason")

## Family results: complete, conservative validity, and documented outcome.
results <- read_tsv(paths[["results"]])
check(nrow(results) == 72L && identical(results$spec_id, specifications$spec_id),
      "family results are not one-for-one with specifications")
check(all(results$retained_draws == results$expected_draws) &&
        all(results$valid_draws + results$invalid_draws == results$retained_draws),
      "family draw denominators are incomplete")
check(sum(results$family_valid) == 64L && sum(!results$family_valid) == 8L,
      "documented 64/72 valid-family result changed")
check(all(results$failure_reason[!results$family_valid] ==
            "empty_community_undefined_for_rclr"),
      "retained family failures are not the documented empty-community rCLR limit")
check(all(results$richness_contrast_lower_minus_upper > 0) &&
        all(draw_results$richness_contrast_lower_minus_upper > 0),
      "a declared family or one-library draw reversed the richness direction")
finite_family_p <- results$exact_p[is.finite(results$exact_p)]
check(all(abs(finite_family_p * 28 - round(finite_family_p * 28)) < 1e-8),
      "a family exact p-value is not on the 1/28 support")

## Regression against the established balanced full-effort Layer-E reference.
legacy <- read_tsv(paths[["legacy_permanova"]])
legacy <- legacy[legacy$grouping == "section", , drop = FALSE]
distance_map <- c(
  bray_curtis_proportions = "Bray-Curtis (rel. abund.)",
  jaccard_occurrence = "Jaccard (presence/absence)",
  robust_aitchison_rclr = "robust Aitchison (rclr)"
)
reference <- results[
  results$panel == "balanced_site_profile" &
    results$taxonomy == "analysis_label" &
    results$library_detection_threshold_reads == 1L &
    results$effort == "all_surviving", , drop = FALSE]
for (i in seq_len(nrow(reference))) {
  old <- legacy[legacy$distance == distance_map[[reference$distance[[i]]]], , drop = FALSE]
  check(nrow(old) == 1L &&
          abs(reference$pseudo_R2[[i]] - old$R2) < 5e-5 &&
          abs(reference$pseudo_F[[i]] - old$F) < 5e-4 &&
          abs(reference$exact_p[[i]] - old$p_value) < 1e-10,
        paste("balanced reference regression failed for", reference$distance[[i]]))
}
check(all(abs(reference$exact_p - 2 / 28) < 1e-10),
      "established balanced-site p=0.0714 boundary changed")

## Exactly 75 valid reference leave-outs and no sign reversal.
influence <- read_tsv(paths[["influence"]])
check(nrow(influence) == 75L &&
        all(table(influence$distance) == 25L) &&
        all(table(influence$omission_type) == c(site = 42L, timepoint = 33L)),
      "influence table is not 3 distances x (14 sites + 11 timepoints)")
check(length(unique(influence$omitted[influence$omission_type == "site"])) == 14L &&
        length(unique(influence$omitted[influence$omission_type == "timepoint"])) == 11L,
      "influence omissions do not cover every reference site/timepoint")
check(all(influence$valid) && !any(influence$sign_reversal) &&
        all(influence$n_spatial_nonidentity_permutations == 27L),
      "a reference influence row failed, reversed direction, or changed resolution")

## Gate arithmetic, claim boundary, and deterministic byte receipts.
gate <- read_tsv(paths[["gate"]])
gate_value <- function(id, field) gate[gate$criterion == id, field][[1]]
check(abs(gate_value("valid_specification_families", "observed_value") - 64 / 72) < 1e-10 &&
        !gate_value("valid_specification_families", "pass"),
      "valid-family gate arithmetic/status changed")
check(gate_value("positive_valid_family_contrasts", "pass") &&
        gate_value("positive_one_library_draw_contrasts", "pass") &&
        gate_value("reference_leaveout_sign_reversals", "pass") &&
        !gate_value("overall_gate", "pass"),
      "documented failed overall gate/components changed")
check(all(grepl("no abundance or confirmed centroid-shift claim",
                results$claim_scope, fixed = TRUE)),
      "result table lost the preregistered claim boundary")

determinism <- read_tsv(paths[["determinism"]])
check(nrow(determinism) == 8L && all(determinism$status == "PASS") &&
        all(determinism$byte_identical_to_preexisting) &&
        all(determinism$preexisting_sha256 == determinism$rebuilt_sha256),
      "second-build determinism receipt is not eight byte-identical PASS rows")
check(all(vapply(seq_len(nrow(determinism)), function(i) {
  identical(sha256(determinism$artifact[[i]]), determinism$rebuilt_sha256[[i]])
}, logical(1))), "current canonical hashes differ from determinism receipt")
check(file.info(paths[["figure_png"]])$size > 100000L &&
        file.info(paths[["figure_svg"]])$size > 100000L,
      "multiverse figures are unexpectedly small")

cat(sprintf("Layer-E multiverse tests PASS (%d unit/invariant/regression checks).\n",
            test_count))
