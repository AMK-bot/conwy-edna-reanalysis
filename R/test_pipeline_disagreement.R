#!/usr/bin/env Rscript

# Independent denominator, bootstrap and claim-boundary checks for R/15.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(root)
options(stringsAsFactors = FALSE, OutDec = ".")
assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
read_tsv <- function(path, ...) read.delim(path, sep = "\t", quote = "", comment.char = "",
                                           check.names = FALSE, ...)

manifest <- jsonlite::fromJSON("data/metadata/impact_sprint_manifest.json", simplifyVector = TRUE)
cfg <- manifest$pipeline_disagreement
cells <- read_tsv("data/metadata/impact_pipeline_disagreement_cells.tsv")
model <- read_tsv("data/derived/impact_pipeline_model.tsv")
boot <- read_tsv("data/derived/impact_pipeline_bootstrap.tsv")
loo <- read_tsv("data/derived/impact_pipeline_leaveout.tsv")
rates <- read_tsv("data/derived/impact_pipeline_rates.tsv")
summary <- read_tsv("data/derived/impact_pipeline_summary.tsv")
receipt <- read_tsv("data/metadata/impact_pipeline_receipt.tsv")

assert(nrow(cells) == as.integer(cfg$source_positive_cells), "Source-positive cell denominator drift")
assert(length(unique(cells$ID)) == as.integer(cfg$shared_libraries), "Shared-library denominator drift")
assert(length(unique(cells$taxon)) == as.integer(cfg$common_source_positive_labels),
       "Source-positive label denominator drift")
assert(sum(cells$miss) == as.integer(cfg$nonrecovered_cells), "Non-recovered cell denominator drift")
assert(sum(cells$miss == 0) == as.integer(cfg$recovered_cells), "Recovered cell denominator drift")
assert(all(cells$source_relative_abundance > 0), "Source-negative cell entered conditional analysis")
assert(all(cells$reprocessed_raw_name_depth > 0), "Invalid raw-name depth")
required_depth_fields <- c(
  "reprocessed_assigned_fish_read_total",
  "reprocessed_assigned_fish_reads_excluding_focal_label"
)
assert(all(required_depth_fields %in% names(cells)),
       "Portable input lacks the focal-excluded assigned-fish depth sensitivity fields")
assert(all(cells$reprocessed_assigned_fish_read_total >= 0) &&
         all(cells$reprocessed_assigned_fish_reads_excluding_focal_label >= 0),
       "Assigned-fish read totals must be non-negative")
assert(all(cells$reprocessed_assigned_fish_reads_excluding_focal_label +
             cells$reprocessed_count == cells$reprocessed_assigned_fish_read_total),
       "Focal-excluded total does not exactly remove the focal harmonised-label count")
raw_counts <- readRDS("data/derived/species_x_sample_counts_rawnames.rds")
assert(ncol(raw_counts) == 35L && "Unassigned" %in% colnames(raw_counts),
       "Raw-name count universe changed")
assigned_columns <- setdiff(colnames(raw_counts), "Unassigned")
expected_assigned_total <- rowSums(raw_counts[cells$ID, assigned_columns, drop = FALSE])
assert(all(cells$reprocessed_assigned_fish_read_total == expected_assigned_total),
       "Assigned-fish total is not the exact 34-column sum excluding Unassigned")
assert(all(cells$reprocessed_raw_name_depth == rowSums(raw_counts[cells$ID, , drop = FALSE])),
       "Preregistered primary total is not the exact 35-column sum")
assert(sum(cells$reprocessed_assigned_fish_reads_excluding_focal_label == 0) == 224L,
       "Focal-excluded zero-cell denominator drift")
assert(all(cells$miss == as.integer(cells$reprocessed_count == 0)), "Miss outcome/count mismatch")
assert(!any(grepl("false-negative|ground truth", cells$evidence_scope, ignore.case = TRUE)),
       "Portable input crosses the comparator-relative claim boundary")

assert(length(unique(boot$draw)) == as.integer(cfg$bootstrap_draws), "Bootstrap draw count drift")
assert(all(table(boot$term) == as.integer(cfg$bootstrap_draws)), "Incomplete coefficient-by-draw grid")
required_boot_sensitivity <- c(
  "focal_excluded_model_term", "focal_excluded_estimate",
  "focal_excluded_fit_ok", "focal_excluded_failure_reason"
)
assert(all(required_boot_sensitivity %in% names(boot)),
       "Bootstrap table lacks the paired focal-excluded sensitivity")
assert(all(table(boot$focal_excluded_model_term) == as.integer(cfg$bootstrap_draws)),
       "Incomplete focal-excluded coefficient-by-draw grid")
assert(all(boot$focal_excluded_fit_ok), "A focal-excluded bootstrap fit failed")
assert(length(unique(loo$group_id[loo$group_type == "taxon"])) == as.integer(cfg$leaveout_taxa),
       "Taxon leave-out denominator drift")
assert(length(unique(loo$group_id[loo$group_type == "site"])) == as.integer(cfg$leaveout_sites),
       "Site leave-out denominator drift")
assert(length(unique(loo$group_id[loo$group_type == "timepoint"])) == as.integer(cfg$leaveout_timepoints),
       "Timepoint leave-out denominator drift")
required_loo_sensitivity <- c(
  "focal_excluded_model_term", "focal_excluded_estimate",
  "focal_excluded_fit_ok", "focal_excluded_failure_reason"
)
assert(all(required_loo_sensitivity %in% names(loo)),
       "Leave-out table lacks the paired focal-excluded sensitivity")
assert(all(loo$focal_excluded_fit_ok), "A focal-excluded leave-out fit failed")
assert(nrow(summary) == 1L && summary$source_positive_cells[[1]] == nrow(cells),
       "Summary denominator mismatch")
assert(summary$nonrecovered_cells[[1]] == sum(cells$miss), "Summary miss mismatch")
assert(receipt$seed[[1]] == as.integer(manifest$seeds$pipeline_bootstrap), "Bootstrap seed drift")
assert(receipt$bootstrap_draws[[1]] == as.integer(cfg$bootstrap_draws), "Receipt draw count drift")

primary_model_name <- if (summary$collinearity_gate_pass[[1]]) {
  "full_primary"
} else {
  "abundance_depth_primary"
}
effect_terms <- intersect(
  c("z_log_source_ra", "z_log_depth", "z_label_prevalence"),
  model$term[model$model == primary_model_name & !is.na(model$bootstrap_finite_fraction)])
for (term in effect_terms) {
  row <- model[model$model == primary_model_name & model$term == term &
                 !is.na(model$bootstrap_finite_fraction), , drop = FALSE]
  assert(nrow(row) == 1L, paste("Missing unique gated effect row for", term))
  expected_gate <- row$bootstrap_finite_fraction >= as.numeric(cfg$min_finite_bootstrap_fraction) &&
    (row$bootstrap_q025 > 0 || row$bootstrap_q975 < 0) &&
    row$leaveout_sign_stability >= as.numeric(cfg$min_leaveout_sign_stability) &&
    row$no_leaveout_reversal
  assert(identical(expected_gate, row$gate_pass), paste("Gate does not reproduce for", term))
}

sensitivity_term <- "z_log_focal_excluded_assigned_fish_reads"
sensitivity_model <- "focal_excluded_assigned_fish_sensitivity"
sensitivity_row <- model[
  model$model == sensitivity_model & model$term == sensitivity_term &
    !is.na(model$bootstrap_finite_fraction), , drop = FALSE]
assert(nrow(sensitivity_row) == 1L,
       "Missing unique machine-readable focal-excluded depth sensitivity row")
expected_sensitivity_gate <-
  sensitivity_row$bootstrap_finite_fraction >= as.numeric(cfg$min_finite_bootstrap_fraction) &&
  (sensitivity_row$bootstrap_q025 > 0 || sensitivity_row$bootstrap_q975 < 0) &&
  sensitivity_row$leaveout_sign_stability >= as.numeric(cfg$min_leaveout_sign_stability) &&
  sensitivity_row$no_leaveout_reversal
assert(identical(expected_sensitivity_gate, sensitivity_row$gate_pass),
       "Focal-excluded depth gate does not reproduce")
assert(summary$sensitivity_bootstrap_draws[[1]] == as.integer(cfg$bootstrap_draws) &&
         summary$sensitivity_leaveout_groups[[1]] == 43L,
       "Summary does not preserve the 999-draw/43-leave-out sensitivity design")
assert(identical(summary$depth_association_survives_focal_exclusion[[1]],
                 sensitivity_row$gate_pass[[1]]),
       "Summary depth-survival status disagrees with the focal-excluded gate")
expected_depth_status <- if (sensitivity_row$gate_pass[[1]]) {
  "SURVIVES_FOCAL_EXCLUSION"
} else {
  "DOES_NOT_SURVIVE_FOCAL_EXCLUSION"
}
assert(summary$depth_association_sensitivity_status[[1]] == expected_depth_status &&
         receipt$depth_association_sensitivity_status[[1]] == expected_depth_status,
       "Machine-readable depth-sensitivity status does not reproduce the gate")
assert(receipt$sensitivity_bootstrap_draws[[1]] == as.integer(cfg$bootstrap_draws) &&
         receipt$sensitivity_leaveout_groups[[1]] == 43L,
       "Receipt does not record the paired focal-excluded sensitivity denominators")

for (kind in c("section", "season", "taxon")) {
  z <- rates[rates$group_type == kind, ]
  assert(sum(z$source_positive_cells) == nrow(cells), paste("Rate denominator mismatch for", kind))
  assert(sum(z$nonrecovered_cells) == sum(cells$miss), paste("Rate miss mismatch for", kind))
}
assert(!grepl("false-negative|ground truth", summary$evidence_scope[[1]], ignore.case = TRUE),
       "Summary crosses the comparator-relative claim boundary")
assert(file.info("figures/impact_pipeline_disagreement.png")$size > 50000,
       "Pipeline-disagreement PNG is unexpectedly small")
assert(file.info("figures/impact_pipeline_disagreement.svg")$size > 5000,
       "Pipeline-disagreement SVG is unexpectedly small")

cat(jsonlite::toJSON(list(
  status = "PASS", checks = list(cells = nrow(cells), nonrecovered = sum(cells$miss),
                                  bootstrap_draws = length(unique(boot$draw)),
                                  leaveout_groups = length(unique(paste(loo$group_type, loo$group_id))))),
  auto_unbox = TRUE, pretty = TRUE), "\n")
