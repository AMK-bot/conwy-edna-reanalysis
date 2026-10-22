#!/usr/bin/env Rscript

# Independent invariants/regression checks for R/14_trait_null_model.R outputs.

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
cfg <- manifest$trait_null
summary <- read_tsv("data/derived/impact_trait_null_summary.tsv")
draws <- read_tsv("data/derived/impact_trait_null_draws.tsv")
gower <- read_tsv("data/derived/impact_trait_gower.tsv")
units <- read_tsv("data/derived/impact_trait_unit_rao.tsv")
receipt <- read_tsv("data/metadata/impact_trait_null_receipt.tsv")
groups <- read_tsv("data/metadata/fish_trait_analysis_groups.tsv", na.strings = character())

assert(identical(summary$analysis, c("primary_12_labels", "strict_5_labels")),
       "Trait-null summary analysis order changed")
assert(nrow(draws) == 2L * as.integer(cfg$null_draws), "Trait-null draw denominator changed")
assert(nrow(units) == 2L * as.integer(cfg$balanced_units), "Trait-null unit denominator changed")
assert(all(table(draws$analysis) == as.integer(cfg$null_draws)), "Unequal null-draw denominators")
assert(all(draws$margin_pass), "Persisted null draw failed its fixed-margin gate")
assert(all(draws$swap_target_pass), "A movable timepoint block missed its declared swap target")
assert(all(draws$accepted_swaps == as.integer(cfg$accepted_swaps_between_draws_per_timepoint) *
             (as.integer(cfg$balanced_timepoints) - draws$structural_fixed_timepoints)),
       "Accepted swap count does not match the movable-block denominator")

primary_labels <- groups$analysis_label[groups$include_in_descriptive_summary == "yes"]
strict_labels <- as.character(cfg$strict_labels)
assert(length(primary_labels) == 12L && length(strict_labels) == 5L, "Trait panel sizes changed")
assert(nrow(gower[gower$analysis == "primary_12_labels", ]) == 12L^2,
       "Primary Gower matrix is incomplete")
assert(nrow(gower[gower$analysis == "strict_5_labels", ]) == 5L^2,
       "Strict Gower matrix is incomplete")

for (analysis in unique(gower$analysis)) {
  z <- gower[gower$analysis == analysis, ]
  labels <- unique(z$label_a)
  m <- matrix(NA_real_, length(labels), length(labels), dimnames = list(labels, labels))
  for (i in seq_len(nrow(z))) m[z$label_a[[i]], z$label_b[[i]]] <- z$distance[[i]]
  assert(!anyNA(m), paste("Incomplete Gower matrix for", analysis))
  assert(max(abs(m - t(m))) < 1e-12, paste("Asymmetric Gower matrix for", analysis))
  assert(max(abs(diag(m))) < 1e-12, paste("Nonzero Gower diagonal for", analysis))
  assert(all(m >= 0 & m <= 1), paste("Gower distance outside [0,1] for", analysis))
}

for (i in seq_len(nrow(summary))) {
  z <- draws[draws$analysis == summary$analysis[[i]], "contrast_lower_minus_upper"]
  expected_p <- (1 + sum(z >= summary$observed_contrast[[i]])) / (length(z) + 1)
  assert(abs(expected_p - summary$monte_carlo_p[[i]]) < 1e-12,
         paste("Monte Carlo p-value does not reproduce for", summary$analysis[[i]]))
  expected_gate <- summary$observed_contrast[[i]] > 0 &&
    expected_p <= as.numeric(cfg$one_sided_alpha) &&
    summary$unique_null_matrices[[i]] >= as.integer(cfg$min_unique_null_matrices) &&
    summary$null_sd[[i]] > 0 && summary$defined_sections[[i]] &&
    summary$every_movable_block_met_swap_target[[i]]
  assert(identical(expected_gate, summary$gate_pass[[i]]),
         paste("Gate does not reproduce for", summary$analysis[[i]]))
}

assert(all(units$richness >= 0), "Negative detected-label richness")
assert(all(units$rao_q[units$richness == 1] == 0), "Singleton Rao Q must be zero")
assert(all(is.na(units$rao_q[units$richness == 0])), "Zero-richness Rao Q must be missing")
assert(all(units$rao_q[is.finite(units$rao_q)] >= 0 &
             units$rao_q[is.finite(units$rao_q)] <= 1), "Rao Q outside [0,1]")
assert(receipt$seed[[1]] == as.integer(manifest$seeds$trait_null), "Trait seed drift")
assert(!grepl("abundance|biomass|diet|realised habitat", receipt$allowed_claim[[1]], ignore.case = TRUE),
       "Trait receipt crosses a prohibited claim boundary")
assert(file.info("figures/impact_trait_null_model.png")$size > 50000,
       "Trait-null PNG is unexpectedly small")
assert(file.info("figures/impact_trait_null_model.svg")$size > 5000,
       "Trait-null SVG is unexpectedly small")

cat(jsonlite::toJSON(list(
  status = "PASS", checks = list(summary_rows = nrow(summary), null_draw_rows = nrow(draws),
                                  unit_rows = nrow(units), gower_rows = nrow(gower))),
  auto_unbox = TRUE, pretty = TRUE), "\n")
