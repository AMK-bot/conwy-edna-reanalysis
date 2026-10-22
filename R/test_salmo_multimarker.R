#!/usr/bin/env Rscript

# Unit, invariant and deterministic-regression tests for the offline Salmo
# multi-marker benchmark.

options(stringsAsFactors = FALSE, digits = 15, scipen = 999)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(script_arg) == 1L)
project_root <- normalizePath(
  file.path(dirname(sub("^--file=", "", script_arg[1])), ".."),
  mustWork = TRUE
)
Sys.setenv(SALMO_BENCHMARK_SOURCE_ONLY = "1")
source(file.path(project_root, "R", "12_salmo_multimarker_benchmark.R"))

stopifnot(requireNamespace("digest", quietly = TRUE))

result <- calculate_benchmark(project_root)

# Frozen input/provenance regression: accession changes require a deliberate
# cache refresh and review, never a silent benchmark drift.
expected_accessions <- c(
  "AF133701.1", "AF154851.1",
  "DQ288268.1", "DQ288269.1", "DQ288270.1", "DQ288271.1",
  "JQ390055.1", "JQ390056.1", "JQ390057.1",
  "KF792729.1", "KP013084.1", "LC012541.1",
  "LT617630.1", "LT617631.1", "LT617632.1", "LT617633.1",
  "MF621740.1", "MF621741.1", "MF621743.1", "MF621760.1",
  "MN530962.1"
)
stopifnot(identical(sort(result$inputs$selected$accession_version), expected_accessions))
expected_query_counts <- c(30L, 21L, 0L, 0L, 0L, 0L, 1L, 111L, 41L)
stopifnot(identical(as.integer(result$inputs$queries$ncbi_total_count), expected_query_counts))
stopifnot(nrow(result$inputs$selected) == 21L)
stopifnot(nrow(result$inputs$marker_manifest) == 63L)
stopifnot(all(table(result$inputs$marker_manifest$marker) == 21L))
stopifnot(all(result$inputs$marker_manifest$sequence_length[
  result$inputs$marker_manifest$marker == "mifish_12S"
] == 170L))
for (label in result$inputs$queries$submitted_label) {
  candidates <- result$inputs$accessions[
    result$inputs$accessions$submitted_label == label, , drop = FALSE
  ]
  if (nrow(candidates) > 1L) {
    stopifnot(identical(candidates$accession_version, sort(candidates$accession_version)))
  }
}

# Complete table dimensions and core metric invariants.
stopifnot(nrow(result$coverage) == 27L)
stopifnot(nrow(result$pairwise) == 3L * choose(21L, 2L))
stopifnot(nrow(result$barcode_gaps) == 6L)
stopifnot(nrow(result$loo) == 60L)
stopifnot(nrow(result$gates) == 3L)
pair_keys <- paste(
  result$pairwise$marker,
  result$pairwise$sequence_id_1,
  result$pairwise$sequence_id_2,
  sep = "\r"
)
stopifnot(!anyDuplicated(pair_keys))
stopifnot(all(result$pairwise$edit_distance >= 0L))
stopifnot(all(result$pairwise$normalization_length %in% c(170L, 648L, 1140L)))
stopifnot(isTRUE(all.equal(
  result$pairwise$normalized_edit_distance,
  result$pairwise$edit_distance / result$pairwise$normalization_length,
  tolerance = 1e-14
)))
stopifnot(all((result$pairwise$edit_distance == 0L) ==
                (result$pairwise$exact_haplotype == "yes")))
stopifnot(all(result$barcode_gaps$barcode_gap ==
                result$barcode_gaps$nearest_between_label_distance -
                result$barcode_gaps$within_max_distance))
stopifnot(all(result$barcode_gaps$positive_barcode_gap == "yes"))

# Coverage, tie, gap and genus-boundary gates.  The narrowed two-label
# diagnostics cannot be promoted through the frozen four-label species gate.
coverage_counts <- result$coverage[
  result$coverage$marker == "mifish_12S",
  c("submitted_label", "passing_accessions"),
  drop = FALSE
]
expected_passing <- c(5L, 5L, 0L, 0L, 0L, 0L, 1L, 5L, 5L)
stopifnot(identical(as.integer(coverage_counts$passing_accessions), expected_passing))
stopifnot(all(result$gates$candidate_salmo_labels == 7L))
stopifnot(all(result$gates$eligible_salmo_labels == 2L))
stopifnot(all(result$gates$required_salmo_labels_for_species == 4L))
stopifnot(all(grepl(
  "2 of 7 candidate Salmo labels eligible; minimum 4 required",
  result$gates$species_gate_reason,
  fixed = TRUE
)))
stopifnot(all(result$gates$species_panel_coverage_adequate == "no"))
stopifnot(all(result$gates$species_gate_status == "coverage-inadequate"))
stopifnot(all(result$gates$species_loo_queries == 10L))
stopifnot(all(result$gates$species_loo_all_ties_accuracy == 1))
stopifnot(all(result$gates$all_loo_ties_same_submitted_species == "yes"))
stopifnot(all(result$gates$cross_label_exact_haplotypes == 0L))
stopifnot(all(result$gates$all_eligible_labels_positive_barcode_gap == "yes"))
stopifnot(all(result$gates$eligible_outgroup_labels == 2L))
stopifnot(all(result$gates$genus_panel_coverage_adequate == "yes"))
stopifnot(all(result$gates$genus_loo_queries == 10L))
stopifnot(all(result$gates$salmo_queries_nearer_than_each_outgroup == 10L))
stopifnot(all(result$gates$all_salmo_queries_nearer_than_each_outgroup == "yes"))
stopifnot(all(result$gates$genus_gate_status == "pass"))
stopifnot(all(result$gates$finest_consistently_supported_rank == "genus"))
stopifnot(all(grepl("ASV0001 species identity", result$gates$claim_boundary, fixed = TRUE)))
stopifnot(!any(result$loo$query_submitted_label == "Salmo obtusirostris"))

# All-tie rule is retained explicitly rather than collapsed to a first hit.
for (i in seq_len(nrow(result$loo))) {
  tied_labels <- strsplit(result$loo$tied_submitted_labels[i], ";", fixed = TRUE)[[1]]
  tied_genera <- strsplit(result$loo$tied_genera[i], ";", fixed = TRUE)[[1]]
  stopifnot((length(tied_labels) == 1L) ==
              (result$loo$all_ties_same_species[i] == "yes"))
  stopifnot((length(tied_genera) == 1L) ==
              (result$loo$all_ties_same_genus[i] == "yes"))
}

# The analysis script is offline: official URLs are provenance strings, never
# opened by R.
r_code <- readLines(file.path(project_root, "R", "12_salmo_multimarker_benchmark.R"), warn = FALSE)
forbidden_network_calls <- c("download.file(", "url(", "curl::", "httr::", "rentrez::")
stopifnot(!any(vapply(forbidden_network_calls, function(token) {
  any(grepl(token, r_code, fixed = TRUE))
}, logical(1))))

# Deterministic double-build regression, including both rendered formats.
build_one <- file.path(tempdir(), "salmo_multimarker_build_one")
build_two <- file.path(tempdir(), "salmo_multimarker_build_two")
dir.create(build_one, recursive = TRUE, showWarnings = FALSE)
dir.create(build_two, recursive = TRUE, showWarnings = FALSE)
write_benchmark(result, build_one)
write_benchmark(result, build_two)
relative_outputs <- c(
  file.path("data", "metadata", "impact_salmo_panel_coverage.tsv"),
  file.path("data", "derived", "impact_salmo_pairwise_distances.tsv"),
  file.path("data", "derived", "impact_salmo_exact_haplotypes.tsv"),
  file.path("data", "derived", "impact_salmo_barcode_gaps.tsv"),
  file.path("data", "derived", "impact_salmo_loo_classification.tsv"),
  file.path("data", "derived", "impact_salmo_marker_gates.tsv"),
  file.path("figures", "impact_salmo_multimarker.png"),
  file.path("figures", "impact_salmo_multimarker.svg")
)
hash_file <- function(path) digest::digest(file = path, algo = "sha256")
hash_one <- vapply(file.path(build_one, relative_outputs), hash_file, character(1))
hash_two <- vapply(file.path(build_two, relative_outputs), hash_file, character(1))
hash_release <- vapply(file.path(project_root, relative_outputs), hash_file, character(1))
stopifnot(identical(unname(hash_one), unname(hash_two)))
stopifnot(identical(unname(hash_one), unname(hash_release)))

cat("PASS: NCBI receipt, accession, marker and sequence-hash invariants\n")
cat("PASS: distance, exact-haplotype, barcode-gap and all-ties LOO invariants\n")
cat("PASS: species coverage withheld and genus boundary gate enforced\n")
cat("PASS: offline double build is byte-identical for 6 TSVs + PNG + SVG\n")
