#!/usr/bin/env Rscript

# Export the Phase-5 impact analyses' small, portable Layer-E input.
#
# This is the only Phase-5 script that reads the two RDS count matrices.  Every
# downstream impact analysis reads the committed TSV instead, so a clean release
# does not depend on R's serialization format.  Structural absence, analytical
# exclusion, and a genuine observed zero are deliberately separate fields.

full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)
options(stringsAsFactors = FALSE, OutDec = ".", digits = 15, scipen = 999)

assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, ...)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
}
sha256 <- function(path) {
  assert(requireNamespace("digest", quietly = TRUE),
         "The digest package is required to make SHA-256 receipts")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

paths <- list(
  species = "data/derived/species_x_sample_counts_clean.rds",
  genus = "data/derived/genus_x_sample_counts_clean.rds",
  covariates = "data/metadata/conwy_covariates.tsv",
  occasions = "data/metadata/occasions.tsv",
  output = "data/derived/impact_library_counts.tsv",
  receipt = "data/metadata/impact_library_export_receipt.tsv"
)
invisible(lapply(paths[c("species", "genus", "covariates", "occasions")], function(path) {
  assert(file.exists(path), paste("Missing required input:", path))
}))

species <- readRDS(paths$species)
genus <- readRDS(paths$genus)
assert(is.matrix(species) && identical(dim(species), c(647L, 14L)),
       "Species RDS must be exactly 647 x 14")
assert(is.matrix(genus) && identical(dim(genus), c(647L, 12L)),
       "Genus RDS must be exactly 647 x 12")
assert(identical(rownames(species), rownames(genus)),
       "Species and genus RDS row names/order differ")
assert(!is.null(rownames(species)) && !anyDuplicated(rownames(species)),
       "Library IDs must be present and unique")
assert(!is.null(colnames(species)) && !anyDuplicated(colnames(species)),
       "Analysis-label names must be present and unique")
assert(!is.null(colnames(genus)) && !anyDuplicated(colnames(genus)),
       "Genus names must be present and unique")
assert(all(is.finite(species)) && all(species >= 0) && all(species == round(species)),
       "Species counts must be finite, non-negative whole reads")
assert(all(is.finite(genus)) && all(genus >= 0) && all(genus == round(genus)),
       "Genus counts must be finite, non-negative whole reads")
assert(all(species <= .Machine$integer.max) && all(genus <= .Machine$integer.max),
       "A count exceeds the portable 32-bit integer range")

species <- round(species)
genus <- round(genus)
storage.mode(species) <- "integer"
storage.mode(genus) <- "integer"
assert(identical(as.integer(rowSums(species)), as.integer(rowSums(genus))),
       "Species and genus fish-read depths are not identical")

# The 12-genus table must be an exact deterministic collapse of the 14 labels.
derived_genus <- sub("_.*$", "", colnames(species))
assert(setequal(unique(derived_genus), colnames(genus)),
       "Analysis-label prefixes do not resolve to the 12 genus columns")
collapsed <- vapply(colnames(genus), function(taxon) {
  rowSums(species[, derived_genus == taxon, drop = FALSE])
}, numeric(nrow(species)))
assert(identical(dim(genus), dim(collapsed)) && all(genus == collapsed),
       "The genus RDS is not an exact collapse of the species RDS")

ids <- rownames(species)
assert(all(grepl("^E[0-9]{2}_R[1-3]_T[0-9]{2}$", ids)),
       "Unexpected library ID; expected E##_R#_T##")
parts <- do.call(rbind, strsplit(ids, "_", fixed = TRUE))
site <- parts[, 1]
replicate <- parts[, 2]
timepoint <- parts[, 3]
unit_id <- paste(site, timepoint, sep = "_")

covariates <- read_tsv(paths$covariates)
occasions <- read_tsv(paths$occasions)
assert("ID" %in% names(covariates) && !anyDuplicated(covariates$ID),
       "Covariate IDs must be unique")
assert("ID" %in% names(occasions) && !anyDuplicated(occasions$ID),
       "Occasion IDs must be unique")
assert(all(ids %in% covariates$ID) && all(ids %in% occasions$ID),
       "Every surviving library must resolve in both metadata tables")
covariates <- covariates[match(ids, covariates$ID), , drop = FALSE]
occasions <- occasions[match(ids, occasions$ID), , drop = FALSE]
assert(identical(covariates$ID, ids) && identical(occasions$ID, ids),
       "Metadata join changed library order")
assert(identical(occasions$site, site) &&
       identical(occasions$timepoint, timepoint) &&
       identical(occasions$rep, replicate),
       "Occasion metadata disagrees with parsed library IDs")
assert(all(occasions$surveyed == 1L) && all(occasions$usable_timepoint == 1L),
       "A reprocessed survivor falls outside the surveyed, usable observation mask")

section <- ifelse(site %in% sprintf("E%02d", 1:5), "upper",
                  ifelse(site %in% sprintf("E%02d", 6:10), "middle", "lower"))
assert(identical(section, covariates$river_section),
       "Derived river sections disagree with conwy_covariates.tsv")
assert(all(covariates$season %in% c("Spring", "Summer", "Autumn", "Winter")),
       "Unexpected or missing season")

unit_n <- table(unit_id)
assert(length(unit_n) == 234L && all(unit_n %in% 1:3),
       "Expected 234 realized site x timepoint units with 1-3 surviving libraries")
unit_frame <- unique(data.frame(site = site, timepoint = timepoint,
                                stringsAsFactors = FALSE))
grid <- xtabs(~ site + timepoint, data = unit_frame) > 0
common_timepoints <- colnames(grid)[colSums(grid) == nrow(grid)]
assert(nrow(grid) == 14L && length(common_timepoints) == 11L,
       "Expected 14 sites and 11 timepoints represented at every site")

fish_depth <- as.integer(rowSums(species))
observed_zero <- fish_depth == 0L
assert(sum(observed_zero) == 18L,
       "Expected exactly 18 observed all-target-zero libraries")

# A portable observed-library table.  Prefixes keep the two namespaces explicit
# and avoid relying on a sidecar mapping when a file is read outside R.
output <- data.frame(
  library_id = ids,
  site = site,
  timepoint = timepoint,
  replicate = replicate,
  unit_id = unit_id,
  section = section,
  season = covariates$season,
  surveyed = TRUE,
  usable_timepoint = TRUE,
  reprocessed_observed = TRUE,
  structural_absence = FALSE,
  analytical_exclusion = FALSE,
  observed_zero = observed_zero,
  observation_state = ifelse(observed_zero, "observed_zero", "observed_positive"),
  surviving_libraries_in_unit = as.integer(unit_n[unit_id]),
  fish_read_depth = fish_depth,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
species_out <- as.data.frame(species, check.names = FALSE)
names(species_out) <- paste0("analysis_label__", names(species_out))
genus_out <- as.data.frame(genus, check.names = FALSE)
names(genus_out) <- paste0("genus__", names(genus_out))
output <- cbind(output, species_out, genus_out)

ord <- order(as.integer(sub("E", "", output$site, fixed = TRUE)),
             as.integer(sub("T", "", output$timepoint, fixed = TRUE)),
             as.integer(sub("R", "", output$replicate, fixed = TRUE)))
output <- output[ord, , drop = FALSE]
rownames(output) <- NULL
assert(nrow(output) == 647L && ncol(output) == 42L,
       "Portable impact table must be 647 rows x 42 columns")
assert(sum(startsWith(names(output), "analysis_label__")) == 14L &&
       sum(startsWith(names(output), "genus__")) == 12L,
       "Portable count namespaces have unexpected dimensions")

dir.create(dirname(paths$output), recursive = TRUE, showWarnings = FALSE)
write_tsv(output, paths$output)

receipt <- data.frame(
  check_id = c(
    "species_dimensions", "genus_dimensions", "library_id_alignment",
    "whole_nonnegative_counts", "genus_exact_collapse", "survey_mask",
    "realized_units", "balanced_panel", "observed_zero_libraries",
    "species_genus_depth_identity", "portable_dimensions",
    "species_input_sha256", "genus_input_sha256", "covariates_input_sha256",
    "occasions_input_sha256", "portable_output_sha256"
  ),
  status = "PASS",
  observed = c(
    paste(dim(species), collapse = "x"), paste(dim(genus), collapse = "x"),
    "identical 647-row IDs", "finite; non-negative; integer",
    "12 exact genus sums", "647/647 surveyed and usable",
    paste0(length(unit_n), "; libraries/unit=", paste(range(unit_n), collapse = "-")),
    paste0(nrow(grid), " sites x ", length(common_timepoints), " common timepoints"),
    as.character(sum(observed_zero)), "647/647 equal row depths",
    paste(dim(output), collapse = "x"), sha256(paths$species), sha256(paths$genus),
    sha256(paths$covariates), sha256(paths$occasions), sha256(paths$output)
  ),
  expected = c(
    "647x14", "647x12", "identical 647-row IDs", "finite; non-negative; integer",
    "12 exact genus sums", "647/647 surveyed and usable",
    "234; libraries/unit=1-3", "14 sites x 11 common timepoints", "18",
    "647/647 equal row depths", "647x42", "recorded digest", "recorded digest",
    "recorded digest", "recorded digest", "recorded digest"
  ),
  stringsAsFactors = FALSE
)
write_tsv(receipt, paths$receipt)

cat(sprintf(
  "Exported %s: %d observed libraries, %d labels + %d genera, %d units, %d observed zeros.\n",
  paths$output, nrow(output), ncol(species), ncol(genus), length(unit_n),
  sum(observed_zero)))
cat("All dimension, mask, metadata, whole-count, collapse, and depth gates passed.\n")
