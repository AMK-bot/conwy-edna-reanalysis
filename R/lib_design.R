# lib_design.R -- canonical Conwy occasion design for the occupancy scripts.
#
# Single source of truth for the "18 usable timepoints; T16 = unsurveyed (NA,
# NOT 0)" rule. Both 05_occupancy_msocc.R and 06_occupancy_occumb.R source this
# instead of re-deriving the site×time×rep grid from memory (where T16 tends to
# get silently encoded as all-zero non-detections — the ψ-biasing error this
# project exists to correct).
#
# Reads data/metadata/occasions.tsv (emitted by data/build_occasions.py).
# Usage:
#   source("R/lib_design.R")
#   occ <- conwy_occasions()                 # the full 798-row design tibble
#   usable_timepoints                        # 18 timepoints (T16 dropped)
#   obs <- observable_mask("published")      # site×tp×rep logical array; set y[!obs] <- NA

suppressPackageStartupMessages(library(readr))

.occasions_path <- function() {
  fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  here <- if (length(fa)) dirname(normalizePath(fa[1])) else "R"
  cand <- c(file.path(here, "..", "data", "metadata", "occasions.tsv"),
            "data/metadata/occasions.tsv")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("occasions.tsv not found — run: python3 data/build_occasions.py")
  normalizePath(hit[1])
}

conwy_occasions <- function() {
  readr::read_tsv(.occasions_path(), col_types = "cccciii", progress = FALSE)
}

# IDs that survived Phase-2 reprocessing into the decontaminated count table (647
# of 755). Keyed off the corrected seqtab_clean.rds rownames (NOT the design grid),
# so the Phase-3 occumb/Layer-E detection histories are built on libraries that
# actually carry data — the 108 dropped-or-empty libraries become NA, not false 0.
survived_reprocessing <- function() {
  cand <- c("data/derived/seqtab_clean.rds",
            "data/derived/species_x_sample_counts_clean.rds")
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("no reprocessed count table in data/derived/ — stage the Phase-2 handback")
  rn <- rownames(readRDS(hit[1]))
  if (!all(grepl("^E[0-9]{2}_R[0-9]_T[0-9]{2}$", rn)))
    stop("unexpected rownames in the reprocessed count table (expect E##_R#_T##)")
  sort(unique(rn))
}

# 18 usable sampling occasions (T16 excluded — unsurveyed for all 14 sites).
usable_timepoints <- local({
  occ <- conwy_occasions()
  sort(unique(occ$timepoint[occ$usable_timepoint == 1L]))
})

# Logical array [site, timepoint, rep] marking cells where a detection CAN be
# observed; everywhere else the detection history must be NA (not 0).
#   resolution = "published": observable = survived QC into the 665-row matrix
#                (the fast path; the 90 QC-dropped + 43 unsurveyed are NA).
#   resolution = "sequenced": observable = a library exists (deposited in ENA;
#                only the 43 unsurveyed design cells are NA).
#   resolution = "reprocessed": observable = the library survived Phase-2 into the
#                647-row decontaminated count table (the Phase-3 input). Everywhere
#                else — the 43 unsurveyed + the 108 dropped/empty libraries — is NA.
observable_mask <- function(resolution = c("published", "sequenced", "reprocessed")) {
  resolution <- match.arg(resolution)
  occ <- conwy_occasions()
  sites <- sort(unique(occ$site)); tps <- sort(unique(occ$timepoint)); reps <- sort(unique(occ$rep))
  if (resolution == "reprocessed") {
    surv <- survived_reprocessing()
    stopifnot(all(surv %in% occ$ID[occ$surveyed == 1L]))   # GATE: survivors ⊆ sequenced design cells
    flag <- as.integer(occ$ID %in% surv)
  } else {
    flag <- if (resolution == "published") occ$in_published_12S else occ$surveyed
  }
  arr <- array(FALSE, dim = c(length(sites), length(tps), length(reps)),
               dimnames = list(site = sites, timepoint = tps, rep = reps))
  idx <- cbind(match(occ$site, sites), match(occ$timepoint, tps), match(occ$rep, reps))
  arr[idx] <- flag == 1L
  arr
}
