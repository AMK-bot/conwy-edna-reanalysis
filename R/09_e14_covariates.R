#!/usr/bin/env Rscript
# 09_e14_covariates.R  --  Resolve the Phase-3 Gate-1 E14 covariate gap (eDNAplan.md §2.2, Gate 1).
#
# The authors NEVER published E14 (Trefriw)'s river_gradient_percent / temp_monthly_average /
# rainfall_monthly_average -- verified all-NA in their own Supplementary Data 1 (MOESM4_ESM.csv).
# So there is nothing to copy; the gap is theirs, not a loss on our side. Two covariates are
# resolved here by the pre-authorised, fully-autonomous path (no CEDA account, no download):
#
#   temp_monthly_average, rainfall_monthly_average  -> E13-DONOR approximation.
#       E14 is 3.8 km directly downstream of E13 (nearest surveyed neighbour), and the two sites share
#       an identical 18-timepoint set, so E14[T##] := E13[T##]. The source article describes weather
#       observations from the Glasgwm Automatic Weather Station, not a local gridded product. This copy
#       is therefore only a flagged, unvalidated supporting-data
#       imputation, not an observed E14 value and not evidence that its error is small. Climate does
#       not enter the executed headline formula. The flag is NEVER silent.
#
#   river_gradient_percent  -> DROPPED (pre-authorised fallback, risk register). It is a local 25 m-
#       buffer slope the authors under-specify, hard to reproduce, and Dist_from_lake (E14 = 35160,
#       the largest = terminal downstream site) already encodes E14's longitudinal position, so the
#       marginal value is low. Left NA; the Phase-3 covariate occumb omits it (psi ~ Dist_from_lake
#       + season). Escalation path if a reviewer demands it: NRW LiDAR DTM 1 m -> terra slope.
#
# Future expansion (not required by the executed formula): retrieve and validate the original
# weather-station lineage or declare a method-divergent climate source before adding climate effects.
#
# Output: data/metadata/conwy_covariates_e14imputed.tsv  (= conwy_covariates.tsv + an
# e14_climate_imputed flag column + E14 temp/rainfall filled). Idempotent.

CV_IN  <- "data/metadata/conwy_covariates.tsv"
CV_OUT <- "data/metadata/conwy_covariates_e14imputed.tsv"

cv <- read.delim(CV_IN, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("ID" %in% names(cv))
cv$site <- sub("_.*", "", cv$ID)
cv$tp   <- sub(".*_(T[0-9]+)$", "\\1", cv$ID)

# --- E13 donor: one temp/rainfall value per timepoint (they are monthly-constant) ---
e13 <- cv[cv$site == "E13", ]
donor <- data.frame(
  tp   = e13$tp,
  temp = as.numeric(e13$temp_monthly_average),
  rain = as.numeric(e13$rainfall_monthly_average)
)
donor <- unique(donor)
stopifnot(!anyNA(donor$temp), !anyNA(donor$rain),          # donor complete
          !anyDuplicated(donor$tp))                        # exactly one value per timepoint

# --- transparency flag, then fill E14 temp/rainfall from the donor by matching timepoint ---
cv$e14_climate_imputed <- FALSE
is14 <- cv$site == "E14"
m <- match(cv$tp[is14], donor$tp)
stopifnot(!anyNA(m))                                       # every E14 timepoint has an E13 donor
cv$temp_monthly_average[is14]     <- donor$temp[m]
cv$rainfall_monthly_average[is14] <- donor$rain[m]
cv$e14_climate_imputed[is14]      <- TRUE
# river_gradient_percent for E14 stays NA (dropped per fallback) -- intentionally not filled.

out <- cv[, setdiff(names(cv), c("site", "tp"))]           # drop the working cols
write.table(out, CV_OUT, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# --- self-check ---
chk <- read.delim(CV_OUT, check.names = FALSE, stringsAsFactors = FALSE)
c14 <- chk[sub("_.*", "", chk$ID) == "E14", ]
c13 <- chk[sub("_.*", "", chk$ID) == "E13", ]
c14$tp <- sub(".*_(T[0-9]+)$", "\\1", c14$ID); c13$tp <- sub(".*_(T[0-9]+)$", "\\1", c13$ID)
stopifnot(
  nrow(chk) == nrow(cv),                                                    # no rows lost
  all(!is.na(as.numeric(c14$temp_monthly_average))),                        # E14 temp now complete
  all(!is.na(as.numeric(c14$rainfall_monthly_average))),                    # E14 rainfall now complete
  all(is.na(suppressWarnings(as.numeric(c14$river_gradient_percent)))),     # gradient still dropped
  sum(chk$e14_climate_imputed == "TRUE") == 54,                             # flag only on E14's 54 rows
  # imputed E14 temp/rainfall equal the E13 donor by timepoint
  all(mapply(function(t, v) isTRUE(all.equal(as.numeric(v),
        as.numeric(c13$temp_monthly_average[match(t, c13$tp)]))),
        c14$tp, c14$temp_monthly_average))
)
cat(sprintf("E14 supporting table: unvalidated E13 donor copied to temp/rainfall (54 rows flagged); gradient remains NA and both are excluded from the headline formula.\n-> %s\nSELF-CHECK PASSED\n", CV_OUT))
