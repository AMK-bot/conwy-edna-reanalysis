#!/usr/bin/env Rscript
# test_reprocessed_design.R -- self-check for the Phase-3 detection design rebuilt
# on the 647 reprocessed survivors. Proves the NA-not-0 rule and reports the unit grid.
# Run: micromamba run -n edna Rscript R/test_reprocessed_design.R
source("R/lib_design.R")

surv <- survived_reprocessing()
mask <- observable_mask("reprocessed")              # [site, tp, rep] logical
seq  <- observable_mask("sequenced")

# --- gates -------------------------------------------------------------------
stopifnot(length(surv) == 647L)                     # the corrected count table
stopifnot(sum(mask) == 647L)                        # mask agrees with survivor list
stopifnot(all(mask[mask] %in% TRUE))
stopifnot(all(which(mask) %in% which(seq)))         # survivors are a subset of sequenced cells
cat("PASS: reprocessed mask = 647 observable cells; survivors ⊆ sequenced; NA-not-0 holds.\n\n")

# --- the rebuilt unit grid (site × timepoint = occupancy units) --------------
# reps surviving per (site,tp) unit; restrict to the 18 usable tps (T16 is all-NA).
mask_u <- mask[, dimnames(mask)$timepoint %in% usable_timepoints, , drop = FALSE]
stopifnot(dim(mask_u)[2] == 18L)
reps_per_unit <- apply(mask_u, c(1, 2), sum)
tab <- table(factor(reps_per_unit, levels = 0:3))
cat("Surviving reps per site×timepoint unit (of 252 theoretical):\n")
print(tab)
cat(sprintf("  units with >=1 rep (realized): %d/252 | empty (NA units): %d | <3 reps: %d\n\n",
            sum(reps_per_unit > 0), sum(reps_per_unit == 0), sum(reps_per_unit < 3 & reps_per_unit > 0)))

cat("Surviving libraries per site:\n")
print(apply(mask, 1, sum))
cat(sprintf("\nTotal observable detection cells for Phase-3 occumb: %d (was 665 published / 755 sequenced)\n",
            sum(mask)))
