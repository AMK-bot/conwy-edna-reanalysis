#!/usr/bin/env Rscript
# 06_occupancy_occumb.R -- Phase 3 flagship: community occupancy (occumb). TWO models:
#   (1) INTERCEPT-ONLY  psi ~ 1, phi ~ 1            (the E14-free baseline; kept unchanged)
#   (2) COVARIATE       psi ~ scale(Dist_from_lake) + season   <-- EXPLORATORY, CLAIM-GATED
# Plan eDNAplan.md §2.2 + §3 Phase 3; decision record docs/prereg.md §1. The executed fit is
# retained as a stress test, not a biological result, because recovery did not converge.
#
# COVARIATE MODEL (exploratory implementation). Occupancy is parameterised by downstream distance + season.
#   river_gradient_percent is DELIBERATELY DROPPED (E14 is NA -> pre-authorised fallback; do
#   NOT use gradient). Covariate source is data/metadata/conwy_covariates_e14imputed.tsv; its flagged,
#   unvalidated E14 climate donor fields are not used by this formula. Dist_from_lake is source-complete
#   across all 14 sites and is z-scored (scale()); season
#   is a factor with reference = Spring. Covariates live at the UNIT (site x timepoint) level and
#   enter occumbData(site_cov = ...) as length-J vectors aligned 1:1 to the cube's unit (dim-2)
#   order. Detection KEEPS the mandated depth safeguard (formula_theta_shared = ~ logdepth) and
#   dominance stays intercept-only (formula_phi = ~ 1). occumb parameter `gamma` = species effects
#   on psi; `psi` = per-(species,unit) occupancy prob.
#
# WHAT THIS FITS
#   The Fukaya et al. (2022) multispecies eDNA site-occupancy model via the `occumb`
#   package (JAGS backend). Model hierarchy (from ?occumb):
#       y[i,j,k] ~ Multinomial(pi[i,j,k], N[j,k])      # read counts; N = per-rep DEPTH
#       pi[i,j,k] = u*r / sum(u*r)
#       r[i,j,k]  ~ Gamma(phi[i,j,k], 1)               # sequence relative dominance
#       u[i,j,k]  ~ Bernoulli(z[i,j] * theta[i,j,k])   # per-replicate CAPTURE (detection)
#       z[i,j]    ~ Bernoulli(psi[i,j])                # site occupancy
#   Occupancy unit j = SITE x TIMEPOINT; replicate k = the 3 water samples.
#   INTERCEPT-ONLY here means psi ~ 1 and phi ~ 1 (NO occupancy/site covariates -- those
#   wait on E14's missing terrain/climate). Detection theta carries the MANDATED depth
#   covariate only (see below); that needs nothing from E14.
#
# DEPTH SAFEGUARD (eDNAplan §2.2, MANDATORY -- non-negotiable). The ~16% detection-cell
#   recall deficit is associated with the 35-column reprocessed raw-name-table count
#   total (not whole-library sequencing depth) and source-study relative abundance in
#   an unclustered diagnostic and can confound occupancy with detection;
#   its effect on psi magnitude or direction is not established. Depth enters detection via
#   TWO distinct, non-interchangeable channels (the plan flags this explicitly):
#     (A) NATIVE multinomial exposure. occumb sets N[j,k] = sum_i y[i,j,k] internally
#         (occumb:::set_const) = the total FISH (in-community) reads per replicate, used as
#         the multinomial size. This is automatic and models relative read partitioning.
#     (B) ADDED logit term on theta. The FULL library sequencing depth (rowSums of the whole
#         1348-ASV seqtab_clean, i.e. ALL reads incl. non-fish/unassigned) is z-scored log10
#         and passed as a REPLICATE covariate on theta via formula_theta_shared = ~ logdepth.
#         logit(theta[i,j,k]) = beta[i,1] + beta_shared[1] * logdepth[j,k]. A shallow library
#         has lower capture prob -> the model inflates psi to explain the missed detection
#         instead of treating a low-depth non-detection as a true absence.
#   (A) and (B) are NOT the same quantity: N[j,k] is in-community fish depth; logdepth is the
#   whole-library depth. Both are used.
#
# STRUCTURAL NA. occumb forbids NA in y and instead reads a replicate with total reads == 0 as
#   MISSING (occumb:::set_const converts any N==0 cell to y<-NA, N<-1 internally -> that cell
#   contributes NO detection likelihood). So the structural-NA map from observable_mask()
#   ("reprocessed") is encoded as ZERO-VECTORS: unobserved (site,tp,rep) cells are left as all
#   zeros and occumb auto-treats them as missing. This is the correct occumb idiom; do NOT
#   insert NA (occumbData rejects it).
#
# MIN READ-DEPTH FILTER (safeguard 3, a researcher degree-of-freedom -> planned sensitivity sweep).
#   Replicate cells whose FULL library depth < MIN_DEPTH are zeroed (-> occumb-missing) so
#   failed/near-failed libraries do not masquerade as confident non-detections. Primary
#   threshold documented below; the production run sweeps it and reports upstream-psi movement.
#
# INPUTS (committed/derived):
#   data/derived/species_x_sample_counts_clean.rds  647 x 14 (double -> cast integer)
#   data/derived/seqtab_clean.rds                    647 x 1348 (rowSums = full depth)
#   R/lib_design.R :: observable_mask("reprocessed") 14 x 19 x 3 structural map (647 TRUE)
# OUTPUTS (all "phase3_"/"occumb_" prefixed -- no collisions):
#   -- intercept model --
#   data/derived/occumb_intercept_estimates.tsv    per-species psi + detection + false-neg rate
#   data/derived/occumb_intercept_fit.rds          the occumbFit object
#   figures/phase3_occupancy_psi_forest.png        per-species psi (95% CrI)
#   figures/phase3_occupancy_detection_fnr.png     per-species detection + cumulative FN rate
#   -- covariate model (headline) --
#   data/derived/occumb_covariate_psi_effects.tsv       per-species Dist_from_lake/season coeffs
#   data/derived/occumb_covariate_corrected_vs_raw.tsv  corrected psi vs raw naive occupancy
#   data/derived/occumb_covariate_design_opt.tsv        eval_util_L design-optimization grid
#   figures/phase3_occupancy_cov_dist_effect.png        per-species downstream-distance coefficient
#   figures/phase3_occupancy_cov_corrected_vs_raw.png   raw vs corrected occupancy shift
#   figures/phase3_occupancy_cov_design.png             expected utility vs replicates/depth
#
# RUN:  micromamba run -n edna Rscript R/06_occupancy_occumb.R
#   Env knobs (optional): OCCUMB_SMOKE=0 for production settings; OCCUMB_MIN_DEPTH=<int>;
#     OCCUMB_SKIP_INTERCEPT=1 to fit ONLY the headline covariate model (dev speed);
#     OCCUMB_SAVE_COV_FIT=1 to persist the ~350MB covariate fit RDS (default: skip).
# Idempotent: overwrites its own outputs only. Sources R/lib_design.R (never re-derives grid).

suppressWarnings(suppressPackageStartupMessages({
  library(occumb)
}))

set.seed(20260705)

## ---- paths (robust to cwd) -------------------------------------------------
args_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
ROOT <- if (length(args_file)) {
  normalizePath(file.path(dirname(args_file[1]), ".."))
} else {
  normalizePath(getwd())
}
p <- function(...) file.path(ROOT, ...)
source(p("R", "lib_design.R"))

DERIVED <- p("data", "derived")
FIGDIR  <- p("figures")
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---- knobs -----------------------------------------------------------------
SMOKE     <- !identical(Sys.getenv("OCCUMB_SMOKE"), "0")   # default: bounded smoke
MIN_DEPTH <- as.integer(Sys.getenv("OCCUMB_MIN_DEPTH", "1000"))  # primary full-depth filter
SKIP_INTERCEPT <- identical(Sys.getenv("OCCUMB_SKIP_INTERCEPT"), "1")  # fit only headline cov model
SAVE_COV_FIT   <- identical(Sys.getenv("OCCUMB_SAVE_COV_FIT"),   "1")  # persist ~350MB cov fit
# Planned production sensitivity sweep (documented; not run): 0, 100, 250, 500, 1000, 2000.

# MCMC settings. Smoke = COMPLETES fast to prove the model fits & yields psi/detection.
# Production = long enough for Rhat < 1.01 & healthy n.eff (documented, run with OCCUMB_SMOKE=0).
if (SMOKE) {
  MCMC <- list(n.chains = 2, n.burnin = 1000, n.iter = 2000, n.thin = 1, parallel = TRUE)
} else {
  MCMC <- list(n.chains = 4, n.burnin = 20000, n.iter = 60000, n.thin = 20, parallel = TRUE)
}

message(sprintf("[06_occumb] mode=%s  MIN_DEPTH=%d  MCMC: %d chains, %d iter (%d burn), thin %d",
                if (SMOKE) "SMOKE" else "PRODUCTION", MIN_DEPTH,
                MCMC$n.chains, MCMC$n.iter, MCMC$n.burnin, MCMC$n.thin))

## ---- load & cast -----------------------------------------------------------
sp_dbl <- readRDS(file.path(DERIVED, "species_x_sample_counts_clean.rds"))  # 647 x 14 double
seqtab <- readRDS(file.path(DERIVED, "seqtab_clean.rds"))                    # 647 x 1348
stopifnot(is.matrix(sp_dbl), ncol(sp_dbl) == 14L)
sp_cnt <- matrix(as.integer(round(sp_dbl)), nrow(sp_dbl), ncol(sp_dbl),
                 dimnames = dimnames(sp_dbl))
species <- colnames(sp_cnt)
lib_ids <- rownames(sp_cnt)
stopifnot(all(grepl("^E[0-9]{2}_R[0-9]_T[0-9]{2}$", lib_ids)))

# FULL library sequencing depth = total reads over ALL 1348 ASVs (safeguard channel B).
full_depth <- rowSums(seqtab)[lib_ids]
stopifnot(!anyNA(full_depth))

## ---- parse IDs -> site / timepoint / rep -----------------------------------
site <- substr(lib_ids, 1, 3)      # E##
rep  <- substr(lib_ids, 5, 6)      # R#
tp   <- substr(lib_ids, 8, 10)     # T##
unit_key <- paste(site, tp, sep = "_")            # site x timepoint = occupancy unit

units <- sort(unique(unit_key))                    # J occupancy units (>=1 observed rep)
reps  <- sort(unique(rep))                          # R1 R2 R3
I <- length(species); J <- length(units); K <- length(reps)
message(sprintf("[06_occumb] I=%d species, J=%d site x timepoint units, K=%d replicates", I, J, K))

j_idx <- match(unit_key, units)
k_idx <- match(rep, reps)

## ---- build y[species, unit, replicate] -------------------------------------
# Structural NA idiom: start ALL zero; fill only the 647 observed libraries. occumb reads any
# remaining all-zero cell as missing (N==0 -> internal NA). We track structural placement to
# assert exactly 647 observed cells (independent of the depth filter).
y <- array(0L, dim = c(I, J, K), dimnames = list(species = species, unit = units, rep = reps))
obs_cell <- matrix(FALSE, J, K, dimnames = list(units, reps))  # a library was placed here
depth_cell <- matrix(NA_real_, J, K)                            # full depth per placed library

for (n in seq_along(lib_ids)) {
  jj <- j_idx[n]; kk <- k_idx[n]
  y[, jj, kk]        <- sp_cnt[n, ]
  obs_cell[jj, kk]   <- TRUE
  depth_cell[jj, kk] <- full_depth[n]
}
n_observed <- sum(obs_cell)
message(sprintf("[06_occumb] observed (site,tp,rep) cells placed: %d", n_observed))

# Cross-check placement against the canonical structural map (both must be 647).
mask <- observable_mask("reprocessed")
message(sprintf("[06_occumb] observable_mask('reprocessed') TRUE cells: %d", sum(mask)))

## ---- min read-depth filter (safeguard 3) -----------------------------------
# Zero out cells whose FULL library depth < MIN_DEPTH -> occumb treats them as missing.
filt <- obs_cell & (!is.na(depth_cell)) & (depth_cell < MIN_DEPTH)
n_filtered <- sum(filt)
for (jj in seq_len(J)) for (kk in seq_len(K)) if (filt[jj, kk]) y[, jj, kk] <- 0L
message(sprintf("[06_occumb] min-depth filter (< %d reads) zeroed %d of %d observed cells",
                MIN_DEPTH, n_filtered, n_observed))

# Informative cells = observed, survived filter, AND carry >0 in-community (fish) reads.
# (occumb's native N==0 rule also drops observed libraries whose 14-species total is 0.)
fish_tot   <- apply(y, c(2, 3), sum)
informative <- (fish_tot > 0)
message(sprintf("[06_occumb] informative detection cells (N>0 after filter): %d", sum(informative)))

## ---- depth covariate matrix (safeguard channel B) --------------------------
# z-scored log10(full depth) over informative cells; non-informative cells -> 0 (the mean;
# harmless because occumb drops them from the theta likelihood). NO NA allowed by occumbData.
logd <- matrix(0, J, K, dimnames = list(units, reps))
info_idx <- which(informative, arr.ind = TRUE)
ld_vals <- log10(depth_cell[informative])
ld_z <- as.numeric(scale(ld_vals))                 # mean 0, sd 1 over informative cells
logd[informative] <- ld_z
stopifnot(!anyNA(logd))

## ---- assemble occumbData ---------------------------------------------------
dat <- occumbData(y = y, repl_cov = list(logdepth = logd))
cat("\n================ occumbData summary ================\n")
summary(dat)

## ---- self-check on the data cube -------------------------------------------
stopifnot(
  "observed cells must equal 647" = n_observed == 647L,
  "mask TRUE cells must equal 647" = sum(mask) == 647L,
  "y must be integer" = is.integer(y),
  "no NA in y" = !anyNA(y)
)

## ---- FIT: intercept-only occupancy + depth on detection --------------------
# Wrapped so OCCUMB_SKIP_INTERCEPT=1 can jump straight to the headline covariate model.
if (!SKIP_INTERCEPT) {
cat("\n================ fitting occumb (intercept-only psi/phi; theta ~ logdepth) ================\n")
t0 <- Sys.time()
fit <- occumb(
  formula_psi          = ~ 1,          # INTERCEPT-ONLY occupancy (no covariates -> no E14 gate)
  formula_phi          = ~ 1,          # intercept-only relative dominance
  formula_theta_shared = ~ logdepth,   # MANDATED depth safeguard on detection (shared slope)
  data     = dat,
  n.chains = MCMC$n.chains,
  n.burnin = MCMC$n.burnin,
  n.iter   = MCMC$n.iter,
  n.thin   = MCMC$n.thin,
  parallel = MCMC$parallel,
  engine   = "JAGS"
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
message(sprintf("[06_occumb] fit complete in %.1f min", elapsed))

## ---- extract per-species psi / detection / false-negative ------------------
psi_s <- get_post_summary(fit, "psi")     # rows psi[i]  (species; intercept-only collapses site)
psi_s <- as.matrix(psi_s)

# theta at REFERENCE (mean) depth = plogis(beta[i,1]); use raw draws for CI + FN posterior.
beta_draws <- get_post_samples(fit, "beta")            # [draws, species, effects]
d <- dim(beta_draws)
beta_int <- if (length(d) == 3) beta_draws[, , 1] else beta_draws  # species intercept draws
theta_ref <- plogis(beta_int)                          # [draws, species] detection @ mean depth
qs <- function(m, f) apply(m, 2, f)
sumr <- function(m) cbind(
  mean = colMeans(m),
  l95  = qs(m, function(x) quantile(x, 0.025)),
  med  = qs(m, function(x) median(x)),
  u95  = qs(m, function(x) quantile(x, 0.975))
)
theta_tab <- sumr(theta_ref)

# False-negative rate given occupancy: per replicate = (1 - theta); across K=3 reps = (1-theta)^3.
fn_rep <- 1 - theta_ref
fn_3   <- (1 - theta_ref)^3
fn_rep_tab <- sumr(fn_rep)
fn_3_tab   <- sumr(fn_3)

# shared depth slope on theta (one value -> get_post_summary returns a NAMED VECTOR, not a matrix).
bshared_raw <- get_post_summary(fit, "beta_shared")
pick <- function(x, nm) if (is.null(dim(x))) x[[nm]] else x[1, nm]
bshared <- c(mean = pick(bshared_raw, "mean"),
             l95  = pick(bshared_raw, "2.5%"),
             u95  = pick(bshared_raw, "97.5%"),
             Rhat = pick(bshared_raw, "Rhat"))

est <- data.frame(
  species              = species,
  psi_mean             = psi_s[, "mean"],
  psi_l95              = psi_s[, "2.5%"],
  psi_med              = psi_s[, "50%"],
  psi_u95              = psi_s[, "97.5%"],
  psi_Rhat             = psi_s[, "Rhat"],
  psi_neff             = psi_s[, "n.eff"],
  detection_perrep_mean = theta_tab[, "mean"],   # theta @ mean sequencing depth
  detection_perrep_l95  = theta_tab[, "l95"],
  detection_perrep_u95  = theta_tab[, "u95"],
  false_neg_perrep_mean = fn_rep_tab[, "mean"],  # 1 - theta
  false_neg_3rep_mean   = fn_3_tab[, "mean"],    # (1 - theta)^3  (missed at a 3-rep unit)
  false_neg_3rep_u95    = fn_3_tab[, "u95"],
  row.names = NULL, check.names = FALSE
)
est <- est[order(-est$psi_mean), ]
write.table(est, file.path(DERIVED, "occumb_intercept_estimates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[06_occumb] wrote data/derived/occumb_intercept_estimates.tsv")

## ---- convergence overview --------------------------------------------------
cat("\n================ convergence (occumb summary) ================\n")
summary(fit)
conv_note <- sprintf(
  "psi Rhat: max %.3f (min n.eff %d) | theta_intercept(beta) draws summarized | depth slope beta_shared mean %.3f [%.3f, %.3f]",
  max(psi_s[, "Rhat"], na.rm = TRUE), as.integer(min(psi_s[, "n.eff"], na.rm = TRUE)),
  bshared[["mean"]], bshared[["l95"]], bshared[["u95"]]
)
message("[06_occumb] ", conv_note)

## ---- save fit --------------------------------------------------------------
saveRDS(fit, file.path(DERIVED, "occumb_intercept_fit.rds"))
message("[06_occumb] saved data/derived/occumb_intercept_fit.rds")

## ---- figures ---------------------------------------------------------------
# 1) per-species psi forest plot (95% CrI), sorted.
ord <- order(est$psi_mean)
png(file.path(FIGDIR, "phase3_occupancy_psi_forest.png"), width = 1100, height = 900, res = 130)
op <- par(mar = c(5, 14, 4, 2))
yy <- seq_len(nrow(est))
plot(est$psi_mean[ord], yy, xlim = c(0, 1), pch = 19, col = "#1b6ca8",
     xlab = "Occupancy probability  psi  (intercept-only, 95% CrI)", ylab = "",
     yaxt = "n", main = "occumb intercept-only: per-species occupancy")
segments(est$psi_l95[ord], yy, est$psi_u95[ord], yy, col = "#1b6ca8", lwd = 2)
axis(2, at = yy, labels = gsub("_", " ", est$species[ord]), las = 1, cex.axis = 0.8)
abline(v = c(0.5), lty = 3, col = "grey60")
par(op); dev.off()
message("[06_occumb] wrote figures/phase3_occupancy_psi_forest.png")

# 2) per-species detection (theta @ mean depth) + cumulative false-negative over 3 reps.
png(file.path(FIGDIR, "phase3_occupancy_detection_fnr.png"), width = 1200, height = 900, res = 130)
op <- par(mar = c(5, 14, 4, 2))
ordd <- order(est$detection_perrep_mean)
yy <- seq_len(nrow(est))
plot(est$detection_perrep_mean[ordd], yy, xlim = c(0, 1), pch = 19, col = "#c0392b",
     xlab = "probability", ylab = "", yaxt = "n",
     main = "Detection per replicate (red) vs 3-rep false-negative rate (grey)")
segments(est$detection_perrep_l95[ordd], yy, est$detection_perrep_u95[ordd], yy,
         col = "#c0392b", lwd = 2)
points(est$false_neg_3rep_mean[ordd], yy, pch = 17, col = "grey35")
axis(2, at = yy, labels = gsub("_", " ", est$species[ordd]), las = 1, cex.axis = 0.8)
legend("bottomright", c("detection / replicate (theta)", "false-negative over 3 reps"),
       pch = c(19, 17), col = c("#c0392b", "grey35"), bty = "n", cex = 0.85)
par(op); dev.off()
message("[06_occumb] wrote figures/phase3_occupancy_detection_fnr.png")

## ---- INTERCEPT SELF-CHECK (runnable) ---------------------------------------
cat("\n================ INTERCEPT SELF-CHECK ================\n")
stopifnot(
  "observed cells == 647"        = n_observed == 647L,
  "psi in [0,1]"                 = all(est$psi_mean >= 0 & est$psi_mean <= 1),
  "detection in [0,1]"           = all(est$detection_perrep_mean >= 0 & est$detection_perrep_mean <= 1),
  "false-neg in [0,1]"           = all(est$false_neg_3rep_mean >= 0 & est$false_neg_3rep_mean <= 1),
  "estimates has 14 analysis labels" = nrow(est) == 14L
)
cat("INTERCEPT SELF-CHECK PASSED:\n")
cat(sprintf("  observed cells = %d (== 647)\n", n_observed))
cat(sprintf("  psi range   = [%.3f, %.3f]\n", min(est$psi_mean), max(est$psi_mean)))
cat(sprintf("  detection/rep range = [%.3f, %.3f]\n",
            min(est$detection_perrep_mean), max(est$detection_perrep_mean)))
cat(sprintf("  depth slope on theta (beta_shared) = %.3f [%.3f, %.3f]\n",
            bshared[["mean"]], bshared[["l95"]], bshared[["u95"]]))
cat(sprintf("  min-depth filter = %d reads; %d cells zeroed; %d informative cells\n",
            MIN_DEPTH, n_filtered, sum(informative)))
cat("Intercept model done.\n")
} else {
  message("[06_occumb] OCCUMB_SKIP_INTERCEPT=1 -> skipping intercept model; fitting headline covariate model only")
}  # end if(!SKIP_INTERCEPT)


## ============================================================================
## COVARIATE MODEL (PHASE-3 EXPLORATORY STRESS TEST): psi ~ scale(Dist_from_lake) + season
## ----------------------------------------------------------------------------
## Occupancy structured by downstream distance + season; gradient DROPPED (E14 NA).
## Reuses the SAME filtered detection cube `y` and depth covariate `logd` as the
## intercept model, so the mandated depth safeguard is identical. Covariates are
## joined by ID at the UNIT (site x timepoint) level and passed as length-J
## site_cov vectors aligned 1:1 to the cube's unit dimension order (`units`).
## ============================================================================
cat("\n\n################ COVARIATE MODEL (Phase-3 exploratory; claim-gated) ################\n")

## ---- load imputed covariate table & collapse to per-unit -------------------
COVTSV <- p("data", "metadata", "conwy_covariates_e14imputed.tsv")
stopifnot("imputed covariate table missing" = file.exists(COVTSV))
cov_raw <- read.delim(COVTSV, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot("ID column required"            = "ID" %in% names(cov_raw),
          "Dist_from_lake column required" = "Dist_from_lake" %in% names(cov_raw),
          "season column required"         = "season" %in% names(cov_raw))
# gradient is DELIBERATELY untouched (E14 NA -> pre-authorised fallback).
cov_raw$.site <- substr(cov_raw$ID, 1, 3)      # E##  -> Dist_from_lake (constant within site)
cov_raw$.tp   <- substr(cov_raw$ID, 8, 10)     # T##  -> season         (constant within timepoint)
dist_by_site  <- tapply(cov_raw$Dist_from_lake, cov_raw$.site, function(x) x[1])
seas_by_tp    <- tapply(cov_raw$season,         cov_raw$.tp,   function(x) x[1])

# align to the cube's unit order (`units` = sort(unique("E##_T##"))).
u_site <- substr(units, 1, 3)
u_tp   <- substr(units, 5, 7)
u_dist <- as.numeric(dist_by_site[u_site])
u_seas <- factor(seas_by_tp[u_tp], levels = c("Spring", "Summer", "Autumn", "Winter"))  # ref=Spring
stopifnot("Dist_from_lake join left NA (source field should be complete)" = !anyNA(u_dist),
          "season join left NA"                                 = !anyNA(u_seas),
          "covariate rows must align 1:1 with J units"          = length(u_dist) == J)

# Authoritative effect ordering: occumb builds this same model.matrix internally.
site_cov_df <- data.frame(Dist_from_lake = u_dist, season = u_seas)
mm_psi    <- model.matrix(~ scale(Dist_from_lake) + season, site_cov_df)
eff_names <- colnames(mm_psi)  # (Intercept), scale(Dist_from_lake), seasonSummer/Autumn/Winter
stopifnot("covariate design matrix has NA (formula fields must be complete)" = !anyNA(mm_psi))
message(sprintf("[06_occumb][cov] J=%d units; psi effects: %s", J, paste(eff_names, collapse = " | ")))

## ---- assemble covariate occumbData -----------------------------------------
dat_cov <- occumbData(
  y        = y,                                        # SAME filtered cube as intercept model
  site_cov = list(Dist_from_lake = u_dist, season = u_seas),
  repl_cov = list(logdepth = logd)                     # depth safeguard, unchanged
)
cat("\n---------------- occumbData (covariate) summary ----------------\n")
summary(dat_cov)

## ---- RAW naive occupancy (this reanalysis's own, pre-model) ----------------
# proportion of surveyed site x timepoint units with >=1 detection, per species,
# on the UNFILTERED counts (what a naive analyst would tabulate).
det_raw <- matrix(0L, I, J, dimnames = list(species, units))
for (n in seq_along(lib_ids)) {
  jj <- j_idx[n]
  det_raw[, jj] <- pmax(det_raw[, jj], as.integer(sp_cnt[n, ] > 0L))
}
raw_occ <- rowMeans(det_raw)   # naive occupancy per species (proportion of units detected)

## ---- FIT covariate model ---------------------------------------------------
# Free the (large) intercept occumbFit BEFORE forking the covariate parallel workers. Retaining
# it caused memory pressure that produced an EINTR ("Interrupted system call") at the covariate
# cluster's socket setup in the full default run; the intercept fit is already persisted to RDS.
if (exists("fit")) { rm(fit); message("[06_occumb][cov] released intercept fit from memory pre-fork") }
invisible(gc())
cat("\n======== fitting occumb COVARIATE (psi ~ scale(Dist_from_lake)+season; phi~1; theta~logdepth) ========\n")
t0c <- Sys.time()
fit_cov <- NULL
for (attempt in 1:3L) {
  fit_cov <- tryCatch(
    occumb(
      formula_psi          = ~ scale(Dist_from_lake) + season,  # exploratory occupancy covariates
      formula_phi          = ~ 1,                                # dominance intercept-only
      formula_theta_shared = ~ logdepth,                         # MANDATED depth safeguard (kept)
      data     = dat_cov,
      n.chains = MCMC$n.chains, n.burnin = MCMC$n.burnin,
      n.iter   = MCMC$n.iter,   n.thin  = MCMC$n.thin,
      parallel = MCMC$parallel, engine  = "JAGS"
    ),
    error = function(e) { message("[06_occumb][cov] fit attempt ", attempt, " errored: ", conditionMessage(e)); NULL })
  if (!is.null(fit_cov)) break
  invisible(gc()); Sys.sleep(3)
}
if (is.null(fit_cov)) stop("[06_occumb][cov] covariate occumb fit failed after 3 attempts")
elapsed_c <- as.numeric(difftime(Sys.time(), t0c, units = "mins"))
message(sprintf("[06_occumb][cov] covariate fit complete in %.1f min", elapsed_c))

## ---- (1) per-species psi covariate effects (gamma) -------------------------
# gamma[draws, species, effect]; effect order == eff_names. Summary df carries Rhat/n.eff.
g_draws <- get_post_samples(fit_cov, "gamma")
dg <- dim(g_draws)
stopifnot("unexpected gamma dims" = length(dg) == 3L, dg[2] == I, dg[3] == length(eff_names))
gsd <- get_post_summary(fit_cov, "gamma", output_dataframe = TRUE)
# Species column may be a name (real data has dimnames) or an index; normalise to name.
sp_lab <- gsd$Species
gsd$species_name <- if (all(as.character(sp_lab) %in% species)) as.character(sp_lab) else species[as.integer(sp_lab)]
qL <- function(x) as.numeric(quantile(x, 0.025))
qU <- function(x) as.numeric(quantile(x, 0.975))

gamma_tab <- data.frame(
  species  = gsd$species_name,
  effect   = as.character(gsd$Effects),
  mean     = gsd$mean,
  l95      = gsd[["2.5%"]],
  med      = gsd[["50%"]],
  u95      = gsd[["97.5%"]],
  Rhat     = gsd$Rhat,
  n_eff    = gsd$n.eff,
  prob_neg = NA_real_,   # posterior P(coef < 0), filled from draws below
  cri_excludes_0 = (gsd[["2.5%"]] > 0) | (gsd[["97.5%"]] < 0),
  row.names = NULL, check.names = FALSE
)
# posterior sign mass from draws (map each summary row back to its (species,effect) slot)
for (r in seq_len(nrow(gamma_tab))) {
  ii <- match(gamma_tab$species[r], species)
  ee <- match(gamma_tab$effect[r],  eff_names)
  if (!is.na(ii) && !is.na(ee)) gamma_tab$prob_neg[r] <- mean(g_draws[, ii, ee] < 0)
}
# order: distance effect first (headline), then season, then intercept; species by effect size
eff_rank <- match(gamma_tab$effect, c("scale(Dist_from_lake)", "seasonSummer", "seasonAutumn",
                                      "seasonWinter", "(Intercept)"))
gamma_tab <- gamma_tab[order(eff_rank, gamma_tab$mean), ]
write.table(gamma_tab, file.path(DERIVED, "occumb_covariate_psi_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[06_occumb][cov] wrote data/derived/occumb_covariate_psi_effects.tsv")

# which species are SPATIALLY STRUCTURED by downstream distance (95% CrI excludes 0)?
dist_rows <- gamma_tab[gamma_tab$effect == "scale(Dist_from_lake)", ]
spatial_sp <- dist_rows$species[dist_rows$cri_excludes_0]
message(sprintf("[06_occumb][cov] Dist_from_lake structures psi (95%% CrI excl. 0) for %d/%d species: %s",
                length(spatial_sp), I,
                if (length(spatial_sp)) paste(spatial_sp, collapse = ", ") else "(none at this MCMC length)"))

## ---- (2) CORRECTED (occumb psi) vs RAW (naive) occupancy -------------------
# psi[draws, species, unit]; corrected per species = posterior mean of the unit-average psi.
psi_draws <- get_post_samples(fit_cov, "psi")
stopifnot("unexpected psi dims" = length(dim(psi_draws)) == 3L, dim(psi_draws)[2] == I)
psi_avg <- apply(psi_draws, c(1, 2), mean)          # [draws, species] mean occupancy over units
corr_mean <- colMeans(psi_avg)
corr_l95  <- apply(psi_avg, 2, qL)
corr_u95  <- apply(psi_avg, 2, qU)
names(corr_mean) <- names(corr_l95) <- names(corr_u95) <-
  if (!is.null(dimnames(psi_draws)) && !is.null(dimnames(psi_draws)[[2]])) dimnames(psi_draws)[[2]] else species

cvr <- data.frame(
  species        = species,
  raw_naive_occ  = as.numeric(raw_occ[species]),
  corrected_psi  = as.numeric(corr_mean[species]),
  corrected_l95  = as.numeric(corr_l95[species]),
  corrected_u95  = as.numeric(corr_u95[species]),
  row.names = NULL, check.names = FALSE
)
cvr$shift <- cvr$corrected_psi - cvr$raw_naive_occ           # >0 = occumb corrects UP (false negs)
cvr <- cvr[order(-cvr$shift), ]
write.table(cvr, file.path(DERIVED, "occumb_covariate_corrected_vs_raw.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[06_occumb][cov] wrote data/derived/occumb_covariate_corrected_vs_raw.tsv")
message(sprintf("[06_occumb][cov] corrected-vs-raw shift: median %+.3f, max %+.3f (%s); largest shifts are the rarer/low-detection taxa",
                median(cvr$shift), max(cvr$shift), cvr$species[which.max(cvr$shift)]))

## ---- (3) DESIGN OPTIMIZATION via occumb eval_util_L ------------------------
# TWO occumb gotchas handled here:
#   (i)  theta is replicate-specific (formula_theta_shared=~logdepth), so eval_util_L will NOT
#        accept the fit's theta directly (it assumes replicate-homogeneous theta). We supply
#        theta explicitly at the REFERENCE (mean) sequencing depth [logdepth=0] = plogis(beta
#        intercept), and phi from the fit.
#   (ii) eval_util_L rejects any posterior z with an all-absent site-draw (a fully empty
#        community -> 0/0 in the utility). Our rare-taxa posterior z has such units. So, per the
#        occumb vignette (util5 pattern), we evaluate the utility on a REPRESENTATIVE full-
#        community unit (z=1): U = expected number of the I focal species detected per unit as a
#        function of replicates K and per-replicate sequencing depth N, CONDITIONAL on presence.
#        This is the detection-side design curve (theta/phi posterior, homogeneous across the
#        representative unit) -- exactly the "recommended replicates/depth" question.
design_ok <- FALSE
design_note <- "design optimization not run"
util <- NULL
try({
  b_draws   <- get_post_samples(fit_cov, "beta")            # [draws, species, 1]
  theta_ref <- plogis(if (length(dim(b_draws)) == 3) b_draws[, , 1] else b_draws)  # [draws, species]
  phi_draws <- get_post_samples(fit_cov, "phi")             # [draws, species] (phi~1)
  phi_mat   <- if (length(dim(phi_draws)) == 3) phi_draws[, , 1] else phi_draws
  nd  <- nrow(theta_ref)
  keep <- if (SMOKE) min(300L, nd) else min(1500L, nd)
  idx <- sort(sample.int(nd, keep))
  theta_s <- theta_ref[idx, , drop = FALSE]
  phi_s   <- phi_mat[idx, , drop = FALSE]
  z_s     <- array(1L, dim = c(keep, I, 1))                 # representative full-community unit
  ncore   <- max(1L, min(4L, parallel::detectCores(logical = FALSE)))
  # (a) technical grid: expected #species detected per unit vs replicates K and depth N.
  grid <- expand.grid(K = 1:5, N = c(1e4, 5e4, 1e5, 5e5))
  util <- eval_util_L(grid, fit = NULL, z = z_s, theta = theta_s, phi = phi_s, cores = ncore)
  util$source <- "grid"
  # (b) budget-constrained frontier via list_cond_L. Cost UNITS are illustrative, but the budget
  #     is PEGGED to the realized design (J units, K=3 reps, ~N0 reads/rep) with 20% headroom, so
  #     the frontier explores trading replicates against sequencing depth around the actual study.
  #     lambda1 = cost per read; lambda2 = cost per replicate (library prep).
  lam1 <- 0.01; lam2 <- 5000; Kcur <- 3L; N0 <- 5e4
  # budget generous enough to afford up to K=5 reps at N0 reads across all J units; K auto-enumerated.
  budget <- lam2 * J * 5 + lam1 * J * 5 * N0
  util_b <- tryCatch({
    cond <- list_cond_L(budget = budget, lambda1 = lam1, lambda2 = lam2, fit = fit_cov)
    ub <- eval_util_L(cond, fit = NULL, z = z_s, theta = theta_s, phi = phi_s, cores = ncore)
    ub$source <- sprintf("budget=%.2g (read=%.3g, rep=%g, %d units)", budget, lam1, lam2, J)
    ub
  }, error = function(e) { message("[06_occumb][cov] list_cond_L skipped: ", conditionMessage(e)); NULL })
  design_tab <- if (!is.null(util_b)) {
    common <- intersect(names(util), names(util_b)); rbind(util[common], util_b[common])
  } else util
  write.table(design_tab, file.path(DERIVED, "occumb_covariate_design_opt.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("[06_occumb][cov] wrote data/derived/occumb_covariate_design_opt.tsv")

  # Recommendation from the technical grid at the study's typical depth (1e5): report the
  # DIMINISHING-RETURNS knee via marginal gain per added replicate (utility is # of the I focal
  # species detected/unit; with z=1 it keeps rising, so a saturation-threshold is meaningless --
  # the honest signal is where the marginal gain of one more replicate becomes small).
  Ntarget <- 1e5
  gg <- util[util$N == Ntarget, ]; gg <- gg[order(gg$K), ]
  dU  <- c(NA_real_, diff(gg$Utility))          # marginal utility gain of adding the k-th replicate
  thr <- 0.5                                    # <0.5 extra species per replicate = "small"
  knee <- which(dU < thr)[1]
  k_rec  <- if (!is.na(knee)) gg$K[knee] else max(gg$K)
  u_at3  <- gg$Utility[gg$K == Kcur]
  design_note <- sprintf(
    "eval_util_L: expected # of the %d focal species detected/unit at N=%.0e reads (z=1 representative community). Utility %s; diminishing-returns knee at K=%d reps (next replicate adds < %.1f species). Study's K=%d design -> %.2f species. Utilities by K: %s.",
    I, Ntarget,
    if (is.na(knee)) "still rising across K=1..5 (ultra-rare taxa need more replication)" else "rises then flattens",
    k_rec, thr, Kcur, if (length(u_at3)) u_at3 else NA_real_,
    paste(sprintf("K%d=%.2f", gg$K, gg$Utility), collapse = " "))
  design_ok <- TRUE
  cat("\n---------------- design-optimization utility grid (head) ----------------\n")
  print(utils::head(util, 12))
  message("[06_occumb][cov] DESIGN RECOMMENDATION: ", design_note)
}, silent = FALSE)
if (!design_ok) message("[06_occumb][cov] design optimization did not complete: ", design_note)

## ---- (4) CONVERGENCE (honest) ----------------------------------------------
gamma_rhat <- gamma_tab$Rhat[is.finite(gamma_tab$Rhat)]
gamma_neff <- gamma_tab$n_eff[is.finite(gamma_tab$n_eff)]
conv_note_cov <- sprintf(
  "psi coefficients (gamma): max Rhat %.3f, min n.eff %d over %d species x %d effects.",
  if (length(gamma_rhat)) max(gamma_rhat) else NA_real_,
  if (length(gamma_neff)) as.integer(min(gamma_neff)) else NA_integer_,
  I, length(eff_names))
message("[06_occumb][cov] ", conv_note_cov)
if (SMOKE) message("[06_occumb][cov] NOTE: SMOKE MCMC (", MCMC$n.chains, " chains, ",
                   MCMC$n.iter, " iter) will NOT converge the dominance (alpha/phi) or ",
                   "cross-species correlation (rho) params -- treat smoke gamma/psi as a ",
                   "plumbing check, not an inferential result. Use the production command below.")

## ---- optionally persist the (large) covariate fit --------------------------
if (SAVE_COV_FIT) {
  saveRDS(fit_cov, file.path(DERIVED, "occumb_covariate_fit.rds"))
  message("[06_occumb][cov] saved data/derived/occumb_covariate_fit.rds (large)")
} else {
  message("[06_occumb][cov] skipped saving covariate fit RDS (~350MB); set OCCUMB_SAVE_COV_FIT=1 to persist")
}

## ---- figures ---------------------------------------------------------------
# (A) per-species downstream-distance coefficient (95% CrI) -- suppressed by the release gate on failure.
dr <- dist_rows[order(dist_rows$mean), ]
png(file.path(FIGDIR, "phase3_occupancy_cov_dist_effect.png"), width = 1150, height = 900, res = 130)
op <- par(mar = c(5, 14, 4, 2))
yy <- seq_len(nrow(dr))
xr <- range(c(dr$l95, dr$u95, 0), na.rm = TRUE)
cols <- ifelse(dr$cri_excludes_0, ifelse(dr$mean < 0, "#c0392b", "#1b6ca8"), "grey60")
plot(dr$mean, yy, xlim = xr, pch = 19, col = cols, yaxt = "n", ylab = "",
     xlab = "psi coefficient on z-scored Dist_from_lake (logit; 95% CrI)",
     main = "Occupancy vs downstream distance (per species)")
segments(dr$l95, yy, dr$u95, yy, col = cols, lwd = 2)
axis(2, at = yy, labels = gsub("_", " ", dr$species), las = 1, cex.axis = 0.8)
abline(v = 0, lty = 3, col = "grey40")
legend("topleft", c("CrI excludes 0 (declines downstream)", "CrI excludes 0 (increases)", "CrI spans 0"),
       pch = 19, col = c("#c0392b", "#1b6ca8", "grey60"), bty = "n", cex = 0.8)
par(op); dev.off()
message("[06_occumb][cov] wrote figures/phase3_occupancy_cov_dist_effect.png")

# (B) corrected (occumb psi) vs raw (naive) occupancy, per species.
png(file.path(FIGDIR, "phase3_occupancy_cov_corrected_vs_raw.png"), width = 1200, height = 900, res = 130)
op <- par(mar = c(5, 14, 4, 2))
cc <- cvr[order(cvr$corrected_psi), ]
yy <- seq_len(nrow(cc))
plot(cc$corrected_psi, yy, xlim = c(0, 1), pch = 19, col = "#1b6ca8", yaxt = "n", ylab = "",
     xlab = "occupancy probability", main = "Corrected (occumb psi, 95% CrI) vs raw naive occupancy")
segments(cc$corrected_l95, yy, cc$corrected_u95, yy, col = "#1b6ca8", lwd = 2)
points(cc$raw_naive_occ, yy, pch = 4, col = "grey30", lwd = 2)
segments(cc$raw_naive_occ, yy, cc$corrected_psi, yy, col = "grey70", lty = 3)
axis(2, at = yy, labels = gsub("_", " ", cc$species), las = 1, cex.axis = 0.8)
legend("bottomright", c("corrected psi (occumb)", "raw naive occupancy"),
       pch = c(19, 4), col = c("#1b6ca8", "grey30"), bty = "n", cex = 0.85)
par(op); dev.off()
message("[06_occumb][cov] wrote figures/phase3_occupancy_cov_corrected_vs_raw.png")

# (C) design-optimization utility curves (if computed).
if (design_ok && !is.null(util)) {
  png(file.path(FIGDIR, "phase3_occupancy_cov_design.png"), width = 1150, height = 850, res = 130)
  op <- par(mar = c(5, 5, 4, 2))
  Ns <- sort(unique(util$N)); pal <- hcl.colors(length(Ns), "Viridis")
  plot(NA, xlim = range(util$K), ylim = range(util$Utility),
       xlab = "replicates per unit (K)", ylab = "expected # species detected / unit",
       main = "occumb design optimization (eval_util_L)")
  for (m in seq_along(Ns)) {
    sub <- util[util$N == Ns[m], ]; sub <- sub[order(sub$K), ]
    lines(sub$K, sub$Utility, col = pal[m], lwd = 2); points(sub$K, sub$Utility, col = pal[m], pch = 19)
  }
  legend("bottomright", legend = sprintf("N = %.0e reads", Ns), col = pal, lwd = 2, pch = 19, bty = "n", cex = 0.85)
  par(op); dev.off()
  message("[06_occumb][cov] wrote figures/phase3_occupancy_cov_design.png")
}

## ---- COVARIATE SELF-CHECK --------------------------------------------------
cat("\n================ COVARIATE SELF-CHECK ================\n")
stopifnot(
  "observed cells == 647"                    = n_observed == 647L,
  "covariate fit is an occumbFit object"     = inherits(fit_cov, "occumbFit"),
  "covariate design matrix has no NA"        = !anyNA(mm_psi),
  "no NA in Dist_from_lake / season join"    = !anyNA(u_dist) && !anyNA(u_seas),
  "corrected psi in [0,1]"                   = all(cvr$corrected_psi >= 0 & cvr$corrected_psi <= 1),
  "raw naive occupancy in [0,1]"             = all(cvr$raw_naive_occ >= 0 & cvr$raw_naive_occ <= 1),
  "gamma effects table has I x eff rows"     = nrow(gamma_tab) == I * length(eff_names),
  "corrected-vs-raw has 14 analysis labels"  = nrow(cvr) == 14L
)
cat("COVARIATE SELF-CHECK PASSED\n")
cat(sprintf("  observed cells         = %d (== 647)\n", n_observed))
cat(sprintf("  J units                = %d ; psi effects = %d [%s]\n",
            J, length(eff_names), paste(eff_names, collapse = ", ")))
cat(sprintf("  corrected psi range    = [%.3f, %.3f]\n", min(cvr$corrected_psi), max(cvr$corrected_psi)))
cat(sprintf("  corrected-vs-raw shift = median %+.3f, max %+.3f (%s)\n",
            median(cvr$shift), max(cvr$shift), cvr$species[which.max(cvr$shift)]))
cat(sprintf("  spatially-structured (Dist_from_lake CrI excl. 0): %s\n",
            if (length(spatial_sp)) paste(spatial_sp, collapse = ", ") else "(none at this MCMC length)"))
cat(sprintf("  %s\n", conv_note_cov))
cat(sprintf("  design recommendation  : %s\n", if (design_ok) design_note else "(not computed)"))
cat("\nPRODUCTION run (4 chains / 60000 iter / 20000 burn / thin 20; ~hours, do NOT run in smoke):\n")
cat("  OCCUMB_SMOKE=0 micromamba run -n edna Rscript R/06_occupancy_occumb.R\n")
cat("  (headline covariate model only, faster: add OCCUMB_SKIP_INTERCEPT=1)\n")
cat("Covariate model done.\n")

## ---- RELEASE CLAIM GATE -----------------------------------------------------
# Rerunning this exploratory fit must not resurrect cached numerical displays
# after the known-truth recovery study has failed its convergence gate.
layerA_status_path <- file.path("data", "metadata", "layerA_recovery_status.tsv")
if (file.exists(layerA_status_path)) {
  layerA_status <- read.delim(layerA_status_path, check.names = FALSE,
                              stringsAsFactors = FALSE)
  if (nrow(layerA_status) == 1L &&
      identical(layerA_status$claim_status[[1]], "WITHHELD_UNCONVERGED")) {
    source(file.path("scripts", "apply_layerA_claim_gate.R"))
    write_layerA_withheld_figures(layerA_status)
  }
}
