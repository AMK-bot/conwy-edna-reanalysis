#!/usr/bin/env Rscript
# =============================================================================
# 05_occupancy_msocc.R  --  Phase 1: single-species 2-level occupancy (msocc).
#
# Plan §2.2 step 1 / §3 Phase 1. Detection histories (site×time unit, 3 water-
# sample replicates) -> per-species occupancy ψ vs detection p. A LEARNING +
# de-risking artifact (not a headline result): build intuition, validate the
# histories that Phase 3's occumb will reuse.
#
# Run inside the env:   micromamba run -n edna Rscript R/05_occupancy_msocc.R
# Re-run is cheap (Gibbs); pass FORCE=1 to ignore the cached .rds.
#
# -----------------------------------------------------------------------------
# THE 2-LEVEL → msocc MAPPING (the crux decision; see "deviations" in eDNAplan.md)
#
# msocc is a 3-level model: site (ψ) → sample (θ) → PCR rep (p). Our design has
# only 2 levels: a site×time occupancy unit, and 3 water-sample replicates (PCR
# triplicates were POOLED before indexing, so each water sample yields ONE call).
#
# We map the 3 water samples to msocc's *sample* level, NOT its *rep* level:
#   site   = site×time unit          -> ψ  (occupancy — the quantity of interest)
#   sample = water-sample replicate  -> θ  (per-replicate detection)
#   rep    = 1 dummy PCR per sample  -> p  (rep~1; collapsed)
# Why: with one call per water sample, putting the 3 reps at the PCR level (1
# sample/unit) makes ψ and θ confounded (only ψ·θ identified) — occupancy, the
# whole point, becomes unidentifiable. Putting them at the sample level gives the
# canonical, identifiable single-season occupancy model y_ij ~ Bern(z_i·θ_ij·p),
# with z_i ~ Bern(ψ_i). ψ is then clean; θ and p are not separately identified
# (1 rep), so we report the identified product θ·p as the per-replicate detection.
#
# NA-NOT-0: we include ONLY the 665 published/QC-passed cells as rows. The 90
# QC-dropped + 43 unsurveyed (incl. all of T16) cells are simply absent — they
# contribute no likelihood (NA), never a false 0. This is validated against
# lib_design.R::observable_mask("published"), the single source of truth.
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(tibble)
  library(msocc); library(ggplot2); library(posterior)
})

# --- repo root (works whether run via Rscript R/05.. or sourced from report/) -
find_root <- function() {
  fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  cands <- c(if (length(fa)) dirname(dirname(normalizePath(fa[1]))),
             normalizePath("."), normalizePath(".."), normalizePath("../.."))
  for (d in cands) if (file.exists(file.path(d, "data/metadata/occasions.tsv"))) return(d)
  stop("repo root not found (need data/metadata/occasions.tsv)")
}
ROOT <- find_root()
.oldwd <- getwd(); setwd(ROOT)   # restored at end-of-script (on.exit misfires at top level of a sourced file)
source("R/lib_design.R")

set.seed(1)
NUM_MCMC <- 11000L; BURNIN <- 1000L      # Gibbs; discard first 1000 draws
MIN_UNITS <- 30L                          # species-inclusion rule (decision record §2)
CACHE <- "data/derived/phase1_msocc.rds"

# --- authors' pre-analysis cleaning (Analysing_BLAST_output_WP2_12S.R ~355–368) -
# Drop the positive-control spike + marine/aquaculture strays; lump the 3 DB
# assignments using the source study's relabel/lump rules. Favonigobius_gymnauchen is on the
# authors' list but absent from our 68-col matrix → its drop is a no-op.
DROP_SPECIES <- c(
  "Gadus_macrocephalus", "Hippoglossus_stenolepis", "Ammodytes_personatus",
  "Clupea_harengus", "Clupea_pallasii", "Copadichromis_virginalis",
  "Favonigobius_gymnauchen", "Pholis_crassispina", "Melanogrammus_aeglefinus",
  "Platichthys_stellatus", "Scomber_scombrus")                 # mackerel = pos. control
RELABEL <- c(Salmo_obtusirostris = "Salmo_trutta",             # DB misassignments → lump
             Anguilla_rostrata = "Anguilla_anguilla",
             Salvelinus_fontinalis_x_Salvelinus_malma = "Salvelinus_malma")

# =============================================================================
# 1. Build cleaned, validated detection histories
# =============================================================================
build_histories <- function() {
  sd  <- read_tsv("data/metadata/source_data_12S_normalised.tsv", show_col_types = FALSE)
  cov <- read_tsv("data/metadata/conwy_covariates.tsv", show_col_types = FALSE)
  lib <- read_tsv("data/metadata/conwy_12S_libraries.tsv", show_col_types = FALSE)

  # -- validate the 665 rows ARE exactly the observable ("published") cells -----
  mask <- observable_mask("published")           # logical [site, tp, rep]
  mdf  <- as.data.frame.table(mask, stringsAsFactors = FALSE)   # cols: site, timepoint, rep, Freq
  obs_ids <- sprintf("%s_%s_%s", mdf$site, mdf$rep, mdf$timepoint)[mdf$Freq]
  stopifnot(setequal(sd$ID, obs_ids))            # GATE: source rows == observable mask
  stopifnot(sum(mask) == nrow(sd))               # 665

  # -- clean species: drop, then OR-merge the relabel pairs (detection = >0) ----
  present_drop <- intersect(DROP_SPECIES, names(sd))
  det <- sd %>% select(-all_of(present_drop))
  spp0 <- setdiff(names(det), "ID")
  det <- det %>% mutate(across(all_of(spp0), ~ as.integer(.x > 0)))   # → 0/1 detection
  for (src in names(RELABEL)) {
    tgt <- RELABEL[[src]]
    if (src %in% names(det)) {
      if (tgt %in% names(det)) det[[tgt]] <- pmax(det[[tgt]], det[[src]]) else det[[tgt]] <- det[[src]]
      det[[src]] <- NULL
    }
  }
  species <- setdiff(names(det), "ID")

  # -- parse ID (E##_R#_T##) → unit; attach covariates; assert completeness -----
  det <- det %>%
    mutate(site = substr(ID, 1, 3),
           rep  = sub("^E[0-9]+_(R[0-9])_.*$", "\\1", ID),
           timepoint = sub("^.*_(T[0-9]+)$", "\\1", ID),
           unit = paste(site, timepoint, sep = "_")) %>%
    left_join(cov %>% select(ID, Dist_from_lake, season), by = "ID") %>%
    left_join(lib %>% select(ID, read_count), by = "ID")
  stopifnot(!anyNA(det$Dist_from_lake), !anyNA(det$season), !anyNA(det$read_count))
  stopifnot(!any(det$timepoint == "T16"))        # T16 never in the published matrix

  # covariates are constant within a unit (ψ level) — verify before collapsing
  chk <- det %>% group_by(unit) %>%
    summarise(nd = n_distinct(Dist_from_lake), ns = n_distinct(season), .groups = "drop")
  stopifnot(all(chk$nd == 1L), all(chk$ns == 1L))

  list(det = det, species = species,
       counts = list(observable = nrow(sd), units = n_distinct(det$unit),
                     units_lt3reps = sum(table(det$unit) < 3L),   # ragged: <3 observed reps
                     dropped_qc = sum(observable_mask("sequenced")) - nrow(sd),  # 90
                     unsurveyed = prod(dim(mask)) - sum(observable_mask("sequenced")), # 43
                     n_species_raw = length(spp0) + length(present_drop),
                     n_species_clean = length(species)))
}

# =============================================================================
# 2. Fit one msocc model for one species
#    model ∈ {null, occ, occ_det}; beta_bin=FALSE so ψ/θ/p come uniformly from
#    the logit coefficients (β/α/δ) for every model.
# =============================================================================
fit_species <- function(det, species, model = "null") {
  uu  <- det %>% distinct(unit, Dist_from_lake)          # z-score the site-level
  dmu <- mean(uu$Dist_from_lake); dsd <- sd(uu$Dist_from_lake)  # covariate per UNIT,
  d <- det %>%                                            # not replicate-weighted
    transmute(site = unit, sample = rep, pcr1 = .data[[species]],
              Dist_z = (Dist_from_lake - dmu) / dsd,
              season = factor(season, levels = c("Spring", "Summer", "Autumn", "Winter")),
              depth_z = as.numeric(scale(log(read_count)))) %>%   # depth IS per-replicate
    arrange(site, sample)

  wide_data <- d %>% select(site, sample, pcr1) %>% as.data.frame()
  site_cov  <- d %>% distinct(site, Dist_z, season) %>% arrange(site) %>% as.data.frame()
  samp_cov  <- d %>% select(site, sample, depth_z) %>% as.data.frame()

  psi_f <- switch(model, null = ~1, occ = ~ Dist_z + season, occ_det = ~ Dist_z + season)
  det_f <- switch(model, null = ~1, occ = ~1,                occ_det = ~ depth_z)

  msocc_mod(
    wide_data,
    site   = list(model = psi_f, cov_tbl = site_cov),
    sample = list(model = det_f, cov_tbl = samp_cov),
    rep    = list(model = ~1,    cov_tbl = samp_cov %>% select(site, sample)),
    num.mcmc = NUM_MCMC, progress = FALSE, seed = 1, beta_bin = FALSE)
}

# posterior ψ (per unit), θ and p (per sample) from the logit coefficients
posterior_probs <- function(mod) {
  keep  <- (BURNIN + 1):nrow(mod$beta)
  X <- mod$model.info$X; W <- mod$model.info$W; V <- mod$model.info$V
  list(psi = plogis(mod$beta[keep, , drop = FALSE]  %*% t(X)),    # draws × units
       det = plogis(mod$alpha[keep, , drop = FALSE] %*% t(W)) *    # θ·p = per-rep detection
             plogis(mod$delta[keep, , drop = FALSE] %*% t(V)))
}
qsum <- function(x) c(mean = mean(x), l95 = unname(quantile(x, .025)), u95 = unname(quantile(x, .975)))

# =============================================================================
# 3. Run: fit the progression for each well-detected analysis label
# =============================================================================
run_phase1 <- function() {
  H <- build_histories()
  det <- H$det

  naive <- vapply(H$species, function(s)
    mean(tapply(det[[s]], det$unit, function(v) as.integer(any(v > 0)))), numeric(1))
  units_det <- vapply(H$species, function(s)
    sum(tapply(det[[s]], det$unit, function(v) any(v > 0))), integer(1))
  selected <- names(sort(units_det[units_det >= MIN_UNITS], decreasing = TRUE))
  message(sprintf("Selected %d analysis labels (≥%d units): %s",
                  length(selected), MIN_UNITS, paste(selected, collapse = ", ")))

  models <- c("null", "occ", "occ_det")
  est <- list(); fits <- list()
  for (sp in selected) for (m in models) {
    mod <- fit_species(det, sp, m)
    pp  <- posterior_probs(mod)
    psi_draw <- rowMeans(pp$psi); det_draw <- rowMeans(pp$det)   # population means per draw
    waic <- tryCatch(as.numeric(msocc::waic(mod, type = 2)), error = function(e) NA_real_)
    # convergence: min effective sample size across the (slowest-mixing) ψ coefficients
    ess  <- tryCatch(min(apply(mod$beta[(BURNIN + 1):nrow(mod$beta), , drop = FALSE],
                               2, posterior::ess_basic)), error = function(e) NA_real_)
    est[[length(est) + 1]] <- tibble(
      species = sp, model = m,
      psi = qsum(psi_draw)["mean"], psi_l = qsum(psi_draw)["l95"], psi_u = qsum(psi_draw)["u95"],
      p   = qsum(det_draw)["mean"], p_l   = qsum(det_draw)["l95"], p_u   = qsum(det_draw)["u95"],
      naive_occ = naive[[sp]], units_detected = units_det[[sp]], waic = waic, min_ess_psi = ess)
    fits[[paste(sp, m, sep = "__")]] <- mod
  }
  estimates <- bind_rows(est)
  # P(an occupied unit is missed) AT THE DESIGN EFFORT of 3 replicates — a forward-looking
  # power statement, not the realized rate (units with <3 observed reps miss more often).
  estimates$fn_rate3 <- (1 - estimates$p)^3

  list(counts = H$counts, species = H$species, selected = selected,
       estimates = estimates, fits = fits)
}

# =============================================================================
# 4. Persist + figure (idempotent; FORCE=1 to recompute)
#    Wrapped so the caller's wd is restored even if a gate/fit errors mid-way
#    (on.exit can't do this at the top level of a sourced file; tryCatch can).
#    `res` is assigned in the enclosing (global) env, so the Quarto note sees it.
# =============================================================================
tryCatch({
res <- NULL
if (file.exists(CACHE) && !nzchar(Sys.getenv("FORCE"))) {
  message("Loading cached fits: ", CACHE, "  (FORCE=1 to recompute)")
  res <- readRDS(CACHE)
} else {
  dir.create("data/derived", showWarnings = FALSE, recursive = TRUE)
  res <- run_phase1()
  saveRDS(res, CACHE)
}

# committed, human-readable estimates (null + covariate models)
write_tsv(res$estimates, "data/metadata/phase1_occupancy_estimates.tsv")

# ψ vs detection (null model) with naive occupancy for contrast
plt <- res$estimates %>% filter(model == "null") %>%
  mutate(species = reorder(species, p))
g <- ggplot(plt, aes(p, species)) +
  geom_errorbarh(aes(xmin = p_l, xmax = p_u), height = .15, colour = "grey40") +
  geom_point(aes(colour = "detection p (theta x p)"), size = 3) +
  geom_point(aes(x = psi, colour = "occupancy psi"), size = 3) +
  geom_errorbarh(aes(xmin = psi_l, xmax = psi_u, x = psi), height = .15, colour = "grey70") +
  geom_point(aes(x = naive_occ, colour = "naive occupancy"), shape = 4, size = 3, stroke = 1.2) +
  scale_colour_manual(NULL, values = c("occupancy psi" = "#1b7837",
                      "detection p (theta x p)" = "#762a83", "naive occupancy" = "black")) +
  labs(x = "probability", y = NULL,
       title = "Phase 1: occupancy vs per-replicate detection (null model)",
       subtitle = "msocc, 6 well-detected analysis labels; occupancy approximates naive occurrence") +
  xlim(0, 1) + theme_minimal(base_size = 12) + theme(legend.position = "top")
dir.create("figures", showWarnings = FALSE)
ggsave("figures/phase1_psi_vs_p.png", g, width = 8, height = 4.5, dpi = 150)

# --- console summary ---------------------------------------------------------
cat("\n================ Phase 1 — detection histories ================\n")
with(res$counts, cat(sprintf(
  "  raw species: %d → cleaned: %d   |   units: %d/252   observable cells: %d\n  NA-not-0: %d QC-dropped + %d unsurveyed cells excluded (never coded 0)\n",
  n_species_raw, n_species_clean, units, observable, dropped_qc, unsurveyed)))
cat("\n================ ψ vs detection (null model) ================\n")
print(res$estimates %>% filter(model == "null") %>%
        transmute(species, naive_occ = round(naive_occ, 2), psi = round(psi, 2),
                  p = round(p, 2), `false-neg(3 reps)` = round(fn_rate3, 2)) %>% as.data.frame())
cat("\nWrote: data/metadata/phase1_occupancy_estimates.tsv  ·  figures/phase1_psi_vs_p.png  ·  ", CACHE, "\n")
}, finally = setwd(.oldwd))   # always restore caller's wd (report/ when sourced by the Quarto note)
