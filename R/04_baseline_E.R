#!/usr/bin/env Rscript
# 04_baseline_E.R  --  Phase 3: Layer E, the compositionally-correct community baseline.
# Plan eDNAplan.md §2.1 + prereg §4.  Runs ONCE, on the reprocessed ASV->species COUNTS
# (Phase-2 handback), NOT on the published normalised table (see §2.1: the DA trio
# requires counts). The planned abundance-mode iNEXT route was rejected because reads are
# not ecological individuals; the executed iNEXT analysis uses incidence frequencies.
#
# Inputs  (repo-root relative; run:  micromamba run -n edna Rscript R/04_baseline_E.R):
#   data/derived/species_x_sample_counts_clean.rds   647 x 14  (double; cast to integer)
#   data/derived/genus_x_sample_counts_clean.rds     647 x 12
#   data/metadata/conwy_covariates.tsv               keyed by ID = E##_R#_T##
#   R/lib_design.R                                    shared design (survivor mask sanity-check)
#
# Outputs:
#   data/derived/layerE_*.tsv     alpha, inference design, iNEXT, PERMANOVA/PERMDISP,
#                                site-level DA screen, indicators
#   figures/phase3_layerE_*.png   iNEXT curves, ordinations, betadisper, DA/indicator summary
#
# What runs here (§2.1):
#   * Alpha : per-library richness/Shannon/Simpson/Pielou/Hill(q=0,1,2), explicitly
#             DESCRIPTIVE, plus the same metrics on the biological site x time unit. Coverage-
#             based iNEXT uses site x time incidence (not replicate libraries) on section- and
#             season-pooled assemblages -- NOT naive rarefying.
#             No phylogeny for a 12S fish community -> Faith's PD is
#             SKIPPED by design (no tree; a 12S rRNA gene tree over 14 fish would be a
#             taxonomy proxy, not an independent phylogenetic-diversity axis).
#   * Beta  : Bray-Curtis (relative counts), Jaccard (presence/absence), robust Aitchison
#             (vegan::decostand "rclr" -> Euclidean). Replicate libraries are first summed to
#             site x time. Spatial inference then uses 14 independent SITE profiles over the 11
#             timepoints observed at every site; seasonal inference uses 18 independent TIMEPOINT
#             profiles over the fixed 10-site panel observed at every timepoint.
#   * Tests : adonis2 PERMANOVA paired with betadisper/PERMDISP. Because sites lie in river order
#             and timepoints in chronological order, p-values use the complete set of cyclic shifts
#             plus mirrored series permutations (27 spatial / 35 temporal), not 629 unrestricted
#             replicate-library permutations.
#   * DA    : exploratory site-level consensus across ALDEx2 + ANCOM-BC2 + corncob, gated by a
#             spatial-series restricted CLR permutation q-value (see CONSENSUS RULE below).
#   * Indic : indicspecies::multipatt on site-level occurrence frequency with the same restricted
#             spatial-series permutations.
#
# Framing: reads are NOT abundance for eDNA (§2.1). Everything defaults to occurrence /
# compositional framing; Bray-Curtis is computed on proportions and rclr is the headline
# beta metric; DA is reported as a compositional-shift consensus, not fold-change truth.
#
# ---------------------------------------------------------------------------------------
# CONSENSUS RULE (prereg §4 left this TODO; this credibility correction makes it explicit):
#   A taxon is a supported DA-consensus hit for river section iff its OMNIBUS (across-section)
#   association is significant at FDR/q < 0.05 in >= 2 of the 3 count methods AND its
#   site-level CLR association has series-permutation BH q < 0.05, where the omnibus
#   test per method is:
#       ALDEx2   : aldex(test="kw") one-way GLM on clr, glm.eBH  < 0.05
#       ANCOM-BC2: ancombc2(global=TRUE) group-level global test, q_val < 0.05
#       corncob  : differentialTest LRT full(~section) vs null(~1) on abundance, p_fdr < 0.05
#   All tests use the SAME 14 site aggregates and prevalence-filtered taxon set (the locked
#   >=20-library rule -> 8 testable species). Counts are summed only over the 11 timepoints
#   observed at every site, eliminating unequal temporal coverage. With 4-5 sites per section
#   this remains an EXPLORATORY screen; no method-specific q-value alone is a confirmatory claim.
# ---------------------------------------------------------------------------------------
#
# Reading: iNEXT vignette (Hsieh et al. 2016) - Gloor et al. 2017 (compositional) -
#          Anderson 2001 (PERMANOVA) / Anderson 2006 (PERMDISP) - De Caceres & Legendre 2009.

suppressPackageStartupMessages({
  library(vegan)
  library(iNEXT)
  library(indicspecies)
  library(ALDEx2)
  library(ANCOMBC)
  library(corncob)
  library(phyloseq)
})

set.seed(20260705)
options(stringsAsFactors = FALSE)

## ----- paths & idempotent output dirs -------------------------------------------------
SP_RDS   <- "data/derived/species_x_sample_counts_clean.rds"
GE_RDS   <- "data/derived/genus_x_sample_counts_clean.rds"
COV_TSV  <- "data/metadata/conwy_covariates.tsv"
FIG_DIR  <- "figures"
OUT_DIR  <- "data/derived"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot(file.exists(SP_RDS), file.exists(GE_RDS), file.exists(COV_TSV))

# Shared design (single source of truth). Used only to SANITY-CHECK that the 647 count-table
# libraries are exactly the Phase-2 reprocessing survivors -- Layer E itself is intercept-free
# community ecology and needs no occasion grid, but this guards a stale-input swap.
source("R/lib_design.R")

## ----- load counts --------------------------------------------------------------------
sp  <- readRDS(SP_RDS)                      # 647 x 14, double
ge  <- readRDS(GE_RDS)                      # 647 x 12, double
stopifnot(is.matrix(sp), nrow(sp) == 647, ncol(sp) == 14)
spc <- round(sp); storage.mode(spc) <- "integer"     # cast to integer counts
gec <- round(ge); storage.mode(gec) <- "integer"

surv <- survived_reprocessing()
if (!setequal(rownames(spc), surv))
  warning("species count-table rownames != reprocessing survivor set from lib_design.R")

## ----- metadata + biological-unit aggregation -----------------------------------------
cov <- read.delim(COV_TSV)
stopifnot("ID" %in% colnames(cov), all(rownames(spc) %in% cov$ID))
md <- cov[match(rownames(spc), cov$ID), ]
md$site  <- sub("_.*", "", md$ID)
md$timepoint <- sub("^.*_(T[0-9]+)$", "\\1", md$ID)
md$rep <- sub("^E[0-9]+_(R[0-9])_.*$", "\\1", md$ID)
md$unit <- paste(md$site, md$timepoint, sep = "_")
# River SECTION per plan §2.1: upper E01-E05 / middle E06-E10 / lower E11-E14.
md$section <- factor(ifelse(md$site %in% sprintf("E%02d", 1:5),  "upper",
                     ifelse(md$site %in% sprintf("E%02d", 6:10), "middle", "lower")),
                     levels = c("upper", "middle", "lower"))
# Cross-check against the covariate file's own river_section (must agree exactly).
if (!identical(as.character(md$section), md$river_section))
  warning("derived section disagrees with covariate river_section column")
md$season <- factor(md$season, levels = c("Spring", "Summer", "Autumn", "Winter"))
md$depth  <- rowSums(spc)                   # curated-fish reads per replicate library
stopifnot(nrow(md) == nrow(spc), !anyNA(md$section), !anyNA(md$season))

# Replicate libraries are technical/detection replicates, not independent ecological samples.
# Sum them to the site x timepoint biological unit before any inferential analysis. Ragged units
# retain 1-3 surviving libraries; n_libraries is recorded so the aggregation is auditable.
sp_unit <- rowsum(spc, md$unit, reorder = TRUE)
u_first <- match(rownames(sp_unit), md$unit)
umd <- data.frame(
  unit = rownames(sp_unit), site = md$site[u_first], timepoint = md$timepoint[u_first],
  section = factor(md$section[u_first], levels = levels(md$section)),
  season = factor(md$season[u_first], levels = levels(md$season)),
  n_libraries = as.integer(table(md$unit)[rownames(sp_unit)]),
  depth = rowSums(sp_unit), stringsAsFactors = FALSE)
stopifnot(
  nrow(sp_unit) == 234L,
  all(vapply(split(as.character(md$section), md$unit), function(x) length(unique(x)) == 1L, logical(1))),
  all(vapply(split(as.character(md$season), md$unit),  function(x) length(unique(x)) == 1L, logical(1))),
  all(umd$n_libraries %in% 1:3),
  identical(rownames(sp_unit), umd$unit)
)

# Balanced panels for the two crossed repeated axes. Spatial inference uses every site but only
# the 11 timepoints represented at all 14 sites. Seasonal inference uses every timepoint but only
# the fixed 10-site panel represented at all 18 timepoints. This avoids filling missing units with
# zeros and avoids giving sites/timepoints with more surviving libraries extra inferential weight.
unit_grid <- xtabs(~ site + timepoint, data = umd) > 0
common_timepoints <- colnames(unit_grid)[colSums(unit_grid) == nrow(unit_grid)]
common_sites <- rownames(unit_grid)[rowSums(unit_grid) == ncol(unit_grid)]
stopifnot(length(common_timepoints) == 11L, length(common_sites) == 10L)

spatial_idx <- umd$timepoint %in% common_timepoints
sp_site <- rowsum(sp_unit[spatial_idx, , drop = FALSE], umd$site[spatial_idx], reorder = TRUE)
site_first <- match(rownames(sp_site), umd$site)
site_md <- data.frame(
  site = rownames(sp_site),
  section = factor(umd$section[site_first], levels = levels(umd$section)),
  n_timepoints = as.integer(table(umd$site[spatial_idx])[rownames(sp_site)]),
  stringsAsFactors = FALSE)

temporal_idx <- umd$site %in% common_sites
sp_time <- rowsum(sp_unit[temporal_idx, , drop = FALSE], umd$timepoint[temporal_idx], reorder = TRUE)
time_first <- match(rownames(sp_time), umd$timepoint)
time_md <- data.frame(
  timepoint = rownames(sp_time),
  season = factor(umd$season[time_first], levels = levels(umd$season)),
  n_sites = as.integer(table(umd$timepoint[temporal_idx])[rownames(sp_time)]),
  stringsAsFactors = FALSE)
stopifnot(nrow(sp_site) == 14L, all(site_md$n_timepoints == 11L),
          nrow(sp_time) == 18L, all(time_md$n_sites == 10L))

# Machine-readable methods receipt. Downstream prose can quote this instead of inferring the
# sampling unit from table row counts.
inference_design <- data.frame(
  analysis = c("alpha_library", "alpha_iNEXT", "section_beta_dispersion", "season_beta_dispersion",
               "section_DA", "section_indicator"),
  analysis_unit = c("water-sample library", "site x timepoint", "site", "timepoint", "site", "site"),
  n_units = c(nrow(spc), nrow(sp_unit), nrow(sp_site), nrow(sp_time), nrow(sp_site), nrow(sp_site)),
  balanced_panel = c(
    "all surviving libraries; descriptive only",
    "234 realized units; 1-3 libraries/unit; coverage standardized",
    sprintf("14 sites x %d common timepoints", length(common_timepoints)),
    sprintf("18 timepoints x %d common sites (%s)", length(common_sites), paste(common_sites, collapse = ",")),
    sprintf("14 sites x %d common timepoints", length(common_timepoints)),
    sprintf("site occurrence frequency over %d common timepoints", length(common_timepoints))),
  permutation = c("none", "none", "27 cyclic/mirrored river-order shifts",
                  "35 cyclic/mirrored chronological shifts",
                  "27 cyclic/mirrored river-order shifts gates method consensus",
                  "27 cyclic/mirrored river-order shifts"),
  scope = c("descriptive", "descriptive diversity standardization", "spatial inference",
            "seasonal inference on fixed spatial panel", "exploratory; 4-5 sites/section",
            "exploratory; 4-5 sites/section"),
  stringsAsFactors = FALSE)
write.table(inference_design, file.path(OUT_DIR, "layerE_inference_design.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf(paste0("Loaded %d replicate libraries (%d all-target-zero) -> %d site x time units ",
                   "(%d all-target-zero). Spatial panel: 14 sites x %d timepoints; ",
                   "seasonal panel: 18 timepoints x %d fixed sites.\n"),
            nrow(spc), sum(md$depth == 0), nrow(sp_unit), sum(umd$depth == 0),
            length(common_timepoints), length(common_sites)))

## =====================================================================================
## 1. ALPHA DIVERSITY
## =====================================================================================
# (a) Per-library OBSERVED metrics. These are explicitly descriptive: replicate libraries are
# useful for showing detection variability, but are not independent ecological observations.
rich   <- rowSums(sp > 0)
shan   <- diversity(sp, index = "shannon")
simp   <- diversity(sp, index = "simpson")
invsimp<- diversity(sp, index = "invsimpson")
zero_rich <- rich == 0
shan[zero_rich] <- simp[zero_rich] <- invsimp[zero_rich] <- NA_real_
pielou <- ifelse(rich > 1, shan / log(rich), NA_real_)
alpha_ps <- data.frame(
  ID = rownames(sp), unit = md$unit, site = md$site, timepoint = md$timepoint,
  section = md$section, season = md$season,
  analysis_scope = "descriptive water-sample library; not an independent inferential unit",
  depth = md$depth, richness_q0 = rich,
  shannon = shan, simpson = simp, invsimpson = invsimp, pielou = pielou,
  hill_q0 = rich, hill_q1 = exp(shan), hill_q2 = invsimp)   # Hill: q0=S, q1=exp(H'), q2=1/D
write.table(alpha_ps, file.path(OUT_DIR, "layerE_alpha_persample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Site x timepoint alpha: the ecological unit used downstream. Counts are summed across the 1-3
# surviving replicate libraries; metrics remain descriptive because depth/effort is ragged.
u_rich <- rowSums(sp_unit > 0)
u_shan <- diversity(sp_unit, index = "shannon")
u_simp <- diversity(sp_unit, index = "simpson")
u_inv  <- diversity(sp_unit, index = "invsimpson")
u_zero_rich <- u_rich == 0
u_shan[u_zero_rich] <- u_simp[u_zero_rich] <- u_inv[u_zero_rich] <- NA_real_
u_pielou <- ifelse(u_rich > 1, u_shan / log(u_rich), NA_real_)
alpha_unit <- data.frame(
  umd[, c("unit", "site", "timepoint", "section", "season", "n_libraries", "depth")],
  analysis_scope = "descriptive site x timepoint biological unit",
  richness_q0 = u_rich, shannon = u_shan, simpson = u_simp, invsimpson = u_inv,
  pielou = u_pielou, hill_q0 = u_rich, hill_q1 = exp(u_shan), hill_q2 = u_inv,
  row.names = NULL, check.names = FALSE)
write.table(alpha_unit, file.path(OUT_DIR, "layerE_alpha_site_time.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# (b) COVERAGE-BASED rarefaction/extrapolation with iNEXT (the plan's requirement -- NOT
#     naive rarefying). Datatype = INCIDENCE_FREQ, deliberately, not abundance:
#       (i)  CONCEPTUAL -- an eDNA read is NOT an ecological individual (PCR amplification,
#            biomass/primer bias); the occurrence framing (reads != abundance, §2.1) makes
#            the sampling unit the SITE x TIMEPOINT biological unit and the datum species
#            incidence across units within a group. Replicate libraries are OR-combined.
#       (ii) NUMERICAL -- abundance-mode treats each of the ~34 M reads as an individual, so
#            iNEXT's reference sample size (and its 2x extrapolation endpoint) balloon to tens
#            of millions and the curve computation does not terminate. Incidence caps effort at
#            the number of realized site x time units -> runs in seconds.
#   Assemblage = river section (and season); incidence-frequency vector = c(T, per-species #
#   of units detected in), zero-incidence species dropped per assemblage. ALL realized units
#   count toward T, including the 4 all-target-zero units (they add no incidence).
inc_by <- function(pa, grp) {
  grp <- factor(grp)
  lst <- lapply(levels(grp), function(g) {
    incid <- colSums(pa[grp == g, , drop = FALSE])
    c(sum(grp == g), sort(incid[incid > 0], decreasing = TRUE))   # T, then incidences>0
  })
  setNames(lst, levels(grp))
}
pa_unit <- (sp_unit > 0) + 0
sec_pool <- inc_by(pa_unit, umd$section)
sea_pool <- inc_by(pa_unit, umd$season)

run_inext <- function(pool, tag) {
  out <- iNEXT(pool, q = c(0, 1, 2), datatype = "incidence_freq", nboot = 50, se = TRUE)
  asy <- out$AsyEst                          # asymptotic Hill estimates per assemblage x q
  write.table(asy, file.path(OUT_DIR, sprintf("layerE_inext_asymptotic_%s.tsv", tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  # Coverage-standardized comparison at the largest two-decimal coverage not exceeding the
  # least-complete reference assemblage (98% for both current pools). Supplying the level is
  # important: estimateD's implicit endpoint is numerically undefined when an assemblage has
  # estimated reference coverage exactly 1. A hard failure prevents a stale file surviving.
  info <- DataInfo(pool, datatype = "incidence_freq")
  target_cov <- min(0.99, floor(min(info$SC) * 100) / 100)
  stopifnot(is.finite(target_cov), target_cov > 0, target_cov < 1)
  # Reset the bootstrap stream because some simulated one-species incidence resamples are
  # numerically degenerate in iNEXT; this fixed stream is reproducible and verified to be valid.
  set.seed(20260705)
  covD <- estimateD(pool, q = c(0, 1, 2), datatype = "incidence_freq",
                    base = "coverage", level = target_cov, nboot = 50, conf = 0.95)
  covD$target_SC_requested <- target_cov
  write.table(covD, file.path(OUT_DIR, sprintf("layerE_inext_coverage_%s.tsv", tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  out
}
inext_sec <- run_inext(sec_pool, "section")
inext_sea <- run_inext(sea_pool, "season")

# iNEXT figures (size-based curves, faceted by Hill order q).
save_inext_fig <- function(out, file, title) {
  g <- ggiNEXT(out, type = 1, facet.var = "Order.q") +
       ggplot2::labs(title = title) + ggplot2::theme_bw(base_size = 12)
  ggplot2::ggsave(file, g, width = 10, height = 4, dpi = 150)
}
save_inext_fig(inext_sec, file.path(FIG_DIR, "phase3_layerE_inext_section.png"),
               "Coverage-based Hill diversity by river section (site x time incidence)")
save_inext_fig(inext_sea, file.path(FIG_DIR, "phase3_layerE_inext_season.png"),
               "Coverage-based Hill diversity by season (site x time incidence)")

## =====================================================================================
## 2. BETA DIVERSITY -- independent spatial and temporal profiles
## =====================================================================================
distance_set <- function(counts) {
  stopifnot(is.matrix(counts), all(rowSums(counts) > 0))
  rel  <- decostand(counts, method = "total")
  rclr <- decostand(counts, method = "rclr")
  list(
    "Bray-Curtis (rel. abund.)" = vegdist(rel, method = "bray"),
    "Jaccard (presence/absence)" = vegdist(counts, method = "jaccard", binary = TRUE),
    "robust Aitchison (rclr)" = dist(rclr, method = "euclidean"))
}
site_dists <- distance_set(sp_site)   # 14 sites; same 11 timepoints at each site
time_dists <- distance_set(sp_time)   # 18 timepoints; same 10 sites at each timepoint

# Complete series permutations preserve river adjacency / chronological adjacency. For n ordered
# units there are 2n-1 non-identity cyclic+mirrored permutations: 27 sites, 35 timepoints.
series_control <- function(nperm = 999L) {
  permute::how(nperm = nperm,
               within = permute::Within(type = "series", mirror = TRUE))
}
site_perm_control <- series_control()
time_perm_control <- series_control()
site_perm <- permute::shuffleSet(nrow(sp_site), nset = 999, control = site_perm_control)
time_perm <- permute::shuffleSet(nrow(sp_time), nset = 999, control = time_perm_control)
stopifnot(nrow(site_perm) == 27L, nrow(time_perm) == 35L)

## =====================================================================================
## 3. ORDINATION + TESTS: PERMANOVA paired with PERMDISP on the same independent units
## =====================================================================================
test_specs <- list(
  section = list(dists = site_dists, group = site_md$section, unit = "site",
                 n = nrow(sp_site), panel = sprintf("14 sites x %d common timepoints", length(common_timepoints)),
                 perm = site_perm, scheme = "cyclic+mirrored river-order series"),
  season = list(dists = time_dists, group = time_md$season, unit = "timepoint",
                n = nrow(sp_time), panel = sprintf("18 timepoints x %d common sites", length(common_sites)),
                perm = time_perm, scheme = "cyclic+mirrored chronological series"))

perma_rows <- list(); disp_rows <- list()
for (gv in names(test_specs)) {
  z <- test_specs[[gv]]
  for (dn in names(z$dists)) {
    D <- z$dists[[dn]]; g <- z$group
    a  <- adonis2(D ~ g, permutations = z$perm, by = NULL)
    bd <- betadisper(D, g)
    pt <- permutest(bd, permutations = z$perm)
    key <- paste(dn, gv)
    perma_rows[[key]] <- data.frame(
      distance = dn, grouping = gv, term = "location", analysis_unit = z$unit,
      n_units = z$n, balanced_panel = z$panel, permutation_scheme = z$scheme,
      n_permutations = nrow(z$perm), Df = a$Df[1], R2 = round(a$R2[1], 4),
      F = round(a$F[1], 3), p_value = a$`Pr(>F)`[1], stringsAsFactors = FALSE)
    disp_rows[[key]] <- data.frame(
      distance = dn, grouping = gv, term = "dispersion", analysis_unit = z$unit,
      n_units = z$n, balanced_panel = z$panel, permutation_scheme = z$scheme,
      n_permutations = nrow(z$perm), Df = pt$tab$Df[1], F = round(pt$tab$F[1], 3),
      p_value = pt$tab$`Pr(>F)`[1], stringsAsFactors = FALSE)
  }
}
permanova <- do.call(rbind, perma_rows); rownames(permanova) <- NULL
permdisp  <- do.call(rbind, disp_rows);  rownames(permdisp)  <- NULL
write.table(permanova, file.path(OUT_DIR, "layerE_permanova.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(permdisp, file.path(OUT_DIR, "layerE_betadisper.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Inference-aligned ordinations: spatial panels show sites; seasonal panels show timepoints.
col_sec <- c(upper = "#1b9e77", middle = "#d95f02", lower = "#7570b3")
col_sea <- c(Spring = "#4daf4a", Summer = "#e41a1c", Autumn = "#ff7f00", Winter = "#377eb8")
plot_ord <- function(pts, grp, pal, main, xlab, ylab, cex = 0.9) {
  plot(pts, col = pal[as.character(grp)], pch = 19, cex = cex,
       main = main, xlab = xlab, ylab = ylab)
  ordiellipse(pts, grp, kind = "sd", col = pal, lwd = 2, label = FALSE)
  legend("topright", legend = names(pal), col = pal, pch = 19, bty = "n", cex = 0.8)
}
site_pcoa <- cmdscale(site_dists[["robust Aitchison (rclr)"]], k = 2, eig = TRUE)
time_pcoa <- cmdscale(time_dists[["robust Aitchison (rclr)"]], k = 2, eig = TRUE)
site_nmds <- tryCatch(metaMDS(site_dists[["Bray-Curtis (rel. abund.)"]], k = 2, trymax = 100, trace = 0),
                      error = function(e) NULL)
time_nmds <- tryCatch(metaMDS(time_dists[["Bray-Curtis (rel. abund.)"]], k = 2, trymax = 100, trace = 0),
                      error = function(e) NULL)

png(file.path(FIG_DIR, "phase3_layerE_ordination.png"), width = 1400, height = 1200, res = 150)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
sv <- round(100 * site_pcoa$eig[1:2] / sum(site_pcoa$eig[site_pcoa$eig > 0]), 1)
tv <- round(100 * time_pcoa$eig[1:2] / sum(time_pcoa$eig[time_pcoa$eig > 0]), 1)
plot_ord(site_pcoa$points, site_md$section, col_sec,
         "Spatial PCoA: 14 sites (11 common timepoints)",
         sprintf("PCo1 (%.1f%%)", sv[1]), sprintf("PCo2 (%.1f%%)", sv[2]))
plot_ord(time_pcoa$points, time_md$season, col_sea,
         "Seasonal PCoA: 18 timepoints (fixed 10-site panel)",
         sprintf("PCo1 (%.1f%%)", tv[1]), sprintf("PCo2 (%.1f%%)", tv[2]))
if (!is.null(site_nmds))
  plot_ord(site_nmds$points, site_md$section, col_sec,
           sprintf("Spatial NMDS Bray-Curtis (stress %.3f)", site_nmds$stress), "NMDS1", "NMDS2") else plot.new()
if (!is.null(time_nmds))
  plot_ord(time_nmds$points, time_md$season, col_sea,
           sprintf("Seasonal NMDS Bray-Curtis (stress %.3f)", time_nmds$stress), "NMDS1", "NMDS2") else plot.new()
par(op); dev.off()

# Preserve the replicate-library ordination as an explicitly DESCRIPTIVE diagnostic. It has no
# p-values and is not used for ellipses, inference, or claims.
lib_nz <- md$depth > 0
lib_rclr <- decostand(spc[lib_nz, , drop = FALSE], method = "rclr")
lib_pcoa <- cmdscale(dist(lib_rclr), k = 2, eig = TRUE)
png(file.path(FIG_DIR, "phase3_layerE_ordination_library_descriptive.png"),
    width = 1200, height = 600, res = 150)
op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(lib_pcoa$points, col = col_sec[as.character(md$section[lib_nz])], pch = 19, cex = 0.45,
     xlab = "PCo1", ylab = "PCo2", cex.main = 0.9,
     main = "Replicate libraries by section\n(descriptive only; no inference)")
legend("topright", names(col_sec), col = col_sec, pch = 19, bty = "n", cex = 0.8)
plot(lib_pcoa$points, col = col_sea[as.character(md$season[lib_nz])], pch = 19, cex = 0.45,
     xlab = "PCo1", ylab = "PCo2", cex.main = 0.9,
     main = "Replicate libraries by season\n(descriptive only; no inference)")
legend("topright", names(col_sea), col = col_sea, pch = 19, bty = "n", cex = 0.8)
par(op); dev.off()

# PERMDISP distance-to-centroid plots at the independent site/timepoint units.
site_bd_ai <- betadisper(site_dists[["robust Aitchison (rclr)"]], site_md$section)
time_bd_ai <- betadisper(time_dists[["robust Aitchison (rclr)"]], time_md$season)
png(file.path(FIG_DIR, "phase3_layerE_betadisper.png"), width = 1200, height = 600, res = 150)
op <- par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
boxplot(site_bd_ai$distances ~ site_md$section, col = col_sec,
        main = "Spatial PERMDISP: sites", xlab = "", ylab = "distance to centroid")
boxplot(time_bd_ai$distances ~ time_md$season, col = col_sea,
        main = "Seasonal PERMDISP: timepoints", xlab = "", ylab = "distance to centroid")
par(op); dev.off()

## =====================================================================================
## 4. DIFFERENTIAL ABUNDANCE -- exploratory SITE-level screen by section
## =====================================================================================
# Keep the documented >=20-library prevalence filter, but perform every test on 14 SITE
# aggregates over the same 11 complete timepoints. This removes technical-replicate and repeated-
# time pseudoreplication. The three package tests do not natively honour river-order dependence,
# so their >=2/3 consensus is additionally gated by a CLR omnibus test using the same 27 restricted
# river-series permutations as PERMANOVA. With only 4-5 sites per section, call this exploratory.
keep_taxa <- names(which(colSums(spc > 0) >= 20))
cat(sprintf("Exploratory site-level DA on %d prevalence-filtered species (>=20 libraries): %s\n",
            length(keep_taxa), paste(keep_taxa, collapse = ", ")))
da_counts <- sp_site[, keep_taxa, drop = FALSE]           # 14 sites x taxa; balanced time panel
da_grp    <- droplevels(site_md$section)
feat_mat  <- t(round(da_counts)); storage.mode(feat_mat) <- "integer"   # taxa x 14 sites

# (a) ALDEx2 -- one-way KW/GLM on CLR across the 14 site aggregates.
aldex_q <- setNames(rep(NA_real_, length(keep_taxa)), keep_taxa)
try({
  ax <- aldex(reads = feat_mat, conditions = as.character(da_grp),
              mc.samples = 128, test = "kw", denom = "all", verbose = FALSE)
  aldex_q[rownames(ax)] <- ax$glm.eBH
}, silent = FALSE)

# (b) ANCOM-BC2 -- group global test across sections on the 14 site aggregates.
#   Fed as matrix + meta_data (NOT phyloseq: ANCOMBC's phyloseq path needs the uninstalled
#   `microbiome` pkg; the group column must be CHARACTER or ANCOMBC's internal reassignment
#   corrupts the factor). pseudo_sens=TRUE runs the recommended sensitivity analysis.
#   CAVEAT: ANCOM-BC's bias estimation wants >~50 taxa; with 8 taxa and 14 sites it warns of
#   inconsistency. This is one reason the result remains exploratory and is permutation-gated.
ancom_q  <- setNames(rep(NA_real_, length(keep_taxa)), keep_taxa)
ancom_ss <- setNames(rep(NA,        length(keep_taxa)), keep_taxa)
try({
  meta <- data.frame(sample_id = colnames(feat_mat),
                     river_section = as.character(da_grp), stringsAsFactors = FALSE)
  rownames(meta) <- meta$sample_id
  ab <- ancombc2(data = feat_mat, taxa_are_rows = TRUE, meta_data = meta,
                 fix_formula = "river_section", group = "river_section",
                 global = TRUE, prv_cut = 0, lib_cut = 0, alpha = 0.05,
                 pseudo_sens = TRUE, verbose = FALSE, n_cl = 1)
  g <- ab$res_global
  ancom_q[g$taxon] <- g$q_val
  if ("passed_ss" %in% names(g)) ancom_ss[g$taxon] <- g$passed_ss
}, silent = FALSE)

# (c) corncob -- LRT full(~section) vs null(~1) on the 14 site aggregates.
corncob_q <- setNames(rep(NA_real_, length(keep_taxa)), keep_taxa)
try({
  ps2 <- phyloseq(otu_table(feat_mat, taxa_are_rows = TRUE),
                  sample_data(data.frame(river_section = da_grp,
                                         row.names = colnames(feat_mat))))
  cc <- differentialTest(formula = ~ river_section, phi.formula = ~ river_section,
                         formula_null = ~ 1, phi.formula_null = ~ river_section,
                         test = "LRT", boot = FALSE, data = ps2, fdr = "BH")
  corncob_q[names(cc$p_fdr)] <- cc$p_fdr
}, silent = FALSE)

# (d) Primary credibility gate: CLR omnibus F with complete river-series permutations.
# A 0.5 pseudocount is used only for this gate; the three count methods above retain their native
# zero handling. `site_perm` contains every non-identity cyclic/mirrored ordering of the 14 sites.
site_clr <- log(da_counts + 0.5)
site_clr <- site_clr - rowMeans(site_clr)
omnibus_f <- function(x, g) {
  a <- anova(lm(x ~ g))
  as.numeric(a$`F value`[1])
}
series_p <- vapply(seq_along(keep_taxa), function(j) {
  obs <- omnibus_f(site_clr[, j], da_grp)
  pf <- apply(site_perm, 1, function(idx) omnibus_f(site_clr[idx, j], da_grp))
  (1 + sum(pf >= obs - sqrt(.Machine$double.eps), na.rm = TRUE)) / (1 + length(pf))
}, numeric(1))
names(series_p) <- keep_taxa
series_q <- p.adjust(series_p, method = "BH")

sig <- function(q) !is.na(q) & q < 0.05
da <- data.frame(
  species            = keep_taxa,
  n_detect           = colSums(spc[, keep_taxa, drop = FALSE] > 0), # legacy: libraries detected
  n_site_time_detect = colSums(sp_unit[spatial_idx, keep_taxa, drop = FALSE] > 0),
  n_sites_detected    = colSums(sp_site[, keep_taxa, drop = FALSE] > 0),
  aldex_glm_eBH      = round(aldex_q[keep_taxa], 4),
  ancombc2_qval      = round(ancom_q[keep_taxa], 4),
  ancombc2_passed_ss = ancom_ss[keep_taxa],
  corncob_p_fdr      = round(corncob_q[keep_taxa], 4),
  river_series_p     = round(series_p[keep_taxa], 4),
  river_series_q     = round(series_q[keep_taxa], 4),
  analysis_unit      = "site",
  balanced_panel     = sprintf("14 sites x %d common timepoints", length(common_timepoints)),
  stringsAsFactors = FALSE)
da$n_methods_sig <- sig(da$aldex_glm_eBH) + sig(da$ancombc2_qval) + sig(da$corncob_p_fdr)
da$method_consensus <- da$n_methods_sig >= 2
da$consensus_DA <- da$method_consensus & sig(da$river_series_q)
da$evidence_tier <- ifelse(da$consensus_DA, "restricted-series-supported exploratory hit",
                    ifelse(da$method_consensus, "method consensus only; not spatially supported",
                           "no >=2-method consensus"))
da <- da[order(!da$consensus_DA, -da$n_methods_sig, da$river_series_q, da$species), ]
write.table(da, file.path(OUT_DIR, "layerE_da_consensus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## =====================================================================================
## 5. INDICATOR SPECIES -- SITE-level occurrence frequency by river section
## =====================================================================================
# Replicate calls are OR-combined within site x timepoint; site profiles are the fraction of the
# same 11 complete timepoints in which each species was detected. IndVal.g therefore has 14
# independent site rows and the same river-series restricted permutation design as spatial beta.
site_occurrence <- rowsum((sp_unit[spatial_idx, , drop = FALSE] > 0) + 0,
                          umd$site[spatial_idx], reorder = TRUE) / length(common_timepoints)
stopifnot(identical(rownames(site_occurrence), site_md$site))
mp <- multipatt(as.data.frame(site_occurrence), site_md$section, func = "IndVal.g",
                control = site_perm_control)
ms <- mp$sign
ms$species <- rownames(ms)
grp_cols <- grep("^s\\.", colnames(ms), value = TRUE)
# A ubiquitous species (present across ALL sections, e.g. Salmo_trutta) belongs to the
# all-groups combination, which multipatt does NOT test -> stat/p = NA; label it explicitly
# rather than leaving a bare NA row.
combo <- apply(ms[, grp_cols, drop = FALSE], 1, function(r)
  paste(sub("^s\\.", "", grp_cols)[r == 1], collapse = "+"))
combo[is.na(ms$stat) | combo == ""] <- "none (ubiquitous/untested)"
indic <- data.frame(species = ms$species, section_group = combo,
                    IndVal_stat = round(ms$stat, 3), p_value = ms$p.value,
                    n_sites_detected = colSums(sp_site[, ms$species, drop = FALSE] > 0),
                    analysis_unit = "site occurrence frequency",
                    balanced_panel = sprintf("14 sites x %d common timepoints", length(common_timepoints)),
                    permutation_scheme = "27 cyclic+mirrored river-order series permutations",
                    stringsAsFactors = FALSE)
indic <- indic[order(indic$p_value, -indic$IndVal_stat), ]
write.table(indic, file.path(OUT_DIR, "layerE_indicator_species.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# DA + indicator summary figure. The fourth DA bar is the credibility gate; a package-method
# consensus without river-series q<0.05 is explicitly not a supported hit.
png(file.path(FIG_DIR, "phase3_layerE_da_indicator.png"), width = 1300, height = 650, res = 150)
op <- par(mfrow = c(1, 2), mar = c(9, 4, 3, 1))
qmat <- t(as.matrix(da[, c("aldex_glm_eBH", "ancombc2_qval", "corncob_p_fdr", "river_series_q")]))
qmat[is.na(qmat)] <- 1
barplot(-log10(qmat + 1e-6), beside = TRUE, las = 2, names.arg = da$species,
        col = c("#66c2a5", "#fc8d62", "#8da0cb", "#1f1f1f"), ylab = "-log10(FDR q)",
        main = "Site-level DA screen\n(restricted-series q is required)", cex.main = 0.9)
abline(h = -log10(0.05), lty = 2)
legend("topright", c("ALDEx2", "ANCOM-BC2", "corncob", "river-series permutation"),
       fill = c("#66c2a5", "#fc8d62", "#8da0cb", "#1f1f1f"), bty = "n", cex = 0.72)
sig_ind <- indic[!is.na(indic$p_value) & indic$p_value <= 0.05 &
                   indic$section_group != "", ]
if (nrow(sig_ind) > 0) {
  barplot(sig_ind$IndVal_stat, names.arg = paste(sig_ind$species, sig_ind$section_group),
          horiz = TRUE, las = 1, col = "#7fbf7b", xlim = c(0, 1),
          xlab = "IndVal.g", main = "Site-level section indicators (series p<=0.05)")
} else {
  plot.new()
  title("Site-level section indicators\n(restricted-series p<=0.05)", cex.main = 0.9)
  text(0.5, 0.5, "No supported indicator species", cex = 1.05)
}
par(op); dev.off()

## =====================================================================================
## 6. SELF-CHECK  (stopifnot on dims / non-empty results) + console summary
## =====================================================================================
outputs <- c(
  file.path(OUT_DIR, "layerE_alpha_persample.tsv"),
  file.path(OUT_DIR, "layerE_alpha_site_time.tsv"),
  file.path(OUT_DIR, "layerE_inference_design.tsv"),
  file.path(OUT_DIR, "layerE_inext_asymptotic_section.tsv"),
  file.path(OUT_DIR, "layerE_inext_asymptotic_season.tsv"),
  file.path(OUT_DIR, "layerE_inext_coverage_section.tsv"),
  file.path(OUT_DIR, "layerE_inext_coverage_season.tsv"),
  file.path(OUT_DIR, "layerE_permanova.tsv"),
  file.path(OUT_DIR, "layerE_betadisper.tsv"),
  file.path(OUT_DIR, "layerE_da_consensus.tsv"),
  file.path(OUT_DIR, "layerE_indicator_species.tsv"),
  file.path(FIG_DIR, "phase3_layerE_inext_section.png"),
  file.path(FIG_DIR, "phase3_layerE_inext_season.png"),
  file.path(FIG_DIR, "phase3_layerE_ordination.png"),
  file.path(FIG_DIR, "phase3_layerE_ordination_library_descriptive.png"),
  file.path(FIG_DIR, "phase3_layerE_betadisper.png"),
  file.path(FIG_DIR, "phase3_layerE_da_indicator.png"))

stopifnot(
  nrow(alpha_ps) == 647,                       # every library retained descriptively
  nrow(alpha_unit) == 234,                     # ecological site x time units
  nrow(sp_site) == 14, nrow(sp_time) == 18,   # independent spatial / temporal profiles
  all(site_md$n_timepoints == 11), all(time_md$n_sites == 10),
  nrow(site_perm) == 27, nrow(time_perm) == 35,
  nrow(permanova) == 6, all(!is.na(permanova$p_value)),
  nrow(permdisp)  == 6, all(!is.na(permdisp$p_value)),
  !any(permanova$analysis_unit == "water-sample library"),
  nrow(da) == length(keep_taxa), length(keep_taxa) >= 2,
  all(da$analysis_unit == "site"),
  nrow(indic) == ncol(sp_site),
  all(file.exists(outputs)))

cat("\n==================== Layer E summary ====================\n")
cat("PERMANOVA (location) + PERMDISP (dispersion), independent units + series permutations:\n")
print(merge(permanova[permanova$term == "location",
              c("distance", "grouping", "analysis_unit", "n_units", "n_permutations", "R2", "F", "p_value")],
            permdisp[, c("distance", "grouping", "p_value")],
            by = c("distance", "grouping"), suffixes = c("_adonis", "_permdisp")),
      row.names = FALSE)
cat("\nExploratory site-level differential abundance (>=2 methods AND river-series q<0.05):\n")
print(da, row.names = FALSE)
cat("\nSite-level section indicator species (restricted-series p<=0.05):\n")
print(indic[!is.na(indic$p_value) & indic$p_value <= 0.05, ], row.names = FALSE)
cat(sprintf("\nWrote %d checked outputs. Layer E credibility correction complete.\n", length(outputs)))
