#!/usr/bin/env Rscript
# 08_validate_concordance.R -- Phase-2 Mac-side VALIDATION GATE.
# Does the reprocessed pipeline (MitoHelper Nov-2025) recover the authors' published
# 12S detections (Perry et al. 2024, 68-species normalised table)?
#
# The two DB vintages assign different nearest-neighbour labels to compatible short-marker evidence
# (reprocessed Salmo_macrostigma = harmonised Salmo; reprocessed Phoxinus names = source-study Phoxinus_ujmonensis),
# so a naive exact-name join is invalid. We bridge namespaces with a consensus crosswalk
# (three independent audit passes plus reconciliation; see data/derived/crosswalk.tsv) and report
# species-level AND genus-level recovery + read-weighted agreement over the shared replicates.
#
# Inputs : data/derived/species_x_sample_counts_rawnames.rds  (reprocessed 647 x 35 counts)
#          data/metadata/source_data_12S_normalised.tsv       (authors 665 x 68 normalised)
# Outputs: data/derived/crosswalk.tsv, concordance_summary.tsv, divergence_map.tsv,
#          species_x_sample_counts_clean.rds, genus_x_sample_counts_clean.rds
#          figures/phase2_concordance_occupancy.png, phase2_concordance_braycurtis.png
suppressMessages({library(readr); library(dplyr); library(tidyr); library(tibble); library(ggplot2)})
d <- "data/derived"; figdir <- "figures"; dir.create(figdir, showWarnings = FALSE)

# ============================================================================
# 1. Authors' table + D3 cleaning (68 raw -> 55 cleaned)  [same rules as R/05]
# ============================================================================
A <- read_tsv("data/metadata/source_data_12S_normalised.tsv", show_col_types = FALSE)
DROP_A  <- c("Gadus_macrocephalus","Hippoglossus_stenolepis","Ammodytes_personatus",
             "Clupea_harengus","Clupea_pallasii","Copadichromis_virginalis",
             "Favonigobius_gymnauchen","Pholis_crassispina","Melanogrammus_aeglefinus",
             "Platichthys_stellatus","Scomber_scombrus")
RELAB_A <- c(Salmo_obtusirostris = "Salmo_trutta",
             Anguilla_rostrata = "Anguilla_anguilla",
             Salvelinus_fontinalis_x_Salvelinus_malma = "Salvelinus_malma")
stopifnot(sum(DROP_A %in% names(A)) == 10)     # GATE: 10 of 11 strays present (Favonigobius absent = no-op)
stopifnot(all(names(RELAB_A) %in% names(A)))    # GATE: all 3 relabel sources present
Am <- as.matrix(A[,-1]); rownames(Am) <- A$ID
Am <- Am[, !colnames(Am) %in% DROP_A, drop = FALSE]
for (src in names(RELAB_A)) { tgt <- RELAB_A[[src]]
  if (src %in% colnames(Am)) {
    if (tgt %in% colnames(Am)) Am[, tgt] <- Am[, tgt] + Am[, src]
    else { Am <- cbind(Am, Am[, src]); colnames(Am)[ncol(Am)] <- tgt }
    Am <- Am[, colnames(Am) != src, drop = FALSE]
  }
}
stopifnot(ncol(Am) == 55)                      # GATE: authors cleaned universe = 55

# ============================================================================
# 2. My table: apply the consensus crosswalk -> authors namespace
# ============================================================================
M <- readRDS(file.path(d, "species_x_sample_counts_rawnames.rds"))  # 647 x 35 (incl Unassigned)
MAP <- c(
  Salmo_macrostigma = "Salmo_trutta", Phoxinus_phoxinus = "Phoxinus_ujmonensis",
  Salmo_salar = "Salmo_salar", Salmo_trutta = "Salmo_trutta",
  `Barbatula_aff._barbatula_lineage_I` = "Barbatula_barbatula",
  Anguilla_anguilla = "Anguilla_anguilla", Gasterosteus_aculeatus = "Gasterosteus_aculeatus",
  Pungitius_pungitius = "Pungitius_pungitius",
  Oncorhynchus_mykiss = "Oncorhynchus_mykiss_x_Salmo_salar",
  Pungitius_kaibarae = "Pungitius_pungitius", Perca_fluviatilis = "Perca_fluviatilis",
  Salvelinus_umbla = "Salvelinus_malma", Cottus_gobio = "Cottus_perifretum",
  `Phoxinus_sp._Morella` = "Phoxinus_ujmonensis", Barbatula_barbatula = "Barbatula_barbatula",
  Anguilla_rostrata = "Anguilla_anguilla", Phoxinus_bigerri = "Phoxinus_ujmonensis",
  Abramis_brama = "Abramis_brama", Salvelinus_fontinalis = "Salvelinus_malma",
  Perca_flavescens = "Perca_flavescens", Barbus_barbus = "Barbus_barbus")
DROP_M <- c("Clupea_harengus","Gadus_morhua","Myoxocephalus_thompsonii","Sprattus_sprattus",
            "Microgobius_gulosus","Pholis_gunnellus","Taurulus_bubalis","Ciliata_mustela",
            "Hippoglossus_stenolepis","Rhamdia_cf._jequitinhonha_NS_3056B","Merluccius_merluccius",
            "Melanogrammus_aeglefinus","Hyperoplus_lanceolatus","Unassigned")
stopifnot(setequal(colnames(M), c(names(MAP), DROP_M)))   # GATE: every reprocessed column handled

unassigned_reads <- if ("Unassigned" %in% colnames(M)) sum(M[, "Unassigned"]) else 0
Mk  <- M[, colnames(M) %in% names(MAP), drop = FALSE]
tgt <- MAP[colnames(Mk)]
harm <- unique(unname(tgt))
Mh <- vapply(harm, function(t) rowSums(Mk[, tgt == t, drop = FALSE]), numeric(nrow(Mk)))
rownames(Mh) <- rownames(Mk)                   # 647 x (harmonised taxa)

# persist the crosswalk as a documented artifact
write_tsv(tibble(my_species = names(MAP), authors_equivalent = unname(MAP), action = "map"),
          file.path(d, "crosswalk.tsv"))
write_tsv(tibble(my_species = DROP_M, authors_equivalent = "DROP", action = "drop"),
          file.path(d, "crosswalk.tsv"), append = TRUE)

# ============================================================================
# 3. Restrict to shared replicates (ID intersection)
# ============================================================================
shared <- intersect(rownames(Mh), rownames(Am))
stopifnot(length(shared) == 642)               # ASSERT the join (E##-vs-C## gotcha guard)
MhS <- Mh[shared, , drop = FALSE]; AmS <- Am[shared, , drop = FALSE]
n_rep <- length(shared)

# ============================================================================
# 4. Per-taxon detection counts over the 642 shared libraries (species-label level)
# ============================================================================
occ_me      <- colSums(MhS > 0)                # reprocessed per-taxon library detections
occ_me_str  <- colSums(MhS >= 2)               # strict (>=2 reads/cell)
occ_authors <- colSums(AmS > 0)
all_taxa <- union(colnames(MhS), colnames(AmS))
g <- function(v, taxa) { x <- setNames(rep(0L, length(taxa)), taxa); x[names(v)] <- v; x[taxa] }
om  <- g(occ_me, all_taxa); oms <- g(occ_me_str, all_taxa); oa <- g(occ_authors, all_taxa)

div <- tibble(taxon = all_taxa, occ_authors = oa, occ_me = om, occ_me_strict = oms) %>%
  mutate(in_my_namespace = taxon %in% colnames(MhS),
         detected_authors = occ_authors > 0,
         detected_me      = occ_me > 0,
         status = case_when(
           detected_authors &  detected_me ~ "both",
           detected_authors & !detected_me ~ "authors_only",
          !detected_authors &  detected_me ~ "me_only",
           TRUE ~ "neither")) %>%
  arrange(desc(occ_authors), desc(occ_me))

# recovery (recall): of source-study taxa detected in the 642 libraries, how many are re-detected
authors_detected <- div %>% filter(detected_authors)
rec_perm <- mean(authors_detected$occ_me     > 0)
rec_str  <- mean(authors_detected$occ_me_strict > 0)
jacc <- with(div, sum(status == "both") / sum(status %in% c("both","authors_only","me_only")))

# ============================================================================
# 5. Read-weighted agreement on COMMON taxa (renorm each row to 1)
# ============================================================================
common <- intersect(colnames(MhS), colnames(AmS))
relrow <- function(X) { rs <- rowSums(X); rs[rs == 0] <- 1; X / rs }
Pm <- relrow(MhS[, common, drop = FALSE]); Pa <- relrow(AmS[, common, drop = FALSE])
bray <- 1 - rowSums(abs(Pm - Pa)) / 2                       # Bray-Curtis similarity per rep
cosm <- rowSums(Pm * Pa) / (sqrt(rowSums(Pm^2)) * sqrt(rowSums(Pa^2)) + 1e-12)
# Spearman of per-label total relative abundance (do the same analysis labels dominate?)
sp_rho <- suppressWarnings(cor(colSums(Pm), colSums(Pa), method = "spearman"))

# ---- 5b. PRECISION: cell-level detection confusion on common taxa ----------
# The gate's recall side asks whether the reprocessing re-detects source-study cells;
# this is the reciprocal precision pass, per (library x taxon) cell.
Dm <- MhS[, common, drop = FALSE] > 0; Da <- AmS[, common, drop = FALSE] > 0
both_c <- sum(Dm & Da); me_only_c <- sum(Dm & !Da); au_only_c <- sum(!Dm & Da)
precision_cell <- both_c / (both_c + me_only_c)     # fraction of reprocessed detection cells also in source
recall_cell    <- both_c / (both_c + au_only_c)     # fraction of source detection cells re-detected
f1_cell <- 2 * precision_cell * recall_cell / (precision_cell + recall_cell)
mcnemar_p <- stats::mcnemar.test(matrix(c(both_c, au_only_c, me_only_c,
                                          sum(!Dm & !Da)), 2))$p.value
# per-taxon reprocessed-only cells, to see whether discordance concentrates in one label
fp_by_taxon <- sort(colSums(Dm & !Da), decreasing = TRUE)

# ============================================================================
# 6. GENUS-level robustness pass (the decisive 12S test)
# ============================================================================
genus <- function(x) sub("_.*$", "", x)
poolg <- function(X) { gg <- genus(colnames(X)); t(rowsum(t(X), gg)) }
Mg <- poolg(MhS); Ag <- poolg(AmS)
gtax <- union(colnames(Mg), colnames(Ag))
gom <- g(colSums(Mg > 0), gtax); goa <- g(colSums(Ag > 0), gtax)
gdiv <- tibble(genus = gtax, occ_authors = goa, occ_me = gom) %>%
  mutate(status = case_when(goa>0 & gom>0~"both", goa>0 & gom==0~"authors_only",
                            goa==0 & gom>0~"me_only", TRUE~"neither")) %>% arrange(desc(occ_authors))
grec <- with(filter(gdiv, occ_authors>0), mean(occ_me>0))
gcommon <- intersect(colnames(Mg), colnames(Ag))
Pmg <- relrow(Mg[,gcommon,drop=FALSE]); Pag <- relrow(Ag[,gcommon,drop=FALSE])
gbray <- 1 - rowSums(abs(Pmg - Pag))/2

# ============================================================================
# 7. Write artifacts + figures + Phase-3-ready cleaned count tables
# ============================================================================
write_tsv(div,  file.path(d, "divergence_map.tsv"))
write_tsv(gdiv, file.path(d, "divergence_map_genus.tsv"))
saveRDS(Mh,  file.path(d, "species_x_sample_counts_clean.rds"))     # 647 x harm taxa (Phase-3 input)
saveRDS(poolg(Mh), file.path(d, "genus_x_sample_counts_clean.rds"))

summ <- tibble(
  metric = c("shared_replicates","authors_cleaned_taxa","my_harmonised_taxa","common_taxa",
             "authors_taxa_detected_in642","species_recovery_permissive","species_recovery_strict",
             "jaccard_detected_sets","braycurtis_mean","braycurtis_median","cosine_mean",
             "spearman_taxon_abund","genus_recovery","genus_braycurtis_mean","unassigned_reads_frac",
             "cell_precision","cell_recall","cell_f1","cells_both","cells_me_only","cells_authors_only",
             "mcnemar_p"),
  value = c(n_rep, ncol(Am), ncol(Mh), length(common),
            sum(authors_detected$detected_authors), round(rec_perm,3), round(rec_str,3),
            round(jacc,3), round(mean(bray),3), round(median(bray),3), round(mean(cosm),3),
            round(sp_rho,3), round(grec,3), round(mean(gbray),3),
            round(unassigned_reads/(sum(M)),3),
            round(precision_cell,3), round(recall_cell,3), round(f1_cell,3),
            both_c, me_only_c, au_only_c, signif(mcnemar_p,3)))
write_tsv(summ, file.path(d, "concordance_summary.tsv"))

# occupancy scatter (species, both-detected and divergent)
pdat <- div %>% filter(occ_authors > 0 | occ_me > 0)
p1 <- ggplot(pdat, aes(occ_authors, occ_me, colour = status)) +
  geom_abline(slope = 1, linetype = 2, colour = "grey60") +
  geom_point(size = 2, alpha = .8) +
  geom_text(aes(label = taxon), size = 2.4, vjust = -0.6, check_overlap = TRUE, show.legend = FALSE) +
  labs(title = "Phase-2 concordance: per-taxon detections across 642 shared libraries",
       x = "source-study detections (libraries)", y = "reprocessed detections (libraries)") + theme_bw()
ggsave(file.path(figdir, "phase2_concordance_occupancy.png"), p1, width = 8, height = 6, dpi = 130)
p2 <- ggplot(tibble(bray = bray), aes(bray)) +
  geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = mean(bray), colour = "red") +
  labs(title = sprintf("Read-weighted agreement per replicate (mean BC sim = %.2f)", mean(bray)),
       x = "Bray-Curtis similarity (common taxa, renormalised)", y = "replicates") + theme_bw()
ggsave(file.path(figdir, "phase2_concordance_braycurtis.png"), p2, width = 8, height = 5, dpi = 130)

# ============================================================================
# 8. Console report
# ============================================================================
cat("\n================ PHASE-2 VALIDATION GATE ================\n")
cat(sprintf("shared libraries: %d | source-study cleaned taxa: %d | reprocessed harmonised taxa: %d | common: %d\n",
            n_rep, ncol(Am), ncol(Mh), length(common)))
cat(sprintf("source-study taxa DETECTED in the 642 shared libraries: %d\n", sum(authors_detected$detected_authors)))
cat(sprintf("SPECIES recovery: permissive %.1f%%  | strict(>=2 reads) %.1f%%  | Jaccard %.2f\n",
            100*rec_perm, 100*rec_str, jacc))
cat(sprintf("GENUS  recovery: %.1f%%   (the decisive 12S pass)\n", 100*grec))
cat(sprintf("read-weighted: Bray-Curtis mean %.2f / median %.2f | cosine %.2f | Spearman(abund) %.2f\n",
            mean(bray), median(bray), mean(cosm), sp_rho))
cat(sprintf("genus read-weighted Bray-Curtis mean: %.2f\n", mean(gbray)))
cat(sprintf("Unassigned reads: %.1f%% of reprocessed total\n", 100*unassigned_reads/sum(M)))
cat(sprintf("PRECISION pass (cell-level, common taxa): precision %.2f | recall %.2f | F1 %.2f | McNemar p=%.2g\n",
            precision_cell, recall_cell, f1_cell, mcnemar_p))
cat(sprintf("  cells: both %d | reprocessed-only %d | source-study-only %d\n", both_c, me_only_c, au_only_c))
cat("  reprocessed-only cells by taxon: "); cat(paste(sprintf("%s=%d", names(fp_by_taxon), fp_by_taxon)[fp_by_taxon>0], collapse=", "), "\n\n")
cat("--- SOURCE-STUDY-ONLY detections (potential recovery gaps) ---\n")
print(as.data.frame(div %>% filter(status == "authors_only") %>% select(taxon, occ_authors)), row.names = FALSE)
cat("\n--- REPROCESSED-ONLY detections (extra; should be sparse after drops) ---\n")
print(as.data.frame(div %>% filter(status == "me_only") %>% select(taxon, occ_me)), row.names = FALSE)
cat("\n--- BOTH (the recovered set) ---\n")
print(as.data.frame(div %>% filter(status == "both") %>% select(taxon, occ_authors, occ_me)), row.names = FALSE)
cat("\nwrote: concordance_summary.tsv, divergence_map.tsv(+_genus), crosswalk.tsv,\n")
cat("       species_x_sample_counts_clean.rds, genus_x_sample_counts_clean.rds, 2 figures\n")
