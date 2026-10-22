#!/usr/bin/env Rscript
# 07_reclass_B.R  --  Phase 3, Layer B: reference-database currency / reclassification.
# Plan eDNAplan.md  §2.3  +  Gate 4 (locked)  +  docs/prereg.md §5.
#
# ---------------------------------------------------------------------------
# GATE-4 BASELINE DECISION (locked) -- paste into docs/prereg.md §5.
#   The authors' Sep-2020 MitoFish build is UNARCHIVED, so a faithful 2020 replay
#   (§2.3 option i/ii) is unattainable. Layer B therefore runs as OPTION (iii):
#   2026 concordance, NOT a validated 2020->2026 change. Two axes:
#     - Baseline in hand : MitoHelper 12S *Nov-2025* first-hit assignment
#                          (data/derived/taxonomy.tsv; 886/1443 ASVs assigned).
#     - Re-assignment    : MMseqs2 @ 97% identity (Bayer et al. 2025: MMseqs2 >
#                          BLAST for short 12S) against the *MitoFish 12S Nov-2025*
#                          build (MitoHelper Zenodo 10.5281/zenodo.17602902,
#                          file mitofish.12S.Nov2025.tsv, 64,059 records).
#   HONEST SCOPE CAVEAT: both Part-B assignments use the same MitoFish vintage,
#   so Part B measures classifier/method concordance (MMseqs2 vs first-hit), not
#   database currency. Part C adds MIDORI2 as a separate current reference
#   resource. Agreement or disagreement among these labels is an assignment audit;
#   it is not biological truth, a local-occurrence conclusion, or a novel-taxon result.
#   (This block is also emitted to data/derived/reclassB_prereg5_addendum.md.)
# ---------------------------------------------------------------------------
#
# PART A (real, no downloads): tabulate the in-hand 2026 (MitoHelper Nov-2025)
#   assignment quality -- per-rank assignment rate (% ASVs, % reads), the
#   Unassigned fraction, and a catalogue of assignment-audit hypotheses.
#   Outputs: data/derived/reclassB_assignment_rates.tsv,
#            data/derived/reclassB_drift_catalogue.tsv,
#            figures/reclassB_assignment_and_drift.png
#
# PART B (run or cached receipt): MMseqs2 @97% re-assignment of the ASVs in
#   repseqs.fasta against the MitoFish 12S Nov-2025 build, reporting
#   concordance vs taxonomy.tsv. MMseqs2 lives in a SEPARATE throwaway env
#   'mmseqs2' (never touches 'edna'); the DB is fetched to a scratch dir, never
#   the repo. If the reference is absent but valid result tables are present, the
#   pipeline retains them, records their date/hash, and marks them not recomputed.
#   Otherwise it stays runnable and emits a reproduction command and DB metadata.
#   Outputs: data/derived/reclassB_mmseqs_concordance.tsv,
#            data/derived/reclassB_mmseqs_summary.tsv,
#            data/derived/reclassB_partB_provenance.txt
#
# PART C (cross-reference comparison): assignment against the separate MIDORI2
#   GenBank271 srRNA reference resource (mitochondrial 12S rRNA). MMseqs2 @97%
#   (same low-RAM params) re-assigns
#   the ASVs; concordance vs the in-hand MitoHelper Nov-2025 taxonomy.tsv is scored
#   per-ASV + summary (species/genus agreement, read-weighted, newly-resolved /
#   changed / lost) and audits how the Salmo_macrostigma-labelled ASVs are labelled
#   by the second resource. Same throwaway 'mmseqs2'
#   env; DB fetched (capped) to a scratch dir outside the repo. Gated on the DB
#   being present or RECLASSB_FETCH_INDEP=1, so the script still runs DB-absent.
#   Outputs: data/derived/reclassB_indepDB_concordance.tsv,
#            data/derived/reclassB_indepDB_summary.tsv,
#            data/derived/reclassB_partC_provenance.txt,
#            figures/reclassB_indepDB.png
#
# Reads != abundance: everything below is occurrence/detection framing.
# Idempotent: overwrites its own reclassB_*/figures outputs; scratch stays outside repo.
suppressWarnings(suppressPackageStartupMessages({
  library(utils); library(stats); library(grDevices); library(graphics)
}))

# ---- repo-root resolution (works regardless of cwd; mirrors R/lib_design.R) ----
.script_dir <- function() {
  fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(fa)) dirname(normalizePath(fa[1])) else getwd()
}
ROOT <- local({
  d <- .script_dir()
  cand <- unique(c(normalizePath(file.path(d, ".."), mustWork = FALSE), getwd(),
                   normalizePath(file.path(getwd(), ".."), mustWork = FALSE)))
  hit <- cand[file.exists(file.path(cand, "data", "derived", "taxonomy.tsv"))]
  if (!length(hit)) stop("07_reclass_B.R: cannot locate repo root (data/derived/taxonomy.tsv not found)")
  hit[1]
})
P <- function(...) file.path(ROOT, ...)

# Ground rule: source lib_design.R. Layer B is design-independent (it works on the
# ASV x reference-DB axis, not the site x time occasion grid) so we source it only
# defensively -- a design-file hiccup must not break the DB-currency analysis.
suppressWarnings(tryCatch(source(P("R", "lib_design.R")),
  error = function(e) message("note: lib_design.R not sourced (Layer B does not use the occasion design): ", conditionMessage(e))))

message("== Layer B (reference-DB currency) ==  ROOT = ", ROOT)

# ===========================================================================
# PART A  --  in-hand 2026 (MitoHelper 12S Nov-2025) assignment quality
# ===========================================================================
tax  <- read.delim(P("data", "derived", "taxonomy.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
amap <- read.delim(P("data", "derived", "asv_map.tsv"),  stringsAsFactors = FALSE, check.names = FALSE)
seqtab <- readRDS(P("data", "derived", "seqtab_clean.rds"))

stopifnot(all(c("asv_id", "Class", "Order", "Family", "Genus", "Species") %in% names(tax)),
          all(c("asv_id", "total_reads", "sequence") %in% names(amap)))

# read bases -------------------------------------------------------------------
#  FULL basis  = asv_map$total_reads over the 1,443-ASV taxonomy universe
#                (matches the "ASV0001 = 18.5M reads = 42% of the table" headline).
#  CLEAN basis = colSums(seqtab_clean) over the 1,348 ASVs in the analysis table.
amap <- amap[match(tax$asv_id, amap$asv_id), ]         # align to taxonomy order
tax$total_reads <- amap$total_reads
TOT_FULL  <- sum(amap$total_reads)
seq2reads <- colSums(seqtab)
id2seq    <- setNames(amap$sequence, amap$asv_id)
tax$reads_clean <- unname(seq2reads[id2seq[tax$asv_id]]); tax$reads_clean[is.na(tax$reads_clean)] <- 0
TOT_CLEAN <- sum(seq2reads)
N_ASV     <- nrow(tax)

is_unassigned <- function(x) is.na(x) | x %in% c("", "Unassigned", "NA", "unassigned", "Unclassified", "unclassified")

ranks <- c("Class", "Order", "Family", "Genus", "Species")
rate_rows <- lapply(ranks, function(r) {
  ok <- !is_unassigned(tax[[r]])
  data.frame(
    rank                    = r,
    n_asv_total             = N_ASV,
    n_asv_assigned          = sum(ok),
    n_asv_unassigned        = sum(!ok),
    pct_asv_assigned        = round(100 * sum(ok) / N_ASV, 2),
    reads_assigned_full     = sum(tax$total_reads[ok]),
    pct_reads_assigned_full = round(100 * sum(tax$total_reads[ok]) / TOT_FULL, 2),
    reads_assigned_clean    = sum(tax$reads_clean[ok]),
    pct_reads_assigned_clean= round(100 * sum(tax$reads_clean[ok]) / TOT_CLEAN, 2),
    stringsAsFactors = FALSE)
})
assign_rates <- do.call(rbind, rate_rows)
write.table(assign_rates, P("data", "derived", "reclassB_assignment_rates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sp_ok  <- !is_unassigned(tax$Species)
UNASSIGNED_ASV      <- sum(!sp_ok)
UNASSIGNED_READS_F  <- sum(tax$total_reads[!sp_ok])
UNASSIGNED_PCT_F    <- round(100 * UNASSIGNED_READS_F / TOT_FULL, 2)

message(sprintf("PART A: species-level assigned = %d/%d ASVs (%.1f%%) | %s reads (%.1f%% full basis)",
        sum(sp_ok), N_ASV, 100 * sum(sp_ok) / N_ASV,
        format(sum(tax$total_reads[sp_ok]), big.mark = ","), assign_rates$pct_reads_assigned_full[ranks == "Species"]))
message(sprintf("PART A: Unassigned = %d/%d ASVs | %s reads (%.1f%% full basis)",
        UNASSIGNED_ASV, N_ASV, format(UNASSIGNED_READS_F, big.mark = ","), UNASSIGNED_PCT_F))

# ---- species-level aggregation (assigned) ------------------------------------
asg <- tax[sp_ok, ]
sp_agg <- aggregate(cbind(total_reads, reads_clean) ~ Species + Genus + Family + Order + Class,
                    data = asg, FUN = sum)
n_by_sp <- as.data.frame(table(asg$Species), stringsAsFactors = FALSE); names(n_by_sp) <- c("Species", "n_asv")
sp_agg <- merge(sp_agg, n_by_sp, by = "Species")
sp_agg$pct_reads_full <- round(100 * sp_agg$total_reads / TOT_FULL, 3)
sp_agg <- sp_agg[order(-sp_agg$total_reads), ]

# ---- assignment audit hypotheses (not biological identity conclusions) -------
# Rows key on the EXACT taxonomy.tsv Species string. The audit columns distinguish
# an observed DB label from a hypothesis to review. They do not assert local
# occurrence, native status, provenance, or a replacement species identity.
audit <- read.csv(text = paste(
'db_species|audit_category|audit_hypothesis|evidence_flag|audit_note',
'Salmo macrostigma|species-resolution-review|Salmo sp.; species unresolved|cross-reference-label-conflict|The same ASV set receives Salmo macrostigma in MitoHelper, predominantly Salmo trutta in Part B, and predominantly Salmo carpio in Part C; genus Salmo is stable in these comparisons. None of the labels establishes a local species or native-status identity. 283 ASVs; ASV0001 has 18.5M reads.',
'Phoxinus phoxinus|namespace-review|Phoxinus sp.; retain reported label for audit|same-genus-label-cluster|Multiple Phoxinus species or lineage labels occur in the assignment table. Reference and taxonomy curation are required before harmonisation; no replacement or local identity is inferred.',
'Phoxinus sp. Morella|namespace-review|Phoxinus sp.; retain reported label for audit|same-genus-label-cluster|Multiple Phoxinus species or lineage labels occur in the assignment table. Reference and taxonomy curation are required before harmonisation; no replacement or local identity is inferred.',
'Phoxinus bigerri|namespace-review|Phoxinus sp.; retain reported label for audit|same-genus-label-cluster|Multiple Phoxinus species or lineage labels occur in the assignment table. Reference and taxonomy curation are required before harmonisation; no replacement or local identity is inferred.',
'Barbatula aff. barbatula lineage I|namespace-review|Barbatula sp.; retain reported label for audit|same-genus-label-cluster|The assignment table contains both an affinis-lineage label and a species label in Barbatula. This is a harmonisation review flag, not evidence that either label is the local species identity.',
'Barbatula barbatula|namespace-review|Barbatula sp.; retain reported label for audit|same-genus-label-cluster|The assignment table contains both an affinis-lineage label and a species label in Barbatula. This is a harmonisation review flag, not evidence that either label is the local species identity.',
'Cottus gobio|taxonomy-review|Cottus sp.; retain reported label for audit|database-label-only|The database label is retained as an observed assignment. Layer B does not provide a curated-reference audit sufficient to replace it or establish local species identity.',
'Salmo trutta|species-resolution-review|Salmo sp.; retain reported label for audit|short-marker-species-review|The reported species label is not treated as local identity evidence because the cross-reference comparisons show species-label instability within Salmo.',
'Salmo salar|species-resolution-review|Salmo sp.; retain reported label for audit|short-marker-species-review|The reported species label is not treated as local identity evidence because the cross-reference comparisons show species-label instability within Salmo.',
'Anguilla anguilla|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence, conservation-status, or corrected-identity claim.',
'Anguilla rostrata|species-resolution-review|Anguilla sp.; no replacement species proposed|same-genus-alternative|One ASV carries this label. A same-genus alternative is a review hypothesis only; geography alone is not used to substitute an identity.',
'Gasterosteus aculeatus|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence or native-status claim.',
'Pungitius pungitius|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence or native-status claim.',
'Pungitius kaibarae|geographic-context-review|Review against curated regional Pungitius references; no replacement proposed|geographic-context-only|One ASV carries this label. Geographic context motivates reference review but does not establish a definitive exclusion or prove another Pungitius identity.',
'Perca fluviatilis|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence or native-status claim.',
'Perca flavescens|geographic-context-review|Review against curated Perca references; no replacement proposed|geographic-context-only|One ASV with four reads carries this label. Low read support and geographic context motivate review but do not prove a mis-hit or another Perca identity.',
'Oncorhynchus mykiss|reported-label-review|Retain reported DB label pending independent validation|database-label-only|Five ASVs carry this label. Layer B provides no evidence for a local-origin, stocking, or native-status interpretation.',
'Salvelinus umbla|species-resolution-review|Salvelinus sp.; no replacement species proposed|same-genus-alternative|Multiple Salvelinus species labels occur at low support. The audit retains genus-level uncertainty and makes no local species or status claim.',
'Salvelinus fontinalis|species-resolution-review|Salvelinus sp.; no replacement species proposed|same-genus-alternative|Multiple Salvelinus species labels occur at low support. The audit retains genus-level uncertainty and makes no local species or status claim.',
'Rhamdia cf. jequitinhonha NS 3056B|geographic-context-review|Review against curated references; no replacement proposed|geographic-context-only|Two ASVs carry this label. Geographic context is a reason to review the match, not a basis for definitive exclusion or assignment of another identity.',
'Myoxocephalus thompsonii|geographic-context-review|Review against curated references; no replacement proposed|geographic-context-only|One ASV carries this label. Geographic context is a review hypothesis only; the Layer B comparison does not establish a replacement identity.',
'Hippoglossus stenolepis|geographic-context-review|Review against curated references; no replacement proposed|geographic-context-only|One ASV carries this label. Geographic context is a review hypothesis only; the Layer B comparison does not establish a replacement identity.',
'Microgobius gulosus|geographic-context-review|Review against curated references; no replacement proposed|geographic-context-only|One ASV carries this label. Geographic context is a review hypothesis only; the Layer B comparison does not establish a replacement identity.',
'Scomber scombrus|control-context|Positive-control-associated assignment|positive-control-metadata|The label is evaluated against the recorded positive-control context and is not interpreted as a field identity.',
'Merlangius merlangus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Clupea harengus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Gadus morhua|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Sprattus sprattus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Merluccius merluccius|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Melanogrammus aeglefinus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Pholis gunnellus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Ciliata mustela|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Osmerus eperlanus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Taurulus bubalis|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Hyperoplus lanceolatus|target-scope-review|Retain reported DB label; review separately from freshwater target inference|analysis-scope-only|The assignment is flagged only against the pre-specified analysis target scope; this is not an identity or occurrence conclusion.',
'Abramis brama|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence or native-status claim.',
'Barbus barbus|reported-label-review|Retain reported DB label pending independent validation|database-label-only|The database label is reported without a local-occurrence or native-status claim.',
sep = "\n"), sep = "|", header = TRUE, stringsAsFactors = FALSE, quote = "")

drift <- merge(sp_agg, audit, by.x = "Species", by.y = "db_species", all.x = TRUE)
# Any assigned species without a curated row remains explicitly unreviewed. The
# observed database label is not silently promoted to a verified identity.
miss <- is.na(drift$audit_category)
if (any(miss)) {
  drift$audit_category[miss] <- "unreviewed"
  drift$audit_hypothesis[miss] <- "Retain reported DB label pending review"
  drift$evidence_flag[miss] <- "database-label-only"
  drift$audit_note[miss] <- "Not in the curated audit table; no identity or status inference made."
}
drift <- drift[order(-drift$total_reads),
               c("Species", "Genus", "Family", "audit_category", "audit_hypothesis", "evidence_flag",
                 "n_asv", "total_reads", "pct_reads_full", "reads_clean", "audit_note")]
names(drift)[1] <- "db_species"
write.table(drift, P("data", "derived", "reclassB_drift_catalogue.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Lamprey-class note: the marker output is unresolved below Class; make no
# species, conservation-status, or local-identity inference.
hyp <- tax[!is_unassigned(tax$Class) & tax$Class == "Hyperoartia", ]
if (nrow(hyp)) message(sprintf("PART A: Class Hyperoartia, %d ASV / %s reads; unresolved below Class, with no species/status inference.",
                               nrow(hyp), format(sum(hyp$total_reads), big.mark = ",")))

message(sprintf("PART A: assignment-audit catalogue = %d species rows (%d rows require contextual or reference review)",
                nrow(drift), sum(drift$audit_category != "reported-label-review")))

# ---- figure: assignment rates (per rank) + drift catalogue (top taxa) --------
figp <- P("figures", "reclassB_assignment_and_drift.png")
dir.create(dirname(figp), showWarnings = FALSE, recursive = TRUE)
cat_cols <- c(`reported-label-review` = "#78909C", `species-resolution-review` = "#E65100",
              `namespace-review` = "#7E57C2", `taxonomy-review` = "#00897B",
              `geographic-context-review` = "#C0A000", `target-scope-review` = "#5C6BC0",
              `control-context` = "#616161", unreviewed = "#9E9E9E")
grDevices::png(figp, width = 1500, height = 720, res = 130)
op <- par(mfrow = c(1, 2), mar = c(6.5, 4.2, 3.2, 1.1), xpd = NA)
# Panel A: per-rank assignment rate (% ASV vs % reads)
m <- rbind(`% ASVs` = assign_rates$pct_asv_assigned,
           `% reads` = assign_rates$pct_reads_assigned_full)
colnames(m) <- assign_rates$rank
bp <- barplot(m, beside = TRUE, ylim = c(0, 108), col = c("#455A64", "#26A69A"),
              ylab = "assigned (%)", main = "A. MitoHelper 12S Nov-2025 assignment rate",
              cex.names = 0.95, las = 1)
text(bp[1, ], m[1, ] + 3, sprintf("%.0f%%", m[1, ]), cex = 0.68, col = "#455A64")
text(bp[2, ], m[2, ] + 3, sprintf("%.0f%%", m[2, ]), cex = 0.68, col = "#00695C")
legend("top", c("% ASVs", "% reads"), fill = c("#455A64", "#26A69A"), bty = "n",
       cex = 0.85, horiz = TRUE, inset = -0.02)
mtext(sprintf("Species-level: %d/%d ASVs (%.0f%%) | reads %.0f%% assigned; Unassigned = %d ASVs / %.0f%% reads",
      sum(sp_ok), N_ASV, 100 * sum(sp_ok)/N_ASV, assign_rates$pct_reads_assigned_full[5],
      UNASSIGNED_ASV, UNASSIGNED_PCT_F), side = 1, line = 5, cex = 0.68, col = "#555555")
# Panel B: top species by reads, coloured by neutral audit category
topn <- head(drift[order(-drift$pct_reads_full), ], 14)
topn <- topn[order(topn$pct_reads_full), ]
bcol <- cat_cols[topn$audit_category]; bcol[is.na(bcol)] <- "#9E9E9E"
bp2 <- barplot(topn$pct_reads_full, horiz = TRUE, col = bcol, border = NA,
               names.arg = topn$db_species, las = 1, cex.names = 0.6, xlab = "% of reads (full basis)",
               main = "B. Top DB calls by reads (assignment-audit categories)", cex.main = 0.98,
               xlim = c(0, max(topn$pct_reads_full) * 1.18))
text(topn$pct_reads_full + max(topn$pct_reads_full) * 0.01, bp2,
     sprintf("%.1f%%", topn$pct_reads_full), pos = 4, cex = 0.6, col = "#222222")
legend("bottomright", names(cat_cols), fill = cat_cols, bty = "n", cex = 0.6, title = "audit category")
par(op); grDevices::dev.off()
message("PART A: figure -> ", figp)

# ===========================================================================
# PART B  --  MMseqs2 @97% re-assignment vs MitoFish 12S Nov-2025 (attempt)
# ===========================================================================
# Config (all overridable by env var so the repo script bakes in no session path).
MMSEQS_ENV  <- Sys.getenv("RECLASSB_MMSEQS_ENV", "mmseqs2")
DB_DIR      <- Sys.getenv("RECLASSB_DBDIR", file.path(tempdir(), "reclassB_dbfetch"))
DB_TSV_NAME <- "mitofish.12S.Nov2025.tsv"
DB_URL      <- "https://zenodo.org/api/records/17602902/files/mitofish.12S.Nov2025.tsv/content"
DB_DOI      <- "10.5281/zenodo.17602902"
DB_NAME     <- "MitoFish 12S (via MitoHelper), Nov-2025 build; 64,059 records >50bp"
ACCESS_DATE <- as.character(Sys.Date())
MM_PARAMS   <- c("--search-type", "3", "--min-seq-id", "0.97", "-c", "0.9",
                 "--cov-mode", "2", "-s", "7.5", "-k", "13", "--max-seqs", "25")
CACHED_PARTB_FILES <- c(P("data", "derived", "reclassB_mmseqs_concordance.tsv"),
                        P("data", "derived", "reclassB_mmseqs_summary.tsv"))
HAS_CACHED_PARTB <- all(file.exists(CACHED_PARTB_FILES))
CACHED_PARTB_DATE <- if (HAS_CACHED_PARTB)
  as.character(as.Date(min(file.info(CACHED_PARTB_FILES)$mtime))) else NA_character_
# -k 13 is REQUIRED on low-RAM hosts: the nucleotide k-mer index table is 4^k
# slots (default k=15 -> ~8.6GB, overflows an 8GB Mac). k=13 -> ~0.5GB.

find_mmseqs <- function() {
  cands <- c(Sys.getenv("RECLASSB_MMSEQS_BIN", ""),
             file.path(Sys.getenv("HOME"), "micromamba", "envs", MMSEQS_ENV, "bin", "mmseqs"),
             file.path(Sys.getenv("HOME"), "micromamba", "envs", "mmseqs2", "bin", "mmseqs"),
             Sys.which("mmseqs")[[1]])
  cands <- cands[nzchar(cands)]
  hit <- cands[file.exists(cands)]
  if (length(hit)) return(hit[1]) else return(NA_character_)
}
run <- function(bin, args, ...) suppressWarnings(system2(bin, args, stdout = TRUE, stderr = TRUE, ...))

partB <- list(status = "not-run", reason = "", n_query = NA, n_hit = NA,
              mmseqs_bin = NA, mmseqs_version = NA, db_present = FALSE)
prov <- character(0)
add <- function(...) prov[[length(prov) + 1L]] <<- paste0(...)

add("Layer B PART B provenance -- ", Sys.time())
add("Gate-4 decision: OPTION (iii) 2026 concordance; NOT a 2020 replay (Sep-2020 MitoFish unarchived).")
add("Part-B ref  : ", DB_NAME)
add("DB DOI      : ", DB_DOI, "  file: ", DB_TSV_NAME)
add("DB URL      : ", DB_URL)
add("current run : ", ACCESS_DATE, "  (receipt-generation date; not the cached-result date)")
add("classifier  : MMseqs2 easy-search  params: ", paste(MM_PARAMS, collapse = " "))
add("env         : throwaway micromamba env '", MMSEQS_ENV, "' (install: micromamba create -y -n ", MMSEQS_ENV, " -c conda-forge -c bioconda mmseqs2)")
add("SCOPE       : taxonomy.tsv and this search use the same MitoFish vintage, so Part B is only a classifier/method comparison (MMseqs2 vs first-hit).")
add("PART C      : the separate MIDORI2 comparison and its run/cache status are recorded independently in reclassB_partC_provenance.txt.")
if (HAS_CACHED_PARTB) {
  cache_md5 <- unname(tools::md5sum(CACHED_PARTB_FILES))
  add("cached tables: ", paste(basename(CACHED_PARTB_FILES), collapse = ", "),
      "  table-date: ", CACHED_PARTB_DATE)
  add("cached md5  : ", paste(paste0(basename(CACHED_PARTB_FILES), "=", cache_md5), collapse = "; "))
}
add("REPRODUCTION COMMAND ONLY (current-run scratch path; not original-result provenance):")
add("  mkdir -p ", DB_DIR, " && \\")
add("  curl -sL --max-time 240 --max-filesize 300000000 -o ", file.path(DB_DIR, DB_TSV_NAME), " \\")
add("    '", DB_URL, "'")

mmseqs_bin <- find_mmseqs()
db_tsv     <- file.path(DB_DIR, DB_TSV_NAME)
repseqs    <- P("data", "derived", "repseqs.fasta")
partB$mmseqs_bin <- mmseqs_bin
partB$db_present <- file.exists(db_tsv) && file.info(db_tsv)$size > 1e6

if (!partB$db_present && HAS_CACHED_PARTB) {
  cached_s <- read.delim(CACHED_PARTB_FILES[2], stringsAsFactors = FALSE)
  cached_cc <- read.delim(CACHED_PARTB_FILES[1], stringsAsFactors = FALSE)
  bgv <- function(m) cached_s$value[cached_s$metric == m]
  required_metrics <- c("n_query_ASV", "n_with_mmseqs_hit", "n_both_assigned",
                        "pct_species_agree_of_both", "pct_genus_agree_of_both",
                        "asv0001_taxonomy", "asv0001_mmseqs")
  stopifnot(all(vapply(required_metrics, function(m) length(bgv(m)) == 1L, logical(1))))
  cached_asv1 <- cached_cc[cached_cc$asv_id == "ASV0001", , drop = FALSE]
  stopifnot(nrow(cached_asv1) == 1L,
            identical(cached_asv1$tax_species, bgv("asv0001_taxonomy")),
            identical(cached_asv1$mm_species, bgv("asv0001_mmseqs")))

  partB$status <- "cached-results"
  partB$n_query <- as.numeric(bgv("n_query_ASV"))
  partB$n_hit <- as.numeric(bgv("n_with_mmseqs_hit"))
  partB$reason <- paste0("Part-B reference TSV absent from the current scratch path; retained result tables from ",
                         CACHED_PARTB_DATE, " without recomputation.")
  if (!is.na(mmseqs_bin)) {
    v <- tryCatch(run(mmseqs_bin, "version"), error = function(e) NA)
    partB$mmseqs_version <- if (length(v)) v[1] else NA
    add("current mmseqs: ", mmseqs_bin, "  version: ", partB$mmseqs_version,
        "  (available now; not used to recompute the cached tables)")
  }
  add("RESULT      : ", partB$reason)
  add(sprintf("CACHED RECEIPT: %s/%s ASVs hit; species-agree=%s%%, genus-agree=%s%% of %s both-assigned.",
              bgv("n_with_mmseqs_hit"), bgv("n_query_ASV"),
              bgv("pct_species_agree_of_both"), bgv("pct_genus_agree_of_both"),
              bgv("n_both_assigned")))
  add(sprintf("CACHED ASV0001 LABELS: MitoHelper=%s; Part-B MMseqs2=%s; pident=%s; comparison=%s. These are reference labels, not a local species identity.",
              cached_asv1$tax_species, cached_asv1$mm_species,
              cached_asv1$mm_pident, cached_asv1$category))
} else if (is.na(mmseqs_bin)) {
  partB$status <- "documented-only"; partB$reason <- "mmseqs binary not found (install env 'mmseqs2')."
  add("RESULT      : mmseqs binary NOT found -> pipeline documented, not executed.")
} else {
  v <- tryCatch(run(mmseqs_bin, "version"), error = function(e) NA)
  partB$mmseqs_version <- if (length(v)) v[1] else NA
  add("mmseqs bin  : ", mmseqs_bin, "  version: ", partB$mmseqs_version)

  # fetch DB if absent (capped) -------------------------------------------------
  if (!file.exists(db_tsv)) {
    dir.create(DB_DIR, showWarnings = FALSE, recursive = TRUE)
    message("PART B: DB not present -> attempting capped fetch ...")
    fetch <- tryCatch(run("curl", c("-sL", "--max-time", "240", "--max-filesize", "300000000",
                                    "-o", db_tsv, DB_URL)), error = function(e) conditionMessage(e))
    if (!file.exists(db_tsv) || file.info(db_tsv)$size < 1e6) {
      add("RESULT      : DB fetch failed/too-heavy -> pipeline documented, not executed.")
    }
  }
  partB$db_present <- file.exists(db_tsv) && file.info(db_tsv)$size > 1e6

  if (!partB$db_present) {
    partB$status <- "documented-only"
    partB$reason <- "DB not fetched within cap; exact fetch command recorded in provenance."
    add("RESULT      : DB absent -> documented-only.")
  } else {
    add("DB present  : ", db_tsv, "  (", round(file.info(db_tsv)$size / 1e6, 1), " MB)")
    # build reference FASTA + accession->taxonomy (pure R; portable, no awk dep)
    reffa   <- file.path(DB_DIR, "reclassB_ref12S.fasta")
    acctax  <- file.path(DB_DIR, "reclassB_acc2tax.tsv")
    if (!file.exists(reffa) || !file.exists(acctax)) {
      message("PART B: building reference FASTA + acc->taxonomy from DB tsv ...")
      ref <- read.delim(db_tsv, stringsAsFactors = FALSE, quote = "", check.names = FALSE, comment.char = "")
      ref <- ref[nzchar(ref$Sequence), ]
      writeLines(paste0(">", ref$Accession, "\n", ref$Sequence), reffa)
      write.table(data.frame(acc = ref$Accession, Class = ref$Class, Order = ref$Order,
                             Family = ref$Family, Genus = ref$Genus, Species = ref$Species,
                             stringsAsFactors = FALSE),
                  acctax, sep = "\t", quote = FALSE, row.names = FALSE)
      rm(ref); invisible(gc())
    }
    tmpdir <- file.path(DB_DIR, "mmtmp"); unlink(tmpdir, recursive = TRUE)
    outm8  <- file.path(DB_DIR, "reclassB_hits.m8")
    message("PART B: running MMseqs2 easy-search on the 1,348 ASVs ...")
    es <- run(mmseqs_bin, c("easy-search", repseqs, reffa, outm8, tmpdir,
                            MM_PARAMS, "--threads", "4",
                            "--format-output", "query,target,pident,alnlen,evalue,bits,qcov"))
    if (!file.exists(outm8)) {
      partB$status <- "run-failed"; partB$reason <- paste(tail(es, 2), collapse = " | ")
      add("RESULT      : MMseqs2 run FAILED: ", partB$reason)
    } else {
      hits <- tryCatch(read.delim(outm8, header = FALSE, stringsAsFactors = FALSE,
                col.names = c("query", "target", "pident", "alnlen", "evalue", "bits", "qcov")),
                error = function(e) data.frame())
      a2t <- read.delim(acctax, header = TRUE, stringsAsFactors = FALSE, quote = "")
      if (!nrow(hits)) {
        partB$status <- "run-empty"; partB$reason <- "MMseqs2 returned 0 hits @97%."
        add("RESULT      : MMseqs2 returned 0 hits.")
      } else {
        # best hit per query (max bits)
        hits <- hits[order(hits$query, -hits$bits), ]
        best <- hits[!duplicated(hits$query), ]
        idx  <- match(best$target, a2t$acc)
        best$mm_species <- a2t$Species[idx]; best$mm_genus <- a2t$Genus[idx]
        best$mm_family  <- a2t$Family[idx]

        # queries = ASV ids in repseqs (parse ">ASV....")
        qids <- sub("^>", "", grep("^>", readLines(repseqs), value = TRUE))
        N_Q  <- length(qids)
        cc <- data.frame(asv_id = qids, stringsAsFactors = FALSE)
        cc$tax_species <- tax$Species[match(cc$asv_id, tax$asv_id)]
        cc$tax_genus   <- tax$Genus[match(cc$asv_id, tax$asv_id)]
        bi <- match(cc$asv_id, best$query)
        cc$mm_species <- best$mm_species[bi]; cc$mm_genus <- best$mm_genus[bi]
        cc$mm_pident  <- best$pident[bi];     cc$mm_bits  <- best$bits[bi]

        tax_asg <- !is_unassigned(cc$tax_species); mm_asg <- !is.na(cc$mm_species) & !is_unassigned(cc$mm_species)
        cc$category <- with(cc, ifelse(!tax_asg & mm_asg, "newly_resolved",
                                ifelse(tax_asg & !mm_asg, "lost_no_hit",
                                ifelse(!tax_asg & !mm_asg, "both_unassigned",
                                ifelse(tax_species == mm_species, "agree_species",
                                ifelse(tax_genus == mm_genus, "genus_agree_species_differ", "changed_genus"))))))
        write.table(cc, P("data", "derived", "reclassB_mmseqs_concordance.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        n_both  <- sum(tax_asg & mm_asg)
        sp_agree_pct <- round(100 * sum(cc$category == "agree_species") / max(1, n_both), 1)
        gn_agree_pct <- round(100 * sum(cc$mm_genus == cc$tax_genus, na.rm = TRUE) / max(1, n_both), 1)
        summ <- data.frame(
          metric = c("n_query_ASV", "n_with_mmseqs_hit", "n_taxonomy_assigned",
                     "n_both_assigned", "pct_species_agree_of_both", "pct_genus_agree_of_both",
                     "n_newly_resolved", "n_lost_no_hit", "n_changed_genus",
                     "n_genus_agree_species_differ", "asv0001_taxonomy", "asv0001_mmseqs"),
          value  = c(N_Q, sum(mm_asg), sum(tax_asg), n_both, sp_agree_pct, gn_agree_pct,
                     sum(cc$category == "newly_resolved"), sum(cc$category == "lost_no_hit"),
                     sum(cc$category == "changed_genus"), sum(cc$category == "genus_agree_species_differ"),
                     cc$tax_species[cc$asv_id == "ASV0001"], cc$mm_species[cc$asv_id == "ASV0001"]),
          stringsAsFactors = FALSE)
        write.table(summ, P("data", "derived", "reclassB_mmseqs_summary.tsv"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        partB$status <- "ran"; partB$n_query <- N_Q; partB$n_hit <- sum(mm_asg)
        add(sprintf("RESULT      : RAN. %d/%d ASVs hit @97%%; species-agree=%.1f%% genus-agree=%.1f%% of %d both-assigned; newly_resolved=%d; ASV0001: taxonomy=%s -> MMseqs2=%s",
            sum(mm_asg), N_Q, sp_agree_pct, gn_agree_pct, n_both,
            sum(cc$category == "newly_resolved"),
            cc$tax_species[cc$asv_id == "ASV0001"], cc$mm_species[cc$asv_id == "ASV0001"]))
        message(sprintf("PART B: RAN -- %d/%d hit; species-agree %.1f%%; ASV0001 taxonomy=%s -> mmseqs=%s",
            sum(mm_asg), N_Q, sp_agree_pct, cc$tax_species[cc$asv_id == "ASV0001"], cc$mm_species[cc$asv_id == "ASV0001"]))
      }
    }
  }
}
add("result_status: ", partB$status)
writeLines(unlist(prov), P("data", "derived", "reclassB_partB_provenance.txt"))
message("PART B: status = ", partB$status, if (nzchar(partB$reason)) paste0(" (", partB$reason, ")") else "")

# ===========================================================================
# PART C  --  separate reference-resource comparison (MIDORI2 GB271 srRNA / 12S)
# ===========================================================================
# MIDORI2 is a separate current, GenBank-derived reference resource. Concordance
# here describes label agreement across resources; it does not validate biological
# identity, local occurrence, or change through time. GenBank271 = 2026-04-07.
# Gate: runs when the MIDORI2 FASTA is present, or when RECLASSB_FETCH_INDEP=1
# triggers a capped fetch (.zip ~10 MB -> ~67 MB); otherwise degrades to
# documented-only so the script still runs with the DB absent. Reuses Part B helpers
# (find_mmseqs/run/is_unassigned), the loaded `tax` (with total_reads/reads_clean),
# `repseqs`, MM_PARAMS, mmseqs_bin. Never touches env 'edna'; DB stays outside repo.
INDEP_DIR      <- Sys.getenv("RECLASSB_INDEP_DBDIR", DB_DIR)
INDEP_FASTA    <- Sys.getenv("RECLASSB_INDEP_FASTA",
                    file.path(INDEP_DIR, "MIDORI2_LONGEST_NUC_GB271_srRNA_BLAST.fasta"))
INDEP_ZIP_URL  <- "https://www.reference-midori.info/download/Databases/GenBank271_2026-04-07/BLAST/longest/fasta/MIDORI2_LONGEST_NUC_GB271_srRNA_BLAST.fasta.zip"
INDEP_DB_NAME  <- "MIDORI2 LONGEST, GenBank271 (release 2026-04-07), gene srRNA = mitochondrial 12S rRNA; 76,451 species-representative records; NCBI taxonomy"
INDEP_DB_HOME  <- "https://www.reference-midori.info/  (Leray, Knowlton & Machida 2022, Environmental DNA 4:894-907; GenBank-derived)"
DO_FETCH_INDEP <- Sys.getenv("RECLASSB_FETCH_INDEP", "0") == "1"
strip_taxid    <- function(x) sub("_[0-9]+$", "", x)
CACHED_PARTC_FILES <- c(P("data", "derived", "reclassB_indepDB_concordance.tsv"),
                        P("data", "derived", "reclassB_indepDB_summary.tsv"))
HAS_CACHED_PARTC <- all(file.exists(CACHED_PARTC_FILES))
CACHED_PARTC_DATE <- if (HAS_CACHED_PARTC)
  as.character(as.Date(min(file.info(CACHED_PARTC_FILES)$mtime))) else NA_character_

partC <- list(status = "not-run", reason = "", n_query = NA, n_hit = NA)
provC <- character(0); addC <- function(...) provC[[length(provC) + 1L]] <<- paste0(...)
addC("Layer B PART C provenance (separate reference-resource comparison) -- ", Sys.time())
addC("reference DB   : ", INDEP_DB_NAME)
addC("DB home / cite  : ", INDEP_DB_HOME)
addC("DB file URL     : ", INDEP_ZIP_URL)
addC("access-date     : ", ACCESS_DATE)
addC("classifier      : MMseqs2 easy-search  params: ", paste(MM_PARAMS, collapse = " "))
addC("env             : throwaway micromamba env '", MMSEQS_ENV, "' (NOT 'edna').")
addC("SCOPE LIMIT     : MIDORI2 is a separate current reference resource. Label agreement/disagreement measures assignment sensitivity across resources, not biological truth or local identity.")
addC("EXACT fetch (documented, reproducible):")
addC("  mkdir -p ", INDEP_DIR, " && \\")
addC("  curl -sL --max-time 300 --max-filesize 1500000000 -o ", file.path(INDEP_DIR, basename(INDEP_ZIP_URL)), " \\")
addC("    '", INDEP_ZIP_URL, "' && \\")
addC("  ( cd ", INDEP_DIR, " && unzip -o ", basename(INDEP_ZIP_URL), " )")

# capped fetch only if opted in AND absent -------------------------------------
if (!file.exists(INDEP_FASTA) && DO_FETCH_INDEP && !is.na(mmseqs_bin)) {
  dir.create(INDEP_DIR, showWarnings = FALSE, recursive = TRUE)
  zipf <- file.path(INDEP_DIR, basename(INDEP_ZIP_URL))
  message("PART C: fetching MIDORI2 (capped 1.5GB / 300s) ...")
  tryCatch(run("curl", c("-sL", "--max-time", "300", "--max-filesize", "1500000000", "-o", zipf, INDEP_ZIP_URL)),
           error = function(e) message("curl: ", conditionMessage(e)))
  if (file.exists(zipf) && file.info(zipf)$size > 1e5)
    tryCatch(utils::unzip(zipf, exdir = INDEP_DIR), error = function(e) message("unzip: ", conditionMessage(e)))
}
INDEP_PRESENT <- file.exists(INDEP_FASTA) && file.info(INDEP_FASTA)$size > 1e6

if ((!INDEP_PRESENT || is.na(mmseqs_bin)) && HAS_CACHED_PARTC) {
  partC$status <- "cached-results"
  partC$reason <- paste0("MIDORI2 reference/classifier unavailable; retained result tables from ",
                         CACHED_PARTC_DATE, " without recomputation.")
  addC("RESULT         : ", partC$reason)
} else if (is.na(mmseqs_bin)) {
  partC$status <- "documented-only"; partC$reason <- "mmseqs binary not found."
  addC("RESULT         : mmseqs binary NOT found -> documented, not executed.")
} else if (!INDEP_PRESENT) {
  partC$status <- "documented-only"
  partC$reason <- "MIDORI2 reference absent (set RECLASSB_FETCH_INDEP=1 to auto-fetch)."
  addC("RESULT         : MIDORI2 reference absent -> documented-only. Set RECLASSB_FETCH_INDEP=1 to fetch + run.")
} else {
  addC("DB present     : ", INDEP_FASTA, "  (", round(file.info(INDEP_FASTA)$size / 1e6, 1), " MB)")
  # build ref FASTA (synthetic IDs, keep wrapped seq lines) + id->taxonomy map
  refC <- file.path(INDEP_DIR, "reclassC_indep_ref.fasta")
  mapC <- file.path(INDEP_DIR, "reclassC_indep_map.tsv")
  message("PART C: parsing MIDORI2 FASTA + building acc->taxonomy ...")
  L    <- readLines(INDEP_FASTA)
  hdr  <- startsWith(L, ">")
  orig <- sub("^>", "", L[hdr])
  new_ids <- sprintf("MID%06d", seq_along(orig))
  L[hdr]  <- paste0(">", new_ids)
  writeLines(L, refC); rm(L); invisible(gc())
  # MIDORI2 taxonomy = root;kingdom;phylum;class;order;family;genus;species (taxid-suffixed)
  parts <- strsplit(sub("^.*###", "", orig), ";", fixed = TRUE)
  getf  <- function(i) strip_taxid(vapply(parts, function(z) if (length(z) >= i) z[i] else NA_character_, character(1)))
  imap  <- data.frame(id = new_ids, Class = getf(4), Order = getf(5),
                      Family = getf(6), Genus = getf(7), Species = getf(8), stringsAsFactors = FALSE)
  write.table(imap, mapC, sep = "\t", quote = FALSE, row.names = FALSE)

  tmpC <- file.path(INDEP_DIR, "mmtmpC"); unlink(tmpC, recursive = TRUE)
  outC <- file.path(INDEP_DIR, "reclassC_hits.m8")
  message("PART C: MMseqs2 easy-search vs MIDORI2 (", nrow(imap), " refs) ...")
  esC <- run(mmseqs_bin, c("easy-search", repseqs, refC, outC, tmpC,
                           MM_PARAMS, "--threads", "4",
                           "--format-output", "query,target,pident,alnlen,evalue,bits,qcov"))
  if (!file.exists(outC)) {
    partC$status <- "run-failed"; partC$reason <- paste(tail(esC, 2), collapse = " | ")
    addC("RESULT         : MMseqs2 FAILED: ", partC$reason)
  } else {
    hitsC <- tryCatch(read.delim(outC, header = FALSE, stringsAsFactors = FALSE,
               col.names = c("query", "target", "pident", "alnlen", "evalue", "bits", "qcov")),
               error = function(e) data.frame())
    if (!nrow(hitsC)) {
      partC$status <- "run-empty"; partC$reason <- "MMseqs2 returned 0 hits @97%."
      addC("RESULT         : MMseqs2 returned 0 hits @97%.")
    } else {
      hitsC <- hitsC[order(hitsC$query, -hitsC$bits), ]
      bestC <- hitsC[!duplicated(hitsC$query), ]           # best-bit hit per ASV
      j <- match(bestC$target, imap$id)
      bestC$mm_species <- imap$Species[j]; bestC$mm_genus <- imap$Genus[j]; bestC$mm_family <- imap$Family[j]

      qids <- sub("^>", "", grep("^>", readLines(repseqs), value = TRUE))
      cc <- data.frame(asv_id = qids, stringsAsFactors = FALSE)
      ti <- match(cc$asv_id, tax$asv_id)
      cc$total_reads   <- tax$total_reads[ti]; cc$reads_clean <- tax$reads_clean[ti]
      cc$tax_species   <- gsub(" ", "_", tax$Species[ti]); cc$tax_genus <- tax$Genus[ti]
      bi <- match(cc$asv_id, bestC$query)
      cc$indep_species <- bestC$mm_species[bi]; cc$indep_genus <- bestC$mm_genus[bi]
      cc$indep_family  <- bestC$mm_family[bi]
      cc$indep_pident  <- bestC$pident[bi];     cc$indep_bits <- bestC$bits[bi]

      tax_asg <- !is_unassigned(cc$tax_species)
      ind_asg <- !is.na(cc$indep_species) & !is_unassigned(cc$indep_species)
      cc$category <- ifelse(!tax_asg & ind_asg, "newly_resolved",
                     ifelse(tax_asg & !ind_asg, "lost_no_hit",
                     ifelse(!tax_asg & !ind_asg, "both_unassigned",
                     ifelse(cc$tax_species == cc$indep_species, "agree_species",
                     ifelse(cc$tax_genus == cc$indep_genus, "genus_agree_species_differ", "changed_genus")))))
      cc$total_reads[is.na(cc$total_reads)] <- 0; cc$reads_clean[is.na(cc$reads_clean)] <- 0
      write.table(cc, P("data", "derived", "reclassB_indepDB_concordance.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)

      both   <- tax_asg & ind_asg
      n_both <- sum(both)
      sp_agree_pct <- round(100 * sum(cc$category == "agree_species") / max(1, n_both), 1)
      gn_agree_pct <- round(100 * sum(both & cc$tax_genus == cc$indep_genus, na.rm = TRUE) / max(1, n_both), 1)
      wb  <- cc$total_reads[both]
      wsp <- round(100 * sum(cc$total_reads[both & cc$tax_species == cc$indep_species]) / max(1, sum(wb)), 1)
      wgn <- round(100 * sum(cc$total_reads[both & cc$tax_genus   == cc$indep_genus])   / max(1, sum(wb)), 1)

      # ---- Salmo assignment audit across the separate reference resource ------
      mac     <- cc[cc$tax_species == "Salmo_macrostigma" & !is.na(cc$tax_species), ]
      mac_tab <- sort(table(mac$indep_species[!is.na(mac$indep_species)]), decreasing = TRUE)
      mac_top <- if (length(mac_tab)) names(mac_tab)[1] else NA_character_
      mac_to_ts       <- sum(mac$indep_species %in% c("Salmo_trutta", "Salmo_salar"))
      mac_genus_salmo <- sum(mac$indep_genus == "Salmo", na.rm = TRUE)

      nr <- cc[cc$category == "newly_resolved", ]
      summ <- data.frame(
        metric = c("independent_DB", "n_query_ASV", "n_with_indepDB_hit", "n_taxonomy_assigned",
                   "n_both_assigned", "pct_species_agree_of_both", "pct_genus_agree_of_both",
                   "pct_species_agree_readwt", "pct_genus_agree_readwt",
                   "n_newly_resolved", "n_newly_resolved_distinct_taxa", "reads_newly_resolved",
                   "n_lost_no_hit", "n_changed_genus", "n_genus_agree_species_differ",
                   "asv0001_taxonomy", "asv0001_indepDB", "asv0001_pident",
                   "macrostigma_n_asv", "macrostigma_reads", "macrostigma_top_indepDB_call",
                   "macrostigma_to_trutta_or_salar", "macrostigma_genus_Salmo_kept"),
        value  = c(INDEP_DB_NAME, length(qids), sum(ind_asg), sum(tax_asg),
                   n_both, sp_agree_pct, gn_agree_pct, wsp, wgn,
                   sum(cc$category == "newly_resolved"), length(unique(nr$indep_species)),
                   sum(nr$total_reads),
                   sum(cc$category == "lost_no_hit"), sum(cc$category == "changed_genus"),
                   sum(cc$category == "genus_agree_species_differ"),
                   cc$tax_species[cc$asv_id == "ASV0001"], cc$indep_species[cc$asv_id == "ASV0001"],
                   cc$indep_pident[cc$asv_id == "ASV0001"],
                   nrow(mac), sum(mac$total_reads), mac_top, mac_to_ts, mac_genus_salmo),
        stringsAsFactors = FALSE)
      write.table(summ, P("data", "derived", "reclassB_indepDB_summary.tsv"),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      partC$status <- "ran"; partC$n_query <- length(qids); partC$n_hit <- sum(ind_asg)

      addC(sprintf("RESULT         : RAN. %d/%d ASVs hit @97%%; species-agree=%.1f%% (readwt %.1f%%), genus-agree=%.1f%% (readwt %.1f%%) of %d both-assigned.",
           sum(ind_asg), length(qids), sp_agree_pct, wsp, gn_agree_pct, wgn, n_both))
      addC(sprintf("Salmo audit    : MitoHelper macrostigma-labelled ASVs=%d (%s reads); MIDORI2 top label=%s; trutta/salar labels=%d; genus Salmo label=%d/%d. ASV0001 label=%s (%s%%).",
           nrow(mac), format(sum(mac$total_reads), big.mark = ","), mac_top, mac_to_ts, mac_genus_salmo, nrow(mac),
           cc$indep_species[cc$asv_id == "ASV0001"], cc$indep_pident[cc$asv_id == "ASV0001"]))
      addC(sprintf("Baseline-unassigned labels: %d ASVs / %d distinct MIDORI2 labels / %s reads. These labels require target-scope and identity validation before biological interpretation.",
           nrow(nr), length(unique(nr$indep_species)), format(sum(nr$total_reads), big.mark = ",")))
      addC("INTERPRETATION : The genus Salmo label is stable in these comparisons while the species label differs. This supports species-level assignment uncertainty only; no local species, provenance, or native-status identity is inferred.")

      # ---- FIGURE ----
      figC <- P("figures", "reclassB_indepDB.png")
      dir.create(dirname(figC), showWarnings = FALSE, recursive = TRUE)
      pbsp <- NA; pbgn <- NA
      pbf <- P("data", "derived", "reclassB_mmseqs_summary.tsv")
      if (file.exists(pbf)) {
        pb <- read.delim(pbf, stringsAsFactors = FALSE)
        pbsp <- suppressWarnings(as.numeric(pb$value[pb$metric == "pct_species_agree_of_both"]))
        pbgn <- suppressWarnings(as.numeric(pb$value[pb$metric == "pct_genus_agree_of_both"]))
      }
      grDevices::png(figC, width = 1500, height = 760, res = 130)
      op <- par(mfrow = c(1, 2), mar = c(7.2, 4.4, 3.4, 1.1), xpd = NA)
      # Panel A: species vs genus agreement, Part B vs the separate Part-C resource
      mA <- rbind(`species agree` = c(pbsp, sp_agree_pct),
                  `genus agree`   = c(pbgn, gn_agree_pct))
      colnames(mA) <- c("Part B\nsame-vintage\nMitoFish", "Part C\nseparate resource\nMIDORI2")
      bpA <- barplot(mA, beside = TRUE, ylim = c(0, 116), col = c("#E65100", "#2E7D32"),
                     ylab = "agreement vs MitoHelper (% of both-assigned)", las = 1,
                     main = "A. Cross-reference agreement by taxonomic rank", cex.main = 0.95)
      for (k in 1:2) text(bpA[, k], ifelse(is.na(mA[, k]), 3, mA[, k]) + 3,
                          ifelse(is.na(mA[, k]), "n/a", sprintf("%.0f%%", mA[, k])),
                          cex = 0.72, col = c("#E65100", "#1B5E20"))
      legend("top", c("species agree", "genus agree"), fill = c("#E65100", "#2E7D32"),
             bty = "n", cex = 0.8, horiz = TRUE, inset = -0.02)
      mtext(sprintf("read-wt: species %.0f%%, genus %.0f%%  |  %d/%d both-assigned",
            wsp, wgn, n_both, length(qids)), side = 1, line = 5.0, cex = 0.62, col = "#555555")
      # Panel B: MIDORI2 labels for ASVs carrying the MitoHelper macrostigma label
      mb <- rev(head(mac_tab, 6))
      bcolB <- rep("#5C6BC0", length(mb))
      bpB <- barplot(as.numeric(mb), horiz = TRUE, names.arg = names(mb), las = 1, cex.names = 0.72,
                     col = bcolB, border = NA, xlab = "n macrostigma-labelled ASVs",
                     main = "B. MIDORI2 labels for MitoHelper macrostigma calls", cex.main = 0.9,
                     xlim = c(0, max(as.numeric(mb)) * 1.16))
      text(as.numeric(mb) + max(as.numeric(mb)) * 0.01, bpB, as.numeric(mb), pos = 4, cex = 0.7)
      mtext(sprintf("MitoHelper: macrostigma | MitoFish/MMseqs2 top: trutta | MIDORI2 top: %s", mac_top),
            side = 1, line = 4.6, cex = 0.6, col = "#555555")
      mtext(sprintf("genus Salmo label %d/%d; trutta/salar labels=%d -- no local species or status identity inferred",
            mac_genus_salmo, nrow(mac), mac_to_ts), side = 1, line = 5.5, cex = 0.6, col = "#555555")
      par(op); grDevices::dev.off()
      message("PART C: figure -> ", figC)
      message(sprintf("PART C: RAN -- %d/%d hit; species-agree %.1f%% genus-agree %.1f%%; MIDORI2 top label for macrostigma-labelled ASVs=%s (trutta/salar labels=%d).",
              sum(ind_asg), length(qids), sp_agree_pct, gn_agree_pct, mac_top, mac_to_ts))

      # ---- Part C self-check ----
      stopifnot(
        "PartC: both-assigned > 0"            = n_both > 0,
        "PartC: genus-agree >= species-agree" = gn_agree_pct >= sp_agree_pct,
        "PartC: agreements in [0,100]"        = all(c(sp_agree_pct, gn_agree_pct, wsp, wgn) >= 0 &
                                                    c(sp_agree_pct, gn_agree_pct, wsp, wgn) <= 100),
        "PartC: macrostigma ASVs present"     = nrow(mac) > 0,
        "PartC: indepDB summary written"      = file.exists(P("data", "derived", "reclassB_indepDB_summary.tsv"))
      )
    }
  }
}
addC("result_status: ", partC$status)
writeLines(unlist(provC), P("data", "derived", "reclassB_partC_provenance.txt"))
message("PART C: status = ", partC$status, if (nzchar(partC$reason)) paste0(" (", partC$reason, ")") else "")

# Part-C bullet for the Gate-4 addendum (adapts to run status) -----------------
pc_line <- if (partC$status %in% c("ran", "cached-results")) {
  s  <- read.delim(P("data", "derived", "reclassB_indepDB_summary.tsv"), stringsAsFactors = FALSE)
  gv <- function(m) s$value[s$metric == m]
  result_status <- if (identical(partC$status, "ran")) "RAN" else
    sprintf("cached results from %s; not recomputed %s", CACHED_PARTC_DATE, ACCESS_DATE)
  sprintf("- **Separate reference-resource comparison (Part C — %s):** %s; home %s. MMseqs2 @97%%: species-agree **%s%%**, genus-agree **%s%%** of %s both-assigned (read-wt species %s%%, genus %s%%). ASVs labelled `Salmo_macrostigma` by MitoHelper receive the top MIDORI2 label **%s** (*trutta/salar* labels = %s of %s), while genus *Salmo* is retained. This is evidence of species-label instability, not a corrected local species or native-status identity. Baseline-unassigned ASVs receiving MIDORI2 labels = %s ASVs / %s labels; these require target-scope and identity validation. See `reclassB_indepDB_summary.tsv`, `figures/reclassB_indepDB.png`.",
    result_status, gv("independent_DB"), INDEP_DB_HOME,
    gv("pct_species_agree_of_both"), gv("pct_genus_agree_of_both"),
    gv("n_both_assigned"), gv("pct_species_agree_readwt"), gv("pct_genus_agree_readwt"),
    gv("macrostigma_top_indepDB_call"), gv("macrostigma_to_trutta_or_salar"), gv("macrostigma_n_asv"),
    gv("n_newly_resolved"), gv("n_newly_resolved_distinct_taxa"))
} else {
  sprintf("- **Separate reference-resource comparison (Part C — %s):** %s. Re-run with `RECLASSB_FETCH_INDEP=1` to fetch (%s) and execute. %s",
    partC$status, INDEP_DB_NAME, INDEP_ZIP_URL, partC$reason)
}

# ---- emit Gate-4 addendum for docs/prereg.md §5 (do NOT edit prereg.md here) --
writeLines(c(
"## 5. Reference-DB reclassification (Layer B) -- Gate-4 addendum (auto-emitted by R/07_reclass_B.R)",
"",
sprintf("- **Baseline decision (Gate 4, locked, option iii):** the Sep-2020 MitoFish build is unarchived, so Layer B is **2026 concordance, not a 2020->2026 replay**. Recorded %s.", ACCESS_DATE),
"- **Databases + exact versions + access dates:**",
"  - *Baseline in hand* — MitoHelper **12S Nov-2025** first-hit assignment (`data/derived/taxonomy.tsv`; 886/1443 ASVs assigned at species).",
sprintf("  - *2nd DB* — %s; DOI %s; file `%s`; accessed **%s**.", DB_NAME, DB_DOI, DB_TSV_NAME, ACCESS_DATE),
"- **Classifier + threshold:** MMseqs2 `easy-search` @ **97%** min-seq-id, `--cov-mode 2 -c 0.9 --search-type 3 -k 13` (k=13 = low-RAM nucleotide index).",
"- **Concordance metric:** best-bit-score hit per ASV; species/genus agreement over both-assigned ASVs; counts of newly-resolved / lost / changed. See `reclassB_mmseqs_summary.tsv`.",
"- **Scope caveat (Part B):** both Part-B assignments use the *same MitoFish vintage* → Part B measures **classifier/method sensitivity (MMseqs2 vs first-hit)**, not database currency. Part C adds MIDORI2 as a separate current reference resource; neither agreement nor disagreement validates biological truth or local occurrence.",
pc_line,
"- **Assignment-audit set (Gate 4):** `Salmo_macrostigma`, `Phoxinus`, `Barbatula`, `Cottus`, `Pungitius_kaibarae`, `Perca_flavescens`, etc. — see `reclassB_drift_catalogue.tsv`. These are review hypotheses, not corrected identities. For the dominant *Salmo* ASVs, genus is stable across the comparisons while species labels differ; no newly observed taxon, local occurrence, or native-status conclusion is made."),
  P("data", "derived", "reclassB_prereg5_addendum.md"))

# ===========================================================================
# SELF-CHECK
# ===========================================================================
stopifnot(
  "assignment rows = 5 ranks"          = nrow(assign_rates) == 5,
  "assigned+unassigned == N_ASV"       = all(assign_rates$n_asv_assigned + assign_rates$n_asv_unassigned == N_ASV),
  "n_asv_assigned monotone non-increasing Class->Species" = all(diff(assign_rates$n_asv_assigned) <= 0),
  "species assigned == 886"            = assign_rates$n_asv_assigned[assign_rates$rank == "Species"] == 886,
  "species unassigned == 557"          = UNASSIGNED_ASV == 557,
  "read %s in [0,100]"                 = all(assign_rates$pct_reads_assigned_full >= 0 & assign_rates$pct_reads_assigned_full <= 100),
  "drift catalogue non-empty"          = nrow(drift) > 0,
  "headline macrostigma present"       = "Salmo macrostigma" %in% drift$db_species,
  "macrostigma is the top-reads assignment" = drift$db_species[which.max(drift$total_reads)] == "Salmo macrostigma",
  "Part B status is a known value"     = partB$status %in% c("ran", "cached-results", "documented-only", "run-failed", "run-empty", "not-run"),
  "Part C status is a known value"     = partC$status %in% c("ran", "cached-results", "documented-only", "run-failed", "run-empty", "not-run"),
  "Part C provenance written"          = file.exists(P("data", "derived", "reclassB_partC_provenance.txt"))
)
message("SELF-CHECK PASSED.")
message(sprintf("DONE. Part A rates+audit+figure written; Part B status = %s; Part C (separate reference resource) status = %s.",
                partB$status, partC$status))
