#!/usr/bin/env Rscript

# Deterministically rebuild and validate the cited fish-trait summaries and the
# Salmo molecular-detective figure. Uses base R only.
#
# From the repository root:
#   Rscript scripts/build_biological_evidence.R
#   Rscript scripts/build_biological_evidence.R --validate-only
#   Rscript scripts/build_biological_evidence.R --no-rds   # exercise clean-clone fallback

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) {
  cat(paste0(
    "Usage: Rscript scripts/build_biological_evidence.R [--validate-only] [--no-rds]\n",
    "  default          Rebuild trait summaries when the RDS is present, validate all evidence,\n",
    "                   and regenerate figures/salmo_molecular_detective.{svg,png}.\n",
    "  --validate-only  Perform checks without writing files.\n",
    "  --no-rds         Skip RDS recomputation; validate committed summary structure and\n",
    "                   regenerate the Salmo figure from committed TSV/FASTA inputs.\n"))
  quit(status = 0)
}
unknown <- setdiff(args, c("--validate-only", "--no-rds"))
if (length(unknown)) stop("Unknown argument(s): ", paste(unknown, collapse = ", "))
validate_only <- "--validate-only" %in% args
skip_rds <- "--no-rds" %in% args

# Resolve repository root from the script location, independent of caller cwd.
full_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", full_args, value = TRUE)
if (length(file_arg)) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)
options(stringsAsFactors = FALSE, OutDec = ".")

assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
need <- function(path) assert(file.exists(path), paste("Missing required file:", path))
read_tsv <- function(path, ...) read.delim(path, sep = "\t", quote = "", comment.char = "",
                                           check.names = FALSE, ...)
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE,
              na = "NA", fileEncoding = "UTF-8")
}

paths <- list(
  traits = "data/metadata/fish_traits.tsv",
  trait_sources = "data/metadata/fish_traits_sources.tsv",
  groups = "data/metadata/fish_trait_analysis_groups.tsv",
  species_summary = "data/metadata/fish_species_section_prevalence.tsv",
  guild_summary = "data/metadata/fish_trait_section_summary.tsv",
  counts = "data/derived/species_x_sample_counts_clean.rds",
  covariates = "data/metadata/conwy_covariates.tsv",
  asv_map = "data/derived/asv_map.tsv",
  fasta = "data/metadata/salmo_reference_amplicons.fasta",
  alignment = "data/metadata/salmo_sequence_alignment.tsv",
  reference_provenance = "data/metadata/salmo_reference_provenance.tsv",
  distances = "data/metadata/salmo_edit_distance_matrix.tsv",
  nodes = "data/metadata/salmo_haplotype_nodes.tsv",
  edges = "data/metadata/salmo_haplotype_edges.tsv",
  svg = "figures/salmo_molecular_detective.svg",
  png = "figures/salmo_molecular_detective.png"
)
invisible(lapply(paths[c("traits", "trait_sources", "groups", "species_summary", "guild_summary",
                         "fasta", "alignment", "reference_provenance", "distances", "nodes", "edges")], need))

## ---- trait-table/source integrity ----------------------------------------------------
traits <- read_tsv(paths$traits)
sources <- read_tsv(paths$trait_sources)
groups <- read_tsv(paths$groups, na.strings = character())
assert(nrow(traits) == 14L && length(unique(traits$analysis_label)) == 14L,
       "fish_traits.tsv must contain 14 unique analysis labels")
assert(nrow(groups) == 14L && identical(groups$analysis_label, traits$analysis_label),
       "fish_trait_analysis_groups.tsv must align one-for-one with fish_traits.tsv")
used_source_ids <- unique(unlist(strsplit(traits$source_ids, ";", fixed = TRUE)))
assert(all(used_source_ids %in% sources$source_id),
       paste("Unresolved fish trait source IDs:",
             paste(setdiff(used_source_ids, sources$source_id), collapse = ", ")))
assert(all(sources$accessed == "2026-07-18"), "Trait source access dates are not uniform")

## ---- balanced site x timepoint summaries -------------------------------------------
build_balanced_summaries <- function() {
  counts <- round(readRDS(paths$counts))
  assert(is.matrix(counts), "Species RDS is not a matrix")
  assert(identical(colnames(counts), groups$analysis_label),
         "Species RDS columns do not match the deterministic trait mapping")
  assert(all(as.integer(traits$detected_libraries) == colSums(counts > 0)),
         "fish_traits.tsv detected-library counts do not reproduce from the species RDS")
  cov <- read_tsv(paths$covariates)
  assert("ID" %in% names(cov) && all(rownames(counts) %in% cov$ID),
         "Covariates do not cover every species-matrix library")
  md <- cov[match(rownames(counts), cov$ID), , drop = FALSE]
  md$site <- sub("_.*", "", md$ID)
  md$timepoint <- sub(".*_", "", md$ID)
  md$unit <- paste(md$site, md$timepoint, sep = "_")
  md$section <- factor(
    ifelse(md$site %in% sprintf("E%02d", 1:5), "upper",
      ifelse(md$site %in% sprintf("E%02d", 6:10), "middle", "lower")),
    levels = c("upper", "middle", "lower"))

  unit_counts <- rowsum(counts, md$unit, reorder = TRUE)
  unit_first <- match(rownames(unit_counts), md$unit)
  unit_md <- data.frame(
    site = md$site[unit_first], timepoint = md$timepoint[unit_first],
    section = md$section[unit_first], stringsAsFactors = FALSE)
  grid <- xtabs(~ site + timepoint, data = unit_md) > 0
  common_timepoints <- colnames(grid)[colSums(grid) == nrow(grid)]
  keep <- unit_md$timepoint %in% common_timepoints
  occurrence <- unit_counts[keep, , drop = FALSE] > 0
  section <- factor(unit_md$section[keep], levels = c("upper", "middle", "lower"))
  site <- unit_md$site[keep]
  assert(length(unique(site)) == 14L && length(common_timepoints) == 11L && sum(keep) == 154L,
         "Balanced spatial panel is no longer 14 sites x 11 common timepoints")

  species_rows <- list()
  k <- 1L
  for (i in seq_len(nrow(groups))) {
    taxon <- groups$analysis_label[[i]]
    for (sec in levels(section)) {
      idx <- section == sec
      n_units <- sum(idx)
      n_detect <- sum(occurrence[idx, taxon])
      species_rows[[k]] <- data.frame(
        analysis_label = taxon,
        section = sec,
        n_sites = length(unique(site[idx])),
        n_common_timepoints = length(common_timepoints),
        n_site_time_units = n_units,
        n_units_detected = n_detect,
        unit_prevalence = sprintf("%.4f", n_detect / n_units),
        included_in_trait_summary = groups$include_in_descriptive_summary[[i]],
        interpretation_note = groups$section_prevalence_note[[i]],
        stringsAsFactors = FALSE)
      k <- k + 1L
    }
  }
  species_summary <- do.call(rbind, species_rows)

  group_order <- list(
    flow = c("rheophilic", "eurytopic", "limnophilic"),
    trophic = c("benthic", "small_invertivore", "omnivore_invertivore",
                "ontogenetic_piscivore"))
  metric_text <- c(
    flow_rheophilic = "A unit has the guild when at least one included label in that group was detected; incidences sum label presences and are not reads",
    flow_eurytopic = "A unit has the guild when at least one included label in that group was detected; incidences sum label presences and are not reads",
    flow_limnophilic = "A unit has the guild when at least one included label in that group was detected; incidences sum label presences and are not reads",
    trophic_benthic = "Broad category pools benthic invertivores, omnivores and generalists; unit prevalence is binary and overlapping across trophic groups",
    trophic_small_invertivore = "Unit prevalence is binary; broad trophic groups describe association, not observed diet",
    trophic_omnivore_invertivore = "This group is represented by the Phoxinus complex proxy",
    trophic_ontogenetic_piscivore = "Potential ontogenetic piscivory is a species trait; eDNA does not establish fish size or realised diet")
  guild_rows <- list()
  k <- 1L
  for (axis in names(group_order)) {
    group_column <- paste0(axis, "_summary_group")
    for (guild in group_order[[axis]]) {
      taxon <- groups$analysis_label[
        groups$include_in_descriptive_summary == "yes" & groups[[group_column]] == guild]
      assert(length(taxon) > 0L, paste("No taxa mapped to", axis, guild))
      for (sec in levels(section)) {
        idx <- section == sec
        z <- occurrence[idx, taxon, drop = FALSE]
        present <- rowSums(z) > 0
        key <- paste(axis, guild, sep = "_")
        guild_rows[[k]] <- data.frame(
          trait_axis = axis,
          trait_group = key,
          section = sec,
          n_sites = length(unique(site[idx])),
          n_common_timepoints = length(common_timepoints),
          n_site_time_units = sum(idx),
          n_units_with_guild = sum(present),
          unit_prevalence = sprintf("%.4f", mean(present)),
          taxon_detection_incidences = sum(z),
          metric_definition = unname(metric_text[[key]]),
          stringsAsFactors = FALSE)
        k <- k + 1L
      }
    }
  }
  list(species = species_summary, guild = do.call(rbind, guild_rows),
       common_timepoints = common_timepoints)
}

validate_recomputed_tables <- function(expected) {
  got_species <- read_tsv(paths$species_summary)
  got_guild <- read_tsv(paths$guild_summary)
  assert(nrow(got_species) == nrow(expected$species), "Taxon-section summary row count changed")
  assert(nrow(got_guild) == nrow(expected$guild), "Guild-section summary row count changed")
  numeric_species <- c("n_sites", "n_common_timepoints", "n_site_time_units",
                       "n_units_detected", "unit_prevalence")
  numeric_guild <- c("n_sites", "n_common_timepoints", "n_site_time_units",
                     "n_units_with_guild", "unit_prevalence", "taxon_detection_incidences")
  for (nm in names(expected$species)) {
    if (nm %in% numeric_species) {
      assert(all(abs(as.numeric(got_species[[nm]]) - as.numeric(expected$species[[nm]])) < 5e-5),
             paste("Taxon-section mismatch in", nm))
    } else {
      assert(identical(as.character(got_species[[nm]]), as.character(expected$species[[nm]])),
             paste("Taxon-section mismatch in", nm))
    }
  }
  for (nm in names(expected$guild)) {
    if (nm %in% numeric_guild) {
      assert(all(abs(as.numeric(got_guild[[nm]]) - as.numeric(expected$guild[[nm]])) < 5e-5),
             paste("Guild-section mismatch in", nm))
    } else {
      assert(identical(as.character(got_guild[[nm]]), as.character(expected$guild[[nm]])),
             paste("Guild-section mismatch in", nm))
    }
  }
}

validate_committed_fallback <- function() {
  species <- read_tsv(paths$species_summary)
  guild <- read_tsv(paths$guild_summary)
  assert(nrow(species) == 42L && length(unique(species$analysis_label)) == 14L,
         "Fallback taxon summary is not 14 labels x 3 sections")
  assert(all(abs(species$n_units_detected / species$n_site_time_units -
                   species$unit_prevalence) < 5e-5),
         "Fallback taxon prevalence arithmetic failed")
  assert(nrow(guild) == 21L && all(guild$n_units_with_guild <= guild$n_site_time_units),
         "Fallback guild summary shape/range failed")
  assert(all(abs(guild$n_units_with_guild / guild$n_site_time_units -
                   guild$unit_prevalence) < 5e-5),
         "Fallback guild prevalence arithmetic failed")
  for (i in seq_len(nrow(guild))) {
    axis <- guild$trait_axis[[i]]
    group <- sub(paste0("^", axis, "_"), "", guild$trait_group[[i]])
    column <- paste0(axis, "_summary_group")
    taxa <- groups$analysis_label[
      groups$include_in_descriptive_summary == "yes" & groups[[column]] == group]
    rows <- species$section == guild$section[[i]] & species$analysis_label %in% taxa
    assert(sum(species$n_units_detected[rows]) == guild$taxon_detection_incidences[[i]],
           paste("Fallback guild incidence mismatch:", guild$trait_group[[i]], guild$section[[i]]))
  }
}

can_recompute <- !skip_rds && file.exists(paths$counts) && file.exists(paths$covariates)
if (can_recompute) {
  expected <- build_balanced_summaries()
  if (!validate_only) {
    write_tsv(expected$species, paths$species_summary)
    write_tsv(expected$guild, paths$guild_summary)
  }
  validate_recomputed_tables(expected)
  cat("PASS trait summaries: 42 taxon-section and 21 guild-section rows reproduce from 154 balanced units.\n")
} else {
  validate_committed_fallback()
  cat("PASS trait fallback: committed summary structure, prevalence arithmetic and guild incidences validated (RDS recomputation skipped).\n")
}

## ---- Salmo sequence/distance integrity ---------------------------------------------
read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^>", lines)
  assert(length(starts) > 0L, "No FASTA records found")
  ends <- c(starts[-1L] - 1L, length(lines))
  ids <- sub("^>([^ ]+).*$", "\\1", lines[starts])
  seqs <- vapply(seq_along(starts), function(i) {
    paste(lines[seq.int(starts[[i]] + 1L, ends[[i]])], collapse = "")
  }, character(1))
  setNames(toupper(seqs), ids)
}
hamming <- function(a, b) {
  assert(nchar(a) == nchar(b), "Hamming distance received unequal lengths")
  sum(strsplit(a, "", fixed = TRUE)[[1]] != strsplit(b, "", fixed = TRUE)[[1]])
}

seqs <- read_fasta(paths$fasta)
alignment <- read_tsv(paths$alignment)
reference_provenance <- read_tsv(paths$reference_provenance)
op931658 <- reference_provenance[
  reference_provenance$sequence_id == "OP931658.1", "historical_role_and_limitation"]
assert(length(op931658) == 1L &&
         grepl("record locality differs from study catchment", op931658, fixed = TRUE) &&
         grepl("exact-match label does not establish local identity", op931658, fixed = TRUE) &&
         !grepl("implausible", op931658, ignore.case = TRUE),
       "OP931658 provenance must use neutral locality/identity wording")
assert(length(seqs) == 9L && all(nchar(seqs) == 170L),
       "Expected nine 170-nt Salmo/query FASTA records")
assert(identical(names(seqs), alignment$sequence_id), "FASTA and alignment order/IDs differ")
query <- seqs[["ASV0001"]]
if (file.exists(paths$asv_map)) {
  asv_map <- read_tsv(paths$asv_map)
  asv1 <- asv_map[asv_map$asv_id == "ASV0001", , drop = FALSE]
  assert(nrow(asv1) == 1L && asv1$length[[1]] == 170L && asv1$sequence[[1]] == query,
         "Committed Salmo query does not match ASV0001 in asv_map.tsv")
}
for (i in seq_len(nrow(alignment))) {
  id <- alignment$sequence_id[[i]]
  s <- seqs[[id]]
  qv <- strsplit(query, "", fixed = TRUE)[[1]]
  sv <- strsplit(s, "", fixed = TRUE)[[1]]
  mismatch <- which(qv != sv)
  mismatch_text <- if (!length(mismatch)) "none" else
    paste(sprintf("%d:%s>%s", mismatch, qv[mismatch], sv[mismatch]), collapse = ";")
  assert(identical(s, alignment$aligned_amplicon_5to3[[i]]), paste("Alignment sequence differs:", id))
  assert(length(mismatch) == alignment$edit_distance_to_ASV0001[[i]], paste("Distance differs:", id))
  assert(length(mismatch) == alignment$mismatches[[i]] && 170L - length(mismatch) == alignment$matches_to_ASV0001[[i]],
         paste("Match/mismatch count differs:", id))
  assert(identical(mismatch_text, alignment$mismatch_positions_and_changes_1based[[i]]),
         paste("Mismatch annotation differs:", id))
  assert(abs(100 * (170L - length(mismatch)) / 170 - alignment$identity_percent[[i]]) < 5e-5,
         paste("Identity percentage differs:", id))
}

distance_table <- read_tsv(paths$distances, row.names = 1)
distance_matrix <- as.matrix(distance_table)
storage.mode(distance_matrix) <- "numeric"
expected_distance <- matrix(0, length(seqs), length(seqs), dimnames = list(names(seqs), names(seqs)))
for (i in seq_along(seqs)) for (j in seq_along(seqs))
  expected_distance[i, j] <- hamming(seqs[[i]], seqs[[j]])
assert(identical(dimnames(distance_matrix), dimnames(expected_distance)) &&
         all(distance_matrix == expected_distance), "Salmo edit-distance matrix does not recompute")

nodes <- read_tsv(paths$nodes)
edges <- read_tsv(paths$edges)
assert(nrow(nodes) == 4L && length(unique(nodes$amplicon_5to3)) == 4L,
       "Expected four unique Salmo haplotypes")
node_seq <- setNames(nodes$amplicon_5to3, nodes$haplotype_id)
for (i in seq_len(nrow(edges))) {
  assert(hamming(node_seq[[edges$from_haplotype[[i]]]], node_seq[[edges$to_haplotype[[i]]]]) ==
           edges$edit_distance[[i]], paste("Haplotype edge mismatch on row", i))
}
cat("PASS Salmo evidence: nine 170-nt sequences, symmetric 9x9 distances, four haplotypes and three checked edges.\n")

## ---- deterministic SVG/PNG ----------------------------------------------------------
ellipse <- function(cx, cy, rx, ry, col, border = NA, lwd = 1) {
  theta <- seq(0, 2 * pi, length.out = 241)
  polygon(cx + rx * cos(theta), cy + ry * sin(theta), col = col, border = border, lwd = lwd)
}
short_taxon <- function(x) sub("^Salmo ", "S. ", x)

draw_salmo_figure <- function() {
  navy <- "#14213d"; muted <- "#526174"; grid <- "#dce4ee"; baseline <- "#bcc7d6"
  teal <- "#0b8995"; blue <- "#486f98"; purple <- "#80549f"; red <- "#d1495b"
  par(mar = c(0, 0, 0, 0), xpd = NA, family = "sans")
  plot.new(); plot.window(xlim = c(0, 1600), ylim = c(0, 1000), xaxs = "i", yaxs = "i")
  rect(0, 0, 1600, 1000, col = "white", border = NA)

  text(70, 950, "Salmo molecular detective", adj = c(0, 0.5), cex = 2.15, font = 2, col = navy)
  text(70, 910, "ASV0001: one 170-bp 12S haplotype, five exact public species labels",
       adj = c(0, 0.5), cex = 1.25, col = muted)
  segments(70, 882, 1530, 882, col = "#d5dce6", lwd = 1.5)
  text(70, 842, "A  Substitutions relative to ASV0001", adj = c(0, 0.5), cex = 1.35, font = 2, col = navy)
  text(1035, 842, "B  Haplotype edit-distance network", adj = c(0, 0.5), cex = 1.35, font = 2, col = navy)

  x0 <- 330; x1 <- 950
  pos_x <- function(p) x0 + (p - 1) / 169 * (x1 - x0)
  row_y <- seq(750, 270, length.out = nrow(alignment))
  all_mismatch <- sort(unique(unlist(lapply(seqs, function(s) {
    qv <- strsplit(query, "", fixed = TRUE)[[1]]
    sv <- strsplit(s, "", fixed = TRUE)[[1]]
    which(qv != sv)
  }))))
  segments(x0, 800, x1, 800, col = "#64748b", lwd = 1.5)
  segments(c(x0, x1), 793, c(x0, x1), 807, col = "#64748b", lwd = 1.5)
  text(c(x0, x1), 820, c("1", "170"), cex = 0.82, col = muted)
  text((x0 + x1) / 2, 820, "amplicon position", cex = 0.82, col = muted)
  for (p in all_mismatch) segments(pos_x(p), 245, pos_x(p), 785, col = grid, lwd = 1)
  text(pos_x(all_mismatch), 235, all_mismatch, cex = 0.78, col = muted)

  for (i in seq_len(nrow(alignment))) {
    id <- alignment$sequence_id[[i]]
    label <- if (id == "ASV0001") "ASV0001 query" else
      paste(sub("\\.1$", "", id), short_taxon(alignment$record_taxon[[i]]), sep = " | ")
    text(310, row_y[[i]], label, adj = c(1, 0.5), cex = 0.88, col = "#1f2937")
    segments(x0, row_y[[i]], x1, row_y[[i]], col = baseline, lwd = 3.5)
    d <- alignment$edit_distance_to_ASV0001[[i]]
    if (id != "ASV0001" && d == 0) segments(x0, row_y[[i]], x1, row_y[[i]], col = teal, lwd = 4)
    qv <- strsplit(query, "", fixed = TRUE)[[1]]
    sv <- strsplit(seqs[[id]], "", fixed = TRUE)[[1]]
    mismatch <- which(qv != sv)
    if (length(mismatch)) points(pos_x(mismatch), rep(row_y[[i]], length(mismatch)),
                                 pch = 21, bg = red, col = "white", lwd = 1.2, cex = 1.5)
    right_label <- if (id == "ASV0001") "query" else paste(d, if (d == 1) "edit" else "edits")
    text(970, row_y[[i]], right_label, adj = c(0, 0.5), cex = 0.88,
         col = if (id != "ASV0001" && d == 0) teal else muted)
  }

  node_pos <- rbind(H0 = c(1190, 545), H1 = c(1440, 720), H2 = c(1490, 425), H3 = c(1360, 265))
  for (i in seq_len(nrow(edges))) {
    a <- node_pos[edges$from_haplotype[[i]], ]; b <- node_pos[edges$to_haplotype[[i]], ]
    segments(a[[1]], a[[2]], b[[1]], b[[2]], col = "#64748b", lwd = 3)
  }
  edge_label_pos <- rbind(H0_H1 = c(1335, 655), H0_H2 = c(1370, 505), H1_H3 = c(1415, 485))
  for (i in seq_len(nrow(edges))) {
    key <- paste(edges$from_haplotype[[i]], edges$to_haplotype[[i]], sep = "_")
    xy <- edge_label_pos[key, ]
    rect(xy[[1]] - 30, xy[[2]] - 13, xy[[1]] + 30, xy[[2]] + 13, col = "white", border = NA)
    text(xy[[1]], xy[[2]], paste(edges$edit_distance[[i]],
         if (edges$edit_distance[[i]] == 1) "edit" else "edits"), cex = 0.82, font = 2, col = "#334155")
  }

  ellipse(1190, 545, 145, 125, teal)
  exact_rows <- alignment$sequence_id != "ASV0001" & alignment$edit_distance_to_ASV0001 == 0
  exact_taxa <- unique(short_taxon(alignment$record_taxon[exact_rows]))
  assert(length(exact_taxa) == 5L, "Figure expected five unique exact-match species labels")
  text(1190, 605, "H0 | ASV0001", cex = 1.18, font = 2, col = "white")
  text(1190, 578, "exact public labels:", cex = 0.82, font = 2, col = "white")
  text(1190, 548, paste(exact_taxa[1:2], collapse = " | "), cex = 0.76, col = "white")
  text(1190, 523, paste(exact_taxa[3:4], collapse = " | "), cex = 0.76, col = "white")
  text(1190, 498, exact_taxa[[5]], cex = 0.76, col = "white")
  text(1190, 462, "same 170 bases", cex = 0.82, font = 2, col = "white")

  small_nodes <- c(H1 = blue, H2 = purple, H3 = red)
  for (h in names(small_nodes)) {
    xy <- node_pos[h, ]; ellipse(xy[[1]], xy[[2]], 78, 78, small_nodes[[h]])
    member <- alignment[alignment$aligned_amplicon_5to3 == node_seq[[h]] &
                          alignment$sequence_id != "ASV0001", , drop = FALSE]
    assert(nrow(member) >= 1L, paste("No alignment member for", h))
    text(xy[[1]], xy[[2]] + 25, h, cex = 1.12, font = 2, col = "white")
    text(xy[[1]], xy[[2]] - 2, short_taxon(member$record_taxon[[1]]), cex = 0.74, font = 2, col = "white")
    text(xy[[1]], xy[[2]] - 27, sub("\\.1$", "", member$sequence_id[[1]]), cex = 0.72, font = 2, col = "white")
  }

  rect(70, 65, 1530, 175, col = "#eef7f7", border = "#91c8ce", lwd = 1.5)
  text(100, 145, "Defensible conclusion", adj = c(0, 0.5), cex = 1.27, font = 2, col = navy)
  text(100, 108, "The dominant haplotype is robustly Salmo, but the 170-bp fragment is not species-diagnostic.",
       adj = c(0, 0.5), cex = 0.98, font = 2, col = navy)
  text(100, 82, "Geographic plausibility may guide interpretation; it cannot turn an identical sequence into species-level molecular proof.",
       adj = c(0, 0.5), cex = 0.86, col = muted)
  text(70, 28, "Distances are ungapped substitutions across 170 bp. Network is a minimum-spanning distance diagram, not a phylogeny. NCBI records accessed 2026-07-18; PX646298 and PX646330 postdate the historical MIDORI2 snapshot.",
       adj = c(0, 0.5), cex = 0.70, col = muted)
}

render_salmo <- function() {
  dir.create(dirname(paths$svg), recursive = TRUE, showWarnings = FALSE)
  svg(paths$svg, width = 16, height = 10, pointsize = 12, bg = "white", onefile = FALSE)
  draw_salmo_figure()
  invisible(dev.off())
  svg_lines <- readLines(paths$svg, warn = FALSE, encoding = "UTF-8")
  svg_root <- grep("^<svg[ >]", svg_lines)
  assert(length(svg_root) == 1L, "Rendered SVG root could not be located")
  svg_lines[[svg_root]] <- sub(
    ">$", " role=\"img\" aria-labelledby=\"title desc\">", svg_lines[[svg_root]])
  svg_lines <- append(svg_lines, c(
    "<title id=\"title\">Salmo molecular detective: one 170-base haplotype, five species labels</title>",
    "<desc id=\"desc\">Mismatch map and edit-distance network showing that ASV0001 exactly matches public records carrying five different Salmo species labels.</desc>"),
    after = svg_root)
  writeLines(svg_lines, paths$svg, useBytes = TRUE)
  png_type <- if (capabilities("cairo")) "cairo-png" else getOption("bitmapType")
  png(paths$png, width = 1600, height = 1000, units = "px", res = 100,
      pointsize = 12, bg = "white", type = png_type)
  draw_salmo_figure()
  invisible(dev.off())
  assert(file.info(paths$svg)$size > 1000L && file.info(paths$png)$size > 1000L,
         "Figure render produced an unexpectedly small file")
}

if (!validate_only) {
  render_salmo()
  cat("WROTE figures/salmo_molecular_detective.svg and figures/salmo_molecular_detective.png (1600x1000).\n")
} else {
  need(paths$svg); need(paths$png)
  cat("PASS validate-only: existing Salmo figure files are present; data inputs and summaries passed semantic checks.\n")
}

cat("DONE biological evidence assets.\n")
