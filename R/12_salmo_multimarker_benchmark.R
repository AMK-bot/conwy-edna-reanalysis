#!/usr/bin/env Rscript

# Offline, accession-backed multi-marker Salmo benchmark.
#
# Network acquisition is intentionally isolated in
# scripts/fetch_salmo_multimarker.py.  This script reads only the committed
# NCBI cache and the frozen impact-sprint manifest, calculates every metric and
# gate, and renders deterministic tables plus PNG/SVG figures.

options(stringsAsFactors = FALSE, digits = 15, scipen = 999)

MARKER_ORDER <- c("mifish_12S", "coi_5prime", "cytb")
MARKER_LABELS <- c(
  mifish_12S = "MiFish 12S (170 bp)",
  coi_5prime = "COI 5-prime (648 bp)",
  cytb = "cyt-b (1,140 bp)"
)
CLAIM_BOUNDARY <- paste(
  "Submitted-label discrimination only; public names are annotations, not truth.",
  "No local occurrence, stock, ancestry, hybrid status, ASV0001 species identity,",
  "or biological abundance is inferred."
)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

assert_that <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

yes_no <- function(value) ifelse(value, "yes", "no")

read_tsv_character <- function(path) {
  assert_that(file.exists(path), paste("Missing required TSV:", path))
  utils::read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    colClasses = "character",
    na.strings = NULL
  )
}

read_fasta <- function(path) {
  assert_that(file.exists(path), paste("Missing required FASTA:", path))
  lines <- readLines(path, warn = FALSE)
  header_rows <- grep("^>", lines)
  assert_that(length(header_rows) > 0L, "Marker FASTA has no records")
  ends <- c(header_rows[-1L] - 1L, length(lines))
  ids <- sub("^>([^[:space:]]+).*$", "\\1", lines[header_rows])
  sequences <- vapply(seq_along(header_rows), function(i) {
    paste0(lines[(header_rows[i] + 1L):ends[i]], collapse = "")
  }, character(1))
  names(sequences) <- ids
  sequences <- toupper(sequences)
  assert_that(!anyDuplicated(names(sequences)), "Duplicate FASTA sequence IDs")
  assert_that(all(grepl("^[ACGT]+$", sequences)), "FASTA contains non-ACGT marker bases")
  sequences
}

write_canonical_tsv <- function(table, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    table,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    eol = "\n",
    fileEncoding = "UTF-8"
  )
}

collapse_sorted <- function(values) {
  values <- unique(as.character(values[nzchar(as.character(values))]))
  paste(sort(values), collapse = ";")
}

load_contract <- function(project_root) {
  assert_that(requireNamespace("jsonlite", quietly = TRUE), "Package 'jsonlite' is required")
  path <- file.path(project_root, "data", "metadata", "impact_sprint_manifest.json")
  assert_that(file.exists(path), paste("Missing frozen manifest:", path))
  manifest <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  contract <- manifest$salmo_benchmark
  assert_that(!is.null(contract), "Frozen manifest lacks salmo_benchmark")
  assert_that(identical(as.character(contract$candidate_labels), c(
    "Salmo trutta", "Salmo salar", "Salmo carpio", "Salmo marmoratus",
    "Salmo macedonicus", "Salmo macrostigma", "Salmo obtusirostris"
  )), "Candidate-label contract changed")
  assert_that(identical(as.character(contract$outgroups), c(
    "Salvelinus alpinus", "Oncorhynchus mykiss"
  )), "Outgroup contract changed")
  assert_that(as.integer(contract$max_query_records_per_label) == 50L, "Query cap changed")
  assert_that(as.integer(contract$max_passing_accessions_per_label) == 5L, "Selection cap changed")
  assert_that(as.integer(contract$min_accessions_per_label) == 2L, "Per-label minimum changed")
  assert_that(as.integer(contract$min_salmo_labels_for_species_gate) == 4L, "Species-panel minimum changed")
  assert_that(as.integer(contract$markers$coi_5prime$coding_bases) == 648L, "COI window changed")
  assert_that(as.integer(contract$markers$cytb$coding_bases) == 1140L, "cyt-b window changed")
  assert_that(identical(as.character(contract$classification_tie_rule),
                        "deepest rank shared by every tied nearest neighbour"),
              "Tie rule changed")
  contract
}

validate_and_load_inputs <- function(project_root, contract) {
  assert_that(requireNamespace("digest", quietly = TRUE), "Package 'digest' is required")
  metadata_dir <- file.path(project_root, "data", "metadata")
  public_dir <- file.path(project_root, "data", "public", "salmo_multimarker")
  queries <- read_tsv_character(file.path(metadata_dir, "impact_salmo_queries.tsv"))
  accessions <- read_tsv_character(file.path(metadata_dir, "impact_salmo_accessions.tsv"))
  marker_manifest <- read_tsv_character(file.path(metadata_dir, "impact_salmo_marker_manifest.tsv"))
  sequences <- read_fasta(file.path(public_dir, "marker_sequences.fasta"))

  expected_labels <- c(as.character(contract$candidate_labels), as.character(contract$outgroups))
  assert_that(identical(queries$submitted_label, expected_labels),
              "Query receipt labels/order do not match the frozen contract")
  assert_that(!anyDuplicated(queries$submitted_label), "Duplicate query label")
  assert_that(all(as.integer(queries$returned_count) <= as.integer(contract$max_query_records_per_label)),
              "A query exceeds the 50-record cap")
  assert_that(all(queries$sort_order == "accession"), "Queries were not accession sorted")
  assert_that(all(grepl("eutils.ncbi.nlm.nih.gov", queries$esearch_url, fixed = TRUE)),
              "Query provenance is not official NCBI E-utilities")
  assert_that(length(unique(queries$retrieval_date)) == 1L &&
                grepl("^[0-9]{4}-[0-9]{2}-[0-9]{2}$", queries$retrieval_date[1]),
              "Retrieval date is missing or inconsistent")
  for (i in seq_len(nrow(queries))) {
    receipt_pairs <- list(c("esearch_receipt", "esearch_sha256"))
    if (as.integer(queries$returned_count[i]) > 0L) {
      receipt_pairs <- c(receipt_pairs, list(c("esummary_receipt", "esummary_sha256")))
    } else {
      assert_that(!nzchar(queries$esummary_receipt[i]) && !nzchar(queries$esummary_sha256[i]),
                  "Zero-result query must not claim a synthetic ESummary receipt")
    }
    for (pair in receipt_pairs) {
      receipt <- file.path(project_root, queries[[pair[1]]][i])
      assert_that(file.exists(receipt), paste("Missing query receipt:", receipt))
      observed_hash <- digest::digest(file = receipt, algo = "sha256")
      assert_that(identical(observed_hash, queries[[pair[2]]][i]),
                  paste("Query receipt SHA-256 mismatch:", receipt))
    }
  }

  assert_that(!anyDuplicated(accessions[, c("submitted_label", "candidate_rank")]),
              "Duplicate label/candidate rank")
  for (label in expected_labels) {
    query_row <- queries[queries$submitted_label == label, , drop = FALSE]
    rows <- accessions[accessions$submitted_label == label, , drop = FALSE]
    assert_that(nrow(rows) == as.integer(query_row$returned_count),
                paste("Candidate rows do not equal query return count for", label))
    if (nrow(rows)) {
      assert_that(identical(as.integer(rows$candidate_rank), seq_len(nrow(rows))),
                  paste("Candidate rank is not contiguous for", label))
      selected_ranks <- as.integer(rows$candidate_rank[rows$selected == "yes"])
      assert_that(length(selected_ranks) <= as.integer(contract$max_passing_accessions_per_label),
                  paste("Selection cap exceeded for", label))
      if (length(selected_ranks)) {
        assert_that(!any(rows$assessed == "no" & as.integer(rows$candidate_rank) < max(selected_ranks)),
                    paste("Unassessed candidate precedes a selected accession for", label))
      }
    }
  }
  selected <- accessions[accessions$selected == "yes", , drop = FALSE]
  assert_that(nrow(selected) > 0L, "No accessions passed acquisition")
  assert_that(!anyDuplicated(selected$accession_version), "Selected accession used more than once")
  assert_that(all(selected$submitted_label == selected$organism_annotation),
              "Selected record organism differs from submitted-label panel")
  assert_that(!any(grepl("^(NC_|NM_|NR_|XM_|XR_|NG_)", selected$accession_version)),
              "RefSeq derivative entered the primary-record panel")
  assert_that(!anyDuplicated(selected[, c("submitted_label", "independence_key")]),
              "Duplicate voucher/isolate/independence key in selected panel")
  for (i in seq_len(nrow(selected))) {
    raw_path <- file.path(project_root, selected$raw_genbank_path[i])
    assert_that(file.exists(raw_path), paste("Missing raw GenBank cache:", raw_path))
    observed_hash <- digest::digest(file = raw_path, algo = "sha256")
    assert_that(identical(observed_hash, selected$raw_genbank_sha256[i]),
                paste("Raw GenBank SHA-256 mismatch:", selected$accession_version[i]))
    assert_that(grepl("^https://www.ncbi.nlm.nih.gov/nuccore/", selected$ncbi_record_url[i]),
                "Selected record URL is not official NCBI Nucleotide")
  }

  selected_keys <- paste(selected$submitted_label, selected$accession_version, sep = "\r")
  marker_keys <- paste(marker_manifest$submitted_label, marker_manifest$accession_version, sep = "\r")
  expected_marker_keys <- rep(selected_keys, each = length(MARKER_ORDER))
  assert_that(setequal(marker_keys, expected_marker_keys),
              "Marker manifest does not cover every selected accession")
  assert_that(all(table(marker_keys) == length(MARKER_ORDER)),
              "Selected accession does not have exactly three marker rows")
  assert_that(setequal(marker_manifest$marker, MARKER_ORDER), "Unexpected marker name")
  assert_that(!anyDuplicated(marker_manifest$sequence_id), "Duplicate marker sequence ID")
  assert_that(setequal(marker_manifest$sequence_id, names(sequences)),
              "FASTA/marker-manifest sequence IDs differ")
  marker_manifest$sequence <- unname(sequences[marker_manifest$sequence_id])
  expected_lengths <- c(mifish_12S = NA_integer_, coi_5prime = 648L, cytb = 1140L)
  mifish_min <- as.integer(contract$markers$mifish_12S$min_product_bp)
  mifish_max <- as.integer(contract$markers$mifish_12S$max_product_bp)
  observed_lengths <- nchar(marker_manifest$sequence)
  assert_that(all(observed_lengths == as.integer(marker_manifest$sequence_length)),
              "Marker sequence-length manifest mismatch")
  assert_that(all(observed_lengths[marker_manifest$marker == "mifish_12S"] >= mifish_min &
                    observed_lengths[marker_manifest$marker == "mifish_12S"] <= mifish_max),
              "MiFish product outside frozen 140--250-bp range")
  for (marker in c("coi_5prime", "cytb")) {
    assert_that(all(observed_lengths[marker_manifest$marker == marker] == expected_lengths[[marker]]),
                paste("Wrong marker length for", marker))
  }
  observed_sequence_hashes <- vapply(
    marker_manifest$sequence,
    digest::digest,
    character(1),
    algo = "sha256",
    serialize = FALSE
  )
  assert_that(identical(unname(observed_sequence_hashes), marker_manifest$sequence_sha256),
              "Marker sequence SHA-256 mismatch")
  mifish <- marker_manifest[marker_manifest$marker == "mifish_12S", , drop = FALSE]
  assert_that(all(as.integer(mifish$forward_primer_mismatches) <=
                    as.integer(contract$markers$mifish_12S$max_mismatches_per_primer)),
              "MiFish forward-primer mismatch cap exceeded")
  assert_that(all(as.integer(mifish$reverse_primer_mismatches) <=
                    as.integer(contract$markers$mifish_12S$max_mismatches_per_primer)),
              "MiFish reverse-primer mismatch cap exceeded")

  marker_manifest$sequence_length <- as.integer(marker_manifest$sequence_length)
  marker_manifest <- marker_manifest[order(
    match(marker_manifest$marker, MARKER_ORDER),
    marker_manifest$submitted_label,
    marker_manifest$accession_version
  ), , drop = FALSE]
  rownames(marker_manifest) <- NULL
  list(
    queries = queries,
    accessions = accessions,
    selected = selected,
    marker_manifest = marker_manifest,
    retrieval_date = queries$retrieval_date[1]
  )
}

make_coverage <- function(inputs, contract) {
  rows <- list()
  cursor <- 0L
  min_n <- as.integer(contract$min_accessions_per_label)
  for (marker in MARKER_ORDER) {
    for (i in seq_len(nrow(inputs$queries))) {
      label <- inputs$queries$submitted_label[i]
      role <- inputs$queries$panel_role[i]
      chosen <- inputs$selected[inputs$selected$submitted_label == label, , drop = FALSE]
      marker_rows <- inputs$marker_manifest[
        inputs$marker_manifest$marker == marker & inputs$marker_manifest$submitted_label == label,
        , drop = FALSE
      ]
      n_passing <- nrow(marker_rows)
      independent_n <- length(unique(chosen$independence_key))
      eligible <- n_passing >= min_n && independent_n >= min_n
      status <- if (n_passing == 0L) "no-complete-records" else if (eligible) "eligible" else "below-minimum"
      limitation <- if (status == "eligible") {
        if (all(chosen$independence_documented == "yes")) {
          "At least two distinct voucher/isolate/strain keys documented."
        } else {
          "Distinct primary accessions; some records lack voucher/isolate qualifiers."
        }
      } else if (n_passing == 0L) {
        "No exact-label complete primary mitogenome passed the frozen common-marker screen."
      } else {
        paste0(n_passing, " passing accession(s); at least ", min_n, " required.")
      }
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        marker = marker,
        submitted_label = label,
        submitted_genus = strsplit(label, " ", fixed = TRUE)[[1]][1],
        panel_role = role,
        ncbi_total_count = as.integer(inputs$queries$ncbi_total_count[i]),
        query_returned_count = as.integer(inputs$queries$returned_count[i]),
        passing_accessions = n_passing,
        unique_independence_keys = independent_n,
        minimum_required = min_n,
        eligible_for_rank_gate = yes_no(eligible),
        coverage_status = status,
        selected_accessions = collapse_sorted(marker_rows$accession_version),
        limitation = limitation,
        stringsAsFactors = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  result[order(match(result$marker, MARKER_ORDER), match(result$submitted_label,
                                                        inputs$queries$submitted_label)), , drop = FALSE]
}

make_pairwise <- function(marker_manifest) {
  rows <- list()
  cursor <- 0L
  for (marker in MARKER_ORDER) {
    panel <- marker_manifest[marker_manifest$marker == marker, , drop = FALSE]
    panel <- panel[order(panel$submitted_label, panel$accession_version), , drop = FALSE]
    distance_matrix <- utils::adist(panel$sequence, panel$sequence, partial = FALSE, ignore.case = FALSE)
    pairs <- utils::combn(seq_len(nrow(panel)), 2L)
    for (column in seq_len(ncol(pairs))) {
      left <- pairs[1L, column]
      right <- pairs[2L, column]
      edit <- as.integer(distance_matrix[left, right])
      denominator <- max(panel$sequence_length[left], panel$sequence_length[right])
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        marker = marker,
        sequence_id_1 = panel$sequence_id[left],
        accession_1 = panel$accession_version[left],
        submitted_label_1 = panel$submitted_label[left],
        genus_1 = panel$submitted_genus[left],
        panel_role_1 = panel$panel_role[left],
        sequence_id_2 = panel$sequence_id[right],
        accession_2 = panel$accession_version[right],
        submitted_label_2 = panel$submitted_label[right],
        genus_2 = panel$submitted_genus[right],
        panel_role_2 = panel$panel_role[right],
        edit_distance = edit,
        normalization_length = denominator,
        normalized_edit_distance = edit / denominator,
        exact_haplotype = yes_no(edit == 0L),
        cross_submitted_label = yes_no(panel$submitted_label[left] != panel$submitted_label[right]),
        cross_genus = yes_no(panel$submitted_genus[left] != panel$submitted_genus[right]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_haplotypes <- function(marker_manifest, coverage) {
  rows <- list()
  cursor <- 0L
  for (marker in MARKER_ORDER) {
    panel <- marker_manifest[marker_manifest$marker == marker, , drop = FALSE]
    eligible_salmo <- coverage$submitted_label[
      coverage$marker == marker & coverage$panel_role == "salmo_candidate" &
        coverage$eligible_for_rank_gate == "yes"
    ]
    groups <- split(seq_len(nrow(panel)), panel$sequence)
    groups <- groups[order(names(groups))]
    for (index in seq_along(groups)) {
      members <- panel[groups[[index]], , drop = FALSE]
      labels <- sort(unique(members$submitted_label))
      genera <- sort(unique(members$submitted_genus))
      eligible_member_labels <- intersect(labels, eligible_salmo)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        marker = marker,
        haplotype_id = sprintf("%s_H%03d", marker, index),
        n_accessions = nrow(members),
        n_submitted_labels = length(labels),
        n_genera = length(genera),
        accession_versions = collapse_sorted(members$accession_version),
        submitted_labels = collapse_sorted(labels),
        genera = collapse_sorted(genera),
        cross_submitted_label = yes_no(length(labels) > 1L),
        cross_genus = yes_no(length(genera) > 1L),
        crosses_eligible_salmo_labels = yes_no(length(eligible_member_labels) > 1L),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_barcode_gaps <- function(pairwise, coverage, contract) {
  rows <- list()
  cursor <- 0L
  species_min <- as.integer(contract$min_salmo_labels_for_species_gate)
  for (marker in MARKER_ORDER) {
    eligible <- coverage$submitted_label[
      coverage$marker == marker & coverage$panel_role == "salmo_candidate" &
        coverage$eligible_for_rank_gate == "yes"
    ]
    panel_adequate <- length(eligible) >= species_min
    marker_pairs <- pairwise[pairwise$marker == marker, , drop = FALSE]
    for (label in eligible) {
      within <- marker_pairs[
        marker_pairs$submitted_label_1 == label & marker_pairs$submitted_label_2 == label,
        "normalized_edit_distance"
      ]
      between <- marker_pairs[
        (marker_pairs$submitted_label_1 == label & marker_pairs$submitted_label_2 %in% setdiff(eligible, label)) |
          (marker_pairs$submitted_label_2 == label & marker_pairs$submitted_label_1 %in% setdiff(eligible, label)),
        "normalized_edit_distance"
      ]
      assert_that(length(within) > 0L, paste("Eligible label has no within-label pair:", label))
      assert_that(length(between) > 0L, paste("Eligible label has no between-label pair:", label))
      within_max <- max(within)
      between_min <- min(between)
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(
        marker = marker,
        submitted_label = label,
        accessions = coverage$passing_accessions[
          coverage$marker == marker & coverage$submitted_label == label
        ],
        within_pairs = length(within),
        within_mean_distance = mean(within),
        within_max_distance = within_max,
        between_pairs = length(between),
        between_mean_distance = mean(between),
        nearest_between_label_distance = between_min,
        barcode_gap = between_min - within_max,
        positive_barcode_gap = yes_no(between_min - within_max > 0),
        eligible_salmo_labels_in_panel = length(eligible),
        species_panel_coverage_adequate = yes_no(panel_adequate),
        interpretation = if (panel_adequate) {
          "Eligible species-panel barcode gap."
        } else {
          "Narrowed-panel diagnostic only; four-label species coverage gate is unmet."
        },
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

nearest_row <- function(marker, scope, query, references) {
  distances <- as.integer(utils::adist(query$sequence, references$sequence, partial = FALSE,
                                       ignore.case = FALSE)[1L, ])
  normalized <- distances / pmax(query$sequence_length, references$sequence_length)
  nearest <- min(normalized)
  tied <- which(abs(normalized - nearest) <= 1e-12)
  tied_references <- references[tied, , drop = FALSE]
  tied_labels <- unique(tied_references$submitted_label)
  tied_genera <- unique(tied_references$submitted_genus)
  deepest_rank <- if (length(tied_labels) == 1L) {
    "species"
  } else if (length(tied_genera) == 1L) {
    "genus"
  } else {
    "family-or-unresolved"
  }
  predicted <- if (deepest_rank == "species") {
    tied_labels[1]
  } else if (deepest_rank == "genus") {
    tied_genera[1]
  } else {
    "Salmonidae-or-unresolved"
  }

  min_group <- function(genus, default = NA_real_) {
    values <- normalized[references$submitted_genus == genus]
    if (length(values)) min(values) else default
  }
  nearest_salmo <- min_group("Salmo")
  nearest_salvelinus <- min_group("Salvelinus")
  nearest_oncorhynchus <- min_group("Oncorhynchus")
  genus_separation <- if (scope == "genus_boundary_panel" &&
                          all(is.finite(c(nearest_salmo, nearest_salvelinus, nearest_oncorhynchus)))) {
    nearest_salmo < nearest_salvelinus && nearest_salmo < nearest_oncorhynchus
  } else {
    NA
  }
  data.frame(
    marker = marker,
    evaluation_scope = scope,
    query_sequence_id = query$sequence_id,
    query_accession = query$accession_version,
    query_submitted_label = query$submitted_label,
    query_genus = query$submitted_genus,
    reference_accessions = nrow(references),
    nearest_edit_distance = min(distances[tied]),
    nearest_normalized_distance = nearest,
    tie_count = length(tied),
    tied_accessions = collapse_sorted(tied_references$accession_version),
    tied_submitted_labels = collapse_sorted(tied_labels),
    tied_genera = collapse_sorted(tied_genera),
    deepest_shared_rank = deepest_rank,
    predicted_annotation = predicted,
    all_ties_same_species = yes_no(length(tied_labels) == 1L),
    species_annotation_correct = yes_no(length(tied_labels) == 1L && tied_labels[1] == query$submitted_label),
    all_ties_same_genus = yes_no(length(tied_genera) == 1L),
    genus_annotation_correct = yes_no(length(tied_genera) == 1L && tied_genera[1] == query$submitted_genus),
    nearest_other_salmo_distance = nearest_salmo,
    nearest_salvelinus_distance = nearest_salvelinus,
    nearest_oncorhynchus_distance = nearest_oncorhynchus,
    salmo_nearer_than_each_outgroup = if (is.na(genus_separation)) "" else yes_no(genus_separation),
    stringsAsFactors = FALSE
  )
}

make_loo <- function(marker_manifest, coverage) {
  rows <- list()
  cursor <- 0L
  for (marker in MARKER_ORDER) {
    panel <- marker_manifest[marker_manifest$marker == marker, , drop = FALSE]
    eligible_salmo_labels <- coverage$submitted_label[
      coverage$marker == marker & coverage$panel_role == "salmo_candidate" &
        coverage$eligible_for_rank_gate == "yes"
    ]
    eligible_outgroups <- coverage$submitted_label[
      coverage$marker == marker & coverage$panel_role == "outgroup" &
        coverage$eligible_for_rank_gate == "yes"
    ]
    salmo_panel <- panel[panel$submitted_label %in% eligible_salmo_labels, , drop = FALSE]
    genus_panel <- panel[panel$submitted_label %in% c(eligible_salmo_labels, eligible_outgroups), , drop = FALSE]
    if (nrow(salmo_panel) > 1L) {
      for (i in seq_len(nrow(salmo_panel))) {
        cursor <- cursor + 1L
        rows[[cursor]] <- nearest_row(
          marker,
          "species_label_panel",
          salmo_panel[i, , drop = FALSE],
          salmo_panel[-i, , drop = FALSE]
        )
      }
      if (length(eligible_outgroups) == 2L) {
        for (i in seq_len(nrow(salmo_panel))) {
          query_id <- salmo_panel$sequence_id[i]
          cursor <- cursor + 1L
          rows[[cursor]] <- nearest_row(
            marker,
            "genus_boundary_panel",
            salmo_panel[i, , drop = FALSE],
            genus_panel[genus_panel$sequence_id != query_id, , drop = FALSE]
          )
        }
      }
    }
  }
  assert_that(length(rows) > 0L, "No eligible leave-one-out queries")
  do.call(rbind, rows)
}

make_gates <- function(coverage, haplotypes, barcode_gaps, loo, contract) {
  rows <- list()
  species_min <- as.integer(contract$min_salmo_labels_for_species_gate)
  min_per_label <- as.integer(contract$min_accessions_per_label)
  required_outgroups <- as.character(contract$outgroups)
  required_accuracy <- as.numeric(contract$species_gate_accuracy)
  for (marker in MARKER_ORDER) {
    marker_coverage <- coverage[coverage$marker == marker, , drop = FALSE]
    candidate_salmo <- marker_coverage$submitted_label[
      marker_coverage$panel_role == "salmo_candidate"
    ]
    eligible_salmo <- marker_coverage$submitted_label[
      marker_coverage$panel_role == "salmo_candidate" &
        marker_coverage$eligible_for_rank_gate == "yes"
    ]
    eligible_outgroups <- marker_coverage$submitted_label[
      marker_coverage$panel_role == "outgroup" &
        marker_coverage$eligible_for_rank_gate == "yes"
    ]
    species_coverage <- length(eligible_salmo) >= species_min
    species_loo <- loo[loo$marker == marker & loo$evaluation_scope == "species_label_panel", , drop = FALSE]
    species_accuracy <- if (nrow(species_loo)) {
      mean(species_loo$species_annotation_correct == "yes")
    } else {
      NA_real_
    }
    all_ties_species <- nrow(species_loo) > 0L && all(species_loo$species_annotation_correct == "yes")
    cross_exact <- sum(
      haplotypes$marker == marker & haplotypes$crosses_eligible_salmo_labels == "yes"
    )
    marker_gaps <- barcode_gaps[barcode_gaps$marker == marker, , drop = FALSE]
    all_positive_gaps <- nrow(marker_gaps) == length(eligible_salmo) &&
      nrow(marker_gaps) > 0L && all(marker_gaps$positive_barcode_gap == "yes")
    species_pass <- species_coverage && all_ties_species &&
      is.finite(species_accuracy) && species_accuracy >= required_accuracy &&
      cross_exact == 0L && all_positive_gaps

    outgroup_counts <- marker_coverage$passing_accessions[
      match(required_outgroups, marker_coverage$submitted_label)
    ]
    genus_coverage <- length(eligible_salmo) >= 2L &&
      all(required_outgroups %in% eligible_outgroups) &&
      all(outgroup_counts >= min_per_label)
    genus_loo <- loo[loo$marker == marker & loo$evaluation_scope == "genus_boundary_panel", , drop = FALSE]
    genus_separation <- nrow(genus_loo) > 0L &&
      all(genus_loo$salmo_nearer_than_each_outgroup == "yes")
    genus_pass <- !species_pass && genus_coverage && genus_separation

    finest <- if (species_pass) {
      "species"
    } else if (genus_pass) {
      "genus"
    } else if (!species_coverage && !genus_coverage) {
      "coverage-inadequate"
    } else {
      "family-or-unresolved"
    }
    species_status <- if (!species_coverage) "coverage-inadequate" else if (species_pass) "pass" else "fail"
    genus_status <- if (!genus_coverage) "coverage-inadequate" else if (genus_pass) "pass" else "fail"
    species_reason <- paste0(
      length(eligible_salmo), " of ", length(candidate_salmo),
      " candidate Salmo labels eligible; minimum ", species_min, " required; ",
      if (is.finite(species_accuracy)) sprintf("all-ties LOO accuracy %.3f", species_accuracy) else "LOO undefined",
      "; cross-label exact haplotypes ", cross_exact,
      "; all positive barcode gaps ", yes_no(all_positive_gaps), "."
    )
    genus_reason <- paste0(
      length(eligible_salmo), " eligible Salmo labels; outgroups ",
      paste(required_outgroups, collapse = " + "), "; ",
      sum(genus_loo$salmo_nearer_than_each_outgroup == "yes"), "/", nrow(genus_loo),
      " Salmo queries strictly nearer to Salmo than to each outgroup genus."
    )
    rows[[length(rows) + 1L]] <- data.frame(
      marker = marker,
      marker_label = unname(MARKER_LABELS[[marker]]),
      passing_salmo_accessions = sum(marker_coverage$passing_accessions[
        marker_coverage$panel_role == "salmo_candidate"
      ]),
      candidate_salmo_labels = length(candidate_salmo),
      eligible_salmo_labels = length(eligible_salmo),
      required_salmo_labels_for_species = species_min,
      species_panel_coverage_adequate = yes_no(species_coverage),
      species_loo_queries = nrow(species_loo),
      species_loo_all_ties_accuracy = species_accuracy,
      all_loo_ties_same_submitted_species = yes_no(all_ties_species),
      cross_label_exact_haplotypes = cross_exact,
      all_eligible_labels_positive_barcode_gap = yes_no(all_positive_gaps),
      species_gate_status = species_status,
      species_gate_reason = species_reason,
      eligible_outgroup_labels = length(eligible_outgroups),
      genus_panel_coverage_adequate = yes_no(genus_coverage),
      genus_loo_queries = nrow(genus_loo),
      salmo_queries_nearer_than_each_outgroup = sum(
        genus_loo$salmo_nearer_than_each_outgroup == "yes"
      ),
      all_salmo_queries_nearer_than_each_outgroup = yes_no(genus_separation),
      genus_gate_status = genus_status,
      genus_gate_reason = genus_reason,
      finest_consistently_supported_rank = finest,
      claim_boundary = CLAIM_BOUNDARY,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

short_species <- function(labels) {
  vapply(strsplit(labels, " ", fixed = TRUE), function(parts) {
    paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
  }, character(1))
}

distance_groups <- function(pairwise, coverage, marker) {
  eligible_salmo <- coverage$submitted_label[
    coverage$marker == marker & coverage$panel_role == "salmo_candidate" &
      coverage$eligible_for_rank_gate == "yes"
  ]
  p <- pairwise[pairwise$marker == marker, , drop = FALSE]
  within <- p$normalized_edit_distance[
    p$submitted_label_1 == p$submitted_label_2 & p$submitted_label_1 %in% eligible_salmo
  ]
  between <- p$normalized_edit_distance[
    p$submitted_label_1 %in% eligible_salmo & p$submitted_label_2 %in% eligible_salmo &
      p$submitted_label_1 != p$submitted_label_2
  ]
  outgroup <- p$normalized_edit_distance[
    (p$submitted_label_1 %in% eligible_salmo & p$panel_role_2 == "outgroup") |
      (p$submitted_label_2 %in% eligible_salmo & p$panel_role_1 == "outgroup")
  ]
  list(within, between, outgroup)
}

draw_benchmark_figure <- function(path, type, coverage, pairwise, gates) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (type == "png") {
    grDevices::png(path, width = 1900, height = 1280, res = 180, bg = "white", type = "cairo")
  } else if (type == "svg") {
    grDevices::svg(path, width = 10.7, height = 7.15, bg = "white", pointsize = 10,
                   onefile = TRUE, family = "sans")
  } else {
    stop("Unknown figure type", call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::layout(matrix(c(1, 2, 3, 3), nrow = 2L, byrow = TRUE),
                   widths = c(0.95, 1.25), heights = c(1.0, 1.05))
  graphics::par(oma = c(2.7, 0.5, 3.6, 0.5), family = "sans")

  first_marker <- coverage[coverage$marker == MARKER_ORDER[1], , drop = FALSE]
  bar_colors <- ifelse(first_marker$panel_role == "salmo_candidate", "#176B87", "#7A8793")
  graphics::par(mar = c(4.0, 8.0, 2.4, 1.0))
  bar_midpoints <- graphics::barplot(
    rev(first_marker$passing_accessions),
    names.arg = rev(short_species(first_marker$submitted_label)),
    horiz = TRUE,
    las = 1,
    col = rev(bar_colors),
    border = NA,
    xlim = c(0, 5.7),
    axes = FALSE,
    cex.names = 0.77,
    main = "A  Common-accession coverage"
  )
  graphics::axis(1, at = 0:5, cex.axis = 0.8)
  graphics::abline(v = 2, lty = 3, col = "#B16A1B", lwd = 1.4)
  graphics::mtext("passing complete mitogenomes", side = 1, line = 2.5, cex = 0.78)
  graphics::legend(2.75, bar_midpoints[7],
                   legend = c("Salmo candidate", "outgroup", "per-label minimum"),
                   fill = c("#176B87", "#7A8793", NA), border = c(NA, NA, NA),
                   lty = c(NA, NA, 3), col = c(NA, NA, "#B16A1B"), bty = "n", cex = 0.7)

  graphics::par(mar = c(5.0, 4.3, 2.4, 0.8))
  group_list <- list()
  positions <- numeric()
  colors <- character()
  marker_colors <- c(mifish_12S = "#3B82B6", coi_5prime = "#D9772A", cytb = "#258C6D")
  for (category in seq_len(3L)) {
    for (j in seq_along(MARKER_ORDER)) {
      marker <- MARKER_ORDER[j]
      group_list[[length(group_list) + 1L]] <- distance_groups(pairwise, coverage, marker)[[category]]
      positions <- c(positions, category + c(-0.23, 0, 0.23)[j])
      colors <- c(colors, marker_colors[[marker]])
    }
  }
  graphics::boxplot(group_list, at = positions, boxwex = 0.18, xaxt = "n", outline = FALSE,
                    col = colors, border = colors, ylab = "normalised global edit distance",
                    main = "B  Separation in the eligible panel", ylim = c(0, max(unlist(group_list)) * 1.08))
  graphics::axis(1, at = 1:3,
                 labels = c("within submitted\nSalmo label", "between eligible\nSalmo labels",
                            "Salmo to\noutgroups"), cex.axis = 0.73)
  graphics::abline(h = 0, col = "#4B5563", lwd = 0.7)
  graphics::legend("topleft", legend = unname(MARKER_LABELS), fill = unname(marker_colors),
                   border = NA, bty = "n", cex = 0.66)

  graphics::par(mar = c(2.0, 0.7, 4.0, 1.0))
  row_labels <- c(
    "Species coverage\n(eligible of 7; minimum 4)",
    "All nearest-neighbour ties retain species label",
    "No exact haplotype crosses eligible labels",
    "Positive barcode gap for every eligible label",
    "Both outgroup genera have\nat least 2 accessions",
    "Every Salmo query is closer than each outgroup"
  )
  state <- matrix("pass", nrow = length(row_labels), ncol = length(MARKER_ORDER))
  state[1, ] <- ifelse(gates$species_panel_coverage_adequate == "yes", "pass", "inadequate")
  state[2, ] <- ifelse(gates$all_loo_ties_same_submitted_species == "yes", "pass", "fail")
  state[3, ] <- ifelse(gates$cross_label_exact_haplotypes == 0L, "pass", "fail")
  state[4, ] <- ifelse(gates$all_eligible_labels_positive_barcode_gap == "yes", "pass", "fail")
  state[5, ] <- ifelse(gates$genus_panel_coverage_adequate == "yes", "pass", "inadequate")
  state[6, ] <- ifelse(gates$all_salmo_queries_nearer_than_each_outgroup == "yes", "pass", "fail")
  tile_colors <- c(pass = "#D5F0E6", fail = "#F7D1D1", inadequate = "#F8E4B5")
  graphics::plot.new()
  graphics::plot.window(xlim = c(-1.0, 3.5), ylim = c(6.5, 0.5), xaxs = "i", yaxs = "i")
  for (row in seq_len(nrow(state))) {
    for (column in seq_len(ncol(state))) {
      graphics::rect(column - 0.46, row - 0.43, column + 0.46, row + 0.43,
                     col = tile_colors[[state[row, column]]], border = "white", lwd = 2)
    }
  }
  text_matrix <- rbind(
    paste0(
      gates$eligible_salmo_labels, " of ", gates$candidate_salmo_labels,
      "\n(min ", gates$required_salmo_labels_for_species, ")"
    ),
    sprintf("%.0f%%", 100 * gates$species_loo_all_ties_accuracy),
    paste0(gates$cross_label_exact_haplotypes, " shared"),
    toupper(gates$all_eligible_labels_positive_barcode_gap),
    paste0(gates$eligible_outgroup_labels, "/2"),
    paste0(gates$salmo_queries_nearer_than_each_outgroup, "/", gates$genus_loo_queries)
  )
  for (row in seq_len(nrow(text_matrix))) {
    for (column in seq_len(ncol(text_matrix))) {
      graphics::text(column, row, text_matrix[row, column], cex = 0.76,
                     font = if (row %in% c(1, 6)) 2 else 1, col = "#24313A")
    }
  }
  graphics::text(0.42, seq_along(row_labels), labels = row_labels, adj = c(1, 0.5),
                 cex = 0.62, col = "#1F2937")
  graphics::axis(3, at = seq_along(MARKER_ORDER), labels = unname(MARKER_LABELS),
                 tick = FALSE, cex.axis = 0.73, font.axis = 2, line = 0.1)
  graphics::rect(0.5, 0.5, 3.5, 6.5, border = "#111827", lwd = 1.0)
  graphics::title(main = "C  Frozen claim gates: narrowed species evidence, tested genus boundary",
                  cex.main = 0.92, line = 2.5)

  graphics::mtext("Accession-backed multi-marker resolution benchmark", side = 3, outer = TRUE,
                  line = 2.0, cex = 1.35, font = 2, col = "#143642")
  graphics::mtext(
    "Complete primary mitochondrial records; exact submitted labels; all nearest-neighbour ties retained",
    side = 3, outer = TRUE, line = 0.55, cex = 0.83, col = "#475569"
  )
  rank_text <- paste0(
    paste(unname(MARKER_LABELS), gates$finest_consistently_supported_rank, sep = ": ", collapse = "   |   "),
    "    Submitted-label benchmark only; no local species, stock or occurrence inference."
  )
  graphics::mtext(rank_text, side = 1, outer = TRUE, line = 1.2, cex = 0.72, col = "#334155")
  invisible(path)
}

calculate_benchmark <- function(project_root) {
  contract <- load_contract(project_root)
  inputs <- validate_and_load_inputs(project_root, contract)
  coverage <- make_coverage(inputs, contract)
  pairwise <- make_pairwise(inputs$marker_manifest)
  haplotypes <- make_haplotypes(inputs$marker_manifest, coverage)
  barcode_gaps <- make_barcode_gaps(pairwise, coverage, contract)
  loo <- make_loo(inputs$marker_manifest, coverage)
  gates <- make_gates(coverage, haplotypes, barcode_gaps, loo, contract)

  expected_pairs <- sum(vapply(MARKER_ORDER, function(marker) {
    n <- sum(inputs$marker_manifest$marker == marker)
    choose(n, 2L)
  }, numeric(1)))
  assert_that(nrow(pairwise) == expected_pairs, "Pairwise-distance denominator mismatch")
  assert_that(all(pairwise$normalized_edit_distance >= 0 & pairwise$normalized_edit_distance <= 1),
              "Normalised edit distance outside [0,1]")
  assert_that(all((pairwise$edit_distance == 0L) == (pairwise$exact_haplotype == "yes")),
              "Exact-haplotype flag disagrees with edit distance")
  assert_that(nrow(gates) == length(MARKER_ORDER), "Expected one gate row per marker")
  assert_that(all(gates$claim_boundary == CLAIM_BOUNDARY), "Claim boundary missing from gate table")

  list(
    contract = contract,
    inputs = inputs,
    coverage = coverage,
    pairwise = pairwise,
    haplotypes = haplotypes,
    barcode_gaps = barcode_gaps,
    loo = loo,
    gates = gates
  )
}

write_benchmark <- function(result, output_root) {
  metadata_dir <- file.path(output_root, "data", "metadata")
  derived_dir <- file.path(output_root, "data", "derived")
  figure_dir <- file.path(output_root, "figures")
  write_canonical_tsv(result$coverage,
                      file.path(metadata_dir, "impact_salmo_panel_coverage.tsv"))
  write_canonical_tsv(result$pairwise,
                      file.path(derived_dir, "impact_salmo_pairwise_distances.tsv"))
  write_canonical_tsv(result$haplotypes,
                      file.path(derived_dir, "impact_salmo_exact_haplotypes.tsv"))
  write_canonical_tsv(result$barcode_gaps,
                      file.path(derived_dir, "impact_salmo_barcode_gaps.tsv"))
  write_canonical_tsv(result$loo,
                      file.path(derived_dir, "impact_salmo_loo_classification.tsv"))
  write_canonical_tsv(result$gates,
                      file.path(derived_dir, "impact_salmo_marker_gates.tsv"))
  draw_benchmark_figure(
    file.path(figure_dir, "impact_salmo_multimarker.png"), "png",
    result$coverage, result$pairwise, result$gates
  )
  draw_benchmark_figure(
    file.path(figure_dir, "impact_salmo_multimarker.svg"), "svg",
    result$coverage, result$pairwise, result$gates
  )
  invisible(result)
}

parse_cli <- function(args) {
  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  inferred_root <- if (length(script_arg)) {
    normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."), mustWork = TRUE)
  } else {
    normalizePath(".", mustWork = TRUE)
  }
  parsed <- list(project_root = inferred_root, output_root = inferred_root, validate_only = FALSE)
  i <- 1L
  while (i <= length(args)) {
    if (args[i] == "--project-root") {
      i <- i + 1L
      assert_that(i <= length(args), "--project-root requires a value")
      parsed$project_root <- normalizePath(args[i], mustWork = TRUE)
    } else if (args[i] == "--output-root") {
      i <- i + 1L
      assert_that(i <= length(args), "--output-root requires a value")
      dir.create(args[i], recursive = TRUE, showWarnings = FALSE)
      parsed$output_root <- normalizePath(args[i], mustWork = TRUE)
    } else if (args[i] == "--validate-only") {
      parsed$validate_only <- TRUE
    } else {
      stop(paste("Unknown argument:", args[i]), call. = FALSE)
    }
    i <- i + 1L
  }
  parsed
}

main <- function() {
  started <- proc.time()[["elapsed"]]
  args <- parse_cli(commandArgs(trailingOnly = TRUE))
  result <- calculate_benchmark(args$project_root)
  if (args$validate_only) {
    required <- c(
      file.path(args$output_root, "data", "metadata", "impact_salmo_panel_coverage.tsv"),
      file.path(args$output_root, "data", "derived", "impact_salmo_pairwise_distances.tsv"),
      file.path(args$output_root, "data", "derived", "impact_salmo_exact_haplotypes.tsv"),
      file.path(args$output_root, "data", "derived", "impact_salmo_barcode_gaps.tsv"),
      file.path(args$output_root, "data", "derived", "impact_salmo_loo_classification.tsv"),
      file.path(args$output_root, "data", "derived", "impact_salmo_marker_gates.tsv"),
      file.path(args$output_root, "figures", "impact_salmo_multimarker.png"),
      file.path(args$output_root, "figures", "impact_salmo_multimarker.svg")
    )
    assert_that(all(file.exists(required)), "One or more benchmark outputs are missing")
    cat("Salmo multi-marker benchmark validation PASS\n")
  } else {
    write_benchmark(result, args$output_root)
    cat("Salmo multi-marker benchmark complete\n")
  }
  cat(sprintf("Runtime: %.3f seconds (telemetry only; not written to canonical outputs)\n",
              proc.time()[["elapsed"]] - started))
  invisible(result)
}

if (!identical(Sys.getenv("SALMO_BENCHMARK_SOURCE_ONLY"), "1")) {
  main()
}
