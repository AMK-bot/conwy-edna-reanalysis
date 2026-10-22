#!/usr/bin/env Rscript

# Phase-5 Layer-E robustness multiverse.
#
# Contract: docs/impact-sprint-preregistration.md section 5 and
# data/metadata/impact_sprint_manifest.json.  This script never reads an RDS;
# scripts/export_impact_inputs.R is the audited bridge to the committed portable
# library table.  Failed specifications and failed draws remain explicit rows.

options(stringsAsFactors = FALSE, OutDec = ".", digits = 15, scipen = 999)

assert <- function(ok, message) {
  if (length(ok) != 1L || is.na(ok) || !ok) stop(message, call. = FALSE)
}
read_tsv <- function(path, ...) {
  read.delim(path, sep = "\t", quote = "", comment.char = "",
             check.names = FALSE, ...)
}
write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "NA", fileEncoding = "UTF-8")
}
sha256 <- function(path) {
  assert(requireNamespace("digest", quietly = TRUE),
         "The digest package is required for deterministic SHA-256 receipts")
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
normalise_numeric <- function(x, digits = 10L) {
  for (nm in names(x)) {
    if (is.numeric(x[[nm]]) && !is.integer(x[[nm]])) {
      x[[nm]] <- round(x[[nm]], digits)
    }
  }
  x
}

layerE_repo_root <- function() {
  full_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", full_args, value = TRUE)
  if (length(file_arg)) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
    normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
  } else {
    normalizePath(getwd(), mustWork = TRUE)
  }
}

# The complete dihedral series set for n ordered river sites: n-1 non-identity
# rotations plus n rotations of the mirrored order.  Identity is not present.
series_permutation_indices <- function(n = 14L) {
  assert(n >= 3L, "Series permutations require at least three ordered units")
  rotate <- function(x, shift) {
    if (shift == 0L) return(x)
    c(x[(shift + 1L):length(x)], x[seq_len(shift)])
  }
  forward <- do.call(rbind, lapply(seq_len(n - 1L), function(k) rotate(seq_len(n), k)))
  mirrored <- do.call(rbind, lapply(0:(n - 1L), function(k) rotate(rev(seq_len(n)), k)))
  out <- rbind(forward, mirrored)
  storage.mode(out) <- "integer"
  assert(nrow(out) == 2L * n - 1L && ncol(out) == n,
         "Unexpected complete-series permutation dimensions")
  assert(!any(apply(out, 1L, identical, y = seq_len(n))),
         "Identity leaked into non-identity permutation set")
  assert(!anyDuplicated(apply(out, 1L, paste, collapse = ",")),
         "Series permutation set contains a duplicate")
  out
}

make_specifications <- function(manifest) {
  cfg <- manifest$layerE_multiverse
  rows <- list()
  k <- 1L
  for (panel in cfg$panels) {
    for (taxonomy in cfg$taxonomy) {
      for (threshold in cfg$library_detection_threshold_reads) {
        for (effort in cfg$effort) {
          for (distance in cfg$distances) {
            rows[[k]] <- data.frame(
              spec_id = sprintf("E%03d", k),
              panel = panel,
              taxonomy = taxonomy,
              library_detection_threshold_reads = as.integer(threshold),
              library_detection_rule = if (threshold == 1L) ">0 reads" else ">=2 reads",
              effort = effort,
              distance = distance,
              analysis_unit = if (panel == "balanced_site_time")
                "site_x_timepoint_with_whole_site_permutation_blocks" else "site_profile",
              expected_analysis_units = if (panel == "balanced_site_time") 154L else 14L,
              expected_draws = if (effort == "one_library_100_draw_ensemble")
                as.integer(cfg$one_library_draws) else 1L,
              spatial_nonidentity_permutations = as.integer(cfg$spatial_nonidentity_permutations),
              exact_p_denominator_with_observed =
                as.integer(cfg$spatial_nonidentity_permutations + 1L),
              claim_scope = paste(
                "descriptive robustness; replicate libraries are observations within",
                "site x timepoint units; no abundance or confirmed centroid-shift claim"),
              stringsAsFactors = FALSE
            )
            k <- k + 1L
          }
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  assert(nrow(out) == as.integer(cfg$expected_specification_families),
         "Specification grid does not match manifest expected count")
  assert(!anyDuplicated(out[, c("panel", "taxonomy",
                                "library_detection_threshold_reads", "effort", "distance")]),
         "Specification grid contains a duplicate family")
  out
}

load_impact_library <- function(path) {
  assert(file.exists(path), paste("Missing portable input:", path,
                                  "-- run scripts/export_impact_inputs.R"))
  x <- read_tsv(path)
  label_cols <- grep("^analysis_label__", names(x), value = TRUE)
  genus_cols <- grep("^genus__", names(x), value = TRUE)
  assert(nrow(x) == 647L && length(label_cols) == 14L && length(genus_cols) == 12L,
         "Portable impact input must contain 647 libraries, 14 labels, and 12 genera")
  required <- c("library_id", "site", "timepoint", "replicate", "unit_id",
                "section", "season", "surveyed", "usable_timepoint",
                "reprocessed_observed", "structural_absence", "analytical_exclusion",
                "observed_zero", "fish_read_depth")
  assert(all(required %in% names(x)), "Portable impact input is missing required metadata fields")
  assert(!anyDuplicated(x$library_id) &&
         all(x$surveyed) && all(x$usable_timepoint) && all(x$reprocessed_observed),
         "Portable library mask is not 647 unique observed, usable libraries")
  assert(!any(x$structural_absence) && !any(x$analytical_exclusion),
         "Observed portable rows cannot be structural absences or analytical exclusions")
  assert(sum(x$observed_zero) == 18L &&
         identical(as.logical(x$observed_zero), x$fish_read_depth == 0L),
         "Observed-zero flag/depth mismatch")
  for (cols in list(label_cols, genus_cols)) {
    vals <- as.matrix(x[, cols, drop = FALSE])
    storage.mode(vals) <- "numeric"
    assert(all(is.finite(vals)) && all(vals >= 0) && all(vals == round(vals)),
           "Portable counts are not finite non-negative whole reads")
  }
  label_depth <- rowSums(as.matrix(x[, label_cols, drop = FALSE]))
  genus_depth <- rowSums(as.matrix(x[, genus_cols, drop = FALSE]))
  assert(identical(as.numeric(x$fish_read_depth), as.numeric(label_depth)) &&
         identical(as.numeric(label_depth), as.numeric(genus_depth)),
         "Portable label/genus depth identity failed")
  x
}

design_from_library <- function(lib) {
  unit_md <- unique(lib[, c("unit_id", "site", "timepoint", "section"), drop = FALSE])
  unit_md <- unit_md[order(as.integer(sub("E", "", unit_md$site, fixed = TRUE)),
                           as.integer(sub("T", "", unit_md$timepoint, fixed = TRUE))), , drop = FALSE]
  rownames(unit_md) <- NULL
  assert(nrow(unit_md) == 234L && !anyDuplicated(unit_md$unit_id),
         "Expected 234 unique realized site x timepoint units")
  grid <- xtabs(~ site + timepoint, data = unit_md) > 0
  site_order <- sprintf("E%02d", 1:14)
  common_timepoints <- colnames(grid)[colSums(grid) == nrow(grid)]
  assert(identical(rownames(grid), site_order) && length(common_timepoints) == 11L,
         "Expected ordered E01-E14 sites and 11 common timepoints")
  site_section <- setNames(ifelse(seq_along(site_order) <= 5L, "upper",
                                  ifelse(seq_along(site_order) <= 10L,
                                         "middle", "lower")), site_order)
  assert(all(site_section[unit_md$site] == unit_md$section),
         "Portable sections disagree with the frozen 5/5/4 river design")
  list(unit_md = unit_md, site_order = site_order,
       site_section = site_section, common_timepoints = common_timepoints)
}

make_library_draws <- function(lib, design, base_seed, n_draws) {
  split_rows <- split(seq_len(nrow(lib)), factor(lib$unit_id,
                                                 levels = design$unit_md$unit_id))
  assert(length(split_rows) == 234L && all(lengths(split_rows) %in% 1:3),
         "Library draw universe must be 234 units with 1-3 choices")
  rows <- vector("list", n_draws)
  for (draw_id in seq_len(n_draws)) {
    seed <- as.integer(base_seed + draw_id)
    set.seed(seed)
    selected <- vapply(split_rows, function(idx) idx[sample.int(length(idx), 1L)], integer(1))
    rows[[draw_id]] <- data.frame(
      draw_id = as.integer(draw_id), seed = seed,
      unit_id = lib$unit_id[selected], site = lib$site[selected],
      timepoint = lib$timepoint[selected],
      n_available_libraries = as.integer(lengths(split_rows)),
      selected_library_id = lib$library_id[selected],
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  assert(nrow(out) == 234L * n_draws &&
         all(table(out$draw_id) == 234L) &&
         !anyDuplicated(out[, c("draw_id", "unit_id")]),
         "One-library selection table is incomplete")
  out
}

taxonomy_matrix <- function(lib, taxonomy) {
  prefix <- if (taxonomy == "analysis_label") "analysis_label__" else
    if (taxonomy == "genus") "genus__" else stop("Unknown taxonomy: ", taxonomy)
  cols <- grep(paste0("^", prefix), names(lib), value = TRUE)
  x <- as.matrix(lib[, cols, drop = FALSE])
  storage.mode(x) <- "numeric"
  colnames(x) <- sub(paste0("^", prefix), "", cols)
  x
}

aggregate_selected_units <- function(lib, design, counts, selected_library_ids, threshold) {
  idx <- match(selected_library_ids, lib$library_id)
  assert(!anyNA(idx) && !anyDuplicated(idx), "Selected library IDs are absent or duplicated")
  x <- counts[idx, , drop = FALSE]
  x[x < threshold] <- 0
  unit <- lib$unit_id[idx]
  aggregated <- rowsum(x, unit, reorder = FALSE)
  assert(setequal(rownames(aggregated), design$unit_md$unit_id),
         "Selected libraries do not cover every realized unit")
  aggregated <- aggregated[design$unit_md$unit_id, , drop = FALSE]
  assert(identical(rownames(aggregated), design$unit_md$unit_id),
         "Unit aggregation order mismatch")
  aggregated
}

site_rowsum <- function(x, site, site_order) {
  z <- rowsum(x, factor(site, levels = site_order), reorder = FALSE)
  z[site_order, , drop = FALSE]
}

make_panel <- function(unit_counts, design, panel) {
  md <- design$unit_md
  if (panel == "balanced_site_time") {
    keep <- md$timepoint %in% design$common_timepoints
    composition <- unit_counts[keep, , drop = FALSE]
    occurrence <- (composition > 0) * 1
    richness <- rowSums(occurrence)
    row_sites <- md$site[keep]
    row_groups <- unname(design$site_section[row_sites])
    assert(nrow(composition) == 154L && all(table(row_sites) == 11L),
           "Balanced site x time panel is not 14 x 11")
    return(list(
      composition = composition, occurrence = occurrence,
      richness = richness, row_sites = row_sites, row_groups = row_groups,
      n_timepoints = 11L, analysis_unit = "site_x_timepoint",
      profile_definition = "154 observed units; whole-site permutation blocks"
    ))
  }

  if (panel == "balanced_site_profile") {
    keep <- md$timepoint %in% design$common_timepoints
    n_by_site <- table(factor(md$site[keep], levels = design$site_order))
    assert(all(n_by_site == 11L), "Balanced site profiles require 11 units/site")
    composition <- site_rowsum(unit_counts[keep, , drop = FALSE], md$site[keep],
                               design$site_order)
    occurrence <- (composition > 0) * 1
    unit_richness <- rowSums(unit_counts[keep, , drop = FALSE] > 0)
    richness <- as.numeric(rowsum(unit_richness, factor(md$site[keep],
                                                        levels = design$site_order),
                                  reorder = FALSE)[, 1] / as.numeric(n_by_site))
    names(richness) <- design$site_order
    return(list(
      composition = composition, occurrence = occurrence,
      richness = richness, row_sites = design$site_order,
      row_groups = unname(design$site_section[design$site_order]),
      n_timepoints = 11L, analysis_unit = "site_profile",
      profile_definition = "14 site profiles summed over 11 common timepoints"
    ))
  }

  if (panel == "all_site_profile") {
    n_by_site <- table(factor(md$site, levels = design$site_order))
    summed <- site_rowsum(unit_counts, md$site, design$site_order)
    composition <- sweep(summed, 1L, as.numeric(n_by_site), "/")
    occurrence <- sweep(site_rowsum((unit_counts > 0) * 1, md$site,
                                    design$site_order), 1L,
                        as.numeric(n_by_site), "/")
    unit_richness <- rowSums(unit_counts > 0)
    richness <- as.numeric(rowsum(unit_richness, factor(md$site,
                                                        levels = design$site_order),
                                  reorder = FALSE)[, 1] / as.numeric(n_by_site))
    names(richness) <- design$site_order
    return(list(
      composition = composition, occurrence = occurrence,
      richness = richness, row_sites = design$site_order,
      row_groups = unname(design$site_section[design$site_order]),
      n_timepoints = as.integer(max(n_by_site)), analysis_unit = "site_profile",
      profile_definition = paste(
        "14 site profiles over all 234 realized units; mean counts and",
        "occurrence frequencies give each site one row")
    ))
  }
  stop("Unknown panel: ", panel)
}

distance_matrix <- function(panel, distance_name) {
  if (distance_name == "bray_curtis_proportions") {
    x <- panel$composition
    totals <- rowSums(x)
    rel <- x
    nz <- totals > 0
    rel[nz, ] <- sweep(rel[nz, , drop = FALSE], 1L, totals[nz], "/")
    rel[!nz, ] <- 0
    manhattan <- as.matrix(stats::dist(rel, method = "manhattan"))
    denominator <- outer(rowSums(rel), rowSums(rel), "+")
    out <- manhattan / denominator
    out[denominator == 0] <- 0
    diag(out) <- 0
    attr(out, "zero_convention") <-
      "two empty observed communities have distance 0; empty versus detected has distance 1"
    return(out)
  }
  if (distance_name == "jaccard_occurrence") {
    x <- panel$occurrence
    manhattan <- as.matrix(stats::dist(x, method = "manhattan"))
    denominator <- (outer(rowSums(x), rowSums(x), "+") + manhattan) / 2
    out <- manhattan / denominator
    out[denominator == 0] <- 0
    diag(out) <- 0
    attr(out, "zero_convention") <-
      "two empty observed occurrence profiles have distance 0; empty versus detected has distance 1"
    return(out)
  }
  if (distance_name == "robust_aitchison_rclr") {
    x <- panel$composition
    if (any(rowSums(x) == 0)) {
      stop("empty_community_undefined_for_rclr")
    }
    assert(requireNamespace("vegan", quietly = TRUE),
           "The vegan package is required for robust Aitchison rCLR")
    rclr <- vegan::decostand(x, method = "rclr")
    out <- as.matrix(stats::dist(rclr, method = "euclidean"))
    if (any(!is.finite(out))) stop("nonfinite_robust_aitchison_distance")
    return(out)
  }
  stop("Unknown distance: ", distance_name)
}

pseudo_f_stat <- function(distance, group) {
  group <- factor(group, levels = c("upper", "middle", "lower"))
  assert(!anyNA(group) && all(table(group) > 0),
         "Every river section must contain at least one analysis unit")
  n <- nrow(distance)
  k <- nlevels(group)
  assert(n == ncol(distance) && n == length(group) && n > k,
         "Distance/group dimensions are invalid")
  A <- -0.5 * distance^2
  G <- sweep(sweep(A, 1L, rowMeans(A), "-"), 2L, colMeans(A), "-") + mean(A)
  total_ss <- sum(diag(G))
  model_ss <- sum(vapply(levels(group), function(level) {
    idx <- which(group == level)
    sum(G[idx, idx, drop = FALSE]) / length(idx)
  }, numeric(1)))
  residual_ss <- total_ss - model_ss
  df_model <- k - 1L
  df_residual <- n - k
  if (!is.finite(total_ss) || total_ss <= 0 || !is.finite(residual_ss) ||
      residual_ss <= 0) stop("invalid_permanova_sum_of_squares")
  list(
    R2 = model_ss / total_ss,
    F = (model_ss / df_model) / (residual_ss / df_residual),
    df_model = df_model, df_residual = df_residual
  )
}

permutation_groups <- function(row_sites, design, permutations = series_permutation_indices(14L)) {
  site_position <- match(row_sites, design$site_order)
  assert(!anyNA(site_position), "Panel contains an unknown site")
  out <- t(vapply(seq_len(nrow(permutations)), function(i) {
    permuted_sections <- unname(design$site_section[design$site_order])[permutations[i, ]]
    permuted_sections[site_position]
  }, character(length(row_sites))))
  colnames(out) <- NULL
  out
}

exact_permanova <- function(distance, observed_group, permuted_groups) {
  observed <- pseudo_f_stat(distance, observed_group)
  perm_f <- apply(permuted_groups, 1L, function(group) pseudo_f_stat(distance, group)$F)
  assert(length(perm_f) == 27L && all(is.finite(perm_f)),
         "All 27 restricted permutation statistics must be finite")
  exceedances <- sum(perm_f >= observed$F - 1e-12)
  list(R2 = observed$R2, F = observed$F,
       exact_p = (1 + exceedances) / (1 + length(perm_f)),
       permutation_exceedances = as.integer(exceedances),
       permutation_F_min = min(perm_f), permutation_F_max = max(perm_f))
}

dispersion_stats <- function(distance, group) {
  assert(requireNamespace("vegan", quietly = TRUE),
         "The vegan package is required for PERMDISP sensitivity")
  group <- factor(group, levels = c("upper", "middle", "lower"))
  bd <- suppressWarnings(vegan::betadisper(stats::as.dist(distance), group,
                                           type = "centroid", bias.adjust = FALSE))
  fit <- stats::anova(bd)
  means <- tapply(bd$distances, group, mean)
  if (any(!is.finite(bd$distances)) || !is.finite(fit$`F value`[1])) {
    stop("nonfinite_dispersion_statistic")
  }
  list(F = unname(fit$`F value`[1]),
       upper = unname(means[["upper"]]),
       middle = unname(means[["middle"]]),
       lower = unname(means[["lower"]]),
       lower_minus_upper = unname(means[["lower"]] - means[["upper"]]))
}

analyse_panel <- function(panel, distance_name, design, permuted_groups = NULL) {
  section <- factor(panel$row_groups, levels = c("upper", "middle", "lower"))
  richness_means <- tapply(panel$richness, section, mean)
  base <- list(
    n_analysis_units = nrow(panel$composition),
    n_nonzero_community_rows = sum(rowSums(panel$composition) > 0),
    richness_upper_mean = unname(richness_means[["upper"]]),
    richness_middle_mean = unname(richness_means[["middle"]]),
    richness_lower_mean = unname(richness_means[["lower"]]),
    richness_contrast_lower_minus_upper =
      unname(richness_means[["lower"]] - richness_means[["upper"]]),
    pseudo_R2 = NA_real_, pseudo_F = NA_real_, exact_p = NA_real_,
    permutation_exceedances = NA_integer_, permutation_F_min = NA_real_,
    permutation_F_max = NA_real_, dispersion_F = NA_real_,
    dispersion_upper_mean = NA_real_, dispersion_middle_mean = NA_real_,
    dispersion_lower_mean = NA_real_, dispersion_lower_minus_upper = NA_real_,
    valid = FALSE, failure_reason = ""
  )
  tryCatch({
    distance <- distance_matrix(panel, distance_name)
    if (any(!is.finite(distance)) || !isTRUE(all.equal(distance, t(distance),
                                                        tolerance = 1e-10))) {
      stop("nonfinite_or_asymmetric_distance_matrix")
    }
    if (is.null(permuted_groups)) {
      permuted_groups <- permutation_groups(panel$row_sites, design)
    }
    exact <- exact_permanova(distance, section, permuted_groups)
    disp <- dispersion_stats(distance, section)
    base$pseudo_R2 <- exact$R2
    base$pseudo_F <- exact$F
    base$exact_p <- exact$exact_p
    base$permutation_exceedances <- exact$permutation_exceedances
    base$permutation_F_min <- exact$permutation_F_min
    base$permutation_F_max <- exact$permutation_F_max
    base$dispersion_F <- disp$F
    base$dispersion_upper_mean <- disp$upper
    base$dispersion_middle_mean <- disp$middle
    base$dispersion_lower_mean <- disp$lower
    base$dispersion_lower_minus_upper <- disp$lower_minus_upper
    base$valid <- TRUE
    base
  }, error = function(e) {
    base$failure_reason <- conditionMessage(e)
    base
  })
}

analysis_result_row <- function(spec, draw_id, seed, selected_count, panel, result) {
  data.frame(
    spec_id = spec$spec_id,
    draw_id = as.integer(draw_id), seed = as.integer(seed),
    panel = spec$panel, taxonomy = spec$taxonomy,
    library_detection_threshold_reads =
      as.integer(spec$library_detection_threshold_reads),
    effort = spec$effort, distance = spec$distance,
    selected_library_count = as.integer(selected_count),
    analysis_unit = panel$analysis_unit,
    n_analysis_units = as.integer(result$n_analysis_units),
    n_nonzero_community_rows = as.integer(result$n_nonzero_community_rows),
    n_spatial_nonidentity_permutations = 27L,
    exact_p_resolution = 1 / 28,
    richness_upper_mean = result$richness_upper_mean,
    richness_middle_mean = result$richness_middle_mean,
    richness_lower_mean = result$richness_lower_mean,
    richness_contrast_lower_minus_upper = result$richness_contrast_lower_minus_upper,
    pseudo_R2 = result$pseudo_R2, pseudo_F = result$pseudo_F,
    exact_p = result$exact_p,
    permutation_exceedances = result$permutation_exceedances,
    permutation_F_min = result$permutation_F_min,
    permutation_F_max = result$permutation_F_max,
    dispersion_F = result$dispersion_F,
    dispersion_upper_mean = result$dispersion_upper_mean,
    dispersion_middle_mean = result$dispersion_middle_mean,
    dispersion_lower_mean = result$dispersion_lower_mean,
    dispersion_lower_minus_upper = result$dispersion_lower_minus_upper,
    valid = result$valid, failure_reason = result$failure_reason,
    stringsAsFactors = FALSE
  )
}

summarise_families <- function(specifications, direct_results, draw_results) {
  metric_names <- c(
    "richness_upper_mean", "richness_middle_mean", "richness_lower_mean",
    "richness_contrast_lower_minus_upper", "pseudo_R2", "pseudo_F", "exact_p",
    "dispersion_F", "dispersion_upper_mean", "dispersion_middle_mean",
    "dispersion_lower_mean", "dispersion_lower_minus_upper"
  )
  quant <- function(x, probability) {
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    unname(stats::quantile(x, probability, names = FALSE, type = 7))
  }
  rows <- vector("list", nrow(specifications))
  for (i in seq_len(nrow(specifications))) {
    spec <- specifications[i, , drop = FALSE]
    z <- if (spec$effort == "all_surviving") {
      direct_results[direct_results$spec_id == spec$spec_id, , drop = FALSE]
    } else {
      draw_results[draw_results$spec_id == spec$spec_id, , drop = FALSE]
    }
    assert(nrow(z) == spec$expected_draws,
           paste("Incomplete result rows for", spec$spec_id))
    failure <- unique(z$failure_reason[!z$valid & nzchar(z$failure_reason)])
    row <- data.frame(
      spec_id = spec$spec_id, panel = spec$panel, taxonomy = spec$taxonomy,
      library_detection_threshold_reads = spec$library_detection_threshold_reads,
      effort = spec$effort, distance = spec$distance,
      claim_scope = spec$claim_scope,
      result_basis = if (spec$effort == "all_surviving")
        "single full-effort result" else "median and 5th-95th percentiles of 100 retained draws",
      expected_draws = as.integer(spec$expected_draws), retained_draws = nrow(z),
      valid_draws = sum(z$valid), invalid_draws = sum(!z$valid),
      family_valid = all(z$valid),
      failure_reason = if (length(failure)) paste(failure, collapse = ";") else "",
      positive_richness_draws = sum(z$richness_contrast_lower_minus_upper > 0,
                                    na.rm = TRUE),
      positive_richness_draw_fraction =
        mean(z$richness_contrast_lower_minus_upper > 0, na.rm = TRUE),
      n_analysis_units_min = min(z$n_analysis_units),
      n_analysis_units_max = max(z$n_analysis_units),
      n_nonzero_community_rows_min = min(z$n_nonzero_community_rows),
      n_nonzero_community_rows_max = max(z$n_nonzero_community_rows),
      n_spatial_nonidentity_permutations = 27L,
      exact_p_resolution = 1 / 28,
      stringsAsFactors = FALSE
    )
    for (metric in metric_names) {
      values <- z[[metric]]
      row[[metric]] <- quant(values, 0.5)
      row[[paste0(metric, "_q05")]] <- quant(values, 0.05)
      row[[paste0(metric, "_q95")]] <- quant(values, 0.95)
    }
    rows[[i]] <- row
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

build_influence <- function(lib, design, specifications, counts) {
  threshold <- 1L
  unit_counts <- aggregate_selected_units(lib, design, counts, lib$library_id, threshold)
  omissions <- rbind(
    data.frame(omission_type = "site", omitted = design$site_order,
               stringsAsFactors = FALSE),
    data.frame(omission_type = "timepoint", omitted = design$common_timepoints,
               stringsAsFactors = FALSE)
  )
  reference_specs <- specifications[
    specifications$panel == "balanced_site_profile" &
      specifications$taxonomy == "analysis_label" &
      specifications$library_detection_threshold_reads == 1L &
      specifications$effort == "all_surviving", , drop = FALSE]
  assert(nrow(reference_specs) == 3L, "Expected three reference distance specifications")
  full_permutations <- series_permutation_indices(14L)
  rows <- list()
  k <- 1L
  for (i in seq_len(nrow(omissions))) {
    omission_type <- omissions$omission_type[[i]]
    omitted <- omissions$omitted[[i]]
    keep_unit <- design$unit_md$timepoint %in% design$common_timepoints
    if (omission_type == "site") keep_unit <- keep_unit & design$unit_md$site != omitted
    if (omission_type == "timepoint") keep_unit <- keep_unit & design$unit_md$timepoint != omitted
    md <- design$unit_md[keep_unit, , drop = FALSE]
    site_order <- if (omission_type == "site")
      setdiff(design$site_order, omitted) else design$site_order
    n_by_site <- table(factor(md$site, levels = site_order))
    composition <- rowsum(unit_counts[keep_unit, , drop = FALSE],
                          factor(md$site, levels = site_order), reorder = FALSE)
    composition <- composition[site_order, , drop = FALSE]
    occurrence <- (composition > 0) * 1
    unit_richness <- rowSums(unit_counts[keep_unit, , drop = FALSE] > 0)
    richness <- as.numeric(rowsum(unit_richness, factor(md$site, levels = site_order),
                                  reorder = FALSE)[, 1] / as.numeric(n_by_site))
    panel <- list(
      composition = composition, occurrence = occurrence, richness = richness,
      row_sites = site_order,
      row_groups = unname(design$site_section[site_order]),
      analysis_unit = "leave_one_out_site_profile"
    )
    # Retain the frozen 27 full-river arrangements even when one physical site is
    # omitted: section labels are permuted over E01-E14 first, then that fixed site
    # is removed.  This preserves the declared 1/28 resolution for all 75 rows.
    perm_groups <- permutation_groups(site_order, design, full_permutations)
    for (j in seq_len(nrow(reference_specs))) {
      spec <- reference_specs[j, , drop = FALSE]
      result <- analyse_panel(panel, spec$distance, design, perm_groups)
      rows[[k]] <- data.frame(
        influence_id = sprintf("I%03d", k), reference_spec_id = spec$spec_id,
        omission_type = omission_type, omitted = omitted,
        distance = spec$distance, analysis_unit = "site_profile",
        remaining_sites = length(site_order),
        remaining_common_timepoints = if (omission_type == "timepoint") 10L else 11L,
        n_spatial_nonidentity_permutations = 27L, exact_p_resolution = 1 / 28,
        richness_upper_mean = result$richness_upper_mean,
        richness_middle_mean = result$richness_middle_mean,
        richness_lower_mean = result$richness_lower_mean,
        richness_contrast_lower_minus_upper =
          result$richness_contrast_lower_minus_upper,
        pseudo_R2 = result$pseudo_R2, pseudo_F = result$pseudo_F,
        exact_p = result$exact_p,
        dispersion_F = result$dispersion_F,
        dispersion_lower_minus_upper = result$dispersion_lower_minus_upper,
        valid = result$valid, failure_reason = result$failure_reason,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  assert(nrow(out) == 75L, "Reference influence table must contain exactly 75 rows")
  out
}

evaluate_gate <- function(results, draw_results, influence, manifest) {
  gate <- manifest$layerE_multiverse$gate
  expected_families <- as.integer(manifest$layerE_multiverse$expected_specification_families)
  valid_families <- sum(results$family_valid)
  positive_valid <- sum(results$family_valid &
                          results$richness_contrast_lower_minus_upper > 0)
  valid_fraction <- valid_families / expected_families
  positive_family_fraction <- if (valid_families) positive_valid / valid_families else 0
  positive_draws <- sum(draw_results$richness_contrast_lower_minus_upper > 0,
                        na.rm = TRUE)
  positive_draw_fraction <- positive_draws / nrow(draw_results)

  reference <- results[
    results$panel == "balanced_site_profile" &
      results$taxonomy == "analysis_label" &
      results$library_detection_threshold_reads == 1L &
      results$effort == "all_surviving" &
      results$distance == "bray_curtis_proportions", , drop = FALSE]
  assert(nrow(reference) == 1L &&
         is.finite(reference$richness_contrast_lower_minus_upper),
         "Reference richness contrast is unavailable")
  ref_sign <- sign(reference$richness_contrast_lower_minus_upper)
  influence$sign_reversal <- influence$valid &
    sign(influence$richness_contrast_lower_minus_upper) == -ref_sign
  reversal_count <- sum(influence$sign_reversal)
  influence_pass <- all(influence$valid) && reversal_count == 0L

  rows <- data.frame(
    criterion = c("valid_specification_families", "positive_valid_family_contrasts",
                  "positive_one_library_draw_contrasts",
                  "reference_leaveout_sign_reversals", "overall_gate"),
    numerator = c(valid_families, positive_valid, positive_draws,
                  reversal_count, NA_integer_),
    denominator = c(expected_families, valid_families, nrow(draw_results),
                    nrow(influence), NA_integer_),
    observed_value = c(valid_fraction, positive_family_fraction,
                       positive_draw_fraction, reversal_count, NA_real_),
    comparison = c(">=", ">=", ">=", "==", "all component gates"),
    threshold = c(as.numeric(gate$min_valid_family_fraction),
                  as.numeric(gate$min_positive_family_fraction),
                  as.numeric(gate$min_positive_one_library_draw_fraction), 0, NA_real_),
    pass = c(valid_fraction >= gate$min_valid_family_fraction,
             positive_family_fraction >= gate$min_positive_family_fraction,
             positive_draw_fraction >= gate$min_positive_one_library_draw_fraction,
             influence_pass, FALSE),
    interpretation = c(
      "A family is valid only when every declared draw has finite distance, PERMANOVA, and dispersion results",
      "Direction is the median lower-minus-upper detected-label richness contrast among valid families",
      "All 3,600 declared one-library family-draw rows remain in the denominator, including distance failures",
      "All 75 reference leave-outs must be valid and none may reverse the reference richness sign",
      "Robust descriptive broadening requires every preregistered component to pass"),
    stringsAsFactors = FALSE
  )
  rows$pass[5] <- all(rows$pass[1:4])
  rows$observed_value[5] <- as.numeric(rows$pass[5])
  list(gate = rows, influence = influence,
       status = if (rows$pass[5]) "PASS: robust descriptive broadening" else
         "FAILED_GATE: robustness wording withheld")
}

make_multiverse_figure <- function(results, influence, gate_status, png_path, svg_path) {
  assert(requireNamespace("ggplot2", quietly = TRUE) &&
         requireNamespace("gridExtra", quietly = TRUE),
         "ggplot2 and gridExtra are required to render the multiverse figure")
  gg <- asNamespace("ggplot2")
  results$panel_label <- factor(results$panel,
    levels = c("balanced_site_profile", "all_site_profile", "balanced_site_time"),
    labels = c("14-site profiles / 11 common times",
               "14-site profiles / all usable times",
               "154 site x time units / site blocks"))
  results$family_status <- ifelse(results$family_valid, "valid", "retained failure")
  results$spec_index <- seq_len(nrow(results))
  results$distance_short <- factor(results$distance,
    levels = c("bray_curtis_proportions", "jaccard_occurrence",
               "robust_aitchison_rclr"),
    labels = c("Bray", "Jaccard", "rCLR"))
  results$config <- paste(
    ifelse(results$taxonomy == "analysis_label", "14 labels", "12 genera"),
    paste0(">=", results$library_detection_threshold_reads, " read"),
    ifelse(results$effort == "all_surviving", "all libraries", "one-library ensemble"),
    sep = " | ")
  results$config <- factor(results$config, levels = rev(unique(results$config)))

  p1 <- gg$ggplot(results, gg$aes(x = spec_index,
                                  y = richness_contrast_lower_minus_upper,
                                  ymin = richness_contrast_lower_minus_upper_q05,
                                  ymax = richness_contrast_lower_minus_upper_q95,
                                  colour = family_status, shape = effort)) +
    gg$geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
    gg$geom_linerange(linewidth = 0.45, alpha = 0.8) +
    gg$geom_point(size = 1.8) +
    gg$facet_wrap(~ panel_label, scales = "free_x", nrow = 1) +
    gg$scale_colour_manual(values = c("valid" = "#176B87", "retained failure" = "#A33B20")) +
    gg$scale_shape_manual(values = c("all_surviving" = 16,
                                     "one_library_100_draw_ensemble" = 17)) +
    gg$labs(
      title = "A. Richness direction across every declared specification family",
      subtitle = "Points are direct estimates or 100-draw medians; bars are 5th-95th percentiles",
      x = "Specification order", y = "Lower - upper mean detected-label richness",
      colour = "Distance-family status", shape = "Library effort") +
    gg$theme_bw(base_size = 10) +
    gg$theme(legend.position = "bottom", panel.grid.minor = gg$element_blank(),
             strip.text = gg$element_text(face = "bold"))

  p2 <- gg$ggplot(results, gg$aes(x = distance_short, y = config,
                                  fill = pseudo_R2)) +
    gg$geom_tile(colour = "white", linewidth = 0.5) +
    gg$geom_text(gg$aes(label = ifelse(family_valid,
      sprintf("R2 %.2f\np %.3f", pseudo_R2, exact_p), "retained\nfailure")),
      size = 2.6, lineheight = 0.92) +
    gg$facet_wrap(~ panel_label, nrow = 1, scales = "free_y") +
    gg$scale_fill_viridis_c(option = "C", limits = c(0, 1), na.value = "grey82") +
    gg$labs(
      title = "B. Restricted-series community result matrix",
      subtitle = "Cell text reports family median pseudo-R2 and exact p; failures remain visible",
      x = NULL, y = NULL, fill = "Pseudo-R2") +
    gg$theme_bw(base_size = 9) +
    gg$theme(panel.grid = gg$element_blank(), legend.position = "bottom",
             axis.text.y = gg$element_text(size = 7.5),
             strip.text = gg$element_text(face = "bold"))

  influence$distance_short <- factor(influence$distance,
    levels = c("bray_curtis_proportions", "jaccard_occurrence",
               "robust_aitchison_rclr"),
    labels = c("Bray", "Jaccard", "rCLR"))
  influence$omission_label <- factor(influence$omission_type,
    levels = c("site", "timepoint"), labels = c("Leave one site out", "Leave one timepoint out"))
  p3 <- gg$ggplot(influence, gg$aes(x = distance_short,
                                    y = richness_contrast_lower_minus_upper,
                                    colour = omission_label)) +
    gg$geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.35) +
    gg$geom_point(position = gg$position_jitter(width = 0.14, height = 0),
                  alpha = 0.75, size = 1.8) +
    gg$scale_colour_manual(values = c("Leave one site out" = "#6A3D9A",
                                     "Leave one timepoint out" = "#E17C05")) +
    gg$labs(
      title = "C. Reference influence analysis (75 retained leave-outs)",
      subtitle = gate_status,
      x = NULL, y = "Lower - upper richness contrast", colour = NULL) +
    gg$theme_bw(base_size = 10) +
    gg$theme(legend.position = "bottom", panel.grid.minor = gg$element_blank())

  footer <- grid::textGrob(paste(
    "Effect sizes and exact 1/28 permutation resolution are emphasized;",
    "the established balanced-site centroid claim remains withheld."),
    gp = grid::gpar(fontsize = 9, col = "grey25"))
  grob <- gridExtra::arrangeGrob(p1, p2, p3, footer, ncol = 1,
                                 heights = c(1.05, 1.45, 0.85, 0.1))
  grDevices::png(png_path, width = 2400, height = 2200, res = 170)
  grid::grid.draw(grob)
  grDevices::dev.off()
  grDevices::svg(svg_path, width = 14.2, height = 13.0, pointsize = 10,
                 onefile = TRUE)
  grid::grid.draw(grob)
  grDevices::dev.off()
}

layerE_main <- function() {
  started <- proc.time()[["elapsed"]]
  repo_root <- layerE_repo_root()
  setwd(repo_root)
  assert(requireNamespace("jsonlite", quietly = TRUE),
         "The jsonlite package is required to read the frozen manifest")

  paths <- list(
    manifest = "data/metadata/impact_sprint_manifest.json",
    input = "data/derived/impact_library_counts.tsv",
    specifications = "data/derived/impact_layerE_specifications.tsv",
    library_draws = "data/derived/impact_layerE_library_draws.tsv",
    draw_results = "data/derived/impact_layerE_draw_results.tsv",
    results = "data/derived/impact_layerE_results.tsv",
    influence = "data/derived/impact_layerE_influence.tsv",
    gate = "data/derived/impact_layerE_gate.tsv",
    figure_png = "figures/impact_layerE_multiverse.png",
    figure_svg = "figures/impact_layerE_multiverse.svg",
    determinism = "data/metadata/impact_layerE_determinism_receipt.tsv",
    run_receipt = "data/metadata/impact_layerE_run_receipt.tsv"
  )
  assert(file.exists(paths$manifest), "Missing frozen impact manifest")
  manifest <- jsonlite::fromJSON(paths$manifest, simplifyVector = TRUE)
  cfg <- manifest$layerE_multiverse
  assert(identical(as.integer(cfg$expected_specification_families), 72L) &&
         identical(as.integer(cfg$one_library_draws), 100L) &&
         identical(as.integer(cfg$spatial_nonidentity_permutations), 27L) &&
         identical(as.integer(cfg$expected_reference_influence_rows), 75L),
         "Frozen Layer-E manifest counts changed")

  canonical_paths <- unlist(paths[c("specifications", "library_draws", "draw_results",
                                    "results", "influence", "gate",
                                    "figure_png", "figure_svg")], use.names = TRUE)
  before_hash <- setNames(vapply(canonical_paths, function(path) {
    if (file.exists(path)) sha256(path) else NA_character_
  }, character(1)), names(canonical_paths))

  lib <- load_impact_library(paths$input)
  design <- design_from_library(lib)
  specifications <- make_specifications(manifest)
  base_seed <- as.integer(manifest$seeds$layerE_multiverse)
  library_draws <- make_library_draws(lib, design, base_seed,
                                      as.integer(cfg$one_library_draws))
  draw_ids <- split(library_draws$selected_library_id, library_draws$draw_id)

  direct_rows <- list()
  draw_rows <- list()
  kd <- 1L
  kr <- 1L
  # Compute one panel per panel x taxonomy x threshold x effort x draw, then
  # evaluate all three distances.  This keeps the 3,600 rows complete without
  # repeating unit aggregation for every distance.
  for (panel_name in cfg$panels) {
    for (taxonomy in cfg$taxonomy) {
      counts <- taxonomy_matrix(lib, taxonomy)
      for (threshold in cfg$library_detection_threshold_reads) {
        for (effort in cfg$effort) {
          selections <- if (effort == "all_surviving") {
            list(`0` = lib$library_id)
          } else draw_ids
          for (draw_name in names(selections)) {
            draw_id <- as.integer(draw_name)
            selected_ids <- selections[[draw_name]]
            unit_counts <- aggregate_selected_units(lib, design, counts, selected_ids,
                                                     as.integer(threshold))
            panel <- make_panel(unit_counts, design, panel_name)
            for (distance_name in cfg$distances) {
              spec <- specifications[
                specifications$panel == panel_name &
                  specifications$taxonomy == taxonomy &
                  specifications$library_detection_threshold_reads == threshold &
                  specifications$effort == effort &
                  specifications$distance == distance_name, , drop = FALSE]
              assert(nrow(spec) == 1L, "Could not resolve specification family")
              result <- analyse_panel(panel, distance_name, design)
              seed <- if (effort == "all_surviving") base_seed else base_seed + draw_id
              row <- analysis_result_row(spec, draw_id, seed,
                                          length(selected_ids), panel, result)
              if (effort == "all_surviving") {
                direct_rows[[kd]] <- row
                kd <- kd + 1L
              } else {
                draw_rows[[kr]] <- row
                kr <- kr + 1L
              }
            }
          }
        }
      }
    }
  }
  direct_results <- do.call(rbind, direct_rows)
  draw_results <- do.call(rbind, draw_rows)
  rownames(direct_results) <- rownames(draw_results) <- NULL
  assert(nrow(direct_results) == 36L && nrow(draw_results) == 3600L,
         "Expected 36 full-effort results and 3,600 one-library draw results")
  results <- summarise_families(specifications, direct_results, draw_results)
  assert(nrow(results) == 72L, "Family result table must contain all 72 specifications")

  label_counts <- taxonomy_matrix(lib, "analysis_label")
  influence <- build_influence(lib, design, specifications, label_counts)
  gate_eval <- evaluate_gate(results, draw_results, influence, manifest)
  influence <- gate_eval$influence
  gate <- gate_eval$gate

  dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
  dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
  dir.create("figures", recursive = TRUE, showWarnings = FALSE)
  specifications <- normalise_numeric(specifications)
  library_draws <- normalise_numeric(library_draws)
  draw_results <- normalise_numeric(draw_results)
  results <- normalise_numeric(results)
  influence <- normalise_numeric(influence)
  gate <- normalise_numeric(gate)
  write_tsv(specifications, paths$specifications)
  write_tsv(library_draws, paths$library_draws)
  write_tsv(draw_results, paths$draw_results)
  write_tsv(results, paths$results)
  write_tsv(influence, paths$influence)
  write_tsv(gate, paths$gate)
  make_multiverse_figure(results, influence, gate_eval$status,
                         paths$figure_png, paths$figure_svg)

  after_hash <- setNames(vapply(canonical_paths, sha256, character(1)),
                         names(canonical_paths))
  determinism <- data.frame(
    artifact = unname(canonical_paths),
    preexisting_sha256 = unname(before_hash),
    rebuilt_sha256 = unname(after_hash),
    byte_identical_to_preexisting = !is.na(before_hash) & before_hash == after_hash,
    status = ifelse(is.na(before_hash), "NOT_EVALUATED_FIRST_BUILD",
                    ifelse(before_hash == after_hash, "PASS", "FAIL")),
    stringsAsFactors = FALSE
  )
  write_tsv(determinism, paths$determinism)

  elapsed <- proc.time()[["elapsed"]] - started
  run_receipt <- data.frame(
    field = c("analysis", "manifest_seed", "runtime_seconds", "R_version",
              "vegan_version", "jsonlite_version", "digest_version",
              "specification_families", "one_library_draw_result_rows",
              "library_selection_rows", "reference_influence_rows",
              "nonidentity_permutations", "exact_p_resolution", "gate_status"),
    value = c("Layer-E robustness multiverse", base_seed, sprintf("%.3f", elapsed),
              R.version.string, as.character(utils::packageVersion("vegan")),
              as.character(utils::packageVersion("jsonlite")),
              as.character(utils::packageVersion("digest")), nrow(results),
              nrow(draw_results), nrow(library_draws), nrow(influence), 27,
              sprintf("%.10f", 1 / 28), gate_eval$status),
    stringsAsFactors = FALSE
  )
  write_tsv(run_receipt, paths$run_receipt)

  cat(sprintf(
    paste0("Layer-E multiverse complete: %d families, %d one-library draw rows, ",
           "%d library selections, %d influence rows, %d/%d valid families.\n"),
    nrow(results), nrow(draw_results), nrow(library_draws), nrow(influence),
    sum(results$family_valid), nrow(results)))
  cat(gate_eval$status, "\n")
  cat(sprintf("Runtime: %.3f seconds. Determinism receipt: %s\n",
              elapsed, paths$determinism))
  invisible(list(specifications = specifications, draws = draw_results,
                 results = results, influence = influence, gate = gate))
}

if (sys.nframe() == 0L) layerE_main()
