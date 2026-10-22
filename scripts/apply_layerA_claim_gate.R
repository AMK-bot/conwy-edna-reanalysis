#!/usr/bin/env Rscript

# Replace smoke-fit Layer-A figures with explicit withheld-status cards whenever
# the known-truth recovery gate fails. This prevents cached exploratory numbers
# from being mistaken for biological results after a failed convergence check.

write_layerA_withheld_figures <- function(status,
                                          figure_dir = "figures",
                                          receipt_path = "data/metadata/layerA_figure_gate.tsv") {
  stopifnot(is.data.frame(status), nrow(status) == 1L)
  required <- c("max_psi_Rhat", "convergence_threshold", "claim_status", "allowed_claim")
  stopifnot(all(required %in% names(status)))

  panels <- c(
    phase3_occupancy_psi_forest.png = "Smoke-fit per-species occupancy forest",
    phase3_occupancy_detection_fnr.png = "Smoke-fit detection and false-negative panel",
    phase3_occupancy_cov_dist_effect.png = "Smoke-fit downstream-distance coefficients",
    phase3_occupancy_cov_corrected_vs_raw.png = "Smoke-fit corrected-versus-naive occupancy shifts",
    phase3_occupancy_cov_design.png = "Smoke-fit survey-design utility curve"
  )

  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  outputs <- file.path(figure_dir, names(panels))
  for (i in seq_along(outputs)) {
    png(outputs[[i]], width = 1500, height = 650, res = 150)
    par(mar = c(0, 0, 0, 0), family = "sans")
    plot.new()
    rect(0, 0, 1, 1, col = "#f7f5f0", border = NA)
    rect(0.05, 0.12, 0.95, 0.88, col = "#edf4f5", border = "#8abfc7", lwd = 2)
    rect(0.05, 0.78, 0.95, 0.88, col = "#a33d46", border = NA)
    text(0.08, 0.83, status$claim_status, adj = c(0, 0.5), col = "white",
         cex = 1.35, font = 2)
    text(0.08, 0.68, "Layer A numerical display withheld", adj = c(0, 0.5),
         col = "#14243d", cex = 2.0, font = 2)
    text(0.08, 0.56,
         sprintf("Recovery gate failed: max psi Rhat = %.2f; required <= %.2f.",
                 status$max_psi_Rhat, status$convergence_threshold),
         adj = c(0, 0.5), col = "#14243d", cex = 1.35)
    text(0.08, 0.44, paste("Suppressed:", panels[[i]]), adj = c(0, 0.5),
         col = "#526176", cex = 1.25)
    text(0.08, 0.34,
         "No rare-species magnitude, direction, or coefficient is an allowed finding.",
         adj = c(0, 0.5), col = "#a33d46", cex = 1.12, font = 2)
    text(0.08, 0.27,
         "No survey recommendation is supported by this smoke fit.",
         adj = c(0, 0.5), col = "#a33d46", cex = 1.12, font = 2)
    text(0.08, 0.17,
         "See data/metadata/layerA_recovery_status.tsv and docs/phase3-claim-status.md.",
         adj = c(0, 0.5), col = "#526176", cex = 1.0)
    dev.off()
  }
  dir.create(dirname(receipt_path), recursive = TRUE, showWarnings = FALSE)
  receipt <- data.frame(
    figure = outputs,
    suppressed_display = unname(panels),
    claim_status = status$claim_status[[1]],
    max_psi_Rhat = status$max_psi_Rhat[[1]],
    convergence_threshold = status$convergence_threshold[[1]],
    md5 = unname(tools::md5sum(outputs)),
    generated_by = "scripts/apply_layerA_claim_gate.R",
    stringsAsFactors = FALSE
  )
  write.table(receipt, receipt_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[Layer-A gate] replaced ", length(outputs), " numerical figures with withheld-status cards")
  message("[Layer-A gate] wrote ", receipt_path)
  invisible(outputs)
}

if (sys.nframe() == 0L) {
  status_path <- "data/metadata/layerA_recovery_status.tsv"
  stopifnot(file.exists(status_path))
  status <- read.delim(status_path, check.names = FALSE, stringsAsFactors = FALSE)
  if (nrow(status) != 1L) stop("Layer-A status receipt must contain exactly one row")
  if (identical(status$claim_status[[1]], "WITHHELD_UNCONVERGED")) {
    write_layerA_withheld_figures(status)
  } else {
    message("[Layer-A gate] convergence receipt does not require suppression; no figures changed")
  }
}
