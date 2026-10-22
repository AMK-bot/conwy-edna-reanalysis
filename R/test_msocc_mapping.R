#!/usr/bin/env Rscript
# test_msocc_mapping.R -- proves the Phase-1 2-level → msocc mapping (deviation D1).
#
# Simulate the true 2-level process y_ij ~ Bern(z_i · d), z_i ~ Bern(ψ), encode it
# the way 05_occupancy_msocc.R does (3 water samples at msocc's SAMPLE level, 1
# dummy PCR/sample, rep ~ 1), fit, and assert msocc recovers ψ and the per-sample
# detection d = θ·p — and that ψ̂ lands ABOVE naïve occupancy (the imperfect-
# detection correction). If the 3 reps were (wrongly) placed at the PCR level
# instead, ψ would confound with θ and this recovery would fail.
#
#   micromamba run -n edna Rscript R/test_msocc_mapping.R

suppressPackageStartupMessages(library(msocc))
set.seed(7)

M <- 250L; J <- 3L                 # 250 occupancy units, 3 water samples each
psi_true <- 0.60; d_true <- 0.40   # occupancy; per-sample detection (= θ·p)

z      <- rbinom(M, 1, psi_true)
units  <- rep(seq_len(M), each = J)
y      <- rbinom(M * J, 1, z[units] * d_true)

wide  <- data.frame(site = sprintf("U%03d", units), sample = rep(seq_len(J), M), pcr1 = y)
scov  <- data.frame(site = unique(wide$site))
mcov  <- wide[, c("site", "sample")]

fit <- msocc_mod(wide,
  site   = list(model = ~1, cov_tbl = scov),
  sample = list(model = ~1, cov_tbl = mcov),
  rep    = list(model = ~1, cov_tbl = mcov),
  num.mcmc = 6000L, beta_bin = FALSE, seed = 1L, progress = FALSE)

keep    <- 1001:6000
psi_hat <- mean(plogis(fit$beta[keep, 1]))
det_hat <- mean(plogis(fit$alpha[keep, 1]) * plogis(fit$delta[keep, 1]))
naive   <- mean(tapply(y, units, function(v) any(v > 0)))   # ≈ ψ·(1−(1−d)^3) < ψ

cat(sprintf("ψ:   true %.2f  est %.3f\n", psi_true, psi_hat))
cat(sprintf("θ·p: true %.2f  est %.3f\n", d_true, det_hat))
cat(sprintf("naïve occupancy %.3f  (must be < true ψ — the false negatives)\n", naive))

stopifnot(abs(psi_hat - psi_true) < 0.10)   # ψ recovered (not confounded)
stopifnot(abs(det_hat - d_true)  < 0.10)    # detection recovered
stopifnot(psi_hat > naive + 0.05)           # model corrects upward for missed detections
cat("\nPASS — sample-level mapping recovers ψ and θ·p; ψ̂ sits above naïve occupancy.\n")
