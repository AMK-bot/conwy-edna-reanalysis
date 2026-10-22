#!/usr/bin/env Rscript
# smoke_test_jags.R -- Phase 0 gate for the FLAGSHIP Bayesian engine.
#
# occumb (Phase-3 community occupancy) runs on JAGS via:
#     occumb -> jagsUI -> rjags -> JAGS
# The env notes record a real rjags<->JAGS "Runtime link error" that had to be
# worked around with conda r-rjags + conda jags. That chain was previously only
# *load*-checked -- while the mere *stretch* engine (Stan) got a full
# compile+sample proof. This test fixes that asymmetry. Two parts:
#
#   A. occumb itself (the exact API 06_occupancy_occumb.R uses): prove the chain
#      compiles+samples+returns a fit with diagnostics computed. NOTE: occumb's
#      hierarchical model (rho/sigma/shared effects, prior_ulim=1e4) is weakly
#      identified, so a degenerate toy legitimately does NOT reach Rhat<1.1 --
#      that is a property of the toy, not the engine. So Part A asserts ENGINE
#      facts (fit returned, samples produced, Rhat computed), not convergence.
#   B. a well-identified single-species occupancy via jagsUI directly: this one
#      DOES converge, so it carries the meaningful Rhat<1.1 + parameter-recovery
#      proof that the sampler produces correct inference.
#
# Run:  micromamba run -n edna Rscript R/smoke_test_jags.R
# Exits non-zero (and writes no PASS receipt) on any failure.

t0 <- Sys.time()
stopifnot("occumb not installed"  = requireNamespace("occumb",  quietly = TRUE),
          "jagsUI not installed"  = requireNamespace("jagsUI",  quietly = TRUE))
suppressPackageStartupMessages({ library(occumb); library(jagsUI) })

ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "??")
jags_ver <- tryCatch(as.character(rjags::jags.version()), error = function(e) "??")
cat("occumb :", ver("occumb"), "\n");  cat("jagsUI :", ver("jagsUI"), "\n")
cat("rjags  :", ver("rjags"),  "\n");  cat("JAGS   :", jags_ver, "\n")
cat("R arch :", R.version$arch, "\n\n")

link_err <- function(e) {
  msg <- conditionMessage(e)
  cat("\nERROR:\n  ", msg, "\n")
  if (grepl("Runtime link error|Cannot load module|failed to load", msg, ignore.case = TRUE))
    cat("\n>>> Known rjags<->JAGS link failure: ensure BOTH r-rjags AND jags come from conda.\n")
  quit(status = 1)
}

# --- Part A: occumb runs end-to-end through JAGS ------------------------------
set.seed(1)
I <- 2; J <- 2; K <- 2
dat <- occumbData(y = array(sample.int(I * J * K), dim = c(I, J, K)))
cat("[A] Fitting trivial occumb null model on JAGS...\n")
fitA <- tryCatch(
  occumb(data = dat, engine = "JAGS",
         n.chains = 2, n.adapt = 500, n.burnin = 500, n.thin = 1, n.iter = 1500),
  error = link_err)
rhatA <- fitA@fit$summary[, "Rhat"]; rhatA <- rhatA[is.finite(rhatA)]
okA <- inherits(fitA, "occumbFit") &&
       length(fitA@fit$samples) == 2 &&            # 2 chains of samples produced
       length(rhatA) > 0                           # Rhat diagnostic computed
cat(sprintf("[A] occumbFit returned, %d chains, %d Rhat values computed -> %s\n",
            length(fitA@fit$samples), length(rhatA), if (okA) "OK" else "FAIL"))

# --- Part B: well-identified single-species occupancy converges ---------------
set.seed(7)
Jb <- 50; Kb <- 3; psi_true <- 0.6; p_true <- 0.5
z <- rbinom(Jb, 1, psi_true)
y <- matrix(rbinom(Jb * Kb, 1, rep(z, Kb) * p_true), Jb, Kb)
model_str <- "model {
  psi ~ dunif(0,1)
  p   ~ dunif(0,1)
  for (j in 1:J) {
    z[j] ~ dbern(psi)
    for (k in 1:K) { y[j,k] ~ dbern(z[j]*p) }
  }
}"
mf <- tempfile(fileext = ".jags"); writeLines(model_str, mf)
cat("[B] Fitting single-species occupancy via jagsUI...\n")
fitB <- tryCatch(
  jagsUI::jags(data = list(y = y, J = Jb, K = Kb),
               inits = function() list(z = rep(1, Jb)),
               parameters.to.save = c("psi", "p"), model.file = mf,
               n.chains = 3, n.adapt = 500, n.burnin = 1000, n.iter = 3000,
               n.thin = 1, verbose = FALSE),
  error = link_err)
rhatB <- c(fitB$Rhat$psi, fitB$Rhat$p)
psi_hat <- fitB$mean$psi; p_hat <- fitB$mean$p
okB <- all(is.finite(rhatB)) && max(rhatB) < 1.1 &&
       is.finite(psi_hat) && psi_hat > 0.3 && psi_hat < 0.9   # recovers psi_true=0.6
cat(sprintf("[B] psi_hat=%.2f (true 0.6)  p_hat=%.2f (true 0.5)  max Rhat=%.3f -> %s\n",
            psi_hat, p_hat, max(rhatB), if (okB) "OK" else "FAIL"))

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Elapsed: %.1fs\n", elapsed))

# One-line auditable PASS receipt (upsert by engine key).
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
.docs <- normalizePath(file.path(if (length(.self)) dirname(.self) else "R", "..", "docs"),
                       mustWork = FALSE)
write_receipt <- function(engine, line) {
  f <- file.path(.docs, "smoke_test_receipts.md")
  hdr <- c("# Smoke-test receipts", "",
           "Auditable PASS lines written by `R/smoke_test_*.R` (latest run per engine).", "")
  old <- if (file.exists(f)) readLines(f, warn = FALSE) else hdr
  old <- old[!grepl(paste0("^- \\*\\*", engine, "\\*\\*"), old)]
  writeLines(c(old, line), f)
}

if (okA && okB) {
  cat("\noccumb/JAGS SMOKE TEST PASSED ✓  -- flagship engine works AND converges.\n")
  write_receipt("JAGS", sprintf(
    "- **JAGS** (occumb %s / jagsUI %s / rjags %s / JAGS %s): PASS ✓ — occumb chain runs; jagsUI occupancy converges (psi_hat=%.2f, max Rhat %.3f), %.1fs",
    ver("occumb"), ver("jagsUI"), ver("rjags"), jags_ver, psi_hat, max(rhatB), elapsed))
} else {
  cat("\noccumb/JAGS SMOKE TEST FAILED ✗ -- inspect Part A/B output above.\n")
  quit(status = 1)
}
