#!/usr/bin/env Rscript
# smoke_test_cmdstanr.R -- Phase 0 gate for the project's #1 setup risk.
#
# Compiles and samples a trivial Stan model to prove the full cmdstanr -> CmdStan ->
# C++ toolchain works natively on this machine BEFORE we depend on it for the
# occumb community-occupancy MCMC in Phase 3.
#
# Run:  micromamba run -n edna Rscript R/smoke_test_cmdstanr.R
# Expects the conda-provided CmdStan (env var CMDSTAN, set by the env's activate hook).

t0 <- Sys.time()
stopifnot("cmdstanr not installed" = requireNamespace("cmdstanr", quietly = TRUE))
library(cmdstanr)

# Point cmdstanr at the conda CmdStan (prebuilt arm64); fall back to the env default.
cmdstan <- Sys.getenv("CMDSTAN", unset = "")
if (nzchar(cmdstan) && dir.exists(cmdstan)) set_cmdstan_path(cmdstan)
cat("cmdstanr  :", as.character(packageVersion("cmdstanr")), "\n")
cat("CmdStan   :", cmdstan_version(), "\n")
cat("CmdStan @ :", cmdstan_path(), "\n")
cat("R arch    :", R.version$arch, "\n\n")

# A trivial Bernoulli model: estimate theta from 10 coin flips (8 heads).
model_code <- "
data { int<lower=0> N; array[N] int<lower=0,upper=1> y; }
parameters { real<lower=0,upper=1> theta; }
model { theta ~ beta(1, 1); y ~ bernoulli(theta); }
"
stan_file <- write_stan_file(model_code)

cat("Compiling trivial model (first compile is the real toolchain test)...\n")
mod <- cmdstan_model(stan_file)

cat("Sampling...\n")
fit <- mod$sample(
  data = list(N = 10, y = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0)),
  seed = 1,                       # deterministic: keeps the committed theta/Rhat receipt stable
  chains = 2, parallel_chains = 2,
  iter_warmup = 500, iter_sampling = 500, refresh = 0, show_messages = FALSE
)

est <- fit$summary("theta")
cat(sprintf("\nRecovered theta = %.3f (true posterior mean ~0.75 for 8/10)\n", est$mean))
ok <- is.finite(est$mean) && est$mean > 0.5 && est$mean < 0.95 && est$rhat < 1.05
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Elapsed: %.1fs\n", elapsed))

# One-line auditable PASS receipt (upsert by engine key).
# ponytail: this 11-line block is duplicated in smoke_test_jags.R; a shared lib
# file + sourcing is more machinery than the duplication saves for 2 call sites.
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

if (ok) {
  cat("\ncmdstanr SMOKE TEST PASSED ✓  -- native Stan MCMC works on this machine.\n")
  write_receipt("Stan", sprintf(
    "- **Stan** (cmdstanr %s / CmdStan %s): PASS ✓ — Bernoulli compile+sample, theta=%.3f, Rhat %.3f, %.1fs",
    as.character(packageVersion("cmdstanr")), cmdstan_version(), est$mean, est$rhat, elapsed))
} else {
  cat("\ncmdstanr SMOKE TEST FAILED ✗ -- inspect compile/sample output above.\n")
  quit(status = 1)
}
