#!/usr/bin/env Rscript
# =============================================================================
# 00_setup.R  --  R-layer setup on top of the conda `edna` env (environment.yml).
#
# Run INSIDE the env:   micromamba run -n edna Rscript R/00_setup.R
#
# Idempotent: only installs what's missing or broken. Does three things the
# conda layer can't:
#   (1) CRAN-only packages (corncob, jagsUI, occumb, msocc + its deps),
#   (2) source rebuilds of dada2 & iNEXT (their osx-arm64 conda builds are broken),
#   (3) the cmdstanr smoke test (the project's #1 setup risk).
#
# Isolation note: this writes ONLY into the active env library
# (~/micromamba/envs/edna/lib/R/library) — never into the congen/CG env.
# =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L))
Sys.unsetenv("GITHUB_PAT")   # avoid a stale/invalid PAT in the credential store

loadable <- function(p) isTRUE(tryCatch({ suppressPackageStartupMessages(
  library(p, character.only = TRUE)); TRUE }, error = function(e) FALSE))

ensure_cran <- function(pkgs) {
  for (p in pkgs) if (!loadable(p)) {
    message("Installing (CRAN, source): ", p)
    install.packages(p, type = "source")
  }
}

# --- (1) CRAN-only packages -------------------------------------------------
ensure_cran(c("corncob", "jagsUI", "occumb",      # occumb needs jagsUI -> rjags(conda) -> JAGS(conda)
              "pgdraw", "mvtnorm"))               # msocc deps

# msocc was archived from CRAN -> install from its public GitHub source tarball.
# PINNED to a commit SHA (was bare `master`, a moving target with no recorded
# RemoteSha). This is the one input renv was supposed to nail and couldn't.
MSOCC_SHA <- "0180cb0489da46208079ae31aac46f6d83c3e058"  # StrattonCh/msocc master @ 2026-06
if (!loadable("msocc")) {
  message("Installing msocc from pinned tarball (StrattonCh/msocc@", substr(MSOCC_SHA, 1, 7), ")")
  tf <- tempfile(fileext = ".tar.gz")
  download.file(sprintf("https://github.com/StrattonCh/msocc/archive/%s.tar.gz", MSOCC_SHA),
                tf, quiet = TRUE)
  td <- tempfile(); untar(tf, exdir = td)
  pkgdir <- list.dirs(td, recursive = FALSE)[1]   # msocc-<sha>; robust to the dir name
  install.packages(pkgdir, repos = NULL, type = "source")
}

# --- (2) Source rebuilds of the broken conda binaries -----------------------
if (!loadable("iNEXT")) { message("Rebuilding iNEXT from source"); install.packages("iNEXT", type = "source") }
if (!loadable("dada2")) {
  message("Rebuilding dada2 from Bioconductor source")
  BiocManager::install("dada2", update = FALSE, ask = FALSE, force = TRUE, type = "source")
}

# --- Roll-call --------------------------------------------------------------
pk <- c("dada2","phyloseq","decontam","ALDEx2","ANCOMBC","corncob","vegan","iNEXT","indicspecies",
        "zCompositions","msocc","occumb","cmdstanr","posterior","bayesplot","shiny","quarto",
        "renv","tidyverse","rjags","jagsUI")
ip <- rownames(installed.packages())
miss <- pk[!pk %in% ip]
cat("\n", length(pk) - length(miss), "/", length(pk), "target packages present.\n", sep = "")
if (length(miss)) { cat("MISSING:", paste(miss, collapse = ", "), "\n"); quit(status = 1) }

# --- Version manifest: the reviewer-facing, preprint-methods version record. --
# conda's environment.lock.yml pins the transitive closure (build-string exact);
# this is the human-readable roll-up of the 21 R targets + the msocc SHA + Bioc
# release, so a clean setup is version-deterministic and citable.
.args0 <- commandArgs(FALSE); .self0 <- sub("^--file=", "", grep("^--file=", .args0, value = TRUE))
.root <- if (length(.self0)) dirname(dirname(normalizePath(.self0[1]))) else "."
bioc <- tryCatch(as.character(BiocManager::version()), error = function(e) NA_character_)
bioc_pkgs <- c("dada2", "phyloseq", "decontam", "ALDEx2", "ANCOMBC")
manifest <- data.frame(
  package = pk,
  version = vapply(pk, function(p) as.character(packageVersion(p)), ""),
  source  = ifelse(pk == "msocc", paste0("github:StrattonCh/msocc@", substr(MSOCC_SHA, 1, 12)),
            ifelse(pk %in% bioc_pkgs, paste0("Bioconductor-", bioc), "CRAN/conda")),
  row.names = NULL)
write.table(manifest, file.path(.root, "environment.R-packages.lock.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote environment.R-packages.lock.tsv (", nrow(manifest), " packages, Bioc ", bioc, ")\n", sep = "")

# --- (3) Bayesian-backend smoke tests: prove BOTH engines end-to-end --------
#   Stan (stretch engine) AND JAGS (flagship occumb engine). Each writes a PASS
#   receipt to docs/smoke_test_receipts.md so the "smoke test passes" claim is
#   auditable, not ephemeral stdout. Sourced in fresh envs so commandArgs/quit
#   inside them won't kill this parent script.
.args <- commandArgs(FALSE); .self <- sub("^--file=", "", .args[grep("^--file=", .args)])
.here <- if (length(.self)) dirname(normalizePath(.self)) else "R"
.rscript <- file.path(R.home("bin"), "Rscript")   # the env's own Rscript, not PATH-dependent
run_smoke <- function(script) {
  cat("\nRunning", script, "...\n")
  rc <- system2(.rscript, file.path(.here, script))
  if (rc != 0L) { cat("SMOKE TEST FAILED:", script, "\n"); quit(status = 1) }
}
run_smoke("smoke_test_cmdstanr.R")   # Stan (stretch)
run_smoke("smoke_test_jags.R")       # JAGS/occumb (flagship)
cat("\nBoth Bayesian backends proven. See docs/smoke_test_receipts.md\n")
