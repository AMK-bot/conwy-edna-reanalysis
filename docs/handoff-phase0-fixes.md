# Handoff — Conwy eDNA reanalysis: Phase-0 hardening before Phase 1

> **Historical audit record — all P1–P6 items are resolved.** This document preserves the
> questions raised at the end of Phase 0 and therefore uses some then-current or hypothetical
> wording. It is not the current status document. For current evidence boundaries, completion,
> and commands, use `eDNAplan.md`, `README.md`, and `START_HERE.md`.

---

## 0. Orientation (read first)

- **Project:** A reanalysis (no wet lab) of the River Conwy 12S fish-eDNA dataset
  (Perry et al. 2024, *Nat. Commun.* 15:4372, doi:10.1038/s41467-024-48640-3, CC BY 4.0).
  Its current layers are **A** = a community-occupancy stress test whose rare-label estimates
  are withheld after a failed recovery/convergence gate, **B** = current-reference concordance
  (not a 2020→2026 replay), and **E** = a compositionally aware community baseline with
  site/time-restricted inference. The finished deliverable is a local website release; public
  preprint publication and Shiny deployment are optional owner-authorized actions.
- **Full plan:** `eDNAplan.md` in the repository root is the canonical current status and plan.
- **Repo:** `~/Desktop/MMP/eDNA/conwy-edna-reanalysis/` is the git repository. Run git there,
  not in its parent directory.

### How to run anything (non-negotiable environment rule)
R lives in an **isolated, native-arm64 conda env named `edna`** built with micromamba.
```bash
micromamba run -n edna Rscript R/<script>.R        # micromamba is at ~/.local/bin/micromamba
micromamba run -n edna quarto render report/conwy-reanalysis.qmd
```
- R is **4.4.3, aarch64**. The env root is `~/micromamba/envs/edna` (separate from the user's
  x86_64 `congen`/anaconda3 env).
- **ISOLATION IS A HARD USER REQUIREMENT.** Never install eDNA packages into `congen`/anaconda3.
  Only ever touch the `edna` env library.
- Setup is reproducible & idempotent: `bash setup_env.sh` →
  `environment.yml` (conda layer) then `R/00_setup.R` (CRAN top-ups + source rebuilds + smoke test).
- **Stale `GITHUB_PAT` gotcha:** there is an invalid PAT in the credential store. Any
  `remotes::install_github` returns HTTP 401. Always `Sys.unsetenv("GITHUB_PAT")` and install
  GitHub packages from a source **tarball** (this is how `00_setup.R` installs `msocc`).

### Hard data constraints you must never violate
1. **`E##`, NOT `C##`.** Conwy = `sample_alias` matching `^E[0-9][0-9]_R[0-9]_T[0-9][0-9]$`
   (sites E01–E14, timepoints T01–T19, reps R1–R3). `C##` is a **different river** in the same
   ENA archive. The filter is enforced in `data/download_ena.sh` and `data/verify_ena_join.py`.
   Never relax it.
2. **T16 is unsurveyed for all 14 sites** (+ `E02_T17_R3` missing) → **18 usable timepoints,
   not 19**. In occupancy detection histories, T16 cells must be **NA (unsurveyed)**, never `0`
   (encoding them as non-detections would bias occupancy ψ downward — the exact error this
   project exists to correct).
3. **Do not commit raw data** (the ~71 GB FASTQs) or re-fetchable publisher files. **Do commit**
   the small derived metadata tables. See `.gitignore`.
4. **Reads ≠ abundance.** Default to occurrence/detection framing everywhere.
5. **Do not download the 71 GB FASTQs yet** — that's Phase 2. P1–P6 below need none of it.

---

## 1. Current infrastructure — what exists and what's verified

**Repo tree (matches plan §4):**
```
conwy-edna-reanalysis/
├── README.md  .gitignore  conwy-edna-reanalysis.Rproj
├── environment.yml  environment.lock.yml  setup_env.sh
├── R/   00_setup.R  smoke_test_cmdstanr.R  01_qc_cutadapt.R 02_dada2.R 03_decontam.R
│        04_baseline_E.R 05_occupancy_msocc.R 06_occupancy_occumb.R 07_reclass_B.R
├── data/  download_ena.sh  download_published.sh  verify_ena_join.py
│          ena/PRJEB48362_read_run.tsv        (committed; 1372 study runs)
│          metadata/  conwy_12S_libraries.tsv  design_gaps.tsv  ena_join_verification.md
│                     data_dictionary.md  supplementary_data_1_covariates.csv (git-ignored)
│          published/  source_data_species_x_replicate.xlsx (git-ignored)
│                      supplementary_code_1/...  (git-ignored)
├── docs/  execution-plan.md  prereg.md  supplementary_information.pdf
├── report/  conwy-reanalysis.qmd  _quarto.yml
├── shiny/  app.R          stan/ (.gitkeep)        figures/ (.gitkeep)
```

**Verified working (do not redo):**
- Env builds; **all 21 target R packages load** (dada2, phyloseq, decontam, ALDEx2, ANCOMBC,
  corncob, vegan, iNEXT, indicspecies, zCompositions, msocc, occumb, cmdstanr, posterior,
  bayesplot, shiny, quarto, renv, tidyverse, rjags, jagsUI).
- **cmdstanr smoke test passes** (`R/smoke_test_cmdstanr.R`, compile + sample ~17 s). The plan's
  named #1 risk (Stan on Apple Silicon) is closed — **for Stan**. (See P2: this is the *stretch*
  engine, not the flagship.)
- **ENA join verified & committed** (`data/verify_ena_join.py` → `ena_join_verification.md`):
  exactly **755** Conwy 12S runs (of 1372), strict 1:1 ID↔run, all paired-end, **43 absent design
  cells** (T16×42 + E02_T17_R3). Footprint: 1510 FASTQs, ~70.9 GB, 495,812,380 reads.
- **Covariate column names verified** against the real CSV header — every column the prereg locks
  exists (see P3 for the BOM caveat).
- **Package install gotchas resolved & documented:** dada2 + iNEXT built from source (broken
  osx-arm64 conda builds); msocc from GitHub tarball (CRAN-archived); occumb via conda `r-rjags`
  + `jags` (source rjags fails to link); corncob from CRAN.

**Status:** R analysis scripts `01_`–`07_` are **annotated stubs that `stop()`** (spec + reading
inline, no implementation) — appropriate for Phase 0. Raw FASTQs not downloaded (correct).
Published supplementary files are present on disk (git-ignored).

---

## 2. Problems to address — priority order

> Severity key: 🔴 blocks/derails Phase 1 · 🟠 correctness/credibility, fix before relying on it ·
> 🟡 cheap hardening, do opportunistically. **P1 gates everything — do it first.**

---

### 🔴 P1 — Confirm the published Source Data actually contains a replicate-level 12S matrix
**This is the gating unknown. The entire "staged / fast-path" plan (Decision 4: published
species×replicate table → quick E + first occupancy, *before* bioinformatics) assumes a clean,
analysis-ready, per-replicate detection table. Nobody has opened the file to confirm it exists.**

**Evidence.** `data/published/source_data_species_x_replicate.xlsx` has **only six figure-keyed
sheets**, none named for a species×replicate matrix:
```
'Figure 1b, 2a&d, 3b&e, 4a&b'   'Figure 5a,d&g, 6a,d&g'
'Figure 1c, 2b&e, 3c&f, 4c&d'   'Figure 5b,e&h, 6b,e&h'
'Figure 1d, 2c&f, 3d&g, 4e&f'   'Figure 5c,f&i, 6c,f&i'
```
The then-current plan §10 asserted "12S Source Data is label-collapsed (68 source-study species-label columns × **665** Conwy replicate
rows)", but `data/metadata/data_dictionary.md` §5 hedges ("the matrix is *extracted in Phase 1*"),
and Phase 0 never actually located or extracted it.

**Two specific risks that would break Phase 1 if true:**
- **(a) Resolution.** If the figure data is collapsed to site×time (or to presence/absence across
  reps), then "first occupancy from the published table" via `msocc` is **impossible** — occupancy
  detection histories need the 3 replicates distinguished. The whole staged rationale (do A+E
  before the 71 GB reprocessing) would collapse and occupancy would have to wait for Phase 2.
- **(b) Row reconciliation.** **755 deposited libraries vs the plan's claimed 665 rows** = a
  90-row gap that is asserted but never explained. This will surface as a join failure in Phase 1.

**Concrete solution.**
1. Write a small, committed inspection script — Python (preferred; no env needed, matches the
   "Python for glue" convention) at `data/inspect_source_data.py`, or R `readxl::read_excel`
   (readxl ships with the conda `r-tidyverse`; verify with
   `micromamba run -n edna Rscript -e 'library(readxl)'`). It must, for **each sheet**: print
   dimensions, the first ~5 column headers and ~5 row labels, and detect whether any sheet is a
   species (rows/cols) × `E##_R#_T##` (cols/rows) matrix.
2. **Identify the 12S fish detection block** specifically (the workbook also covers CO1/18S and
   other figure panels). Cross-reference `docs/supplementary_information.pdf` and the authors'
   script `data/published/supplementary_code_1/Supplementary Code 1/Analysing_BLAST_output_WP2_12S.R`
   to learn how they constructed the published table and which figure(s) it feeds.
3. **Assert replicate resolution:** confirm rows/columns carry distinct `R1/R2/R3` per site×time.
   If only site×time (no reps) → STOP and escalate (see decision gate).
4. **Reconcile the row count** against `data/metadata/conwy_12S_libraries.tsv` (755 IDs): join the
   published table's sample IDs to that file on `ID`; report matched / published-only / library-only.
   Explain the 755-vs-665(-or-actual-N) gap (expected causes: QC-failed libraries dropped, control
   rows, or a timepoint subset). Document the exact derivation.
5. **Write the finding** into `data_dictionary.md` §5 (replace the hedge with verified facts: which
   sheet, which columns, exact N rows, exact N species, resolution) and into a committed
   `data/metadata/source_data_verification.md`.

**Decision gate (record the outcome explicitly):**
- ✅ *Replicate-level matrix exists* → fast path is viable. Document sheet + column map; Phase 1
  (`04_baseline_E.R`, `05_occupancy_msocc.R`) can proceed.
- ❌ *Collapsed below replicate resolution* → occupancy cannot start from the published table.
  Escalate to the owner: either (i) reorder to do Phase 2 (raw reads → ASVs) before Layer A, or
  (ii) confirm whether replicate-level detections live in a different supplementary file. Update
  the plan's staging assumption and the prereg accordingly. **Do not silently work around it.**

**Definition of done.** A committed inspection script + `source_data_verification.md` that states,
with evidence, exactly which sheet/columns hold the 12S fish detections, at what resolution, with
what N(species) × N(samples), and how those samples reconcile to the 755-row join table — and a
recorded fast-path go/no-go decision.

---

### 🟠 P2 — Add a JAGS/occumb *sample* smoke test (the flagship engine was only load-checked)
**The smoke test you have exercises cmdstanr/Stan — but Stan is only the *stretch* (custom dynamic
occupancy). The *core Phase-3 deliverable*, `occumb`, runs on JAGS** (occumb → jagsUI → rjags →
JAGS), and that chain was verified only by "the package loads." The asymmetry is backwards: the
stretch engine got a full compile+sample proof; the flagship got none — despite the env notes
recording a real rjags↔JAGS *"Runtime link error"* that had to be worked around with conda
`r-rjags` + conda `jags`. A link/runtime failure here would not surface until deep in Phase 3.

**Concrete solution.** Create `R/smoke_test_jags.R`, mirroring `smoke_test_cmdstanr.R`:
1. **Prefer testing the exact API:** run `occumb`'s own example through `occumb()` with tiny
   settings. Check `?occumb` / the package vignette for the bundled example data and the
   `occumbData()` constructor; fit with minimal MCMC (e.g. `n.chains = 2`, a few hundred
   iterations) so it runs in seconds. This proves the full occumb → jagsUI → rjags → JAGS path
   *and* the API you'll use in `06_occupancy_occumb.R`.
2. **Fallback if the occumb example is awkward:** fit a trivial single-species occupancy (or even a
   coin-flip) model directly via `jagsUI::jags(...)` to prove rjags can compile+sample against the
   conda JAGS.
3. Assert: a fit object returns, parameter estimates are finite, Rhat is computed and < ~1.1, and
   **no "Runtime link error"** appears. Print versions (occumb, jagsUI, rjags, JAGS) and elapsed.
4. Wire it into `00_setup.R` after the cmdstanr smoke test (so `setup_env.sh` proves *both*
   Bayesian backends), and have both smoke tests write a one-line **receipt** to
   `docs/smoke_test_receipts.md` (version + elapsed + PASS) so the README's "smoke test passes"
   claim is auditable like the ENA gate — not ephemeral stdout.

**Definition of done.** `R/smoke_test_jags.R` exists, is called by `00_setup.R`, returns a valid
occumb (or jagsUI) fit on trivial data, and writes a committed PASS receipt. The flagship engine is
now proven end-to-end, not assumed.

---

### 🟠 P3 — Make the covariate join a committed, BOM-safe, re-runnable assertion
Three coupled defects around the single most important key in the project (`ID`):

**(a) UTF-8 BOM on the `ID` column — latent join corruptor.** First bytes of
`data/metadata/supplementary_data_1_covariates.csv` are `EF BB BF` then `ID,Bl…`. A naive
`read.csv("...csv")` in R names the first column `ï..ID` / `﻿ID`, so a join on `ID` silently
yields zero matches. Fix the *convention*, not just one call site: standardize all CSV reads on
`readr::read_csv()` (strips the BOM) or `read.csv(..., fileEncoding = "UTF-8-BOM")`, and note the
rule in `data_dictionary.md`.

**(b) The headline "755/755 covariate join" claim has no re-runnable backing.**
`verify_ena_join.py` only checks the ENA TSV against the `E##` regex/design grid — it **never reads
the covariate CSV**. The most-cited integration fact in the docs is reproducible by zero scripts.

**(c) The covariate CSV is git-ignored** (`.gitignore`), so a fresh clone can't even re-derive the
join. The file is small (294 KB) and CC BY 4.0; committing a *derived, Conwy-only, modelling-column*
subset is consistent with the repo already committing the ENA TSV and the derived join tables (the
"don't commit raw data" rule is about the 71 GB FASTQs, not metadata).

**Concrete solution.** Add `data/verify_covariate_join.py` (or extend `verify_ena_join.py`) that:
1. Reads the covariate CSV with `encoding="utf-8-sig"` (drops the BOM).
2. Filters Conwy rows (`River == "Conwy"` and/or `ID` matches the `E##_R#_T##` regex).
3. **Asserts 755/755 join** to `conwy_12S_libraries.tsv` on `ID` (report matched / left-only /
   right-only — must be 755/0/0).
4. **Asserts no NA** in the prereg's locked modelling covariates: `Dist_from_lake`,
   `river_gradient_percent`, `season`, `week_before_flow_by_catchment`, `Water_temp_oC`
   (and `river_section`, `Latitude`, `Longitude`). Report any NAs by column.
5. Writes (i) a committed derived `data/metadata/conwy_covariates.tsv` = `ID` + the modelling
   columns for the 755 rows (the analysis-ready, BOM-free, version-controlled covariate table that
   R scripts consume), and (ii) `data/metadata/covariate_join_verification.md`.
6. Have R scripts read `conwy_covariates.tsv` (committed, clean) rather than the raw publisher CSV.

**Definition of done.** A re-runnable script proves 755/755 join with no NA in the locked
covariates; a committed BOM-free `conwy_covariates.tsv` exists; the BOM-handling rule is documented;
the "755/755" claim is now backed by code, not prose.

---

### 🟡 P4 — Reconcile the plan's `renv` mandate with reality (decision is settled: conda is the backbone, not renv)
The **plan** commits to `renv` as *the* reproducibility backbone in six places (`eDNAplan.md`
lines 33, 78, 116, 128, 141 — *"One environment (`renv`)"* — and 165; copied into
`docs/execution-plan.md`). **The README does NOT mention renv** (it is already conda-honest — do not
"fix" it). Reality: there is no `renv/` and no `renv.lock`; reproducibility is 100% conda. The
`renv` package is installed but never called.

**This was tested, not assumed — and the test settles the choice (do not re-litigate):**
- `environment.lock.yml` is **417 build-string-exact entries** — the full transitive closure of the
  conda layer, including the three things **renv structurally cannot manage**: `r-base=4.4.3`,
  `jags=4.3.2` (the flagship `occumb` engine), `cmdstan=2.39.0` (the stretch engine), plus the
  compiler toolchain. **renv pins R *packages* only**, so it can never be the "one environment" the
  plan names — it cannot reproduce R itself, JAGS, or CmdStan. Conda already does, *and* exact-pins
  19 of the 21 R packages transitively.
- Bioconductor is pinned to **3.20 by construction** (BiocManager keys off R 4.4), so `dada2`/`ALDEx2`
  /`ANCOMBC`/etc. are version-stable as long as `r-base` stays 4.4 — the earlier "dada2 could drift"
  concern was **overstated**. The real residual holes are narrow (see below).
- The one thing renv is *supposed* to do better here — capture the `msocc` GitHub SHA — **does not
  materialize**: `msocc` was installed from a bare `master` tarball, so its DESCRIPTION has no
  `RemoteSha`/`GithubSHA1`. `renv::snapshot()` would record no SHA either. Both approaches need the
  identical fix.

**Objective verdict:** renv cannot serve as the backbone (system-dep gap) and conda already covers
what renv aims at, *more completely*. So **make conda the single source of truth and harden it** —
do **not** adopt renv as the primary mechanism. (renv remains a reasonable *optional supplement*
**only if** the Posit Cloud x86_64 fallback later becomes a primary workflow, since the conda lock
is arm64-only and `renv.lock` is platform-agnostic for R packages. Flag that to the owner; don't
build it pre-emptively.)

**The narrow residual holes to close (this is the actual work):**
1. `msocc` is installed from `refs/heads/master` (a moving target) with no recorded commit. In
   `00_setup.R`, pin the tarball URL to a specific commit SHA. Current master is
   **`0180cb0489da46208079ae31aac46f6d83c3e058`** (the SHA originally installed is unrecoverable;
   pin current master and note that).
2. The source-built `iNEXT` and CRAN top-ups (`corncob`, `occumb`, `jagsUI`, …) are installed
   "latest" — CRAN serves only the newest, so they can drift on a future rebuild. After the
   roll-call in `00_setup.R`, write a committed `environment.R-packages.lock.tsv` capturing
   `installed.packages()[targets, c("Version")]` for all 21 packages + the msocc SHA + Bioc release
   (3.20). This is the human-readable, reviewer-facing version record for the preprint methods.
   *(Optional, only if exact CRAN-version restore is wanted: switch those installs to
   `remotes::install_version(..., version = <pinned>)`.)*
3. **Reconcile the plan** (`eDNAplan.md`): replace the renv mandate with
   "conda env (`environment.lock.yml`) + pinned R-package manifest" so the stated contract matches
   the (better) mechanism. Leave the README as-is.

**Definition of done.** The plan no longer mandates a backbone that can't work here; `msocc` is
pinned to a commit SHA; a committed `environment.R-packages.lock.tsv` records exact versions of all
21 targets + Bioc release; a clean `setup_env.sh` run is version-deterministic. renv is explicitly
recorded as a deferred, optional, cross-platform-only supplement — not adopted.

---

### 🟡 P5 — Encode the T16 / usable-occasions rule as a machine-readable artifact
The "18 usable timepoints; T16 = unsurveyed (NA, not 0)" rule currently lives only in prose +
`design_gaps.tsv` (a raw list of absent cells). When detection-history construction is written
(possibly weeks later in `05_`/`06_`), the rule is likely to be re-derived from memory and the
**structured missingness mishandled** — e.g. T16 encoded as all-zero detections, which manufactures
false non-detections and biases ψ downward (the precise failure this project is built to expose).

**Concrete solution.** Emit a canonical, committed occasion map now, while it's fresh. Minimum: an
`R/lib_design.R` helper (or a `data/metadata/occasions.tsv`) that, from `conwy_12S_libraries.tsv` +
`design_gaps.tsv`, produces the **site × occasion × replicate** design array with status
`surveyed` / `unsurveyed` and **NA (not 0)** at every T16 cell and at `E02_T17_R3`, plus a
`usable_timepoints` vector (T01–T19 minus T16 = 18). Both `05_occupancy_msocc.R` and
`06_occupancy_occumb.R` should consume this single source rather than re-deriving the grid.

**Definition of done.** A committed artifact/helper yields the occasion design with NA (verified
non-zero) at all unsurveyed cells and exactly 18 usable occasions, ready for both occupancy scripts.

---

### 🟡 P6 — Decide the Quarto reproducibility/caching strategy before Phase 2 writes heavy scripts
The plan's contract is "one `quarto render` re-runs everything." That collides with reality once
`02_dada2.R` (71 GB + DADA2) and `06_occupancy_occumb.R` (MCMC) exist — every render would re-run
hours of compute. No caching/checkpointing is set up. The `.qmd` currently keeps Phase-1/3 sources
commented out, so it renders today, but the strategy must be chosen *before* data-producing scripts
land.

**Concrete solution.** Pick and document one: (a) Quarto `freeze: auto` in `_quarto.yml` +
explicit guards on heavy steps (`if (!file.exists(rds)) { … } else readRDS(rds)`); or (b) a
`targets` pipeline (cleaner dependency graph, steeper learning curve — flag the reading for the
owner); or (c) manual `.rds` checkpoints under a git-ignored `data/derived/`. Recommend (a) for the
lowest friction given the owner is new to R. Write the chosen convention into `report/_quarto.yml`
comments and the README.

**Definition of done.** A documented caching convention exists; the "single render reproduces
everything" claim is achievable without re-running multi-hour steps on every render.

---

### Minor cleanup (do opportunistically, not gating)
- **cmdstanr receipt:** fold into P2's `docs/smoke_test_receipts.md`.
- ~~**`docs/execution-plan.md` duplicates `../eDNAplan.md`** → drift risk.~~ **Resolved:**
  `docs/execution-plan.md` was removed; `eDNAplan.md` (now vendored into the repo root) is canonical.
- **Read `Analysing_BLAST_output_WP2_12S.R` before Phase 2/B** and document its expected inputs
  (BLAST output format) and outputs — the "reproduce the authors' 2020 baseline" gate depends on it.
  Not needed for P1–P6, but cheap insurance against a Phase-2 surprise.

---
