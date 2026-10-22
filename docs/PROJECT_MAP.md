# Project map — River Conwy fish-eDNA reanalysis

**For:** whoever picks this up next.
**Map date:** 2026-10-21 · **Repo state:** git clean, HEAD `1a00a9f` "Complete local release and Phase 5 credibility audits".
**One-line:** a reproducible *reanalysis* + credibility audit of the Perry et al. 2024 River Conwy 12S fish-eDNA dataset, packaged as a local static website ("portfolio edition"). Not published, not deployed.

---

## 0. The three things that bite you first

1. **`E##` vs `C##`.** Conwy = site codes **`E01`–`E14`**. The same ENA archive holds `C##` = a *different river*. Pulling `C##` silently gives wrong data. The download filter + `data/verify_ena_join.py` hard-enforce `E##`. Never relax this.
2. **Layer A numbers are diagnostic, not inferential.** Every occupancy/occumb number in the repo came from a 2-chain/2000-iter **smoke** config and/or failed the convergence gate (`WITHHELD_UNCONVERGED`). The gate (`scripts/apply_layerA_claim_gate.R`) replaces stale numeric figures with withheld-status cards and the release verifier fails if an old number reappears. See `docs/phase3-claim-status.md` for what may and may not be claimed.
3. **Two-machine history.** Phase 0–2 raw-read work ran on a WSL2 PC + USB hand-back drive; Phase 3+ audit layers ran on the Mac. The repo here is the Mac side. Production MCMC and full FASTQ reruns are the PC job and are **not** done.

---

## 1. What this project is

The source study (Perry 2024, *Nat. Commun.* 15:4372) supplied Conwy samples, reads, published matrix, covariates, methods. This project adds a reanalysis + audit layer that separates **observed evidence** from **claim decisions**, under three analytical layers:

| Layer | Name | Question |
|---|---|---|
| **E** | Compositional community baseline | What community pattern survives when the biological unit + repeated structure are handled honestly? |
| **A** | Bayesian site-occupancy (msocc / occumb) | Can imperfect-detection correction recover truth when read depth and biological signal are confounded? (Currently: exposes what *cannot* be recovered.) |
| **B** | Reference-DB-currency reclassification | Which taxonomic ranks stay stable when the same sequences are re-assigned against current databases? |

**Headline scale:** 14 sites · 18 usable timepoints · 647 validated libraries · 1,348 ASVs.
**Final claim tally:** **1 ACCEPTED · 1 REJECTED · 3 NON-CONCLUSIVE** (pre-registered decision rules, see §4).

---

## 2. Directory map

```
conwy-edna-reanalysis/
├─ README.md               project entry point + the E##/C## rule
├─ START_HERE.md           5-min reproducible route (website build only)
├─ eDNAplan.md             the full locked plan (phases 0–5)
├─ EDIT_WEBSITE.md/.command no-HTML text editor for the website
├─ CITATION.cff / LICENSE  citation metadata; all rights reserved
│
├─ R/                      the analysis pipeline — numbered 00→15 (see §3)
│   ├─ 00..15_*.R          ordered stages
│   ├─ lib_design.R        shared helpers
│   ├─ smoke_test_*.R      cmdstanr / jags environment smoke checks
│   └─ test_*.R            per-module checks (msocc mapping, atlas, multiverse, …)
│
├─ scripts/                orchestration + build (Python + a little R)
│   ├─ build_website.py      rebuilds website figures + PDFs (small inputs)
│   ├─ package_website.py       stages + verifies the shareable site zip
│   ├─ verify_release.py        release gate (links, checksums, claim gate)
│   ├─ apply_layerA_claim_gate.R enforces "no stale Layer A numbers"
│   ├─ build_biological_evidence.R / build_claim_survival.py
│   ├─ run_impact_sprint.py / export_impact_inputs.R
│   ├─ fetch_salmo_multimarker.py / edit_website.py
│
├─ data/
│   ├─ derived/   (78) analysis outputs: seqtab, summaries, estimates, gate tables
│   ├─ metadata/  (42) site/sample maps, gate status (layerA_*.tsv), traits
│   ├─ published/ (3)  Perry source matrix/covariates
│   ├─ ena/ public/    ENA join + fetched reference material
│   └─ verify_ena_join.py  ← the E##/C## enforcement
│
├─ stan/                   (empty) — production Stan models go here (PC job)
├─ shiny/app.R             exploratory app: corrected-vs-raw, diagnostics
├─ report/                 Quarto technical report (conwy-reanalysis.qmd + _freeze)
│
├─ website/             THE LOCAL WEBSITE (11 HTML pages)
│   ├─ index.html, research_brief.(md|html|css), visual-brief.html
│   ├─ deliverables.html   catalogue linking every named output
│   ├─ claim-survival.html, detection-lab.html, poster.html, editor.html
│   ├─ graphical_abstract.svg, assets/, data/, content-overrides.json
│   └─ build-receipt.json
│
├─ output/
│   ├─ pdf/    RESEARCH_BRIEF.pdf, POSTER.pdf
│   └─ share/  conwy-edna-project-website.zip (+ .sha256, .bundle.json) ← the deliverable
│
├─ docs/                   methods, prereg, claim status, handoff notes (see §5)
├─ figures/                generated figures
├─ environment.yml / .lock.yml / environment.R-packages.lock.tsv / setup_env.sh
└─ tmp/                    scratch (brief-redesign, pdfs, quarto-shim, svg-preview)
```

---

## 3. Analysis data flow (R pipeline)

Numbered stages, run in order. Early stages (00–03) are the raw-read pipeline (PC history); 04+ is the audit layer (Mac).

```
00_setup → 01_qc_cutadapt → 02_dada2 → 03_decontam
                                          │  (→ 647 libraries × 1,348 ASVs: data/derived seqtab)
                                          ▼
                     04_baseline_E  ──────────────► Layer E compositional baseline
                                          │
          ┌───────────────────────────────┼───────────────────────────────┐
          ▼                               ▼                               ▼
   05_occupancy_msocc          06_occupancy_occumb               07_reclass_B
   (2-level occupancy)         (community occupancy; smoke)      (DB-currency reclass)
          │                               │                               │
          ▼                               ▼                               ▼
   08_validate_concordance   09_e14_covariates   10_sim_recovery (→ WITHHELD_UNCONVERGED)
                                    │                    │
                                    ▼                    ▼
   11_occupancy_recoverability_atlas   12_salmo_multimarker_benchmark
   13_layerE_multiverse (72 families)  14_trait_null_model   15_pipeline_disagreement
```

Then: `scripts/build_biological_evidence.R` + `build_claim_survival.py` → claim tables → `build_website.py` → website/PDFs → `package_website.py` + `verify_release.py` → `output/share/`.

**Gotchas inside the pipeline (from memory + docs):**
- `06_occupancy_occumb.R`: has memory/EINTR robustness fixes; `OCCUMB_SAVE_COV_FIT` gate controls saving the covariate fit.
- occumb structurally drops all-zero community vectors → 18 libraries that the naïve comparator counts as non-detections are masked. This observation-mask mismatch is a core reason Layer A shifts aren't inferential.
- `10_sim_recovery.R`: the depth safeguard **over-corrects** rare ψ — the key negative finding. Writes `data/metadata/layerA_recovery_status.tsv`.
- E14 covariates: supporting climate metadata uses E13 as donor (documented Gate-1 fallback).
- `Catchment_flow` is a time series, not a scalar.

---

## 4. Claim outcomes (the scientific bottom line)

| Claim | Verdict | Why |
|---|---|---|
| Pipeline disagreement (`SURVIVES`) | **ACCEPTED** | comparator-relative association met every stability check (non-causal; not a false-negative rate) |
| Trait diversity broadens downstream (`FAILED_GATE`) | **REJECTED** | fixed-margin null: one-sided p = 0.3159 (strict 0.1142) |
| River community section shift (`SURVIVES_CONDITIONAL`) | **NON-CONCLUSIVE** | direction stable but 64/72 multiverse families valid (<90% rule); ref p = 0.0714 |
| Longer markers resolve *Salmo* species (`WITHHELD`) | **NON-CONCLUSIVE** | only 2 of 7 candidate labels eligible; genus stable, species not |
| Imperfect-detection recovery (`WITHHELD`) | **NON-CONCLUSIVE** | atlas maps design limits; real Conwy model did not converge (`WITHHELD_UNCONVERGED`) |

Authoritative sources: `data/metadata/claim_survival.tsv`, `docs/phase3-claim-status.md`. Do not upgrade any of these from the smoke fits.

---

## 5. Key docs to read, in order

1. `README.md` → `eDNAplan.md` — what/why/plan.
2. `docs/phase3-claim-status.md` — **what you are allowed to claim** (read before touching Layer A).
3. `docs/phase3-readiness.md` — historical hand-back audit (GO-WITH-CAVEATS); note the 2026-07-18 current-status banner.
4. `docs/prereg.md` + `docs/impact-sprint-preregistration.md` — the pre-declared decision rules.
5. `docs/impact_*_methods.md` (×5) — methods per audit (atlas, multiverse, pipeline-disagreement, salmo, trait-null).
6. `docs/handoff-phase0-fixes.md`, `docs/smoke_test_receipts.md` — environment + early-stage gotchas.
7. `docs/preprint-draft.md`, `report/conwy-reanalysis.qmd` — the manuscript/report scaffold.

---

## 6. Open items / what's NOT done

**Deferred research work** (deliberately out of scope for the local edition):
- Production-scale MCMC (longer chains) — the **PC job**; lands in `stan/` (currently empty). Note: longer chains fix Monte-Carlo convergence, **not** the structural recoverability problem the simulation exposed.
- Full FASTQ rerun, reference-DB rebuild, hosted Shiny server, peer review / journal submission.

---

## 7. How to run

**5-min reproducible route (no FASTQs, no MCMC):**
```bash
python3 -m venv .venv-website && source .venv-website/bin/activate
python3 -m pip install -r requirements-website.txt
python3 scripts/build_website.py      # rebuilds website figures + PDFs
```

**R environment** (full pipeline): isolated micromamba env `edna` (Mac, native-arm64) / `edna-phase2` (PC). Run R stages via `micromamba run -n edna Rscript R/NN_*.R`. Env pinned in `environment.yml` / `environment.R-packages.lock.tsv`; `setup_env.sh` bootstraps.

**Verify a release:** `python3 scripts/verify_release.py` (checks links, checksums, and the Layer A claim gate).

**Explore interactively:** `shiny/app.R`.
```
```
