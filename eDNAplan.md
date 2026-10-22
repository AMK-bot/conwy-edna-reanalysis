# eDNA Metabarcoding Reanalysis — Full Execution Plan

**Status:** ✅ Locally complete website release and verified Phase-5 scientific-impact extension.
**Derived from:** `eDNA-reanalysis-project-scoping.md` via a full decision-tree grilling.
**Owner profile:** Programmer (Python + bash); learning R, Bayesian stats, and ecology from scratch — *recommended reading is embedded at every gap.*
**Original planning deadline:** January 2027. The locally achievable release scope was completed on 2026-07-18; optional owner-authorized public actions are not on the critical path.
**Last updated:** 2026-07-19 — *Phases 0–5 and the presentation sprint are complete for the local-release scope. Phase 5 adds five deterministic scientific audits plus a claim-survival map, with failed gates and withheld claims retained as results. Layer A remains complete as a failed-claim-gate result rather than a final rare-species estimate.*

### Local-release completion ledger

Percentages use each phase's declared locally achievable scope as its denominator. Phase 5 is a new post-release extension with its own denominator; adding it does not invalidate the verified Phase 0–4 release.

| Phase | Completion | Verified local checkpoint | Missing local subphase |
|---|---:|---|---|
| Phase 0 — setup/provenance | 100% | Environment, strict ENA join, occasions, covariate/data-source audits | None |
| Phase 1 — detection occupancy on-ramp | 100% | Six well-detected analysis-label models, structural missingness, mapping test, bounded interpretation | None |
| Phase 2 — reads-to-ASV handback and validation | 100% | Healthy 647 × 1,348 count table and PASS-WITH-CAVEATS concordance gate | None; a fresh 71 GB clean-clone rerun is optional, not needed for this local checkpoint |
| Phase 3 — credibility corrections (E/A/B) | 100% | Site/time restricted inference; rare-occupancy claim gate; current-resource genus-level concordance | None; production/sensitivity MCMC is optional research and cannot repair the identified mask/confounding issue by itself |
| Phase 4 — local synthesis/release | 100% | Synchronized report/site, documentation, local explorer schema, citation/license notices | None |
| Presentation sprint | 100% | Project website with complete deliverables catalogue, self-contained paper-style research brief, four-page visual PDF with five model architecture diagrams, molecular vignette, ecology table, simulator, poster, storyboard, interview set, <5-minute capsule | None; personal contribution text and voice/video require the owner |
| Phase 5 — scientific-impact extension | 100% | Frozen local specification; recoverability atlas; accession-backed marker benchmark; Layer-E multiverse; trait null; clustered pipeline disagreement; claim-survival map; deterministic release verification | None |

> ### Phase-0 reality — three findings that reshaped Phases 1–4
> Phase 0 was completed **and** hardened (the P1–P6 review in `docs/handoff-phase0-fixes.md` is all resolved). Three verified findings changed the downstream plan:
> 1. **The published 12S table is _normalised relative abundance_ (rows sum to 1), not raw read counts.** Detection (>0 vs 0) survives, so occupancy is unaffected — but the count-based DA models, any abundance-mode iNEXT analysis, and `occumb` (which models read counts) **cannot run on it.** ⟹ **Layer E moves out of Phase 1 and runs once, fully, on the reprocessed counts in Phase 3;** its executed iNEXT analysis uses incidence frequencies across site×time units. The published table is used only for `msocc` detection histories and as the Phase-2 validation target.
> 2. **Compute is now two machines.** The 8 GB Mac can't hold the 71 GB download + DADA2; a **Windows + WSL2 PC (64 GB RAM)** runs Phase 2, the Mac runs everything else, and a USB stick (`PS/eDNA`) carries the code + the small derived outputs back. **Phase 1 (Mac) and Phase 2 (PC) now run in parallel** — Phase 1 needs none of Phase 2's output. Posit Cloud is no longer needed.
> 3. **Covariate gaps (final local disposition):** source E14 climate/gradient values are absent. A flagged, unvalidated E13 donor exists only in supporting metadata; climate, flow, and gradient are excluded from the executed headline formula rather than fabricated.

> ### Phase-2 reality — completion & validation (2026-06-28), and the one caveat that now governs Phase 3
> Phase 2 (raw reads → ASVs on the WSL2 PC) is **DONE and adversarially validated.** Four findings reshape Phase 3:
> 1. **The count table is healthy and already on the Mac.** `data/derived/seqtab_clean.rds` = **647 libraries × 1,348 ASVs, 35.7 M reads** (dims re-verified 2026-07-05), collapsed to `species_x_sample_counts_rawnames.rds` (**35 columns: 34 assigned labels + `Unassigned`**) / `_clean.rds` (14 harmonised labels) + a genus table. *(The 34 field-assigned labels reconcile with the broader taxonomy universe's 37 assigned labels: 3 taxonomy labels have no surviving field ASV; `Unassigned` is the 35th collapse column.)* This is `occumb`'s read-count input — **Blocker 0 is genuinely cleared.** *(An earlier handback shipped a degenerate 375 × 829 / 173 k-read `seqtab_clean.rds` from a **tag-jump cutoff = 1.0** bug; it was regenerated on the PC with a 0.001 floor. ⚠️ The red-flag tables in `docs/phase3-readiness.md` §3–§4 describe that **superseded corrupt state** — do not read them as current.)*
> 2. **Validation gate = PASS-WITH-CAVEATS** (`docs/phase2-validation-gate.md`). `R/08_validate_concordance.R` bridges the 2020↔2025 database namespaces with a consensus **crosswalk** (`data/derived/crosswalk.tsv`), then compares the published and reprocessed tables over **642 shared libraries**: dominant composition agrees well (mean read-weighted Bray–Curtis **similarity** 0.84, cosine similarity 0.90), all 10 source-study taxa detected in ≥3 libraries are recovered, cell precision is 0.99, and the 5.9% Unassigned bin is dominated by non-target vertebrate bycatch. Trace one- or two-library label disagreements remain taxonomically uncertain rather than being treated as confirmed local species.
> 3. **The governing caveat — a ~16% detection-cell recall deficit creates a recoverability risk.** Reprocessed detections run at **~84% of the source-study cell total** over recovered taxa (selected label-specific deficits 22–29%), and misses track source-study relative abundance and library depth along the river gradient. Layer A therefore tested a depth-aware detection model, structural NA, and known-truth recovery. The recovery MCMC failed convergence and an 18-library comparison-mask inconsistency remained for vectors that are zero only within the 14-label harmonised input, so **the direction and magnitude of sparse-label occupancy adjustment are withheld**. The failed gate, not a claimed correction, is the completed result.
> 4. **Layer B's 2020 baseline is constrained and the choice is resolved.** The authors' **MitoFish build (23 Sep 2020) is unarchived**, so a faithful 2020 reproduction is not claimed. The completed layer is current-resource concordance: three comparisons across two current reference resources, supported by a current public-record witness set. The stable assignment is genus *Salmo*; the short fragment does not establish species identity or a 2020→2026 taxonomic change.
>
> **Provenance / data-loss note (updated 2026-07-18):** the Phase-2 scripts and earlier documentation were committed (`08eae8e`), but that commit does **not** cover the count table. `data/derived/*.rds` (including `seqtab_clean.rds`) is git-ignored and regenerable only by re-running the 71 GB PC process. Phase 3 has consumed the local hand-back; keep the USB copy as an ongoing off-Mac backup.

---

## 0. One-paragraph definition

> **"Imperfect detection and reference-database currency in River Conwy fish eDNA: a spatio-temporal reanalysis."**
> Take the *Nature Communications* River Conwy 12S fish-eDNA dataset and audit three distinct inferential layers: **(E)** a compositionally aware community baseline with the biological unit and repeated structure handled correctly; **(A)** occupancy/detection modelling that reports results for six well-detected analysis labels but demotes sparse-label estimates when recovery fails; and **(B)** current-database concordance that identifies the stable taxonomic rank. Ship a transparent local technical record plus a five-minute, production-MCMC-independent portfolio edition. Public release, preprint, and deployment remain optional owner actions.

**Contribution boundary:** The original Conwy analysis did not include occupancy/imperfect-detection modelling. This repository adds that diagnostic layer and a database-concordance audit, but makes no claim about being the only later reanalysis and no claim of publication.

---

## 1. Final decision ledger

| # | Decision | Locked choice | Notes |
|---|----------|---------------|-------|
| — | **Scope** | Full **A + B + E** core, staged, no time-box | Stretch (labeled): dynamic occupancy, D, I. Future notes: community layer |
| — | **Deadline / pace** | Jan 2027, ~8–12 hrs/wk | Comfortable slack; phases below are effort-estimated, not calendar-locked |
| — | **Skills posture** | Learn R + Bayesian + ecology as I go | Reading embedded per phase (§7) |
| 1 | **Dataset** | **River Conwy** — ENA **PRJEB48362** | ⚠️ Conwy = site codes **`E01`–`E14`** (`T01`–`T19`, `R1`–`R3`). **NOT `C##`** (different river in same archive) |
| 2 | **Marker** | **12S fish (MiFish-U)** | Amplicon ~172 bp. COI/18S available for a future multi-marker extension |
| 3 | **Angles** | **A** (occupancy) **+ B** (DB-currency reclassification) **+ E** (correct baseline) | Stretch: **D** (pipeline benchmarking), **I** (river-transport), dynamic occupancy |
| 4 | **Start point** | **Staged, two machines, parallel** | Mac runs **msocc** (Phase 1) on the published *detection* table; PC reprocesses **raw 12S reads** (Phase 2) → ASVs; Phase 3 joins them for **occumb + Layer E + B**. Phase 1 ∥ Phase 2 (different machines, no dependency). *(Published table is normalised relative abundance → detection-only; full E + occumb need the reprocessed counts.)* |
| 5 | **Language / stack** | **R-centric** | DADA2-in-R + phyloseq/vegan/iNEXT + msocc/occumb. Python/bash for ENA download + glue |
| 6 | **Compute** | **Two machines** | **Mac (8 GB, arm64 `edna` env)** = light work (msocc, occumb, Layer E on counts, Shiny, writing) + git repo. **Windows+WSL2 PC (64 GB)** = the RAM-heavy Phase 2 (71 GB download + DADA2), env `edna-phase2`. **USB `PS/eDNA`** carries the Phase-2 code + derived outputs. `cmdstanr` on Apple Silicon (not `rstan`); **Posit Cloud no longer needed** (the 64 GB PC covers DADA2 and overnight MCMC). |
| 7 | **Deliverables** | **Local evidence release + portfolio edition** | Repo, Quarto technical record, exploratory local Shiny app, static website pages/PDFs/simulator. Public preprint, deployment, release, and DOI are optional owner actions |
| — | **Occupancy design** | **Static multiscale *community*** (`occumb`); single-species (`msocc`) to learn | Dynamic (colonization/extinction) = labeled stretch |
| — | **Community contribution** | **Out of core scope** | Preserved as future-expansion notes (§9) |

---

## 2. The scientific design (what we are actually computing)

### 2.1 Layer E — the correct baseline (community ecology, done properly)
> **Runs in Phase 3, on the reprocessed read counts — not on the published table.** The DA trio requires counts; the published Source Data is normalised relative abundance, so feeding it shares is invalid. Layer E is therefore run **once** on the Phase-2 ASV→label counts. The planned abundance-mode iNEXT route was rejected because reads are not ecological individuals; the executed analysis uses incidence-frequency iNEXT across site×time units. For compositional beta diversity, robust Aitchison uses rCLR→Euclidean and handles zeros without `zCompositions` replacement.

Run a community analysis with explicit compositional, replication, and claim-boundary checks:
- **Alpha diversity:** richness, Shannon, Simpson, Pielou, **Hill numbers**, with **coverage-based rarefaction/extrapolation (`iNEXT`)** on site × time incidence — not naïve rarefying or replicate-library pseudoreplication. (No defensible phylogenetic tree for this short 12S panel → skip Faith's PD and say why.)
- **Beta diversity:** Bray–Curtis, **Jaccard** (presence/absence — defensible given read-count caveats), and **robust Aitchison** (rCLR → Euclidean) as the compositional option.
- **Ordination + tests:** PCoA/NMDS; **PERMANOVA paired with PERMDISP**, after replicate libraries are aggregated. Spatial inference uses 14 site profiles over 11 common timepoints; seasonal inference uses 18 timepoint profiles over a fixed 10-site panel. Exact cyclic/mirrored river or chronological series restrict p-values.
- **Differential abundance:** consensus of **≥2–3 FDR methods** (ANCOM-BC2 + ALDEx2 + `corncob`) plus a restricted site-series BH gate. Package-only consensus is never a final hit by itself.
- **Indicator species:** `indicspecies::multipatt` on site occurrence frequency with the restricted spatial series.
- **Compositional hygiene throughout:** rCLR for robust-Aitchison beta diversity, ALDEx2 within the multi-method DA screen, and a documented +0.5 CLR only for the restricted site-series gate; default to **detection/occurrence** framing (reads ≠ abundance).

### 2.2 Layer A — multiscale occupancy (the flagship)
**Structure (verified fittable):** Each **site × timepoint** is an occupancy unit; up to **3 water-sample libraries** were sequenced separately and are the detection replicates (1–3 survive per realized Phase-3 unit). This is a **2-level** multiscale model (site→replicate). *(PCR triplicates were pooled before indexing, so a third detection level is unavailable.)*

**Model progression:**
1. **Learn (Phase 1, Mac):** single-label static occupancy with **`msocc`** (fast Gibbs sampler) on six well-detected analysis labels → build intuition + verify detection histories. Runs on the published *detection* table (only needs presence/absence), so it proceeds in parallel with Phase 2 and needs no read counts. The labels are the source study's harmonised analysis units, not independent species identifications.
2. **Stress test (Phase 3, Mac, on Phase-2 counts):** community static occupancy with `occumb`. The short fit proves the design and code path, but rare-taxon estimates are unconverged and fail a known-truth recovery study under the observed depth × biological-gradient confound. They are retained only as exploratory diagnostics; the portfolio edition does not depend on them.
3. **Stretch (labeled):** **dynamic** occupancy (colonization/local-extinction across the 19 timepoints) — likely **custom Stan** (no canned eDNA-multiscale-dynamic package), attempted only after the static core lands.

**Executed covariate specification (recorded plan plus documented availability amendment):**
- Occupancy ψ ~ `Dist_from_lake` + `season`.
- Detection θ ~ z-scored `log(read depth)`; `occumb` also uses in-community total reads as its native multinomial exposure.
- `river_gradient_percent` is dropped because E14 is missing; no substitute gradient is invented. A flagged, unvalidated E13 donor exists only in supporting metadata; climate and catchment-wide flow are not in the executed headline formula. See `docs/prereg.md` for the planned set and its recorded deviation.

**Reported outputs:** validated detection estimates for six well-detected analysis labels from Phase 1; smoke-scale corrected-vs-raw and utility tables only as diagnostics; convergence and recovery evidence beside every Layer-A display. No count of "recovered absences," final sparse-label occupancy, or sampling recommendation is allowed from the present Layer-A fit.

**Detection-sensitivity safeguard and disposition.** The ~84% detection-cell recall result established a
confounding risk, not the direction of occupancy bias. The executed model included per-library read depth,
the native read-count exposure, and structural missingness for unavailable cells. The project then ran a
known-truth recovery study. That MCMC failed the Rhat gate (maximum ψ Rhat 6.61; threshold 1.05), and 18
observed libraries that are zero within the 14-label harmonised Layer-A input exposed an inconsistent comparison mask. These libraries contain reads under raw taxonomy labels and are not zero across the raw-name fish table. Therefore neither the recovery
direction nor magnitude is reportable. A separate fixed-seed deterministic simulator illustrates a generic
failure mode without being presented as a Conwy estimate.

### 2.3 Layer B — reference-database currency
**Clean, pre-defined target (verified):** the authors assigned 12S taxonomy with **BLAST vs NCBI nt + MitoFish downloaded 23 Sep 2020**. So "reclassify against a current database" is concretely:
- **Executed comparison:** three classification outcomes spanning two current reference resources (MitoFish via two methods and MIDORI2 as the independent resource), plus a direct current-record witness set for ASV0001.
- **Quantified:** assignment agreement by rank and read weighting, with the dominant sequence audited at base-pair level.
- **Result boundary:** genus *Salmo* is stable; species is not identifiable from the 170-bp fragment/reference set. No invasive, threatened, or local species record is inferred solely from a classifier label.

> **Resolved scope:** the unarchived MitoFish-2020 build prevents a faithful historical replay, and a full current nt build is unnecessary for the local question. Layer B is explicitly current-resource concordance, not evidence of taxonomic change through time.

---

## 3. Phase-by-phase execution plan

> Each phase ends with a **shippable checkpoint**, so there's always a presentable artifact. Effort is in focused work-weeks at ~8–12 hrs/wk; with 7 months you have ~2× slack.

### Phase 0 — Setup, verification & orientation · ✅ **DONE** (incl. P1–P6 hardening)
**Goal:** a working R environment and the *correct* data on disk, with the join validated.
- [x] Stand up the isolated, native-arm64 conda **`edna`** env (`environment.yml` → `environment.lock.yml`) and a **git** repo. *(Conda — not `renv` — is the reproducibility backbone: it pins R itself, JAGS, and CmdStan, which `renv` structurally cannot. See `README` + `environment.R-packages.lock.tsv`.)*
- [x] Install core packages — **all 21 load** (dada2, phyloseq, vegan, iNEXT, decontam, ALDEx2, ANCOMBC, corncob, zCompositions, indicspecies, msocc, occumb, cmdstanr, shiny, quarto, …); two built from source (dada2, iNEXT — broken osx-arm64 conda builds); `msocc` pinned to a commit SHA.
- [x] **Smoke-test both Bayesian engines:** `cmdstanr`/Stan (stretch) **and** the flagship JAGS→`occumb` chain — auditable PASS receipts in `docs/smoke_test_receipts.md`.
- [x] **Verified the right data:** ENA join is strict **1:1, 755 Conwy 12S runs** (E## hard-filter; C## decoy excluded); **43 unsequenced design cells** (T16 ×42 + `E02_T17_R3`) → **18 usable timepoints**, encoded NA-not-0 in `occasions.tsv` + `R/lib_design.R`. Covariate join **755/755**, BOM-safe (`conwy_covariates.tsv`). Published supplementary fetched (raw FASTQs deliberately **not** — that's Phase 2 on the PC).
- [x] **P1 fast-path gate:** published 12S Source Data is replicate-level **but normalised relative abundance** (665 reps × 68 source-study species-label columns), not counts → **detection-only** (see the reality callout up top). 665 = 755 − 90 QC-dropped (<1000 reads); the 90 dropped are NA-not-0 too.
- **Checkpoint:** ✅ repo + reproducible env + validated joined dataset + data dictionary + the P1–P6 hardening (`docs/handoff-phase0-fixes.md`).
- **Reading (done):** [Happy Git with R](https://happygitwithr.com/) · *R for Data Science* (2e) [r4ds.hadley.nz](https://r4ds.hadley.nz/) (you know Python — skim for R idioms) · [ENA file-download guide](https://ena-docs.readthedocs.io/). *(renv reading dropped — conda is the backbone; skim renv only if a cross-platform R-package lock is ever wanted.)*

### Phase 1 — Bayesian on-ramp: single-species occupancy (`msocc`) · ✅ **DONE** (2026-06-24) · *Mac · runs in parallel with Phase 2 on the PC*
**Goal:** start the Bayesian climb on clean detection data — no bioinformatics, no read counts, nothing from Phase 2. *(Layer E is **not** here anymore: the published table is normalised, so the count-based baseline runs once on real counts in Phase 3.)*
- [x] **Build detection histories** from `source_data_12S_normalised.tsv` (detection = value > 0): site×time = occupancy unit, 3 replicates = detections. **Consume `R/lib_design.R` → `observable_mask("published")`** so the 43 unsurveyed cells (T16) **and** the 90 QC-dropped libraries are **NA, never 0** (validated by an equality gate — see deviation D5). 665 cells → 238 units; 68 raw source-study labels → 55 cleaned analysis labels.
- [x] **Apply the authors' pre-analysis cleaning** the published matrix still lacks: drop the `Scomber_scombrus` positive-control spike + the marine/aquaculture strays they remove (`Analysing_BLAST_output_WP2_12S.R`, ~lines 355–368). *(Also applied the 3 relabel/lump fixes — D3.)*
- [x] **Fit `msocc`** (fast Gibbs) to 6 **well-detected analysis labels** (≥30 units); interpret occupancy ψ vs detection θ·p; sanity-checked the detection histories. **Covariates complete across all 14 sites:** ψ ~ `Dist_from_lake` + `season`, detection ~ read depth. Flow was dropped (no spatial signal — D2); `river_gradient_percent` was deferred historically and ultimately excluded because the source E14 value remained absent.
- **Checkpoint:** ✅ `R/05_occupancy_msocc.R` (+ `R/test_msocc_mapping.R` proving the mapping) → `report/conwy-reanalysis.qmd` (renders to `_site/`), `data/metadata/phase1_occupancy_estimates.tsv`, `figures/phase1_psi_vs_p.png`. **Honest result:** for these six well-detected analysis labels ψ ≈ naïve occurrence and detection is high (θ·p 0.79–1.00), so imperfect-detection corrections are small *for this selected label set*. Sparse-label adjustment and design optimisation were deferred to the Phase-3 community stress test, then ultimately withheld after the recovery gate failed. *(Learning + de-risking artifact, not a headline result.)*
- **Reading:** [msocc paper & vignette](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13442) · MacKenzie et al., *Occupancy Estimation and Modeling* (ch. 1–4) · **McElreath, *Statistical Rethinking*** (the gentlest serious Bayesian on-ramp — ch. 1–4, 8).

### Phase 2 — Raw-reads pipeline → validated ASVs · ✅ **DONE & validated (2026-06-28)** · *PC (Windows/WSL2, 64 GB) · ran in parallel with Phase 1*
**Goal:** your own **ASV count table** (taxonomy-agnostic), validated against the authors', so B, D, Layer E, and `occumb`'s count input are unlocked.
> **Outcome (2026-06-28, verified):** `seqtab_clean.rds` = **647 field libraries × 1,348 ASVs, 35.7 M reads**. The broader **1,443-ASV field+control classification universe** contains 37 assigned taxonomy labels (886/1,443 ASVs assigned; 61% by ASV / 93% read-weighted); the field collapse has 35 columns (34 assigned labels + `Unassigned`). The validation gate is **PASS-WITH-CAVEATS** (all 10 source-study taxa detected in ≥3 shared libraries were recovered, mean Bray–Curtis similarity 0.84, precision 0.99). The one finding carried into Phase 3 is the **~16% detection-cell recall deficit and resulting occupancy-recoverability risk** (Phase-2 reality callout #3, §2.2 safeguard); no ψ-bias direction is established. Execution deltas vs the checklist below are catalogued in **Deviations from Phase 2** at the end of this doc.

> ### Phase-1 → Phase-2 adaptations (folded in 2026-06-25, post-grilling)
> Five changes from the deployed state of Phase 1 + a verification of the archive (the deviation catalogue at the end of this doc is the lens):
> 1. **The archive PRJEB48362 is 12S-only** — the `E##_R#_T##` filter returns *exactly* 755 (one marker's worth, not 3×), so COI/18S live in separate studies and **every control in this archive is a 12S control.** No marker-filtering of controls is needed.
> 2. **The download script excludes all controls** (`download_ena.sh` hard-filters to `E##_R#_T##`), so `decontam` and the tag-jump step have nothing to run on. **Fetching + identifying the controls is now Phase 2's first task** (see the new "Controls first" step).
> 3. **Hand-off is ASV-level, not species-level.** Layer B re-assigns taxonomy downstream, so collapsing ASVs→species in Phase 2 would bake in a DB that B immediately undoes. `seqtab_clean.rds` = **ASV×sample counts**; the species collapse moves to Phase 3, once per database. *(Corrects the `04_`/`06_` stubs' "ASV→species" wording.)*
> 4. **The validation gate targets the 68 *published* species by concordance**, then applies the **D3** cleaning (drops + 3 relabels) to reach the same **55** — keeping pipeline reproduction separable from the authors' taxonomic cleaning. Exact ASV/species-count match is *not* the gate (DADA2 1.14.1 → env 1.34.0 drift guarantees a mismatch).
> 5. **decontam thresholds are recorded pre-run in `docs/prereg.md` §3** (prevalence P\*=0.1 reported + 0.5 sensitivity; tag-jump cutoff calibrated to the PC spike, swept 0 / 0.1% / 1%). Unlike the Phase-1/Layer-E choices, this record appears in repository history before the Phase-2 run; none of it was externally registered.

- [x] **Run on the PC, not the Mac.** This is the only RAM/disk-heavy phase. Set up WSL2 + the lean `edna-phase2` conda env and download the 71 GB to the **PC's disk** (`export EDNA_RAW=/mnt/d/…`). Everything Phase-2 lives in the USB kit **`PS/eDNA/`** (`download_ena.sh`, `environment.phase2.yml`, `README_PC.md`); the `01_`–`03_` R scripts are authored on the Mac when this phase starts (DADA2 trim lengths need the real QC curves) and copied to the USB then. Raw FASTQs stay on the PC disk — never on the USB, never in git.
- [x] **Controls first — fetch + identify the Conwy 12S controls** *(new; `decontam` is blocked without them).* Extend `download_ena.sh` with a **second filter** for the control families, and add `data/verify_controls.py` (a Phase-0-style alias audit) to pin exactly which libraries are Conwy 12S **negatives** — `ENC_R#_T##` (45) + `EBlank_R#_T##` (9) + the shared WP2 lab/PCR blanks (`BLANK_WP2C_plate*`, `BLANK-<index>`) — vs **positives** — the `PC_*` mackerel-spike libraries — while **excluding the other rivers' controls** (`TWNC/SKNC/GWNC/GLNC`, `Blank_*_Season`, anything `C##_*`). Reconcile the downloaded set against the paper's stated **66 field + 68 lab** controls and document any difference (it may pool markers/rivers).
- [x] **QC:** `FastQC` + `MultiQC` (read the per-base quality curves; high duplication is normal for amplicons).
- [x] **Primer removal:** `cutadapt` with **MiFish-U** primers (F `GTCGGTAAAACTCGTGCCAGC`, R `CATAGTGGGGTATCTAATCCCAGTTTG`) as **linked adapters**, `--discard-untrimmed`, IUPAC degeneracy passed directly.
- [x] **DADA2:** `filterAndTrim` (set `truncLen`/`maxEE` from the QC curves; keep ≥12 bp merge overlap), `learnErrors` **per sequencing run**, `dada(pool="pseudo")`, `mergePairs`, `removeBimeraDenovo(method="consensus")`. (Authors used DADA2 v1.14.1 — note version drift to env 1.34.0.) **Process the control libraries through the same pipeline** so they share the ASV space `decontam` needs.
- [x] **Validation gate (concordance, not exact match):** use the current MitoHelper/MitoFish approximation to the authors' classification workflow, collapse ASVs to labels, and compare species×sample **detections** with the published 68-label table — report **concordance** (% recovered, read-weighted agreement, and *where* they diverge). The unarchived 23-Sep-2020 MitoFish build prevents a faithful historical replay, and DADA2/database-version drift makes the ~2726-ASV / 68-label / >20%-length-filter numbers loose **sanity targets, not a pass/fail.** Apply the documented D3 drops and relabels to the comparison table; do not present current labels as a reproduction of the 2020 database.
- [x] **Decontamination (rules locked in `docs/prereg.md` §3):** `decontam` **prevalence** method (no DNA-concentration data → the frequency method is N/A) using the negative controls, reported at threshold **0.1** with **0.5** as a stated **sensitivity** check; an **index-hopping / tag-jump** abundance cutoff **calibrated to the `PC` mackerel-spike bleed** into real samples, **swept over 0 / 0.1% / 1%** relative abundance with the detection changes reported; `metabaR` optional for tag-jump modelling.
- [x] **Hand back to the Mac (via USB):** copy the small derived outputs only — `seqtab_clean.rds` (the decontaminated **ASV×sample count** table), `repseqs.fasta`, the MultiQC report, and the **controls + decontam audit** — into the Mac repo's git-ignored `data/derived/`. The ASV→species collapse happens in Phase 3, per database; Phase 3 consumes these.
- [x] **Close with a 3-lens adversarial review** (statistical · code · fidelity/honesty), exactly as Phase 1 did (deviation **D8**) — record any corrections before the checkpoint. *(Delivered Mac-side as the validation-gate review — four lenses: recompute · crosswalk-audit · detection-offset · honesty — in `docs/phase2-validation-gate.md` + the six-audit synthesis in `docs/phase3-readiness.md`, rather than as a PC-handback artefact. See Deviation P2-D5.)*
- **Checkpoint:** ✅ validated, decontaminated **ASV count table** (`seqtab_clean.rds` 647 × 1,348, 35.7 M reads) + `repseqs.fasta` + `taxonomy.tsv` (MitoHelper 12S Nov-2025) + species/genus collapses + a documented **concordance** reproduction (`R/08_validate_concordance.R` → `docs/phase2-validation-gate.md`, PASS-WITH-CAVEATS) + control/decontam audits — all on the Mac. The former 79-vs-77 control discrepancy is reconciled: 79 negatives are in the manifest, while 77 enter decontam because `BLANK_D1_A1` and `BLANK_C3_A1` had no ASV-table row.
- **Reading:** [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial.html) · [cutadapt linked-adapters](https://cutadapt.readthedocs.io/) · Callahan et al. 2017, *Exact sequence variants (ASVs)* · [decontam vignette](https://benjjneb.github.io/decontam/) · **MIEM** reporting guidelines.

### Phase 3 — Baseline (E), occupancy audit (A), database concordance (B) · ✅ **LOCALLY COMPLETE WITH A CLAIM DEMOTION**

**Goal:** finish every defensible local analysis on the Phase-2 read counts and prevent a failed recovery gate from becoming a headline.

- [x] **E14 gate resolved for the executed formula:** source climate/gradient gaps remain explicit; a flagged, unvalidated E13 donor exists only in supporting metadata; climate, flow, and gradient are excluded. `Dist_from_lake + season` remains the occupancy specification.
- [x] **Layer E corrected:** 647 libraries → 234 site×time units. Spatial inference uses 14 sites × 11 common timepoints with 27 cyclic/mirrored river-series permutations; seasonal inference uses 18 timepoints × a fixed 10-site panel with 35 chronological-series permutations. Section PERMANOVA p=0.0714 for every distance; season p=0.333–0.778. The two unadjusted restricted-series section PERMDISP sensitivities for Bray–Curtis and rCLR each give p=0.0357; this is an exploratory dispersion signal, not a confirmed centroid shift. No restricted-series DA or indicator result passes 0.05/FDR.
- [x] **Layer A executed and bounded:** structural-NA design and smoke-scale `occumb` run end-to-end; the known-truth recovery MCMC fails its convergence gate, and 18 observed libraries that are zero within the 14-label harmonised Layer-A input expose an observation-mask/comparator mismatch. They contain reads under raw taxonomy labels and are not zero across the raw-name fish table. Sparse-label magnitudes, directions, coefficients, and utility curves are therefore withheld/diagnostic only. A separate deterministic lab demonstrates a generic over-correction counterexample. Production MCMC is not required for the website release and would not by itself repair structural non-identifiability.
- [x] **Layer B completed:** current-database concordance is reported honestly, not as a faithful 2020→2026 replay. Genus *Salmo* is stable across the ASV0001 comparisons; the species label is not.
- [x] **Control provenance reconciled:** 79 negative controls are listed in the manifest; 77 enter decontam because two well blanks have no ASV-table row.
- [x] **Optional post-release research:** the bounded, ordinary-laptop simulation, marker-resolution, robustness, trait and pipeline-disagreement work is completed in Phase 5 below. Production MCMC, dynamic occupancy, river-transport modelling, full FASTQ clean-clone reconstruction, and large reference-database rebuilds remain optional and do not block either the completed website release or the Phase-5 low-hardware scope.

**Checkpoint:** Layer E and B are locally complete; Layer A has a scientifically useful negative result and an explicit claim gate rather than an unconverged final estimate.

### Phase 4 — Local synthesis and website shipping · ✅ **LOCALLY COMPLETE**

**Goal:** one coherent local release that an unfamiliar reader can understand in five minutes.

- [x] **Technical synthesis:** report and manuscript evidence draft state the same Layer E, Layer A, and Layer B boundaries; neither claims publication or complete raw-read clean-clone reproduction.
- [x] **Local explorer:** Shiny is retained as an exploratory local app with schema validation and diagnostic Layer-A banners; no deployment is implied.
- [x] **Portfolio edition and project website:** static narrative landing page, explained catalogue linking every release deliverable, self-contained A4 research brief, exactly four-page visual PDF with five model architecture diagrams, standalone synthetic simulator, one-page poster, three hero figures, recording storyboard, and ten interview answers.
- [x] **Portable website handoff:** `python3 scripts/package_website.py` creates one sendable ZIP containing every local deliverable directly linked by the website plus its static dependencies, adds an obvious offline entry point and checksums, then verifies a clean extraction and all local HTML links.
- [x] **Laptop capsule:** `START_HERE.md` plus `python3 scripts/build_website.py` rebuild the website figures and PDFs from small included local inputs in under five minutes without MCMC, FASTQs, a server, or specialist hardware.
- [x] **Repo polish:** accurate README, license-status notice, citation metadata, scripted figures, and local links.

**Checkpoint:** release-ready local edition. A public repo, DOI, live URL, or preprint is optional and requires explicit owner action.

### Phase 5 — Post-release scientific-impact sprint · ✅ **LOCALLY COMPLETE · 100%**

**Goal:** turn the project's existing limitations into generalisable, testable results using existing, public-cached or synthetic data on an ordinary laptop. This phase adds scientific evidence; it must not weaken the claim boundaries that made the Phase 0–4 release defensible.

**Priority rule:** complete the recoverability atlas, multi-marker benchmark and Layer-E robustness multiverse first. Trait, pipeline-disagreement and interactive synthesis work follow. A failed predeclared gate is a valid negative result only when the affected numerical claim is withheld.

#### P5.0 — Audit and preregister the extension *(P0)*

- [x] Re-audit the real repository state and read the governing claim, inference, taxonomy and trait evidence before changing analysis code.
- [x] Add a concise local preregistration defining each question, estimand, biological unit, simulation grid, robustness specification, exclusion rule, randomisation unit and pass/demote gate (`docs/impact-sprint-preregistration.md` plus its machine-readable manifest and SHA-256 freeze receipt).
- [x] Record the timing limitation honestly: this is a repository-local analysis specification, not an externally registered protocol unless the owner later performs that action.
- [x] Establish one deterministic Phase-5 entry command, expected ordinary-laptop runtime, fixed seeds, small inputs and output manifest.

**Gate:** no Phase-5 result is promoted until its decision rule exists in the local specification and its denominator/unit has an automated check.

#### P5.1 — Occupancy recoverability atlas *(P1; best single addition)*

- [x] Expand the fixed-seed teaching scenario into a factorial known-truth experiment varying true occupancy, per-replicate detection, taxon prevalence, one-to-three replicates, sequencing-depth variation, depth–river-gradient coupling, structured missingness and observation-mask error.
- [x] Compare a small set of defensible fast estimators, such as naïve occurrence, replication-aware correction and a depth-aware model; do not require production MCMC.
- [x] Report bias, uncertainty calibration or interval coverage where available, failure rate and claim-gate status in a machine-readable result grid.
- [x] Produce a high-information **recoverable / unreliable / non-identifiable** phase diagram and a concise methods/limitations note.
- [x] Prove in tests and captions that these known-truth results are independent of the unconverged real Conwy Layer-A fit.

**Gate:** an estimator/scenario may be called recoverable only under a predeclared error/coverage rule. Failed fits and non-identifiable cells remain visible and cannot be averaged away.

#### P5.2 — In-silico multi-marker *Salmo* resolution benchmark *(P1)*

- [x] Assemble an accession-backed, locally cached sequence and provenance manifest for relevant *Salmo* labels using bundled or small public inputs.
- [x] Where sequence and primer coverage are adequate, compare the MiFish 12S fragment with suitable COI and cyt-b regions.
- [x] Measure exact-haplotype sharing, within- and between-label distances, barcode-gap behaviour and leave-one-accession-out nearest-neighbour classification.
- [x] Determine the finest taxonomic rank each marker consistently supports and visualise the resolution boundary.
- [x] Record missing taxa, uncertain labels, accession conflicts, retrieval dates and primer/amplicon definitions; narrow the benchmark rather than filling gaps with invented data.

**Gate:** this benchmark may establish sequence discrimination within the audited public records only. It cannot establish local occurrence, stock identity, database-label truth or ASV0001's species identity.

#### P5.3 — Layer-E robustness multiverse *(P1)*

- [x] Predeclare a finite specification set covering balanced versus all usable timepoints, site×time versus site-profile summaries, defensible taxonomic aggregation, reasonable detection thresholds, effort equalisation, supported distance measures and leave-one-site/timepoint-out checks.
- [x] Preserve the correct biological unit and repeated structure in every inferential specification; library-level tests remain descriptive only.
- [x] Store every specification, denominator, effect size, direction, uncertainty measure, permutation resolution and gate status in one result table.
- [x] Produce a specification curve or robustness matrix that shows the full result family rather than selecting favourable analyses.

**Gate:** the downstream pattern may be called robust only if it passes a predeclared stability rule. No collection of exploratory specifications may be used to override the restricted-series p=0.0714 claim boundary without a genuinely new, defensible estimand.

#### P5.4 — Functional-trait broadening null model *(P2)*

- [x] Use the existing cited trait table to construct a documented mixed-trait distance matrix, such as Gower distance.
- [x] Select and justify an occurrence-based conditional functional-diversity statistic for the balanced site×time panel.
- [x] Build a null model that preserves the relevant detected-label richness and prevalence structure.
- [x] Test whether downstream functional broadening exceeds the expectation created merely by detecting more analysis labels.
- [x] Repeat the result after excluding uncertain taxonomic or trait proxies and publish the sensitivity beside the primary result.

**Gate:** conclusions remain conditional on analysis-label assignments and cited trait proxies. They cannot claim abundance, diet, realised habitat use, causal environmental filtering or independently verified species occurrence.

#### P5.5 — Cluster-aware pipeline-disagreement analysis *(P2)*

- [x] Analyse source-positive/reprocessed-negative disagreements over the verified shared-library comparison using already justified predictors such as source relative abundance, reprocessed depth, label prevalence, site and time.
- [x] Respect dependence among taxon, library, site and time observations using a defensible clustered, hierarchical or block-bootstrap design.
- [x] Add leave-one-taxon-out and influential-group sensitivity checks where feasible.
- [x] Produce a machine-readable coefficient/effect table, calibration or diagnostic evidence, and a figure explaining where the pipelines disagree.

**Gate:** the outcome is pipeline disagreement or non-recovery relative to the source table. Neither pipeline is biological ground truth, so the result cannot be called a biological false-negative rate.

#### P5.6 — Claim-survival map and reader-facing synthesis *(P3)*

- [x] Build a lightweight offline view tracing **initial claim → credibility audit → evidence → limitation → surviving wording** for river community, detection/recoverability, pipeline disagreement, molecular resolution and functional traits.
- [x] Link every surviving statement to its machine-readable evidence, method and limitation.
- [x] Update the technical report, manuscript evidence draft, README, research brief, website, catalogue, interview set and recording storyboard only after the corresponding new result passes its gate.
- [x] Keep the page and completed release usable without a server, persistent network connection or specialist software.

#### P5.7 — Reproducibility, verification and packaging *(P0 throughout; final gate)*

- [x] Add unit, invariant and regression tests for seeds, masks, aggregation units, denominators, randomisation rules, accessions, provenance and claim gates.
- [x] Make every public input small, cached and provenance-traceable so the completed analysis works offline after acquisition.
- [x] Run the Phase-5 command twice and verify deterministic outputs where the method is designed to be deterministic.
- [x] Run the existing website build and full release verifier; extend the verifier for every new required artefact and scientific boundary.
- [x] Render and inspect every affected report, page, figure and PDF for overflow, legibility, misleading scales and unsupported captions.
- [x] Rebuild the portable ZIP, extract it into a clean temporary directory whose path contains spaces, and pass the bundled integrity/link verifier with zero failures.

#### Phase-5 completion condition

Mark Phase 5 complete only when:

- [x] all six impact additions have either passed their predeclared gate and produced a deterministic script, machine-readable result, figure, methods record, limitation statement and verified integration, **or** failed a predeclared gate with the numerical claim withheld and the required next evidence documented;
- [x] the three P1 additions—recoverability atlas, multi-marker benchmark and Layer-E multiverse—are implemented rather than merely scoped;
- [x] no headline depends on pseudoreplication, an unconverged model, an unresolved mask, unsupported taxonomy, unrecorded provenance or selective specification reporting;
- [x] the report, manuscript draft, README, plan, research brief, website and machine-readable status agree;
- [x] the full release verifier and clean-extraction bundle verifier pass with zero failures;
- [x] the Phase-5 workflow runs on an ordinary laptop without fieldwork, wet-lab work, specialist hardware, production MCMC, a raw FASTQ rerun, a massive database build or a persistent network service; and
- [x] an unfamiliar biosciences reader can understand the expanded question, at least three new defensible findings, their limitations and the scientific reasoning within five minutes.

#### Phase-5 non-blocker rule

Difficulty, local runtime, insufficient public sequence coverage, an inconclusive analysis or a failed scientific gate are not blockers: reduce the design defensibly, narrow the claim or record the failed gate.

**Checkpoint:** a verified, offline-capable scientific-impact extension in which the new evidence is more impressive because its failure conditions and inferential limits are visible.

---

## 4. Reproducibility & repo structure
```
conwy-edna-reanalysis/
├── README.md                  # what, why, how-to-run, links to preprint + tool
├── environment.lock.yml       # pinned conda env — MAC (arm64): R + JAGS + CmdStan; 21/21 required packages load
├── environment.R-packages.lock.tsv  # human-readable R-package version manifest
├── data/
│   ├── download_ena.sh        # pulls E## 12S runs (+ Phase-2 controls filter) from PRJEB48362 → $EDNA_RAW (PC disk); raw NOT committed
│   ├── derived/               # .rds checkpoints (DADA2 seqtab, MCMC) — git-ignored, regenerable (P6)
│   ├── metadata/              # covariates + occasions.tsv + data dictionary + verification reports
│   └── README.md              # exact accessions, versions, access dates
├── R/
│   ├── 00_setup.R  lib_design.R (canonical occasion/NA map)  01_qc_cutadapt.R  02_dada2.R  03_decontam.R
│   ├── 04_baseline_E.R (Phase 3, on counts)  05_occupancy_msocc.R (Phase 1)  06_occupancy_occumb.R  07_reclass_B.R
├── stan/                      # reserved for optional dynamic occupancy; no model implemented
├── report/                    # Quarto technical report: renders from verified local checkpoints; not submission-ready
├── shiny/                     # occupancy explorer app
└── figures/

  (off-repo) USB PS/eDNA/      # Phase-2 field kit for the PC: download_ena.sh + environment.phase2.yml
                               #   + ena/PRJEB48362_read_run.tsv + README_PC.md; also carries outputs back
```
- **Two execution environments were used:** `environment.lock.yml` records the Mac analysis backbone. The x86_64 WSL2 Phase-2 environment and full raw-read kit remain on the off-repo USB hand-back and are not claimed to rebuild from this repository. A local Quarto render reads the current derived checkpoints; the separate quick-start capsule rebuilds only the small presentation layer from included inputs.
- **Report DB name + exact version + access date** for every assignment; report assignment rate (% ASVs, % reads) per rank.

---

## 5. Rigor & honest-reporting checklist (build in from day 1)
- [x] Reads ≠ abundance → occurrence/detection framing.
- [x] Compositionality handled with relative-count, presence/absence, and robust-CLR distances.
- [x] Contamination + tag-jump checked; the manifest/decontam control denominators are reconciled and thresholds stated.
- [x] Dispersion checked (`betadisper`) alongside restricted PERMANOVA.
- [x] DA uses multi-method consensus plus a restricted-series gate, not a single liberal method.
- [x] **River transport:** detection location ≠ organism location; the authors' ≤5 km experiment is retained as a bound.
- [x] Reference-resource incompleteness and the missing MIDORI top-hit accession are reported.
- [x] Dedicated limitations state what eDNA cannot say and withhold failed Layer-A claims.

---

## 6. Analysis-decision history (with timestamp limitation)
The project intended to commit covariate sets, species inclusion, decontamination thresholds, and the DA rule before analysis. Git history co-commits the first Phase-1 record with its first fit and the Layer-E/DA additions with their outputs, so those sections are **not independently timestamped preregistration**. Phase-2 thresholds were present before the Phase-2 run. `docs/prereg.md` is retained as a transparent decision/amendment record, not described as externally registered.

---

## 7. Learning curriculum (sequenced — the steepest part is Bayesian/occupancy)
1. **R for a Python user:** *R for Data Science* (2e) — skim; focus on tidyverse + `ggplot2`.
2. **Reproducibility:** conda env-lock + R-package manifest (backbone), Quarto, Happy Git with R; the two-machine Mac↔USB↔PC workflow (`environment.phase2.yml` on the PC). *(`renv` optional — skim only if a cross-platform R-package lock is ever wanted.)*
3. **Amplicon bioinformatics:** DADA2 tutorial; Callahan 2017 (ASVs).
4. **Community ecology stats:** phyloseq + vegan tutorials; `microeco`; Gloor 2017 (compositional); iNEXT.
5. **Bayesian thinking (steepest):** McElreath, *Statistical Rethinking* (ch. 1–4, 8) — the on-ramp; budget 1–2 weeks.
6. **Occupancy modeling:** MacKenzie et al. (concepts) → Kéry & Royle AHM vol. 1 → `msocc` → `occumb`; foundational eDNA-occupancy papers: Schmidt 2013, Dorazio & Erickson 2018, Stratton 2020, **Fukaya et al. 2022**.
7. **eDNA rigor:** MIEM guidelines; `decontam`/`metabaR`; Bayer 2025 benchmark; Skelton 2023 (reads≠abundance).

*Claude Code removes the mechanical friction (scaffolding code, explaining errors); these readings build the conceptual judgment the code can't supply.*

---

## 8. Risk register (updated with verified findings)
| Risk | Likelihood | Mitigation |
|------|:----------:|------------|
| Pull the **wrong river** (`C##` vs `E##`) | ~~Med~~ **Closed** | Phase 0 verified strict 1:1, 755 E## runs; hard-filter in `download_ena.sh` + `verify_ena_join.py` |
| `cmdstanr`/Stan won't build on Apple Silicon | ~~Med~~ **Closed** | Phase-0 smoke test passes (receipt logged); JAGS/`occumb` chain also proven |
| **Published table is normalised, not counts** | **Closed** (designed-around) | Layer E + `occumb` moved to Phase-3 counts; published table = detection-only (msocc + validation) |
| DADA2 / 71 GB too big for the 8 GB Mac | ~~Med~~ **Mitigated** | Two-machine split: 64 GB PC runs Phase 2 via the USB kit; Posit Cloud no longer needed |
| **Two-machine sync friction** (env drift, lost outputs) | Low–Med | Separate per-platform envs by design; git repo single-homed on the Mac; USB carries only small derived `.rds` back; raw stays on PC |
| **E14 source covariates unavailable** | **Closed for executed formula** | A flagged E13 donor fills only a supporting climate table and is unvalidated; climate, flow, and `river_gradient_percent` are excluded from the headline formula. `Dist_from_lake + season` remains the executed specification |
| Community occupancy MCMC slow or unstable | **Managed by claim gate** | Smoke-scale work ran locally; the rare-species recovery fit failed convergence, so its numerical result is withheld. Production MCMC is optional future research, not an website-release dependency |
| ASV reproduction doesn't match authors' species table | ~~Med~~ **Closed** | Validation gate ran (`R/08`): **PASS-WITH-CAVEATS** — all ≥3-rep taxa recovered, mean Bray–Curtis similarity 0.84, precision 0.99; the one-library *Leuciscus* alternative remains unresolved. Concordance, not exact match, as designed |
| **Tag-jump miscalibration corrupts the count table** | ~~High~~ **Closed** | A cutoff = 1.0 bug zeroed 354 field libs (0.35 % reads) in the first handback; **regenerated on the PC with a 0.001 floor** → healthy 647 × 1,348 / 35.7 M-read table verified on the Mac |
| **Structured detection-cell recall deficit can confound ψ** | **High → bounded** | Across recovered taxa and 642 shared libraries, the pooled reprocessed/source-study detection-cell ratio is ~84% (selected label differences 22–29%) and varies with depth/relative-abundance signals. Structural NA and a depth term were tested, but recovery failed convergence and a comparison-mask mismatch remained for 18 libraries that are zero only within the 14-label harmonised input. Sparse-label magnitude and direction are therefore demoted rather than “corrected” |
| **Derived count tables live only on the Mac working tree + USB (data-loss)** | **Med** | Code/docs/plan committed 2026-07-05 (`08eae8e`), but `data/derived/*.rds` (incl. `seqtab_clean.rds`) is git-ignored and regenerable only by the 71 GB PC rerun → keep the USB copy as an ongoing off-Mac backup. The decision-record timing limitation is documented separately |
| **DB vintage can manufacture apparent Layer-B novelty; biological change is not established** | Med | Sep-2020 MitoFish unarchived → can't faithfully reproduce; frame as current-resource concordance and report genus stability/species uncertainty rather than temporal novelty (Gate 4) |
| **Controls excluded by the `E##` download filter** | ~~Med~~ **Closed** | Phase 2 fetched and processed the Conwy 12S controls; `decontam` (prevalence P\*=0.1 → 13 contaminant ASVs) ran on them. The denominator is reconciled: 79 negatives in the manifest, 77 in decontam because two well blanks have no ASV-table row |
| 5.9% Unassigned bin | **Characterised** | The audited bin is dominated by non-target vertebrate bycatch; it is reported with that scope. A curated regional database remains optional future work |
| Bayesian/occupancy conceptually hard | Med | Staged msocc→occumb; §7 curriculum; budget the time — you have slack |
| Unsupported uniqueness claim | **Closed by scope** | The source analysis did not report occupancy modelling; this project makes no claim that no other reanalysis exists or that its contribution is unique |
| Scope creep into stretch angles | Low–Med | D/I/dynamic are *labeled* stretch with explicit gates after the core lands |

---

## 9. Future-expansion notes (preserved, out of core scope)
- **Community contribution:** optionally build a **CRABS-curated UK fish 12S reference DB** and release it only after ownership/licensing review. Any future GBIF/OBIS submission must use curated detections and may be called occupancy-corrected only if a later model passes convergence and recovery gates.
- **Multi-marker extension:** Conwy also has **COI (invertebrates) + 18S** — a cross-marker community reanalysis.
- **Temporal extension:** richer **dynamic** (colonization/extinction) and **river-network transport (`eDITH`)** models using the ≤5 km empirical prior.
- **Method extension:** full pipeline-benchmarking (angle D) as a standalone "researcher degrees of freedom" study.

---

## 10. Verified-facts appendix (so nothing has to be re-derived)
- **Data:** ENA **PRJEB48362** (`ERP132733`). Conwy = sites **`E01`–`E14`**, timepoints **`T01`–`T19`**, replicates **`R1`–`R3`**; **755** deposited 12S replicate-libraries of 798 theoretical cells. The 43-cell difference is absent/unsequenced design space (including all T16), not a QC drop. **⚠️ `C##` = a different river** in the same archive.
- **Design:** 14 sites × 19 timepoints × 3 water-sample replicates; water-sample reps **sequenced separately** (2-level occupancy ✓); **PCR triplicates pooled** before indexing (no 3rd level).
- **Controls (Phase-2-verified):** the paper states **66 field + 68 lab** negatives across its reporting scope. The local 12S manifest identifies 79 negative controls; 77 enter decontam because `BLANK_D1_A1` and `BLANK_C3_A1` have no ASV-table row. Conwy-tagged families include `ENC_R#_T##`, `EBlank_R#_T##`, and shared WP2 lab/PCR blanks; `PC_*` are positive mackerel-spike controls. Other-river controls are excluded. See `data/metadata/control_reconciliation.md` for the exact denominator boundary.
- **Marker / primers:** MiFish-U 12S — F `5′-GTCGGTAAAACTCGTGCCAGC-3′`, R `5′-CATAGTGGGGTATCTAATCCCAGTTTG-3′`; amplicon ~172 bp.
- **Original pipeline:** 2-step PCR (Bista 2017), 384 dual indexes; demux by index then gene region via **Cutadapt 2.3**; **DADA2 v1.14.1**; ASVs >20% longer than 172 bp discarded; ~2726 12S ASVs pre-curation.
- **Original 12S taxonomy:** BLAST+ v2.10.0 — broad vs **NCBI nt** (≥70% id, ≥80% qcov, e≤1; keep Actinopteri), then specific vs **MitoFish (downloaded 23 Sep 2020)** (≥90% id, ≥90% qcov, e≤0.001), first hit.
- **Published tables (P1-verified):** 12S Source Data is **label-collapsed AND normalised relative abundance** — **665 Conwy replicate rows × 68 source-study species-label columns**, each row summing to 1 (*not* raw counts, *not* ASVs, and not independent species validation). Detection (>0 vs 0) is preserved, so it feeds msocc + the Phase-2 validation target, but **count-based methods (`occumb`, ALDEx2/ANCOM-BC2/corncob) require reprocessing**. **665 = 755 deposited − 90 QC-dropped** (authors drop libraries with <1000 reads); the 90 dropped (concentrated upstream E01–E04 + summer T07–T10) and T16 are **NA, not 0**. Supplementary Data 1 = per-sample covariates joined on `ID` (no ENA-accession column). E14 source climate/gradient values remain absent; the executed formula excludes climate, flow, and gradient. A flagged E13 donor is supporting-only and unvalidated.
- **Phase-2 outputs (2026-06-28, dims re-verified 2026-07-05):** reprocessing yielded **`seqtab_clean.rds` = 647 field libraries × 1,348 ASVs, 35.7 M reads** (regenerated after a tag-jump = 1.0 bug; 0.001 floor). The broader **1,443-ASV field+control classification universe** has 37 assigned taxonomy labels (886/1,443 ASVs assigned, 61% by ASV / 93% read-weighted); the field count collapse has 35 columns (34 assigned labels + `Unassigned`). The top 20 ASVs within the 5.9% Unassigned bin were non-target vertebrate labels; the smaller tail is unexamined. Taxonomy uses **MitoHelper 12S Nov-2025**. Validation (`R/08` + `crosswalk.tsv`) is **PASS-WITH-CAVEATS** over 642 shared libraries — mean Bray–Curtis similarity 0.84, cosine similarity 0.90, precision 0.99, all ≥3-library taxa recovered. **Governing caveat:** the pooled reprocessed/source-study detection-cell ratio across recovered taxa is ~84%, structured by depth and source-study relative-abundance signals. ASV0001 is stable at genus *Salmo* but not species; no local species identity is inferred from this fragment alone.
- **Transport:** authors' mackerel-release experiment bounds eDNA transport at **≤5 km**.
- **Tooling:** occupancy — `msocc` (Stratton 2020), `occumb` (Fukaya et al. 2022), `eDNAoccupancy` (Dorazio 2018). Classifier — MMseqs2@97% > BLAST for 12S (Bayer 2025).

---

*This plan is the committed answer to every open decision in `eDNA-reanalysis-project-scoping.md` §12. Stretch items are labeled and gated; the core (A + B + E → repo + preprint + tool) is the contract.*

---

## Deviations from Phase 1
*Catalogue of every change made while executing Phase 1 vs. the plan as written above. Appended 2026-06-24. Each entry: what changed · why · impact. Implemented in `R/05_occupancy_msocc.R`, `docs/prereg.md`, `report/conwy-reanalysis.qmd`.*

### D1 — `msocc` 2-level mapping: the 3 water samples enter at the **sample (θ)** level, not the **rep (p)** level *(method correction)*
- **Plan said:** "the 3 water-sample replicates … are the detection replicates … 2-level multiscale model (site→replicate)" (§2.2), implying the 3 reps sit at `msocc`'s replicate/PCR level.
- **Changed to:** site×time → `msocc` *site* (ψ); the 3 water samples → `msocc` *sample* (θ); a single dummy PCR per sample → *rep* (`rep ~ 1`, collapsed).
- **Why (objective, identifiability):** PCR triplicates were pooled before indexing, so each water sample yields exactly **one** call. Read directly from the `msocc_mod` source: occupancy `z` is "any positive across a site's samples," and a level with one unit confounds with its neighbour. Placing the 3 reps at the PCR level (1 sample/unit) makes only the product **ψ·θ** identifiable → occupancy ψ, the entire point of the project, becomes unrecoverable. Placing them at the sample level gives the canonical, identifiable single-season occupancy model `y_ij ~ Bern(z_i·θ_ij·p)`, `z_i ~ Bern(ψ_i)`: **ψ is clean**. With one call per sample, θ and p are not separately identified, so we report the identified product **θ·p as the per-replicate detection probability**. Verified empirically against the `msocc` package internals + `fung`/`goby` examples before committing.
- **Impact:** none downstream — Phase-3 `occumb` is purpose-built for this 2-level metabarcoding structure and is unaffected. Phase-1 reports ψ and θ·p.

### D2 — Flow **dropped** from the Phase-1 detection covariates *(degenerate covariate)*
- **Plan said:** Phase-1 covariates include "flow-with-`Catchment_flow` substitution" (§3 Phase 1; prereg §1 amendment).
- **Changed to:** Phase-1 detection model uses **read depth only**; flow deferred to Phase 3.
- **Why (objective):** the only flow column complete across all 14 sites, `Catchment_flow_m3_s`, is a **catchment-wide aggregate that is identical across all 14 sites at each timepoint** — it varies only over *time* (18 distinct values, T01–T19), carrying **no spatial signal** and being effectively collinear with timepoint (and hence with the `season` factor already on ψ). So it cannot inform the spatial occupancy story that is Phase 1's focus. The time-varying `*_flow_by_catchment` family that *would* add per-site/temporal signal is **100% NA for E13+E14** (both sites entirely). Net: no flow covariate that is both complete across 14 sites *and* spatially informative exists for Phase 1, and `Catchment_flow` adds only a temporal term redundant with the near-ceiling detection of these selected well-detected labels.
- **Impact (final disposition):** flow did not re-enter the executed Phase-3 headline formula. Source E14 climate/gradient values remain unavailable; a flagged E13 donor exists only in a supporting table and is unvalidated. The executed formula excludes flow, climate, and gradient and uses `Dist_from_lake + season`. Phase-1 detection covariate = per-sample **z-scored `log(read_count)`**.
- **Correction (2026-06-24, post-review):** an earlier draft of this entry (and `docs/prereg.md`) claimed `Catchment_flow_m3_s` was *"a single constant 2.53, zero variance."* That was **wrong** — an artifact of a probe that took the first value per site (all T01 = 2.53). The column does vary over time (variance ≈ 306). The *decision* to drop flow stands (no spatial signal, as above); only the stated reason was corrected.

### D3 — Applied the authors' **3 relabel/lump** operations, not only the drops *(fidelity to the authors' cleaned set)*
- **Plan said:** "drop the `Scomber_scombrus` positive-control spike + the marine/aquaculture strays" (§3 Phase 1) — drops only.
- **Changed to:** also apply the 3 DB-misassignment relabels from the **same script block** (lines 366–368): `Salmo_obtusirostris`→`Salmo_trutta`, `Anguilla_rostrata`→`Anguilla_anguilla`, `Salvelinus_fontinalis_x_Salvelinus_malma`→`Salvelinus_malma`, OR-combining detections (all 3 targets already exist as columns).
- **Why:** these are the source analysis's explicit taxonomic harmonisation operations; applying the drops without the relabels would leave the comparison namespace inconsistent and inflate the label count. They are bookkeeping for concordance, not independent proof of biological identity. **68 raw → 55 cleaned labels.**
- **Impact:** the Phase-1 analysis-label set matches the authors' cleaned 12S set; it is carried forward as the label universe for Phase 3, not as independent species validation.

### D4 — `Favonigobius_gymnauchen` drop is a **no-op** *(documented, not an error)*
- The authors' drop list includes `Favonigobius_gymnauchen`, but it is **absent** from our published 68-column matrix, so its removal does nothing. Recorded so the discrepancy in "11 listed vs 10 removed" is explained, not mistaken for a bug.

### D5 — NA-not-0 enforced by **row inclusion**, validated against `observable_mask` *(equivalent, simpler)*
- **Plan said:** "set `y[!observable_mask] <- NA`" — materialise a full site×tp×rep array and NA-fill the 90 QC-dropped + 43 unsurveyed cells.
- **Changed to:** include **only** the 665 observable cells as `msocc` rows; unobserved cells are simply absent → contribute no likelihood (`msocc` drops all-NA rows anyway). A `stopifnot(setequal(source rows, observable cells))` gate proves the 665 rows are *exactly* `observable_mask("published")`, so the NA rule is verified rather than hand-applied.
- **Why:** identical statistical result (a dropped/unsurveyed cell contributes nothing either way), fewer moving parts, and it keeps `lib_design.R` as the single source of truth via an explicit equality check.
- **Impact:** none — same likelihood; the QC-dropped/unsurveyed cells are never coded as false zeros.

### D6 — Implementation choices (minor) *(no statistical effect on reported quantities)*
- **`beta_bin = FALSE`** for every model so ψ/θ/p are extracted uniformly from the logit coefficients (β/α/δ × design matrices) across the null and covariate models — avoids `msocc`'s dual Beta-Binomial vs logit return paths.
- **Model progression per species** (M0 `ψ~1, det~1` → M1 `ψ~Dist+season` → M2 `+ det~depth`) with **WAIC**, slightly beyond the plan's "fit msocc to several species," to teach model comparison (the Phase-1 learning goal). `num.mcmc = 11000`, 1000 burn-in, `seed = 1`.
- **Analysis-label selection rule:** detected in **≥ 30 of 238 units** → 6 labels spanning the detection gradient (now recorded in `docs/prereg.md` §2).

### D7 — Analysis-decision record finalized; timestamp limitation documented
- `docs/prereg.md` records transformations, species inclusion, the Phase-1 covariate subset, later amendments, and the failed-claim-gate disposition. The first Phase-1 record and first fit share the 2026-06-24 commit, and Layer-E decisions share a commit with their outputs; repository history therefore does not independently establish those choices as pre-result. Phase-2 thresholds were recorded before that run. No external registration is claimed.

### D8 — Adversarial review pass + resulting corrections *(2026-06-24; quality gate)*
After the first working build, a 3-lens independent review (statistical correctness · code correctness · fidelity/honesty) was run against the deliverables. It surfaced — and these were then fixed — :
- **Corrected the D2 flow rationale (was factually wrong).** The original claim that `Catchment_flow_m3_s` is "a single constant 2.53, zero variance" came from a probe that took the first value per site (all T01 = 2.53). The column actually varies over time (18 values, ~2.5–58 m³/s). The *drop decision* stands (it is a catchment aggregate with **no spatial signal**), but the stated reason — in D2, `docs/prereg.md`, `data_dictionary.md`, and `covariate_join_verification.md` — was rewritten to be accurate.
- **Fixed a WAIC bug.** `msocc::waic()` returns a bare numeric, not a list; the original `…$waic` call errored and was silently `NA`-ed (so the committed estimates shipped an all-NA WAIC column). Corrected to `as.numeric(msocc::waic(mod, type = 2))`; the model-comparison table is now populated (and the report cites it via live inline-R, not literals, so it cannot drift).
- **Relabelled the false-negative quantity** `(1−θ·p)³` as a *design-effort* (nominal-3-replicate) power statement, and surfaced that **40 of 238 units have <3 observed replicates** (their realized miss rate is higher) — in both the script and the report.
- **Added a convergence diagnostic** (`posterior::ess_basic`): the null-model minimum ψ ESS is about 212 (the saturated *Salmo trutta* boundary case), while the exploratory occupancy-covariate models reach a minimum of about 117 and remain preliminary. The later Phase-3 recovery fit added multi-chain R̂ and failed its claim gate, so rare-species recovery numbers are not results.
- **Unit-level z-scoring** of the site-level covariate `Dist_from_lake` (was replicate-weighted across the 665 rows — an affine reparam with no effect on reported ψ/WAIC/detection, fixed for cleanliness).
- Reconciled the **252-theoretical vs 238-realized** unit denominators across the report and prereg.
- *Confirmed sound (no change needed):* the D1 sample-level mapping (re-derived as identifiable), the `observable_mask` equality gate, the drop/relabel cleaning vs the authors' script, the ID parsing, and idempotent re-runs.

---

## Deviations from Phase 2
*Catalogue of every material change made while executing Phase 2 vs. the plan as written above. Appended 2026-07-05. Each entry: what changed · why · impact. Sources: `docs/phase2-validation-gate.md`, `docs/phase3-readiness.md`, `R/08_validate_concordance.R`, `R/lib_design.R`.*

### P2-D1 — A **tag-jump = 1.0** calibration bug corrupted the first handback; fixed with a 0.001 floor *(critical, caught by audit)*
- **What happened:** the first `seqtab_clean.rds` was structurally degenerate — 375 × 829, 173 k reads, every ASV in exactly one sample — because the tag-jump cross-talk cutoff was applied at **1.0** (treat 100 % of reads as index-hopping), which zeroed 354 healthy field libraries (0.35 % of field reads retained). Its sibling `species_table.rds` (43 M reads) contradicted it 250×.
- **Fix:** regenerated on the PC with a sane **0.001** floor (sweep-validated: 0.001 → 77.5 % retained). Result is the healthy **647 × 1,348 / 35.7 M-read** table now on the Mac.
- **Impact:** Blocker 0 cleared. ⚠️ The `docs/phase3-readiness.md` §3 red-flag table and §4 reconciled-numbers still quote the **pre-fix** figures — treat them as superseded provenance, not current state.

### P2-D2 — Validation is **crosswalk-bridged concordance**, not a naive name-join or exact-count match *(method correction)*
- **Plan said:** "assign the ASVs with the authors' 2020-era approach, collapse to species, compare … report concordance."
- **Changed to:** the two database vintages give different nearest-neighbour labels to compatible short-marker evidence, so a raw name overlap (15 taxa) manufactures false discordance. A **consensus crosswalk** (`data/derived/crosswalk.tsv`, three independent audit passes plus reconciliation) maps 35 reprocessed raw columns to 14 source-study namespace targets, applies the D3 bookkeeping to both sides (68→55 source-study labels; 35→14 reprocessed labels), and compares **642 shared libraries** with read-conservation gates.
- **Impact:** recovery is measured against the **16 source-study taxa actually detected** in the shared libraries (the correct recall denominator), not 68 or 55. Recovery is **68.8 % (11/16)**; genus recovery 66.7 %; every ≥3-library taxon was recovered, while trace alternatives remain unresolved.

### P2-D3 — The load-bearing result is the **detection-offset (recall deficit)**, not the recovery percentage *(re-framing; reshapes Phase 3)*
- **Emergent finding (not in the plan):** across recovered taxa and 642 shared libraries, the reprocessed/source-study **pooled detection-cell ratio is ~84%** — a 16% deficit (selected label-specific differences reach 22–29%). Miss indicators are associated with library depth and source-study relative abundance in an unclustered descriptive logistic diagnostic. Precision is 0.99; this is detection agreement, not an occupancy estimate or a per-library rate.
- **Impact:** promoted to a **hard modelling constraint** on Layer A (§2.2 safeguard) and a primary risk-register row. The project tested whether depth-aware occupancy adjustment could recover known truth; because that recovery fit failed convergence and retained a comparison-mask problem, the outcome is to bound and demote rare-species claims rather than claim correction.

### P2-D4 — Phase-2 assignment labels were handed to Phase 3 for bounded audit *(scope handoff completed)*
- Phase 2 delivered the ASV table plus a 35-column raw collapse (34 assigned labels + `Unassigned`) and a 14-label harmonised collapse for the gate. Phase 3 did **not** force a `Salmo_macrostigma`→*trutta* replacement: current-resource comparisons support genus *Salmo* but not a local species identity. The 37-label broader taxonomy universe versus 55-label source-study universe is not a biological recovery statement. The validated recall result is the **~16% pooled detection-cell deficit across recovered taxa and 642 shared libraries** (denominator 16 source-study taxa); geography-based artifact/stray conclusions are not inferred from short labels alone.

### P2-D5 — The closing adversarial review was delivered **Mac-side on the validation**, not as a PC-handback artefact *(equivalent, relocated)*
- **Plan said:** close Phase 2 with a 3-lens (statistical/code/fidelity) review in the handback.
- **Changed to:** the review ran on the Mac against the concordance run — **four lenses** (recompute · crosswalk-audit · detection-offset · honesty) in `docs/phase2-validation-gate.md`, plus a six-deep-reader-audit synthesis in `docs/phase3-readiness.md`. Every headline reproduced under independent rewrite.
- **Impact:** the quality gate is met and lives in the Mac documentation rather than the USB handback. The remaining well-blank discrepancy is resolved: 79 negatives are present in the manifest and 77 enter decontam because two blanks have no ASV-table row.

### P2-D6 — `lib_design.R` gained a **`reprocessed`** survivor mask *(fulfils Phase-3 Gate 2 ahead of time)*
- Added `observable_mask("reprocessed")`, keyed off the 647 count-table rownames with an equality gate (survivors ⊆ sequenced design cells) and `R/test_reprocessed_design.R`, so Phase-3 detection histories code dropped / `E02_R3_T17` / NA-trim libraries as **structural NA, never observed-zero** — the coding the §2.2 safeguard depends on.

---

## 11. Portfolio enhancement — low-hardware, self-controlled sprint

**Added:** 2026-07-18

**Purpose:** turn the technically strong reanalysis into a concise scientific story that demonstrates biological insight, quantitative reasoning, reproducibility, and willingness to revise a conclusion when its assumptions fail.

**Scope constraint:** laptop-first; reuse existing and synthetic data; no new sampling, wet-lab work, massive database build, production-scale MCMC, hosted server, peer review, or publication dependency.

### 11.1 Reader-facing thesis

The portfolio should centre on one defensible claim:

> **More data do not automatically produce more biological certainty.** eDNA reprocessing can create structured non-detection; a short reference sequence can receive conflicting species-level labels; and a complex model output must be withheld when its convergence or recovery gate fails. A separate known-truth lab demonstrates why confounded adjustment can mislead.

The strongest version is a scientific detective story, not simply a larger pipeline. It should make four qualities visible to a biosciences reader:

1. **Biological interpretation:** connect molecular detections to fish ecology and the river system.
2. **Quantitative judgement:** distinguish read counts, detection probability, and occupancy.
3. **Scientific self-correction:** show where audits changed the analysis or weakened a claim.
4. **Independent execution:** provide a small, reproducible result that runs without specialist hardware.

### 11.2 Protect credibility before adding polish *(P0)*

- [x] **Correct Layer E's unit of inference:** 234 site×time units feed descriptive diversity; inference uses 14 site or 18 timepoint profiles with restricted series permutations.
- [x] **Demote the smoke-scale sparse-label occupancy result:** recovery MCMC convergence fails and 18 observed libraries that are zero only within the 14-label harmonised input expose a comparison-mask mismatch; magnitude and direction are withheld.
- [x] **Keep the portfolio edition independent of production `occumb`:** it uses Phase-1 detection estimates for six well-detected analysis labels, corrected Layer E, the database-resolution result, and a separate deterministic teaching lab.
- [x] **Reconcile the 79-vs-77 well-blank count:** 79 manifest negatives, 77 decontam inputs; two zero-read well blanks have no ASV-table row.

**Gate:** no website-facing claim may depend on unresolved pseudoreplication, a mismatched observation mask, or unconverged production MCMC.

### 11.3 Core presentation artefacts

#### A. Four-page visual investigation brief · ✅ **complete locally**

Create canonical `website/research_brief.html` in the paper-style format used by the companion website project, retain `visual-brief.html` as a compact compatibility view, and provide a matching four-page visual PDF containing:

- the biological question in plain language;
- five model architecture diagrams that make the multiverse, recoverability simulation, pipeline model, multi-marker benchmark, and trait null visually legible;
- three headline findings, each paired with its evidence and limitation;
- a short **“What surprised me / what changed my mind”** panel;
- the important mistakes or assumptions caught by audit: wrong-river prevention, normalized-vs-raw counts, tag-jump calibration failure, structured recall deficit, unstable *Salmo* species labels, and the risk of depth-based over-correction;
- a clear separation between **established**, **provisional**, and **not claimed**;
- the next experiment that would be run with additional resources.

#### B. Molecular detective vignette: the *Salmo* sequence · ✅ **complete locally**

Build a self-contained figure and explanation titled:

> **How can the same approximately 170-bp sequence receive three different *Salmo* species-level labels?**

Include:

- ASV0001 aligned against 5–10 relevant references, including *Salmo trutta*, *S. macrostigma*, *S. carpio*, and *S. salar*;
- annotated informative and non-informative bases;
- a small edit-distance network or tree;
- a comparison of sequence similarity with geographic plausibility; and
- the biological conclusion: genus-level evidence may be stable while species-level assignment is not identifiable from this marker/reference set.

The figure must avoid implying that a short tree proves species identity.

#### C. Trait-based ecological interpretation · ✅ **complete locally**

Create a cited table for the 14 detected fish covering, where evidence is available:

- migratory strategy;
- rheophilic/lentic association;
- typical river zone or habitat;
- trophic guild;
- native/introduced status; and
- conservation relevance.

Frame one narrow question before presenting the plot (the repository history does not establish this as independently pre-result):

> **Does upstream–downstream community turnover correspond to known habitat or life-history traits?**

Use site-aware incidence summaries and uncertainty appropriate to the sampling design. Treat the result as ecological interpretation, not proof of causation.

#### D. Small detection-and-occupancy laboratory · ✅ **complete locally**

Produce a standalone HTML page using synthetic data only. It should run locally with no Shiny server and expose a small set of controls:

- true occupancy;
- per-replicate detection probability;
- number of replicates;
- sequencing depth; and
- correlation between depth and occupancy/abundance.

It should show:

- naive detection prevalence versus true occupancy;
- probability of missing a present species;
- the benefit of replication versus simply increasing depth; and
- a counterexample in which depth correction over-corrects because depth is entangled with biology.

**Definition of done:** the page teaches why `reads ≠ abundance`, `non-detection ≠ absence`, and `statistical adjustment ≠ automatic truth` in under three minutes.

#### E. Five-minute reproducibility capsule · ✅ **complete locally**

Add one documented command that runs on bundled, small data in under five minutes and recreates three hero figures:

1. the river community gradient;
2. the annotated *Salmo* resolution result; and
3. the occupancy/detection simulation.

Add `START_HERE.md` with prerequisites, inputs, outputs, expected runtime, and the exact command. The capsule must not require private data, specialist hardware, a full FASTQ rerun, or a network service after setup.

#### F. Busy-reader presentation layer · ✅ **complete locally; recording requires owner voice**

Build a static landing page in this order:

1. a 30-second graphical abstract;
2. three conclusions;
3. three hero figures;
4. how the conclusions changed during the investigation;
5. the synthetic simulator;
6. methods and reproducibility; and
7. links to the report, code, and data provenance.

Also prepare:

- a one-page poster;
- a three-minute narrated screen recording; and
- ten short interview answers covering reads versus abundance, occupancy versus detection, conflicting databases, pseudoreplication, negative results, and an improved field design.

The static site can be built with [Quarto](https://quarto.org/docs/output-formats/html-publishing.html) and, if desired, served through [GitHub Pages](https://docs.github.com/en/pages); neither requires a permanent application server.

### 11.4 Optional instant citable release *(P3; not a completion dependency)*

- [ ] Create a versioned GitHub release after the reproducibility capsule passes.
- [ ] Optionally archive that release with Zenodo for a DOI using its [GitHub software archiving workflow](https://help.zenodo.org/docs/github/archive-software/github-upload/).

This is a stable, citable snapshot—not a claim of peer review or formal publication. The portfolio edition is complete without it.

### 11.5 Ordered execution and completion criteria

If time is limited, complete these five items in order:

1. **Correct Layer E replication and demote any unsupported rare-species claim.**
2. **Write the four-page visual investigation brief.**
3. **Complete the *Salmo* molecular detective vignette.**
4. **Build the small synthetic detection/occupancy laboratory.**
5. **Build the local static landing page with the three hero figures.**

The presentation sprint is complete when:

- every headline claim is linked to a figure/table and a stated limitation;
- established and provisional conclusions are visibly separated;
- an unfamiliar reader can understand the biological problem and three findings in five minutes;
- the small reproducibility command succeeds on a clean local setup;
- scientific self-corrections are explicit; and
- none of the artefacts depend on external review, a hosted compute service, new hardware, or unfinished production MCMC.

### 11.6 Explicit non-goals

- no new field sampling or wet-lab experiment;
- no dynamic occupancy model;
- no full FASTQ rerun for presentation purposes;
- no massive reference-database download or rebuild;
- no additional machine-learning model merely for novelty;
- no production MCMC requirement for the portfolio edition;
- no Shiny server or other persistent hosting;
- no full journal submission or publication dependency; and
- no visual polish before the scientific coherence and units of inference are correct.
