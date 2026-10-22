# Analysis-decision record — Conwy 12S occupancy & reclassification reanalysis

*This file was intended as a repository-timestamped preregistration. Git history,
however, first records the Phase-1 decisions in the same commit as the first Phase-1
fit and records the Layer-E/DA additions with their first outputs. It is therefore
**not** an externally registered preregistration and does not independently prove
those choices predated results. Phase-2 decontamination thresholds were present in
the repository before the Phase-2 run, but they were not externally registered.*

**Status:** Historical analysis-decision record plus a dated post-fit disposition in §7. The original Phase-3 set was not executable as written because E14 remained incomplete; its documented amendment and failed-claim-gate outcome are recorded below.

> Sections 1–6 preserve the decisions as originally recorded. They are not a current
> results summary. Where their planned language conflicts with the completed release,
> the §7 post-fit disposition and `docs/phase3-claim-status.md` control.

## 1. Occupancy covariates (historically recorded decision)

- **Occupancy ψ** ~ `Dist_from_lake` + `river_gradient_percent` + `season`  *(position/habitat drivers)*
- **Detection p** ~ `week_before_flow_by_catchment` + per-sample read depth + `temp_monthly_average`  *(eDNA-availability/sampling drivers)*
- **Transformations / standardization (fixed):** z-score continuous covariates (`Dist_from_lake`; read depth entered as z-scored `log(read_count)`); `season` as a factor with reference level **Spring**.
- Random effects / structure: 2-level occupancy (site×time unit; 3 water-sample replicate detections). **T16 = unsurveyed across all sites** (NA occasion, not non-detections). In `msocc`'s 3-level parameterisation the 3 water samples enter at the *sample* (θ) level, not the *rep* (p) level, so occupancy ψ stays identifiable (PCR reps were pooled → 1 call/sample → θ and p are not separately identified; the product θ·p is the reported per-replicate detection). See eDNAplan.md "Deviations from Phase 1."

> **Phase-1 covariate subset (recorded 2026-06-24; first commit coincides with the first fit):** Phase-1 `msocc` uses only covariates complete across all 14 sites, so the progression is ψ ~ `Dist_from_lake` + `season`, detection ~ z-scored `log(read_count)`. **`river_gradient_percent`** is held for Phase 3 (E14 missing). **Flow is dropped from Phase 1**: the only flow column complete across all 14 sites, `Catchment_flow_m3_s`, is a catchment-wide aggregate identical across all 14 sites at each timepoint (varies only over time, 18 distinct values) → **no spatial signal**, collinear with timepoint/`season` → uninformative for the spatial occupancy model; the spatially-resolved `*_flow_by_catchment` family is 100% NA for E13+E14. Re-entry with externally sourced E14 data was planned at the time but is superseded by §7: the executed formula excludes flow, climate, and gradient.

> **Covariate-availability note (dated 2026-06-19 in the project record; repository timing limitation above applies):**
> - **Temperature:** the originally-named `Water_temp_oC` is **100% NA** for Conwy in Supplementary Data 1. The record selects **`temp_monthly_average`** (the column the authors used; complete for 13/14 sites). Git history does not independently establish this note as pre-fit.
> - **Flow for E13/E14:** `week_before_flow_by_catchment` is NA for the two lowest sites (E13, E14); the catchment-wide **`Catchment_flow_m3_s`** (a per-timepoint series, ~2.5–58 m³/s; *not* a constant — an earlier note erroneously called it "= 2.53", the T01 value) is on hand for both and carried in `conwy_covariates.tsv` as the available substitute. *(Phase 1 does not use it — see the covariate-subset note above; it is a catchment aggregate with no spatial signal.)*
> - **Historical E14 decision (superseded by §7):** `river_gradient_percent`, `temp_monthly_average`, and `rainfall_monthly_average` are NA for E14 in the on-hand source table. A prior note incorrectly inferred a gridded-climate provenance; the source Methods instead describe the Glasgwm Automatic Weather Station. The final local analysis excludes climate, flow, and gradient from the headline formula; a separately flagged E13 donor exists only in a supporting table and is not treated as observed or validated E14 climate.

## 2. Species inclusion rule

- **Cleaning (fixed):** reproduce the source study's pre-analysis bookkeeping (`Analysing_BLAST_output_WP2_12S.R` ~355–368) — drop the `Scomber_scombrus` positive-control spike + 9 out-of-target labels present in our matrix (`Favonigobius_gymnauchen` is on their list but absent here → no-op), and apply its 3 harmonisation relabels (`Salmo_obtusirostris`→`Salmo_trutta`; `Anguilla_rostrata`→`Anguilla_anguilla`; `Salvelinus_fontinalis_x_Salvelinus_malma`→`Salvelinus_malma`) by OR-combining detections. These reproduce the source namespace; they are not independent identity validation. 68 raw labels → **55 cleaned labels**.
- **Phase-1 fitting set (recorded with the first fit):** fit single-label `msocc` to source-study analysis labels detected in **≥ 30 of the 238 occupancy units**. This yields 6 well-detected labels spanning the detection gradient (`Salmo_trutta`, `Phoxinus_ujmonensis`, `Salmo_salar`, `Barbatula_barbatula`, `Anguilla_anguilla`, `Gasterosteus_aculeatus`). Repository history does not independently prove the rule preceded inspection; the remaining sparse labels enter only the later claim-gated stress test.

## 3. Decontamination thresholds (Phase 2) — FIXED pre-run (2026-06-25)

- **decontam method:** **prevalence** (presence across negative controls vs samples). The *frequency* method is **N/A** — no per-sample DNA-concentration data exists for this dataset.
- **decontam threshold:** primary **P\* = 0.1** (package default) **reported**; **P\* = 0.5** (more aggressive) run as a stated **sensitivity** check.
- **Negative controls used:** Conwy 12S `ENC_R#_T##` + `EBlank_R#_T##` + shared WP2 lab/PCR blanks (`BLANK_WP2C_plate*`, `BLANK-<index>`); other rivers' controls excluded.
- **Index-hopping / tag-jump:** an abundance cutoff **calibrated to the `PC` mackerel-spike (`Scomber_scombrus`) bleed** into real samples, **swept over 0 / 0.1% / 1%** within-sample relative abundance, **reporting how detections change** across the sweep. `metabaR` optional for explicit tag-jump modelling.
- *Rationale + the controls-download gap (the `E##` filter excludes controls) are documented in eDNAplan.md Phase 2.*

## 4. Differential-abundance consensus rule (Layer E)

- Methods: ANCOM-BC2 + ALDEx2 + corncob (each an omnibus across-section test). **Consensus rule documented 2026-07-05 in the same commit as the first Layer-E outputs:** a taxon is called DA if significant at **FDR/q < 0.05 in ≥ 2 of the 3** methods, computed on the SAME prevalence-filtered set (species detected in ≥ 20 / 647 libraries → 8 testable; the 6 ultra-rare species are reported *not-testable*, not *not-DA*). This is not independently pre-result. The final release adds a restricted-series site-level gate, under which no DA claim passes.

## 5. Reference-DB reclassification (Layer B)

- **Baseline decision (Gate 4, recorded 2026-07-05, option iii):** the authors' Sep-2020 MitoFish build is **unarchived**, so Layer B is current-resource concordance, not a faithful 2020→2026 replay. No externally registered or independently pre-result timing is claimed.
- **Databases + exact versions + access dates:** *local baseline in hand* = MitoHelper **12S Nov-2025** first-hit assignment (`data/derived/taxonomy.tsv`; 886/1,443 ASVs assigned in the broader field+control classification universe). This is not the authors' unarchived 2020 output. *2nd DB* = MitoFish 12S via MitoHelper, **Nov-2025** build, 64,059 records, DOI 10.5281/zenodo.17602902, accessed **2026-07-05**.
- **Classifier + threshold:** MMseqs2 `easy-search` @ **97%** min-seq-id (`--cov-mode 2 -c 0.9 --search-type 3 -k 13`). Concordance metric: best-bit-score hit per ASV; species/genus agreement over both-assigned ASVs; counts of newly-resolved / lost / changed (`data/derived/reclassB_mmseqs_summary.tsv`).
- **Scope caveat + cross-reference comparison:** Part B uses the *same MitoFish vintage* twice, so it measures **classifier/method sensitivity (MMseqs2 vs first-hit BLAST)**, not database currency. Part C (completed 2026-07-05) adds a separate current reference resource — **MIDORI2 LONGEST, GenBank271 (2026-04-07), 12S srRNA; accessed 2026-07-05**: genus agreement is **99.4%** and species agreement **46.7%** (`reclassB_indepDB_summary.tsv`). ASVs labelled *Salmo_macrostigma* by MitoHelper receive a predominant *S. carpio* label under MIDORI2, while genus *Salmo* is retained. This is evidence of species-label instability, not a corrected local species or native-status identity.
- **Assignment-audit set (Gate 4):** `Salmo_macrostigma`, `Phoxinus`, `Barbatula`, `Cottus`, `Pungitius_kaibarae`, `Perca_flavescens`, … (`data/derived/reclassB_drift_catalogue.tsv`). The catalogue records `audit_category`, `audit_hypothesis`, `evidence_flag`, and `audit_note`; its alternatives are review hypotheses, not substitutions or local-occurrence/status conclusions. A newly observed taxon claim requires independent sequence and reference validation beyond these label comparisons.

## 6. Primary outcomes (declared)

- Per-species detection probability; estimated false-negative rate among published absences;
  corrected-vs-raw occupancy shift; occumb design-optimization recommendation;
  per-rank assignment rate (% ASVs, % reads) per database.

## 7. Post-fit disposition and availability amendment (2026-07-18)

This section records outcomes without rewriting the historical choices above.

- **Executed Layer-A formula:** occupancy ψ ~ `Dist_from_lake + season`; detection θ ~
  z-scored `log(read depth)`, in addition to `occumb`'s native read-count exposure.
  `river_gradient_percent` was dropped because E14 is missing; it was not fabricated.
  E14 temperature/rainfall use a documented E13 donor in supporting metadata, but
  temperature, rainfall, and catchment-wide flow were not in the executed headline formula.
- **Claim gate:** the smoke-scale known-truth recovery MCMC failed convergence
  (maximum ψ Rhat 6.61 against a 1.05 threshold), and 18 observed libraries that
  are zero within the 14-label harmonised Layer-A input exposed an inconsistent
  comparison mask. They contain reads under raw taxonomy labels. Corrected-vs-raw magnitude,
  direction, coefficients, and the utility curve remain diagnostics, not primary
  outcomes or survey recommendations. Production MCMC is optional future research.
- **Layer E inference amendment:** replicate libraries were aggregated before
  inference. Spatial tests use 14 site profiles over 11 common timepoints with all
  27 cyclic/mirrored river-order permutations; seasonal tests use 18 timepoint
  profiles over a fixed 10-site panel with 35 chronological-series permutations.
  Package-level DA consensus is additionally gated by restricted site-series BH;
  no taxon passes the final gate.
- **Layer B scope:** the completed audit has three classification comparisons across
  two current reference resources, not three independent databases. A current
  public-record witness set shows that ASV0001's 170-bp haplotype carries at least
  five *Salmo* species labels. Genus is stable; species is not identifiable here.
- **Control denominator:** 79 negatives occur in the manifest and 77 enter decontam;
  two well blanks have no ASV-table row. This is reconciled in
  `data/metadata/control_reconciliation.md`.
