# Imperfect detection and reference-database currency in River Conwy fish eDNA: a spatio-temporal reanalysis

> **STATUS: LOCAL EVIDENCE DRAFT (updated 2026-07-19).** Layers E and B have
> completed local outputs. Layer E now uses aggregated site/time profiles and
> restricted-series inference. Layer A is retained only as an exploratory smoke-fit
> stress test because the rare-species correction is unconverged and structurally
> confounded. A frozen Phase-5 credibility sprint adds synthetic, multiverse,
> accession, fixed-margin-null, and comparator-relative gates without promoting that
> real-data correction. This is **not** a submission-ready manuscript. Nothing
> has been submitted.

**Data/code:** ENA PRJEB48362 (River Conwy 12S); reanalysis repo `conwy-edna-reanalysis` (this repo).
**Original study reanalysed:** Perry et al. 2024, *Nature Communications* **15**:4372 (CC BY 4.0).

---

## Abstract *(local draft; authorship and references pending)*
The River Conwy source study reports eDNA *detections* without an occupancy model that separates
absence from missed detection. Reference-database choice also affects what a short sequence is
called. We reanalyse the River Conwy
12S fish eDNA dataset (14 sites × 18 usable timepoints, with up to 3 water-sample libraries per
site-time) and add three analytical audit layers: a compositionally aware community baseline (E);
a multispecies Bayesian site-occupancy stress test that separates occupancy from detection when the
data support doing so (A); and a re-classification of the same amplicon sequence variants (ASVs) against current
reference databases (B). We reprocessed the raw reads to an ASV count table (647 libraries × 1,348
ASVs) validated against the published result (all source-study taxa detected in at least three shared
libraries recovered; mean read-weighted Bray–Curtis similarity 0.84). Validated single-label models
show that six well-detected analysis labels have per-replicate detection estimates of 0.79–1.00.
A subsequent community-model stress test
produced unstable rare-taxon adjustment diagnostics, while its known-truth
recovery MCMC fails the relevant convergence gate. Neither the rare-species
correction magnitude nor direction is therefore reported as recovered truth. A
separate deterministic teaching lab demonstrates only the generic possibility of
over-correction under confounding.
On **reference-database currency**, the local Nov-2025 *Salmo macrostigma* first-hit label
accounts for 46.0% of reads in the 647-field-library count table (48.8% of assigned field reads), while cross-resource *species* labels disagree
across three comparisons (macrostigma / trutta / carpio) and **genus
is stable**. Short-12S evidence therefore supports *Salmo*, not a unique species.
After aggregating water-sample replicate libraries within each site×time unit, detected-label richness
broadens downstream, but restricted-series PERMANOVA does not cross 0.05. Two
unadjusted distance sensitivities give the same exploratory section-dispersion p-value. We therefore report a
descriptive gradient rather than a confirmed section centroid shift. A frozen credibility sprint
finds 613/1,944 synthetic estimator–scenario rows recoverable, but none at baseline
per-replicate detection `p = 0.20`; the full Layer-E multiverse gate misses its 90%
validity threshold (64/72), and fixed-margin trait nulls do not support a
beyond-richness interpretation (`p = 0.3159` primary; `p = 0.1142` strict). A
three-marker *Salmo* species benchmark is coverage-inadequate, although the genus
boundary passes 10/10 queries for every marker. Cluster-aware comparator analysis
supports higher pipeline non-recovery at lower source signal and a lower focal-excluded
assigned-fish count total, but neither count total is whole-library sequencing depth
and neither pipeline is biological ground truth. We give an honest account of
what eDNA here can and cannot say.

---

## 1. Introduction *(stable — scope and gap framing)*
- eDNA metabarcoding as a biomonitoring tool; the River Conwy dataset (Perry et al. 2024) as a
  spatio-temporal 12S fish survey.
- **Gap 1 — imperfect detection.** The source study reports observed detections without modelling detection probability; non-detection ≠ absence. No
  occupancy / detection-probability modelling was part of the original analysis.
- **Gap 2 — reference-database currency.** 2020-era BLAST vs NCBI nt + MitoFish (23 Sep 2020);
  taxonomy assignments can differ as databases grow. What differs under 2025/26 resources?
- **Gap 3 — compositional correctness.** Reads ≠ abundance; community inference here therefore requires
  compositional handling, independent biological units, dispersion checks, and bounded DA screens.
- **Contribution.** A three-layer reanalysis (E + A + B) on a validated reprocessed count table, with
  a transparent analysis-decision record and a dedicated honest-limitations treatment.

## 2. Methods

### 2.1 Data and reprocessing *(stable)*
- ENA **PRJEB48362** (River Conwy, site codes **E01–E14**, timepoints T01–T19, replicates R1–R3);
  ⚠️ `C##` in the same archive is a **different river** and is excluded. **755** deposited 12S
  replicate-libraries; T16 absent for all sites → **18 usable timepoints** (unsurveyed cells coded
  **NA, not 0**).
- Marker: **MiFish-U 12S** (~172 bp). Primers removed with `cutadapt` (linked adapters); ASVs by
  **DADA2** (per-run error learning, `pool="pseudo"`, consensus bimera removal). Contamination via
  `decontam` (prevalence, P\*=0.1 primary / 0.5 sensitivity) + an index-hop/tag-jump abundance floor
  (0.001), calibrated to the mackerel-spike positive control.
- Result: **`seqtab_clean` = 647 field libraries × 1,348 ASVs, 35.7 M reads** (survival 647/755 vs the
  authors' 665/755). Validation gate (`R/08`): a consensus DB-namespace crosswalk bridges the 2020↔2025
  name vintages, then compares to the published table over 642 shared replicates — **PASS-WITH-CAVEATS**
  (all ≥3-replicate taxa recovered; mean read-weighted Bray–Curtis similarity 0.84, cosine similarity 0.90, cell precision 0.99).
- **Governing caveat (motivates the detection stress test):** the reprocessing runs at ~84% of the source study's
  detection-cell total—a ~16% pooled deficit associated with a 35-column reprocessed
  raw-name-table count total (not whole-library sequencing depth) and low source-study
  relative abundance, which co-vary with the upstream-to-downstream biological-signal gradient. This
  is tested as a recoverability risk rather than assumed correctable (§2.3).

### 2.2 Layer E — compositionally-correct community baseline *(locally complete)*
- Alpha: richness/Shannon/Simpson/Pielou/**Hill q=0,1,2** with **coverage-based iNEXT** in **incidence**
  mode after replicate libraries are combined to **234 site × timepoint biological units**;
  reads ≠ individuals. Faith's PD omitted (no defensible 12S phylogeny).
- Beta: Bray–Curtis, Jaccard (presence/absence), **robust Aitchison** (rCLR→Euclidean); PCoA + NMDS.
- Spatial tests: 14 site profiles over 11 common timepoints; temporal tests: 18 timepoint profiles
  over a fixed 10-site panel. **PERMANOVA** is paired with **PERMDISP** and uses the complete cyclic
  plus mirrored river/chronological series (27/35 permutations).
- Differential abundance: the documented ≥2-of-3 package consensus (ALDEx2 + ANCOM-BC2 +
  corncob) is additionally gated by a site-level restricted-series CLR test. Indicator analysis uses
  the same restricted spatial series. These are conservative exploratory screens at 14 sites.

### 2.3 Layer A — multispecies Bayesian occupancy (`occumb`) *(exploratory stress test)*
- 2-level design: occupancy unit = **site × timepoint** (234 realized of 252); the 3 water samples =
  detection replicates. Read counts as input; structural-NA (never observed-zero) via
  `observable_mask("reprocessed")`; a 1,000-read minimum-depth filter in the smoke run. The planned
  0/100/250/500/1,000/2,000-read threshold sweep did not run, so no threshold-robustness claim is made.
- Occupancy ψ ~ `Dist_from_lake` + `season` (`river_gradient_percent` dropped — see §2.5).
- **Detection safeguard tested:** per-library z-scored log read depth enters the
  detection process because an unclustered diagnostic associates the recall deficit with
  the reprocessed raw-name-table count total and source-study relative abundance; both co-vary with
  the river gradient. The processed total is not whole-library sequencing depth.
  The current outputs come from a short two-chain smoke fit and are diagnostic.
- **Recovery check:** the known-ψ recovery MCMC fails its convergence gate. Its
  direction is withheld, and sparse-label corrected-vs-raw magnitudes, covariate
  coefficients, and design utility are not inferential results. The portfolio
  edition instead uses a separate deterministic teaching counterexample and does
  not depend on production MCMC.

### 2.4 Layer B — reference-database currency *(completed local comparison)*
- Re-assign the same ASVs with **MMseqs2 @ 97%** against current 12S databases and report per-rank
  assignment rate (% ASVs, % reads), concordance across DBs, and taxa that gain resolution / change /
  newly resolve, with DB-drift control.
- **Honest scope (Gate 4):** the authors' Sep-2020 MitoFish build is unarchived, so Layer B is **2026
  concordance, not a faithful 2020→2026 replay.** A same-reference-set classifier comparison measures method sensitivity;
  the independent current comparison uses the GenBank-derived MIDORI2 GenBank271 resource.

### 2.5 Covariates and the E14 gap *(stable)*
- Recorded covariate-availability choices: `Water_temp_oC` (100% NA) → `temp_monthly_average`; catchment flow handled per
  the data dictionary. **E14 (Trefriw):** its `river_gradient_percent`/temp/rainfall were never published
  by the authors (all-NA in their own Supplementary Data 1). A clearly flagged E13 donor was used only to
  complete the supporting climate table; it is an unvalidated spatial imputation, not an observed E14
  measurement. The executed headline occupancy formula excludes climate and gradient, and
  `river_gradient_percent` remains missing rather than being fabricated.

### 2.6 Reproducibility & decision history *(stable)*
- Two conda environments (arm64 Mac analysis + x86_64 WSL2-PC bioinformatics). The local report and
  quick-start capsule rebuild from included small inputs and verified checkpoints; a clean-clone raw-read
  reconstruction still requires the external FASTQs and the Phase-2 PC environment. `docs/prereg.md` records covariate sets, species inclusion, decontamination thresholds,
  the DA consensus rule, and the Layer B decision. Repository history co-commits the Phase-1 and Layer-E
  choices with their first outputs, so this is a transparent decision record, **not** an externally
  registered or independently pre-result preregistration. Phase-2 thresholds do appear before that run.

### 2.7 Phase-5 credibility sprint *(locally frozen before new runs)*
- **Occupancy atlas:** 648 synthetic known-truth scenarios × three estimators × 100
  simulations varied occupancy, prevalence structure, detection, replication, depth
  coupling, missingness, and mask error. Recoverability required the frozen bias,
  RMSE, coverage, and failure-rate gates. It does not estimate Conwy occupancy.
- ***Salmo* benchmark:** accession-backed MiFish 12S, 5-prime COI, and cyt-b windows
  used exact sharing, barcode gaps, all-ties leave-one-out classification, two outgroup
  genera, and a minimum four-label species-coverage gate. Submitted names were treated
  as annotations, not taxonomic truth.
- **Layer-E multiverse and trait null:** 72 finite specification families retained the
  site×time biological unit; 100 one-library draws per applicable family and 75
  reference leave-outs tested direction. Rao's Q was separately compared with 9,999
  fixed-margin checkerboard draws in primary 12-label and strict five-label panels.
- **Pipeline disagreement:** exactly 2,344 source-positive cells were modelled with
  site-within-section and circular time-block resampling, taxon-cluster weights, and
  43 taxon/site/time leave-outs. The response is non-recovery relative to the source
  table, not a biological false negative.
- The [frozen contract](impact-sprint-preregistration.md) and the five
  `impact_*_methods.md` records define these extensions; the
  [machine-readable claim record](../data/metadata/claim_survival.tsv) is the governing
  claim-by-claim synthesis.

## 3. Results *(local evidence draft; apply the claim boundary below)*
- **3.1 Community baseline (E).** Across all 234 realized units, median site-time detected-label richness
  rises from one label at E01–E05 to seven at E14; this descriptive summary retains unequal timepoint
  availability and 1–3 libraries per unit. The balanced 11-common-timepoint website hero preserves
  those endpoints (only E10 changes, from 4.5 to 4.0). In the corrected analysis, section PERMANOVA has p=0.0714 for Bray–Curtis,
  Jaccard, and robust Aitchison; season p-values are 0.333–0.778. Two unadjusted restricted-series
  section PERMDISP sensitivities reach p=0.0357 for Bray–Curtis and robust Aitchison, but Jaccard does not;
  no seasonal dispersion test reaches 0.05. This repeated-distance signal is treated as exploratory.
  Two taxa meet package-only DA consensus, but neither passes the restricted-series BH gate, and no
  indicator reaches p≤0.05 (minimum p=0.077). The result is **descriptive** downstream
  broadening with an exploratory dispersion signal, not a confirmed section centroid shift.
- **3.2 Occupancy & imperfect detection (A) — validated on-ramp plus an exploratory failure signal.**
  **(i) Validated well-detected-label result.** Phase-1 single-label models estimate
  per-replicate detection of 0.79–1.00 for six selected analysis labels; their average
  occupancy is consequently close to naïve occurrence. This result does not
  extend to sparse labels.
  **(ii) Exploratory community-model result.** Numerical rare-taxon adjustments from the smoke occupancy
  fits are withheld: their magnitude and direction cannot be interpreted. The **simulation recovery study** (§2.3 / §Methods), which fits
  `occumb` to data simulated at known ψ, is itself smoke-scale and fails the ψ convergence gate.
  Its apparent recovery direction and magnitude are therefore withheld. In addition, 18 observed
  libraries that are zero within the 14-label harmonised Layer-A input are dropped by the community model while the naïve comparator counts them as
  non-detections. **Rare-species occupancy is not reliably recoverable from the current evidence.**
  The local methodological lesson is the claim-gate decision: the well-detected-label estimates are usable,
  while unconverged and comparison-mismatched rare-taxon corrections are not. A longer run could
  reduce Monte Carlo error but would not by itself remove depth-biological non-identifiability or
  align the observation masks.
- **3.3 Reference-database currency (B).** In the 647-field-library count table, MitoHelper Nov-2025
  assigns 838/1,348 ASVs (62.2%) and 94.1% of reads. The broader 1,443-ASV field+control classification
  universe has 886 assigned ASVs (61.4%) and 93.3% assigned reads. **Headline assignment audit:** in the
  field table, the single most abundant first-hit call, *Salmo macrostigma*, accounts for 46.0% of all
  reads (48.8% of assigned reads; ASV0001 alone = 39.5% of the field table), but
  that label alone does not establish local identity. Cross-resource **species** labels disagree
  across three comparisons spanning two current reference resources — first-hit → *macrostigma*,
  same-vintage MitoFish + MMseqs2 → *S. trutta*, independent MIDORI2 GenBank271 → *S. carpio* —
  while **genus *Salmo* is stable.** Concordance with the independent DB bears this out: **genus agreement 99.4%, species
  agreement only 46.7%** (vs 67.3% for the same-reference-set classifier comparison). **The molecular
  assignment is stable only at genus; locality or native identity is not established by this 170 bp
  fragment, and none of the individual species labels is sufficient to establish a local species
  identity** — a concrete audit example for short-amplicon eDNA. A
  direct current witness alignment shows that the exact 170-bp haplotype carries at least five public
  *Salmo* species labels (`docs/salmo_molecular_detective_methods.md`). The
  all-eukaryote independent DB additionally surfaced labels consistent with potential non-target
  riparian vertebrate bycatch (otter *Lutra lutra*, water vole, roe deer, dipper); these assignments
  are outside the fish occupancy analysis and are not independent organism validation.
- **3.4 Phase-5 claim gates.** The synthetic atlas classifies **613/1,944**
  estimator–scenario rows as recoverable, **732** as unreliable, and **599** as
  non-identifiable. None of the 648 rows at baseline per-replicate detection
  `p = 0.20` is recoverable; at one replicate, both correction estimators are
  structurally non-identifiable in all 216 applicable rows each. These design results
  do not revive the real Layer-A correction. In the *Salmo* benchmark, only **2 of 7**
  candidate labels pass coverage for each of MiFish 12S, COI, and cyt-b (minimum
  **4** required), so
  every species gate is coverage-inadequate; the genus boundary passes **10/10**
  *Salmo* queries against each outgroup for all three markers. The Layer-E full gate
  fails at **64/72 valid families** (88.9% < 90%), despite positive direction in
  **64/64** valid families and **3,600/3,600** one-library draws and **0/75**
  leave-out reversals. The beyond-richness trait interpretation fails in the primary
  (`p = 0.3159`) and strict (`p = 0.1142`) fixed-margin nulls. Finally, the scoped
  pipeline-disagreement gate passes: standardized source-signal coefficient
  **-1.574** (cluster-aware bootstrap 95% interval **-2.950 to -0.807**) and
  focal-excluded assigned-fish count-total coefficient **-1.149** (**-3.086 to
  -0.579**), with all **43/43** leave-outs negative. The preregistered 35-column
  total contains `Unassigned` and the focal count; its -1.112 coefficient is retained
  as an audit, not used alone. These processed totals are not whole-library sequencing
  depth, and the coefficients describe non-causal, source-table-relative pipeline
  non-recovery—not biological false-negative rates.

## 4. Discussion *(stable framing; sharpen with numbers)*
- Why the exploratory occupancy correction is not a recovered biological result;
  distinguish the approximately 16% pipeline recall deficit from ecological false
  negatives in the original.
- DB currency: assignment stability at genus vs species; how DB vintage can produce apparent novelty.
- Credibility gates distinguish stable descriptive direction (Layer E), a supported
  comparator-relative pipeline association, and hypotheses that remain withheld
  (real rare-label correction, broader *Salmo* species resolution, and trait broadening
  beyond richness/prevalence).
- Management relevance is limited to claims supported in the cited fish-trait/ecology note; database
  labels and eDNA transport prevent literal species-location interpretation for ambiguous calls.

## 5. Limitations (MIEM-aligned) *(stable — a strength, stated up front)*
- Reads ≠ abundance (occurrence framing throughout). Non-detection ≠ absence — the quantity Layer A models.
- **Rare-species occupancy is not reliably recoverable here (the key Layer-A limitation):**
  the recovery MCMC fails convergence, and the community-model/naïve comparator masks differ for 18
  observed libraries that are zero within the 14-label harmonised input but contain raw-name fish reads.
  The validated well-detected-label Phase-1 result remains usable; smoke-fit
  rare-taxon directions, magnitudes, and design-utility values are withheld or diagnostic-only.
- The pooled ~16% detection-cell recall deficit; the flagged, unvalidated E14 supporting-table donor
  (excluded from the headline formula); Layer B as current-resource concordance rather than a
  2020 replay; low-resolution restricted permutations; exploratory unadjusted PERMDISP sensitivities; eDNA transport (authors'
  ≤5 km bound) — detection location ≠ organism location; 12S-MiFish misses lamprey/flounder resolution.
- Phase-5 extensions remain conditional on a synthetic occupancy generator, incomplete
  public *Salmo* record coverage, undefined rCLR distances for genuine empty communities,
  trait proxies, and a comparison in which neither pipeline is ground truth.

## 6. Data & code availability
- ENA PRJEB48362; local reanalysis repository with environment lock, scripts,
  Quarto technical report, and a small quick-start capsule. No public repository,
  DOI, or deployment is claimed.
- The offline Phase-5 rebuild is `python3 scripts/run_impact_sprint.py`. The readable
  [decision map](../website/claim-survival.html) and its
  [machine source](../data/metadata/claim_survival.tsv) preserve every gate and claim boundary.

## Key references

- Perry, W. B. et al. (2024). *An integrated spatio-temporal view of riverine biodiversity using environmental DNA metabarcoding.* **Nature Communications 15**, 4372. https://doi.org/10.1038/s41467-024-48640-3
- Fukaya, K., Kondo, N. I., Matsuzaki, S.-I. S. & Kadoya, T. (2022). *Multispecies site occupancy modeling and study design for spatially replicated environmental DNA metabarcoding.* **Methods in Ecology and Evolution 13**, 183–193. https://doi.org/10.1111/2041-210X.13732
- Stratton, C., Sepulveda, A. J. & Hoegh, A. B. (2020). *msocc: Fit and analyse computationally efficient multi-scale occupancy models in R.* **Methods in Ecology and Evolution 11**, 1113–1120. https://doi.org/10.1111/2041-210X.13442
- Bayer, P. E. et al. (2025). *A comprehensive evaluation of taxonomic classifiers in marine vertebrate eDNA studies.* **Molecular Ecology Resources 25**, e14107. https://doi.org/10.1111/1755-0998.14107
- Gloor, G. B., Macklaim, J. M., Pawlowsky-Glahn, V. & Egozcue, J. J. (2017). *Microbiome datasets are compositional: And this is not optional.* **Frontiers in Microbiology 8**, 2224. https://doi.org/10.3389/fmicb.2017.02224
- Hsieh, T. C., Ma, K. H. & Chao, A. (2016). *iNEXT: An R package for rarefaction and extrapolation of species diversity (Hill numbers).* **Methods in Ecology and Evolution 7**, 1451–1456. https://doi.org/10.1111/2041-210X.12613
- Klymus, K. E. et al. (2024). *The MIEM guidelines: Minimum information for reporting of environmental metabarcoding data.* **Metabarcoding and Metagenomics 8**, e128689. https://doi.org/10.3897/mbmg.8.128689
