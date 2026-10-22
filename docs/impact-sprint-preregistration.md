# Phase 5 scientific-impact sprint — local analysis specification

**Specification date:** 2026-07-19  
**Status:** frozen before execution or inspection of any new Phase-5 result  
**Scope:** repository-local decision record, not an externally registered preregistration

## 1. Purpose and inherited claim boundaries

This specification turns limitations of the completed River Conwy fish-eDNA release into six
bounded extensions that can run on an ordinary laptop. It was written after inspecting the
existing release and its already-known results, but before running the new Phase-5 analyses.
It therefore constrains Phase-5 choices without claiming independent or external timing proof.

The following inherited boundaries cannot be relaxed by a Phase-5 result:

- replicate libraries are observations within site x timepoint units, not independent river
  communities;
- reads are not abundance, biomass, diet, realised habitat use or organism location;
- the source-study table is a comparator, not biological truth;
- the unconverged real Conwy Layer-A recovery fit supplies no numerical direction or magnitude;
- ASV0001 is not identified to species by its existing 170-bp 12S fragment;
- the existing restricted-series result does not establish a section centroid shift; and
- public-record organism names are submitted annotations, not unquestionable identities.

The machine-readable companion is `data/metadata/impact_sprint_manifest.json`. If prose and
that manifest disagree on a factor level, seed, denominator or numeric gate, the manifest
controls. Changes after this freeze must be recorded as deviations before rerunning, never
silently substituted after seeing results.

## 2. Shared execution and exclusion rules

- Fixed seeds are derived from the manifest plus integer scenario/draw identifiers; language
  runtime hashes are never used as seeds.
- Structural absence, analytical exclusion and observed zero are distinct fields in every
  portable input. No unavailable library or site x timepoint cell is zero-filled unless the
  atlas is explicitly simulating that error.
- A failed fit, invalid interval, ambiguous accession, zero-variance null, singular design or
  missing specification remains in the output denominator with a failure reason.
- No outlier, taxon, site, timepoint, accession or specification is removed because its result
  is inconvenient. Only the exclusions below are permitted.
- Canonical tables use stable row ordering, decimal formatting and no wall-clock values. Runtime
  is reported separately. A second build must reproduce their SHA-256 digests.
- Reader-facing integration occurs only after the relevant gate is evaluated. A failed gate is
  integrated as a withheld claim with the evidence needed to change that outcome.

## 3. Occupancy recoverability atlas

### Question and estimand

Under which combinations of occupancy, detection, replication, sequencing-depth structure and
observation-mask integrity can marginal occupancy be recovered? The estimand is the mean true
occupancy probability across 120 ordered synthetic site x timepoint units. `prevalence_structure`
describes how the same target mean is distributed along the synthetic river; realised latent
prevalence is an outcome, not a second synonym for true occupancy.

### Known-truth design

The full 648-scenario factorial crosses:

- target mean occupancy: 0.10, 0.35, 0.70;
- prevalence structure: uniform or downstream-gradient (logit slope 2.0, intercept calibrated
  to the target mean);
- baseline per-replicate detection at mean depth: 0.20, 0.50, 0.80;
- planned replicates: 1, 2, 3;
- log-depth standard deviation: 0.15 or 0.75;
- depth-river correlation: 0 or +0.70; and
- mask scenario: complete/correct, 20% depth-linked structural missingness encoded correctly,
  or the same unavailable cells incorrectly encoded as observed zeros.

Each scenario has 100 fixed-seed simulations. Detection follows a logit-linear depth effect of
0.75. The generator persists separate truth, surveyed-mask, analysis-mask and detection fields.

### Estimators and uncertainty

1. **Naive occurrence:** fraction of analysable units with at least one detection; Wilson 95%
   interval, evaluated against true marginal occupancy while retaining its known estimand mismatch.
2. **Replication-aware likelihood:** constant occupancy and constant detection, integrating the
   latent state over binary replicate histories.
3. **Depth-aware likelihood:** occupancy on river position and detection on log depth, using the
   same integrated binary-history likelihood.

The likelihood estimators use deterministic bounded multi-start optimisation. Wald/delta 95%
intervals use the observed-information Hessian. One-replicate correction models are structurally
non-identifiable and must be labelled as such, not regularised into apparent estimates. Boundary,
rank-deficient and invalid-interval fits are failures.

### Gate

For each scenario x estimator, `recoverable` requires absolute bias <= 0.05, RMSE <= 0.10,
empirical 95% coverage in [0.88, 0.99], fit-failure rate <= 0.05 and interval-failure rate <= 0.05.
`non-identifiable` applies to a structural one-replicate correction or a combined fit/boundary/
singularity failure rate >= 0.20. Every other valid cell is `unreliable`. Coverage is conditional
on valid intervals; its denominator and the unconditional failure rates are reported beside it.

The atlas may support general synthetic design conclusions only. It cannot estimate Conwy
occupancy, repair the real mask or revive the withheld Layer-A correction.

## 4. In-silico multi-marker *Salmo* benchmark

### Question and estimand

Within a provenance-audited public-record panel, what is the finest submitted taxonomic rank
consistently discriminated by MiFish 12S, a 5-prime COI barcode-sized window and cyt-b? Estimands
are cross-label exact-haplotype sharing, normalised edit-distance distributions, label-specific
barcode gaps and conservative leave-one-accession-out nearest-neighbour classification.

### Accession panel and extraction

After this specification is frozen, a network-only fetch step may query NCBI Nucleotide for up to
50 complete mitochondrial records per submitted label, cache the raw GenBank records and query
receipts, then select at most five passing accessions per label in accession-version order.
Candidate *Salmo* annotations are `S. trutta`, `S. salar`, `S. carpio`, `S. marmoratus`,
`S. macedonicus`, `S. macrostigma` and `S. obtusirostris`; `Salvelinus alpinus` and
`Oncorhynchus mykiss` are outgroups for a genus boundary. Duplicate voucher/isolate records,
ambiguous marker sequences, missing features and non-unique primer products are excluded with
their reason retained.

Marker definitions are independent of ASV0001:

- MiFish-U 12S: the unique 140-250 bp in-silico product between the documented primers, allowing
  at most three substitutions per primer and recording both mismatch counts;
- COI: the first 648 unambiguous coding bases of the annotated mitochondrial COI/COX1 feature;
- cyt-b: the first 1,140 unambiguous coding bases of the annotated mitochondrial CYTB feature.

Distances use deterministic global edit distance divided by the longer sequence length; exact
haplotypes require identical strings. The same complete-mitogenome accessions are used across
markers whenever possible. A submitted label enters species-level summaries only with at least
two independent passing accessions. The species panel requires at least four eligible *Salmo*
labels; the genus boundary requires at least two eligible *Salmo* labels and both outgroup genera
with at least two accessions each. Inadequate coverage narrows or fails the relevant comparison.

### Gate

For each marker, species rank is supported only if every eligible leave-one-out query has all tied
nearest neighbours under the same submitted species annotation, no exact haplotype crosses labels,
and every eligible label has a positive barcode gap. Genus rank is supported only if species fails
but every eligible *Salmo* query remains nearer to *Salmo* than to either outgroup. Otherwise the
marker is `family-or-unresolved`; insufficient panels are `coverage-inadequate`.

These are submitted-label discrimination results. They cannot establish label truth, local
occurrence, stock, ancestry, hybrid status or a species identity for ASV0001.

## 5. Layer-E robustness multiverse

### Question and estimands

Does the descriptive downstream broadening remain directionally stable under a finite family of
defensible community summaries? The primary estimand is the lower-minus-upper mean detected-label
richness contrast. Supporting estimands are section pseudo-R2, exact restricted-series p-value,
and dispersion sensitivity. The existing balanced site-profile p=0.0714 estimand remains a
withheld centroid claim regardless of the multiverse.

### Specification families

The 72 core families cross:

- panel/summary: 14-site profiles over 11 common timepoints; 14-site profiles over all usable
  timepoints using per-site means/frequencies; or 154 balanced site x timepoint units with whole-site
  permutation blocks;
- taxonomy: 14 harmonised analysis labels or 12 genera;
- library-level detection threshold before aggregation: >0 or >=2 reads;
- effort: every surviving library or a 100-draw fixed-seed ensemble selecting one surviving library
  per site x timepoint unit; and
- distance: Bray-Curtis on proportions, Jaccard on occurrence, or robust Aitchison rCLR-Euclidean.

Section labels are permuted only through the 27 non-identity cyclic/mirrored river arrangements,
giving resolution 1/28 with the observed arrangement. Repeated site x timepoint rows always move as
whole-site blocks. Full one-library draw results are retained; the family row reports their median
and 5th-95th percentile envelope. The reference full-effort, >0, analysis-label balanced panel adds
75 influence rows: each of 14 sites and 11 timepoints omitted for each of three distances.

### Gate

`robust descriptive broadening` requires at least 90% of families to be valid, a positive richness
contrast in at least 80% of valid families, a positive contrast in at least 90% of all one-library
draws, and no sign reversal in any reference leave-one-site or leave-one-timepoint result. Otherwise
the effect-size distribution is reported and robustness wording is withheld. No exploratory
specification can convert the established restricted-series centroid result into confirmation.

## 6. Functional-trait null model

### Question and estimand

Does downstream functional broadening exceed the expectation created by detecting more labels?
Equal-axis Gower distance is computed from the machine-coded flow, trophic and movement factors.
For each of the 154 balanced site x timepoint units, presence/absence-weighted Rao quadratic entropy
is calculated (zero for richness one; missing for richness zero). Values are averaged within site;
the estimand is the lower-minus-upper mean site-level Rao contrast.

The primary panel contains the 12 labels already marked for descriptive inclusion. The strict
sensitivity retains only `Anguilla_anguilla`, `Barbatula_barbatula`,
`Gasterosteus_aculeatus`, `Phoxinus_ujmonensis` and the genus-stable `Salmo_trutta` proxy.

### Null and gate

Within each of the 11 timepoints, fixed-fixed checkerboard swaps preserve every unit's richness and
every label's timepoint prevalence exactly while breaking river-section placement. Seed 2026071903
generates 9,999 recombined null matrices after 5,000 accepted burn-in swaps per timepoint and at
least 100 accepted swaps between draws. Every draw is margin-checked.

The `beyond richness/prevalence` claim requires a positive observed contrast and one-sided Monte
Carlo p <= 0.05 in both primary and strict analyses, at least 500 unique null matrices, non-zero null
variance and defined site means in every section. Failure retains the effect and null interval but
withholds that wording. All interpretation remains conditional on analysis-label and trait proxies.

## 7. Cluster-aware pipeline disagreement

### Question and estimand

Among the 2,344 source-positive cells in the verified 642-library x 14-label shared comparison,
which recorded features are associated with 382 reprocessed-negative outcomes? The response is
pipeline non-recovery relative to the source table. Five source-only labels outside the mapped
namespace and reprocessed-only detections are excluded by design, not treated as errors.

The portable cell table must reproduce 1,962 recovered and 382 non-recovered cells and include
library ID, site, timepoint, replicate, site x timepoint unit, section, season, taxon, source relative
abundance, raw-name reprocessed depth, reprocessed count and source-label prevalence. Depth is never
computed from the 14-label table because 18 observed libraries are zero only in that harmonised
view.

### Model, dependence and sensitivity

The full logistic model uses z-scored log10 source relative abundance, z-scored log10(1 + raw-name
reprocessed depth), z-scored label prevalence, section and season. Absolute predictor correlation
>0.70 or design condition number >30 triggers a predeclared nested model: abundance/depth remain
primary and prevalence is reported separately. No coefficient is silently dropped.

Uncertainty comes from 999 fixed-seed block-bootstrap fits. Sites are resampled with replacement
within the 5/5/4 sections; the 18 ordered timepoints are resampled as circular moving blocks of
length three; taxa are resampled as whole clusters. Their multiplicities form observation weights,
so all rows within selected libraries and taxa travel together. Sensitivities omit each of 11 taxa,
14 sites and 18 timepoints in turn.

### Gate

A predictor association is promoted only when at least 95% of bootstrap fits are finite, its 95%
percentile interval excludes zero, its sign is unchanged in at least 90% of relevant leave-one-group
fits, and no single taxon/site/timepoint reverses the principal conclusion. Failed optimisation,
separation and omitted groups remain in denominators. If no predictor passes, clustered descriptive
rates are still a completed result. Neither outcome is a biological false-negative rate.

## 8. Interactive claim-survival map

The canonical table must contain exactly these claim IDs: `river-community`,
`imperfect-detection`, `pipeline-disagreement`, `taxonomic-resolution` and `trait-diversity`.
Each row records initial claim, audit, evidence, method, limitation, gate status, surviving wording
and next evidence. `admissions/claim-survival.html` must be keyboard usable, work from `file://`,
require no network and link to existing local evidence/method/limitation files. A static fallback is
required. The page passes only when all local links resolve and every wording/status agrees with
the machine-readable claim ledger and reader-facing release.

## 9. Completion and demotion rule

Each extension is complete only with a deterministic script, machine-readable result, figure,
method, limitation, tests and synchronized reader-facing status—or with a failed gate, withheld
claim and explicit next evidence. The impact sprint is complete only after a second deterministic
build, full release verification, rendered-artefact inspection and clean ZIP extraction under a
path containing spaces all pass with zero failures.
