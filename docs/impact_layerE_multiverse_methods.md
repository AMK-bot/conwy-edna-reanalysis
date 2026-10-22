# Layer-E robustness multiverse — methods and gate record

**Specification:** `docs/impact-sprint-preregistration.md`, section 5  
**Machine-readable contract:** `data/metadata/impact_sprint_manifest.json`  
**Seed:** 2026071902  
**Gate status:** **FAILED_GATE — robustness wording withheld**

## Question and claim boundary

This analysis asks whether the observed downstream increase in detected analysis-label
richness remains directionally stable across the finite 72-family multiverse frozen on
2026-07-19. It does not test abundance, biomass, organism location or a causal river-section
effect. Replicate libraries remain observations within site x timepoint biological units;
they are never treated as independent river communities. The established balanced-site
centroid result remains exact `p = 2/28 = 0.0714` and therefore withheld, regardless of any
multiverse result.

## Portable input and observation mask

`scripts/export_impact_inputs.R` is the only impact-sprint bridge that reads the existing
647 x 14 analysis-label and 647 x 12 genus RDS matrices. It writes the committed,
serialization-independent `data/derived/impact_library_counts.tsv`, with one row per
observed reprocessed library and explicit fields for structural absence, analytical
exclusion and an observed all-target zero.

The exporter stops unless all of the following hold:

- both matrices have their frozen dimensions, identical unique library IDs, finite
  non-negative whole counts and no missing values;
- the 12-genus matrix is the exact column collapse of the 14-label matrix;
- analysis-label, genus and stored fish-read depth agree for every library;
- all 647 rows are surveyed, usable observations in the occasion mask;
- there are exactly 234 realized site x timepoint units with 1-3 surviving libraries,
  14 ordered sites, 11 timepoints represented at every site and 18 genuine observed-zero
  libraries; and
- metadata-derived sections reproduce the fixed upper E01-E05, middle E06-E10 and lower
  E11-E14 design.

Input/output SHA-256 values and every passed check are stored in
`data/metadata/impact_library_export_receipt.tsv`.

## The 72 specification families

The complete Cartesian grid is retained in
`data/derived/impact_layerE_specifications.tsv`:

| Dimension | Frozen levels |
|---|---|
| Panel / summary | 14 site profiles over 11 common timepoints; 14 site profiles over all 234 realized units; 154 balanced site x timepoint units |
| Taxonomy | 14 harmonised analysis labels; 12 genera |
| Library detection threshold | at least 1 read; at least 2 reads, applied before aggregation |
| Library effort | every surviving library; 100-draw ensemble selecting one surviving library per realized site x timepoint unit |
| Distance | Bray-Curtis on proportions; Jaccard on occurrence; robust Aitchison rCLR-Euclidean |

For the balanced site-profile panel, filtered library counts are summed within site x
timepoint and then across the 11 shared timepoints, reproducing the established Layer-E
reference. For the all-usable profile, each site contributes one row: counts are averaged
over its realized timepoints and occurrences are represented as timepoint frequencies. For
the 154-row panel, replicate libraries are first aggregated or equalised within each unit,
and all repeated rows travel together as whole-site permutation blocks.

The primary richness value is always based on detected labels at the site x timepoint unit.
Site-profile specifications first average unit richness within each site and then average
sites within section; because the balanced panel has 11 timepoints per site, its site-weighted
and unit-weighted section means coincide. The estimand is lower-section minus upper-section
mean detected-label richness.

Bray-Curtis uses row proportions. Jaccard uses binary occurrence in the balanced and unit
panels and occurrence frequency in the all-usable site profiles. For these two distances,
two genuinely empty observed communities have distance zero and an empty-versus-detected
pair has distance one. rCLR is not defined for an entirely empty community: no pseudocount,
row deletion or invented composition is substituted.

## Equal-effort draws

The seed for draw `d` is `2026071902 + d`. Each of the 100 draws selects uniformly among the
1-3 surviving libraries in every one of the 234 realized units, using sorted, stable unit
order. The resulting 23,400 selections are committed in
`data/derived/impact_layerE_library_draws.tsv`. The same draw is reused across taxonomy,
threshold and distance families, so differences among specifications are not created by
different random library choices.

All 3,600 one-library family-draw results are retained in
`data/derived/impact_layerE_draw_results.tsv`, including 400 rCLR failures. Ensemble family
rows report medians and type-7 5th-95th percentiles; no failed draw is silently removed from
the family-validity denominator.

## Restricted inference and dispersion sensitivity

The exact spatial set contains 13 non-identity cyclic rotations and 14 rotations of the
mirrored E01-E14 order. Together with the observed arrangement this yields 28 arrangements
and exact p-value resolution `1/28 = 0.035714...`. Section labels are permuted at the site
level. In the 154-row panel, all 11 timepoint rows for a site receive the same permuted
section label.

Pseudo-R2 and pseudo-F are calculated from the doubly centred squared-distance matrix. The
exact p-value is `(1 + number of permuted F values at least as large as observed F) / 28`.
Dispersion sensitivity uses distances to section centroids and its observed one-way F;
the multiverse does not promote dispersion p-values as confirmatory evidence.

The reference influence analysis fixes the full-effort, at-least-one-read, 14-label,
balanced site profile and omits each of 14 sites and 11 common timepoints for all three
distances: exactly 75 rows. For a site omission, the frozen 27 full-river section
arrangements are generated first and the fixed physical site is then removed, preserving
the declared 1/28 resolution. All influence results are retained in
`data/derived/impact_layerE_influence.tsv`.

## Validity rule and observed result

A family counts as valid only when every declared evaluation row returns finite distance,
PERMANOVA and dispersion results. This conservative rule operationalises the frozen
requirement that failed fits and invalid specifications remain in the denominator. The
richness estimand is still reported when a supporting distance is undefined.

The observed results are:

- **64/72 families valid (0.8889): failed the required 0.90 validity fraction.** The eight
  retained failures are every rCLR family on the 154-unit panel. The full-effort panel has
  152/154 non-empty communities; equalised draws have 143-151/154. rCLR has no composition
  for the remaining observed-zero rows, so those families were not repaired by deletion or
  pseudocount.
- The lower-minus-upper richness contrast is positive in all 72 family summaries, including
  the eight with an undefined supporting rCLR distance. Among the 64 valid families, 64/64
  are positive, passing the required 0.80 directional fraction. Family medians range from
  3.0773 to 5.0924 detected labels.
- All **3,600/3,600** declared one-library family-draw contrasts are positive, passing the
  required 0.90 draw fraction. Their full range is 2.7364-4.2832 labels (median 3.4362;
  5th-95th percentiles 2.9227-4.1097).
- All 75 influence rows are valid and none reverses the reference sign. Their richness
  contrasts range from 4.6848 to 5.2000 labels.
- The established full-effort, >0-read, 14-label balanced-site result reproduces pseudo-R2
  0.8788 (Bray-Curtis), 0.7532 (Jaccard) and 0.5638 (rCLR), with exact `p = 0.0714` for each.
  Across the 64 valid family summaries, 58 have median exact `p = 0.0714` and six have
  `p = 0.1429`; none reaches the minimum possible 0.0357 threshold.

Because 64/72 = 0.8889 misses the frozen 0.90 validity requirement, the composite
`robust descriptive broadening` wording is withheld even though every defined richness
contrast is positive and no influence analysis reverses direction. The defensible result is
therefore: **downstream detected-label broadening is strongly directionally stable across
the declared richness summaries, but the preregistered multiverse robustness gate fails
because rCLR is undefined for empty observed site x timepoint communities.** This does not
establish a section centroid shift.

The complete gate arithmetic is in `data/derived/impact_layerE_gate.tsv`; the complete
family result table is `data/derived/impact_layerE_results.tsv`.

## Rebuild, tests and deterministic receipts

From the repository root in the documented `edna` environment:

```bash
Rscript scripts/export_impact_inputs.R
Rscript R/13_layerE_multiverse.R
Rscript R/13_layerE_multiverse.R
Rscript R/test_layerE_multiverse.R
```

The last verified rebuild used R 4.4.3, vegan 2.7-5, jsonlite 2.0.0 and digest 0.6.39.
Its R elapsed time was 53.035 seconds; runtime is kept separately from canonical tables in
`data/metadata/impact_layerE_run_receipt.tsv`. The second build reproduced byte-identical
SHA-256 digests for all six canonical analysis tables and both figure formats; those eight
comparisons are in `data/metadata/impact_layerE_determinism_receipt.tsv`. The test script
performs 38 unit, invariant and regression checks, including masks, depths, the complete
factorial, draw denominators, whole-site blocks, exact p support, the established Layer-E
reference, gate arithmetic, influence coverage and current artifact hashes.

The high-information specification curve, result matrix and influence display are in
`figures/impact_layerE_multiverse.png` and `.svg`.
