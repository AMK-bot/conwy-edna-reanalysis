# Cluster-aware pipeline-disagreement analysis

## Question

Among cells that are positive in the source-study table, where does the reprocessed table fail to
recover the corresponding harmonised analysis label? This is a comparator-relative pipeline
question. Neither table is biological ground truth, so the response is not a biological
false-negative rate.

## Denominator and portable input

The analysis reconstructs the verified comparison used by `R/08_validate_concordance.R` and
persists it as a small offline TSV. The exact denominator is **2,344 source-positive cells** across
**642 shared libraries** and **11 source-positive labels** within the 14-label mapped namespace:
1,962 are recovered and **382 are source-positive/reprocessed-negative**. Five source-only labels
outside that mapped namespace and reprocessed-only detections do not enter this conditional model.

Each row carries library, site, timepoint, replicate, site x timepoint unit, section, season, label,
source relative abundance, reprocessed count, source-label prevalence and two explicitly different
read totals. The preregistered `reprocessed_raw_name_depth` is the sum of all 35 raw-name columns,
including `Unassigned` and the focal reprocessed count. Because `miss = 1` exactly when that focal
count is zero, this primary predictor has a direct part-whole coupling to the response. It is
retained to preserve the frozen primary analysis, but cannot by itself establish a depth
association.

The added sensitivity is cell-specific and removes that coupling. First, the assigned-fish total
is summed over the 34 named fish columns, excluding `Unassigned`. Then every raw-name count mapped
to the focal harmonised label is subtracted. Thus
`reprocessed_assigned_fish_reads_excluding_focal_label + reprocessed_count` exactly reproduces the
assigned-fish total in all 2,344 cells. The resulting total is non-negative; 224 cells are zero and
remain defined through `log10(1 + total)`.

## Model and dependence

The predeclared full logistic model considered source relative abundance, the 35-column raw-name total,
source-label prevalence, section and season. Source relative abundance and label prevalence had
absolute correlation **0.735**, exceeding the frozen 0.70 gate, so prevalence was removed from the
primary model and retained in the predeclared nested model. The primary estimand is therefore:

```text
pipeline non-recovery ~ z(log10(source relative abundance))
                      + z(log10(1 + 35-column raw-name total))
                      + section + season
```

The focal-excluded assigned-fish sensitivity replaces only the coupled total:

```text
pipeline non-recovery ~ z(log10(source relative abundance))
                      + z(log10(1 + assigned-fish reads excluding focal label))
                      + section + season
```

Uncertainty uses 999 fixed-seed block-bootstrap fits. Sites are resampled within river section,
timepoints as circular moving blocks of length three, and taxa as whole clusters. All taxon rows in
the selected library/site/time blocks travel together. A bootstrap draw remains valid when the two
primary continuous coefficients are finite; a season or section absent from a resample is retained
as a recorded nuisance alias rather than falsely counted as failure. The primary and focal-excluded
models use the **same 999 resampling weights**, and all 999 fits were finite for both. Both models
also use the same omissions of each of 11 taxa, 14 sites and 18 timepoints: 43 influence fits in
total.

## Result and gate

Both preregistered primary associations pass the clustered interval and influence gates:

- **Source relative abundance:** log-odds coefficient -1.574 per SD; 95% block-bootstrap interval
  **[-2.950, -0.807]**. The corresponding conditional odds ratio is 0.21 per SD.
- **35-column raw-name total:** coefficient -1.112 per SD; interval **[-1.972, -0.719]**;
  conditional odds ratio 0.33 per SD. This total includes both `Unassigned` and the focal count, so
  this primary result is recorded but not used alone for the surviving depth statement.

The focal-excluded model also passes:

- **Assigned-fish reads excluding the focal label:** coefficient **-1.149** per SD; 95%
  block-bootstrap interval **[-3.086, -0.579]**; conditional odds ratio 0.32 per SD. All **999/999**
  bootstrap fits were finite, and all **43/43** leave-one-group estimates were negative (range
  **[-1.738, -1.054]**).
- **Source relative abundance in the same sensitivity model:** coefficient -1.615; interval
  **[-2.987, -0.835]**, also with no leave-out sign reversal.

Therefore the processed count-total direction **survives focal exclusion**. The surviving result is:

> Pipeline disagreement is concentrated in source-positive cells with lower source relative
> abundance and where fewer reprocessed reads are assigned to non-focal fish labels, after section
> and season are included; the directions survive site, timepoint and taxon influence checks.

This quantifies non-recovery relative to a comparator. It does not show that the source call is a
true biological presence, that the reprocessed zero is a biological false negative, or that deeper
sequencing would causally recover a fish occurrence.

The nested prevalence coefficient is reported but is not promoted: its collinearity with source
relative abundance prevented inclusion in the primary clustered claim.

## Limits

- Source relative abundance is itself a processed compositional value and not organism abundance.
- Both read totals are downstream pipeline quantities and may be endogenous to other processing
  decisions. Focal exclusion removes the direct part-whole term but not compositional dependence,
  omitted processing factors or confounding; the association is not a causal depth intervention.
- Only 14 sites and 11 source-positive mapped labels are available, although the block and influence
  design retains their dependence.
- Harmonisation groups database labels for comparison and is not independent taxonomic validation.
- The block bootstrap represents the observed site/time/taxon structure; it does not create new
  rivers or independent catchments.

## Reproducibility and files

```bash
Rscript R/15_pipeline_disagreement.R --rebuild-input
Rscript R/test_pipeline_disagreement.R
```

The paired 999-draw analysis took about eleven seconds locally and uses seed `2026071904`.

- `data/metadata/impact_pipeline_disagreement_cells.tsv` — 2,344-row portable comparator with the
  assigned-fish and focal-excluded per-cell totals
- `data/derived/impact_pipeline_model.tsv` — coefficients, clustered intervals and claim gates
- `data/derived/impact_pipeline_bootstrap.tsv` — paired primary/focal-excluded coefficients for the
  same complete draw grid, including aliases
- `data/derived/impact_pipeline_leaveout.tsv` — paired coefficients for every taxon/site/timepoint
  omission
- `data/derived/impact_pipeline_rates.tsv` — transparent grouped disagreement counts
- `data/derived/impact_pipeline_summary.tsv` and `data/metadata/impact_pipeline_receipt.tsv`
- `figures/impact_pipeline_disagreement.svg` / `.png`
- `R/15_pipeline_disagreement.R` / `R/test_pipeline_disagreement.R`
