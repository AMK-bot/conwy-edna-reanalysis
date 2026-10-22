# Functional-trait broadening against a fixed-margin null

## Question

Does the visually strong downstream expansion of conditional fish-trait space exceed what is
expected merely because downstream site x timepoint units contain more detected analysis labels?

## Design

The analysis uses the 154 balanced units (14 sites x 11 common timepoints) and only occurrence,
never read abundance. Three machine-coded categorical axes from the cited trait table—flow
association, trophic group and movement group—receive equal weight in a Gower distance. A unit's
conditional functional diversity is presence/absence-weighted Rao quadratic entropy; singleton
units have Rao Q = 0 and zero-label units are undefined. Unit values are averaged within site,
then the lower-minus-upper mean site contrast is evaluated.

The primary panel contains the 12 labels previously admitted to descriptive trait summaries. The
strict panel retains five better-supported proxies: `Anguilla_anguilla`, `Barbatula_barbatula`,
`Gasterosteus_aculeatus`, `Phoxinus_ujmonensis` and the genus-stable `Salmo_trutta` proxy.

For each panel, 9,999 fixed-seed null matrices preserve every unit's detected-label richness and
every label's prevalence separately within each timepoint. Checkerboard swaps break the placement
of labels along the river while retaining those margins. Two primary and seven strict timepoint
blocks have only one state under their fixed margins; they remain unchanged as structural
single-state components. Every movable block met the declared burn-in/thinning target, every draw
passed its margins, and the primary and strict runs produced 9,999 and 9,988 unique matrices,
respectively.

## Result and gate

The observed downstream functional contrast is large descriptively: **0.5094 Rao-Q units** in the
12-label panel and **0.5288** in the strict panel. The fixed-margin null, however, expects nearly the
same contrasts because it preserves the strong downstream richness structure: null means are
**0.5081** and **0.5256**. The one-sided Monte Carlo p-values are **0.3159** and **0.1142**.

The preregistered claim gate therefore fails in both panels. The surviving conclusion is:

> Conditional trait space broadens downstream in the detected-label panel. Both observed contrasts
> were slightly above their null means, but neither was unusually large under the fixed-margin
> richness and label-prevalence null; the beyond-margins mechanism is unsupported.

This is a substantive negative result: the earlier ecological description remains useful, while a
stronger claim of functional assembly beyond detected-label richness is withheld.

## Limits

- Traits are conditional proxies for analysis-label assignments, not verified local species traits.
- The statistic does not measure abundance, diet, realised habitat use, environmental filtering or
  causal community assembly.
- Three broad categorical axes deliberately favour a transparent, auditable comparison over a
  high-dimensional trait synthesis.
- The null conditions on the observed richness and timepoint prevalence margins; it asks a narrow
  question and is not a universal ecological null.
- Fixed blocks are a property of the observed margins, not failed randomisation. Their counts are
  retained in the receipt.

## Reproducibility and files

Run:

```bash
Rscript R/14_trait_null_model.R
Rscript R/test_trait_null_model.R
```

The full run uses seed `2026071903`, writes 19,998 null-draw rows, and took about 78 seconds in the
verified local environment. `Rcpp` accelerates the preregistered checkerboard swaps; all analysis
and validation remain offline.

- `data/derived/impact_trait_gower.tsv` — complete primary and strict Gower matrices
- `data/derived/impact_trait_unit_rao.tsv` — observed unit richness and Rao Q
- `data/derived/impact_trait_null_draws.tsv` — all null contrasts and invariant receipts
- `data/derived/impact_trait_null_summary.tsv` — effect, null interval, p-value and gate
- `data/metadata/impact_trait_null_receipt.tsv` — seed, denominators and allowed wording
- `figures/impact_trait_null_model.svg` / `.png` — site pattern and both null distributions
- `R/14_trait_null_model.R` / `R/test_trait_null_model.R` — deterministic analysis and tests
