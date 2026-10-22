# Occupancy recoverability atlas — methods, result and claim boundary

**Status:** completed synthetic known-truth analysis  
**Frozen specification:** docs/impact-sprint-preregistration.md and
data/metadata/impact_sprint_manifest.json  
**Evidence scope:** synthetic_known_truth_only

## Question and estimand

The atlas asks when marginal occupancy is recoverable under different combinations of true
occupancy, prevalence structure, per-replicate detection, replication, log-depth variation,
depth–river coupling and observation-mask integrity. The estimand is the mean true occupancy
probability across 120 ordered synthetic site x timepoint units, not the realised fraction of
latent occupied units in one draw.

This is a general design experiment. It does not estimate River Conwy occupancy, repair the real
Conwy observation mask, recover an effect from the unconverged Layer-A fit, or establish the
direction or magnitude of a real sparse-label correction.

## Frozen design and generator

The complete 648-scenario factorial and 100 simulations per scenario were run without shrinking
the preregistered design: 64,800 scenario simulations and 194,400
scenario x simulation x estimator rows.

Synthetic river position is evenly spaced from -0.5 to +0.5. Uniform scenarios use constant
occupancy. Gradient scenarios use the frozen logit slope 2.0, with the intercept solved so the
mean of 120 occupancy probabilities equals the target (0.10, 0.35 or 0.70) within 1e-10.

For each replicate, centred log depth is constructed from standardised river position and a
random residual made exactly orthogonal to it. The realised sample SD (0.15 or 0.75) and
depth–river correlation (0 or 0.70) therefore equal their factor levels to floating-point
tolerance. Detection follows:

    logit(p_ij) = logit(p_baseline) + 0.75 * log_depth_ij

Structured unavailability removes the lowest-depth 20% of planned cells exactly.
structured_correct retains these as structural NA; structured_zero_fill uses the identical
latent draw, depths, detections and surveyed mask, but incorrectly makes unavailable cells
analysable zeroes. All mask scenarios share a seed within each truth-scenario x simulation pair,
isolating mask encoding as the only contrast.

The generator holds distinct true-probability, latent-state, planned-depth,
detection-probability, generated-detection, surveyed-mask, analysis-mask, observed-detection and
analysis-detection fields. The 1,080-row fixture exposes these fields cell by cell; the
64,800-row simulation audit retains denominators and event counts for every draw.

Seeds are 2026071901 + 1000 × truth_scenario_id + simulation_id, using R
Mersenne-Twister/Inversion/Rejection with RNG version 4.0.0. Mask scenarios share the same
truth_scenario_id and seed.

## Estimators, uncertainty and failures

1. **Naive occurrence** is the fraction of analysable units with at least one detection, with a
   Wilson 95% interval. It is deliberately evaluated against true marginal occupancy despite the
   known estimand mismatch.
2. **Replication-aware likelihood** uses constant occupancy and detection and integrates the
   latent state over each binary replicate history.
3. **Depth-aware likelihood** models occupancy against river position and detection against log
   depth using the same integrated likelihood.

For an all-zero history the integrated likelihood contribution is
(1 - psi) + psi × product(1 - p_j); a positive history requires occupancy and contributes psi
times the corresponding Bernoulli history probability.

Correction estimators use three deterministic starts, L-BFGS-B coefficient bounds [-12, +12],
and the best finite converged likelihood. The analytic score and observed-information Hessian use
the conditional latent-state probability (Louis identity). The 95% interval is a logit-scale
delta/Wald interval for mean predicted occupancy over all 120 positions. One-replicate correction
models are returned directly as structural_nonidentifiable without fitting or regularisation.

A correction fit fails if optimisation fails, a coefficient reaches a bound, or information is
non-positive/rank deficient (eigenvalue tolerance 1e-8; maximum condition number 1e8).
Non-finite, unordered or degenerate intervals are interval failures. Every failed row and reason
remains in the fit table.

Bias and RMSE use valid point estimates. Coverage is conditional on valid intervals and reports
its numerator and denominator. Fit- and interval-failure rates are unconditional over all 100
simulations. Boundary and singularity flags can overlap.

The frozen gate requires absolute bias ≤0.05, RMSE ≤0.10, coverage 0.88–0.99, fit failure ≤0.05
and interval failure ≤0.05. One-replicate corrections and cells with combined identification
failure ≥0.20 are non-identifiable; other failed cells are unreliable.

## Result

Across all 1,944 scenario-estimator cells, **613 were recoverable, 732 unreliable and 599
non-identifiable**. Naive occurrence contributed 173 recoverable and 475 unreliable cells;
replication-aware likelihood contributed 236 recoverable, 114 unreliable and 298
non-identifiable cells; depth-aware likelihood contributed 204, 143 and 301, respectively.

Three bounded design conclusions survive:

- No estimator passed the joint gate in any baseline-detection 0.20 cell, even with three
  replicates. Low detection remained a recoverability barrier in this design.
- Correction models with one replicate were structurally non-identifiable in all 432 applicable
  cells. With three replicates and baseline detection 0.50, replication-aware and depth-aware
  likelihoods were recoverable in 68/72 and 69/72 cells, versus 28/72 for naive occurrence. At
  baseline detection 0.80 and three replicates the counts were 71/72, 60/72 and 66/72.
  Adding depth parameters was not uniformly beneficial: finite-sample confounding control
  sometimes traded against extra estimation instability.
- Recoverable cells declined from 241 under a complete correct mask to 207 under correctly
  encoded structured missingness and 165 when the same unavailable cells were zero-filled.
  Relative to the paired correct structured mask, zero-filling demoted 53 recoverable cells and
  promoted 11: a net loss of 42. This is not an estimate of error in the real Conwy mask.

There were no numerical optimiser failures. Of 64,800 fits per estimator, replication-aware had
25,621 fit/interval failures and depth-aware had 26,212; each includes 21,600 predeclared
one-replicate structural failures. Naive occurrence had no fit or interval failures, but often
failed bias, RMSE or coverage against the occupancy estimand.

## Provenance, determinism and runtime

The analysis source contains no RDS, CSV or TSV input reader. Its runtime audit records only the
two frozen specification files and verifies both SHA-256 hashes. The manifest and receipt record
layerA_real_inputs_read as an empty list, and every portable row is tagged
synthetic_known_truth_only.

On the verified local arm64 R 4.4.3 environment with four workers, the final baseline build took
**89.658 seconds** and the second full build took **80.971 seconds**. SHA-256 digests matched for
all five tables and both figure formats. R/test_impact_atlas.R checks derivatives, fixed RNG
regressions, masks, dimensions, raw-to-summary denominators, gates, provenance, hashes and render
dimensions.

Run:

    Rscript R/11_occupancy_recoverability_atlas.R
    Rscript R/11_occupancy_recoverability_atlas.R --verify-rebuild
    Rscript R/test_impact_atlas.R

## Outputs

- data/derived/impact_occupancy_scenarios.tsv — 648 scenarios
- data/derived/impact_occupancy_simulation_audit.tsv — 64,800 generator audits
- data/derived/impact_occupancy_fit_results.tsv — 194,400 estimator rows
- data/derived/impact_occupancy_atlas.tsv — 1,944 gate cells with denominators
- data/derived/impact_occupancy_generator_fixture.tsv — explicit cell-field audit
- data/derived/impact_occupancy_manifest.json — input, design, hash and scope manifest
- data/derived/impact_occupancy_receipt.json — runtime and rebuild receipt
- data/derived/impact_occupancy_test_receipt.tsv — independent test receipt
- figures/impact_occupancy_atlas.png and .svg — complete phase diagram

Transfer to the River Conwy still requires a resolved real observation mask and an identifiable,
converged real-data analysis. The synthetic atlas substitutes for neither.
