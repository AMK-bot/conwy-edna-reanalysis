# Phase 3 Layer A — claim and provenance status

**Status date:** 2026-07-18  
**Purpose:** define what the Shiny app may safely show before production occupancy MCMC exists.

## Current evidence boundary

- `data/derived/occumb_intercept_estimates.tsv` and the covariate correction/design tables were
  produced by the default **smoke** configuration in `R/06_occupancy_occumb.R`: 2 chains, 2,000
  iterations, 1,000 burn-in. They prove that the data design and model plumbing run; they are not
  publication estimates.
- Some smoke-fit parameters do not meet production convergence targets. The app therefore reports
  Rhat/effective sample size and labels every Layer A number as diagnostic.
- The recovery study (`R/10_sim_recovery.R`) simulated known occupancy under the observed depth ×
  biological-gradient confound, but its smoke-scale MCMC failed the convergence gate. Its tabulated
  direction and magnitude are therefore withheld, not interpreted. The script writes
  `data/metadata/layerA_recovery_status.tsv`, marks numerical rows `WITHHELD_UNCONVERGED`, and
  replaces its figure with the failed-gate status whenever convergence does not pass.
- A failed gate also runs `scripts/apply_layerA_claim_gate.R`, replacing all five cached
  occupancy/correction/coefficient/detection/design figures with explicit withheld-status cards.
  `data/metadata/layerA_figure_gate.tsv` records and hashes those replacements so the release verifier
  fails if an old numerical image reappears.
- The 14-label `occumb` input has 18 observed libraries whose harmonised vector is all zero;
  the raw-name count table contains positive fish reads in every one of those libraries.
  `occumb` structurally drops an all-zero community vector, while the naïve comparator counts those
  observed libraries as non-detections. This observation-mask/comparator mismatch is another reason
  the smoke corrected-vs-raw shifts are not inferential results.
- The standalone teaching simulator is a separate deterministic teaching model. At its fixed
  default it demonstrates that an incomplete depth-only calculation can overshoot known truth; it
  does not estimate or validate a direction for the Conwy Layer A model.

## Claims that are not currently allowed

- A rare species' corrected occupancy magnitude.
- The direction or size of a rare-species correction relative to naïve detection.
- A count or percentage of published absences that have been recovered as true presences.
- A final sampling-effort recommendation from the smoke-fit utility grid.
- “Production MCMC will confirm the numbers.” Longer chains address Monte Carlo convergence, not the
  structural recoverability problem exposed by the simulation.

## Claims the current artifacts can support

- The two-level site × timepoint / water-sample-replicate design is implemented and passes its
  structural-NA schema test.
- The `occumb` model executes end-to-end on the 647-library count table.
- The present rare-occupancy recovery exercise fails its convergence gate; this alone is sufficient
  to demote every recovery magnitude and direction and guide the next model comparison.
- The smoke fit generated heterogeneous adjustment diagnostics, but no cross-taxon magnitude or direction
  is interpretable.
- The present utility grid is withheld and cannot support sampling advice; the standalone fixed-seed
  simulator illustrates a generic confounding failure mode without estimating the Conwy design.

## Gate before any Layer A result becomes final

1. Run the production configuration (4 chains, 60,000 iterations, 20,000 burn-in, thin 20) and retain
   machine-readable run metadata plus the fitted object or sufficient posterior summaries.
2. Run the planned minimum-depth sweep (0, 100, 250, 500, 1,000, 2,000 reads) and compare depth-aware
   with depth-blind or redesigned detection models.
3. Persist Rhat, effective sample size, posterior predictive checks, and recovery diagnostics with the
   result tables.
4. Version and commit lightweight production summaries required by the Shiny app.
5. Update this note and the app's hard-coded evidence label only after those gates pass.

## Local app check

From the repository root:

```bash
EDNA_APP_SCHEMA_CHECK=1 micromamba run -n edna Rscript shiny/app.R
```

The check does not launch a browser. It validates the Layer A input schemas and numeric ranges,
constructs the UI/server, and exercises the species-selection reactive.
