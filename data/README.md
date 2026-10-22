# `data/` — provenance, accessions, versions, access dates

Core source-data accesses described below: **2026-06-18**. Later biological-evidence
access dates are recorded per source in `metadata/fish_traits_sources.tsv` and
`metadata/salmo_reference_provenance.tsv` (including 2026-07-18 accesses). The
accession-backed Phase-5 marker benchmark records its separate 2026-07-19 NCBI
retrieval in `metadata/impact_salmo_queries.tsv`.

## Phase-5 portable evidence

The impact sprint is specified by `metadata/impact_sprint_manifest.json` and its
repository-local freeze receipt. Small `derived/impact_*` outputs are committed
release evidence; analysis telemetry is kept separately from deterministic tables.

- `derived/impact_occupancy_*` — synthetic-only 648-scenario recoverability design,
  194,400 fit rows, result atlas, generator fixture and scope manifest. These files
  never estimate or repair real Conwy occupancy.
- `public/salmo_multimarker/` — official cached NCBI ESearch/ESummary receipts, 29
  assessed primary GenBank records and the extracted 63-sequence marker FASTA.
  Every receipt, raw record and sequence has a checked SHA-256 digest.
- `metadata/impact_salmo_*` and `derived/impact_salmo_*` — accession decisions,
  marker windows, coverage, distances, haplotypes, all-ties classification and
  rank gates. Submitted organism labels are annotations, not biological truth.
- `derived/impact_layerE_*` — the complete 72-family community multiverse, all
  3,600 one-library draw results, 75 influence rows and the failed full gate.
- `derived/impact_trait_*` — equal-axis Gower distances, 154-unit Rao values and
  both 9,999-draw fixed-margin null panels.
- `metadata/impact_pipeline_disagreement_cells.tsv` plus
  `derived/impact_pipeline_*` — the exact 2,344-cell comparator, 999 block-bootstrap
  draws, 43 leave-out groups and scoped association results.
- `metadata/claim_survival.tsv` — the checked five-claim synthesis used by the
  offline reader map. It links each surviving sentence to evidence, method,
  figure, limitation and next evidence.

Run `python3 scripts/run_impact_sprint.py` from the repository root to rebuild and
test the complete offline impact workflow. Running the same command a second time
must reproduce all 47 canonical scientific artifacts byte-for-byte.

## Sequence data — ENA `PRJEB48362` (ERP132733)

- **What we use:** the **755** River Conwy 12S read-runs whose `sample_alias` matches
  `E##_R#_T##` (14 sites E01–E14 × 19 timepoints T01–T19 × 3 reps R1–R3, minus 43 absent).
- `ena/PRJEB48362_read_run.tsv` — full study read_run metadata (1372 runs; **committed**).
- `download_ena.sh` — hard-filters `E##` and pulls the 755 paired FASTQs (~71 GB) into
  `ena/raw/` (**git-ignored**), md5-verified, resumable. `--metadata-only` refreshes the TSV.
- `verify_ena_join.py` — the Phase-0 gate; writes the two derived tables below + a report.

### Derived (committed)
- `metadata/conwy_12S_libraries.tsv` — clean `ID → run_accession → FASTQ` join (755 rows).
- `metadata/design_gaps.tsv` — the 43 absent site×timepoint×rep cells.
- `metadata/ena_join_verification.md` — computed verification report (re-run to refresh).
- `metadata/data_dictionary.md` — the one-page data dictionary.

## Published artifacts — Perry et al. 2024 (Nat. Commun., CC BY 4.0)

doi:10.1038/s41467-024-48640-3. Fetched by `download_published.sh` (**not committed**):

| file | local path | content |
|---|---|---|
| MOESM7 (Source Data, .xlsx) | `published/source_data_species_x_replicate.xlsx` | normalised relative abundance by species label × replicate; detection is preserved (6 figure-keyed sheets) |
| MOESM4 (Supp. Data 1, .csv) | `metadata/supplementary_data_1_covariates.csv` | per-sample metadata + environmental covariates (1141 rows × 40) |
| MOESM5 (Supp. Code 1, .zip) | `published/supplementary_code_1/` | authors' R scripts (incl. `Analysing_BLAST_output_WP2_12S.R`) |
| MOESM1 (Supp. Info, .pdf) | `../docs/supplementary_information.pdf` | methods detail |

## Key facts (verified against the live archive)

- ID ↔ ENA run: **1:1** (755↔755). ENA → covariates: **755/755**, all `River=Conwy`.
- **T16 is absent for every site** (+`E02_T17_R3`) → 18 usable timepoints; treat T16 as *unsurveyed*.
- Original taxonomy: BLAST+ v2.10.0 vs NCBI nt + **MitoFish (downloaded 23 Sep 2020)**.
  That exact MitoFish build is unarchived. Layer B uses the recorded source-study calls
  in current-resource concordance; it does **not** claim a faithful 2020→2026 replay.
