# Data dictionary — Conwy 12S eDNA reanalysis

*Phase 0 deliverable. Companion to the computed `ena_join_verification.md`.*

## 1. Source

| | |
|---|---|
| **Paper** | Perry, W.B. *et al.* (2024) *An integrated spatio-temporal view of riverine biodiversity using environmental DNA metabarcoding.* **Nat. Commun. 15:4372**, doi:[10.1038/s41467-024-48640-3](https://doi.org/10.1038/s41467-024-48640-3) (CC BY 4.0) |
| **Sequence archive** | ENA **PRJEB48362** (ERP132733), project "LOFRESH", Bangor University |
| **Marker** | 12S rRNA fish, **MiFish-U** primers, ~172 bp amplicon |
| **System** | River Conwy, Wales — headwater lake (Llyn Conwy) → estuary |
| **Design** | 14 sites × 19 timepoints × 3 water-sample replicates; replicates sequenced separately (→ 2-level occupancy) |

> ⚠️ **Site-code trap.** Conwy = **`E01`–`E14`**. The same archive contains **`C##` = a DIFFERENT river** (coded by *season*, not timepoint), plus other rivers (`G*`, `TW*`, `SK*`) and markers. **Only `E##_R#_T##` is ours.** The download filter and `verify_ena_join.py` hard-enforce this.

## 2. Verified counts (see `ena_join_verification.md`)

- Total read-runs in PRJEB48362: **1372**. Conwy 12S (`E##_R#_T##`): **755**.
- ID ↔ ENA `run_accession` join: strictly **1:1** (755 ↔ 755). ENA → covariate table: **755/755**, all `River=Conwy`.
- All 755 runs are **paired-end** (1510 FASTQs, ~71 GB).
- **43 design cells absent**, and **not random**: timepoint **T16 missing for all 14 sites** (42) + `E02_T17_R3` (1). → **18 usable timepoints, not 19**; treat T16 as *unsurveyed* (not as non-detections) in occupancy histories.

## 3. `conwy_12S_libraries.tsv` — the clean join table (committed, 755 rows)

| column | meaning |
|---|---|
| `ID` | sample key, `E##_R#_T##` (= ENA `sample_alias`; joins to covariates `ID` and to Source Data) |
| `site` / `timepoint` / `rep` | parsed `E01`–`E14` / `T01`–`T19` / `R1`–`R3` |
| `run_accession` | ENA run (`ERR…`); `experiment_accession`, `sample_accession` also kept |
| `read_count`, `base_count` | per-run yields |
| `fastq_ftp`, `fastq_md5`, `fastq_bytes` | `;`-separated R1/R2 URLs, checksums, sizes (used by `download_ena.sh`) |

## 4. `supplementary_data_1_covariates.csv` — per-sample metadata (1141 rows; 755 are Conwy)

*Verified by `data/verify_covariate_join.py` → `covariate_join_verification.md`. The P3 gate.*

> ⚠️ **UTF-8 BOM convention.** The CSV's `ID` header is prefixed with a BOM (`EF BB BF`). A naive `read.csv("…csv")` in R names the first column `ï..ID`, so a join on `ID` silently returns **zero rows**. **Always** read it with `readr::read_csv()` (strips the BOM) or `read.csv(…, fileEncoding = "UTF-8-BOM")`. Better: R scripts consume the derived **`conwy_covariates.tsv`** below, which is already BOM-free.
>
> ⚠️ **Locked covariates with gaps — final local disposition:** `Water_temp_oC` was **100% NA** → substitute `temp_monthly_average`; the spatially resolved flow family is incomplete for E13/E14, while the on-hand `Catchment_flow_m3_s` is only a per-timepoint catchment series. Source E14 climate and gradient values are absent. A separately flagged E13 donor fills the supporting E14 climate table but is unvalidated; it is not an observed E14 measurement. The executed headline occupancy formula excludes climate, flow, and `river_gradient_percent`, using `Dist_from_lake + season`; gradient remains missing rather than being fabricated. See `covariate_join_verification.md` and `R/09_e14_covariates.R`.

- **`conwy_covariates.tsv`** (committed, 755 rows, BOM-free) — `ID` + the decision record's 8 locked covariates + 8 documented alternatives (incl. `temp_monthly_average`). The analysis-ready table R scripts read instead of the raw publisher CSV. Empty cell = true NA.
- The covariate CSV is a **superset**: it covers the full 798-cell Conwy design grid, of which **43 are unsequenced** (= the ENA gate's absent cells), plus other rivers + controls.

Join key `ID`. Modelling-relevant columns (plan §2.2):

| role | columns |
|---|---|
| **Occupancy ψ** (position/habitat) | `Dist_from_lake`, `river_gradient_percent`, `season`, `river_section` (upper/middle/lower), `Latitude`, `Longitude` |
| **Detection p** (eDNA availability/sampling) | `daily_flow_by_catchment`, `week_before_flow_by_catchment` (also monthly/seasonal), `Water_temp_oC`, + per-run read depth (from §3) |
| **Other environment** | `river_width_m`, `pH`, `conductivity_uS_cm`, `gran_alk_ueq_L`, `rainfall_monthly_average`, `Potential_Evapo_Transpiration`, `temp_monthly_average`, `Catchment_flow_m3_s` |
| **Indexing / structure** | `Date_sampled`, `Week`, `Month`, `Year`, `timepoint`, `Replicate`, `River`, `Region`, `SampleSite`, `Blank`, `WP` |

## 5. `source_data_species_x_replicate.xlsx` — published detections (Phase-1 input) — **VERIFIED**

*Verified by `data/inspect_source_data.py` → `source_data_verification.md`. The P1 fast-path gate.*

6 figure-keyed sheets. **The two 12S fish sheets keep replicates distinct → the fast path is VIABLE.**

| Sheet | Rows × Cols | Marker | Resolution |
|---|---|---|---|
| `Figure 1b, 2a&d, 3b&e, 4a&b` | 790 × 68 | **12S fish ← USE** | **REPLICATE (`E##_R#_T##`)** |
| `Figure 5a,d&g, 6a,d&g` | 790 × 68 | 12S fish (2nd copy) | REPLICATE (`E##_R#_T##`) |
| `Figure 1c, 2b&e, 3c&f, 4c&d` | 287 × 3018 | CO1/18S eukaryote | site×time (collapsed) |
| `Figure 5b,e&h, 6b,e&h` | 287 × 3018 | CO1/18S eukaryote | site×time (collapsed) |
| `Figure 1d, 2c&f, 3d&g, 4e&f` | 1395 × 3 | CO1/18S eukaryote | long (site×time) |
| `Figure 5c,f&i, 6c,f&i` | 1395 × 3 | CO1/18S eukaryote | long (site×time) |

- The 12S sheets are **665 Conwy replicate rows × 68 source-study fish species-label columns**, as *normalised within-sample relative abundance* (each row sums to 1). The labels are analysis units rather than independent species validations. The other 125 rows are controls/other rivers (`BLANK*`, `ENC_*`, `GL/TW/SK*`) — excluded by the `E##_R#_T##` regex. The 68-column tally is confirmed exact.
- **Reads ≠ counts here.** The Source Data is normalised, not raw reads — but detection (value > 0 vs 0) is fully preserved (normalising divides by a positive total; it never flips a zero). This is sufficient for the Phase-1 detection-only occupancy models and validation. Count-based Layer E and the Phase-3 community occupancy stress test require the reprocessed counts. Raw read depth, if wanted as a detection covariate, is in §3 (`read_count`).
- **665 = 755 deposited − 90 QC-dropped.** The authors drop libraries with <1000 reads (`Analysing_BLAST_output_WP2_12S.R`). The join to §3 is a clean **665 matched / 0 published-only / 90 dropped**. The 90 dropped concentrate in upstream low-yield sites (E01–E04) and summer low-flow timepoints (T07–T10) — *not random*.
- ⚠️ **The 90 dropped replicates — and T16 — are NA, not 0** in any detection history (sequenced-but-unusable / unsurveyed ≠ non-detection). Same rule as §2. T16 is absent from the published matrix too (18 published timepoints), consistent with the design grid.
- **Derived, analysis-ready table:** `data/metadata/source_data_12S_normalised.tsv` (665 × 68, BOM-free, committed). Phase-1 R scripts consume this, **not** the 5.8 MB git-ignored xlsx.
- The matrix is *pre-cleaning*: it still contains the positive-control spike (`Scomber_scombrus`) and marine/aquaculture strays the authors later remove (their script lists exactly which, lines ~355–368). The Phase-1 occupancy workflow applies that removal.

## 6. Not used (present in archive, out of scope)

Other rivers (`C##`, `G*`, `TW*`, `SK*`), other markers (COI/18S — see `Analysing_BLAST_output_WP2_{CO1,18S}.R`), and controls (`BLANK*`, `EBlank`, `ENC`, `PC`) — the 66 field + 68 lab negatives are revisited for `decontam` in Phase 2.
