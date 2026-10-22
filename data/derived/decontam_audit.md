# Phase-2 decontam audit

## decontam prevalence
- threshold 0.1 (reported): **13** contaminant ASVs (0.9% of ASVs, 0.66% of reads).
- threshold 0.5 (sensitivity): 33 contaminant ASVs.
- negative set used by decontam = **77** present in this seqtab, of **79** in the controls manifest (2 absent: BLANK_D1_A1, BLANK_C3_A1 -- 0 reads / not in seqtab). PC positive is never a negative.

## Tag-jump (index-hopping) filter
- Scomber-assigned ASVs: 26; of these, present in >=1 PC (positive) well: **0**.
- **positive control carries no usable mackerel spike** (PC reads in Scomber ASVs = 0); cutoff CANNOT be spike-calibrated. Applied the pre-registered default per-ASV index-hop floor **0.0010** (0.1%).
- applied tag-jump cutoff: **0.0010** -> retains **77.54%** of raw detections.

| cutoff | nonzero cells | detections retained |
|---|---|---|
| 0.0000 | 13013 | 100.00% |
| 0.0005 | 11123 | 85.48% |
| 0.0010 | 10090 | 77.54% |
| 0.0050 | 6422 | 49.35% |
| 0.0100 | 4900 | 37.65% |

## Final clean table
- ASV columns dropped: 13 decontam contaminant + 26 Scomber scombrus (mackerel) ASV(s).
  Scomber = the marine spike *target* species, dropped per prereg. NB: the spike produced 0 reads
  in its own (empty) PC wells, so these reads are spike-derived lab cross-talk (e.g. ASV0011 is 100%
  mackerel yet ~99% upstream/freshwater and absent at the well-sequenced marine mouth -- the inverse
  of genuine eDNA). Dropping them is conservative and matches the authors.

### Field-library survival (tracked separately from control-row removal)
- field libraries in seqtab: **755**; non-empty before cleaning: **724**.
- **survive into seqtab_clean.rds: 647 / 755 field libraries (85.7%)** [authors' target 665/755].
- field libraries dropped by cleaning: 108 (of which 31 were already empty pre-clean).
- control/non-sample rows removed from the seqtab axis: 134 (these are NOT field-library losses).
- NOTE: 647 survival reflects the applied 0.001 floor. The brief's '~729' target corresponds to NO
  tag-jump (cutoff 0 -> 724 survive); applying any index-hop floor necessarily trims low-signal libraries.
  0.001 is the pre-registered sweep value; 647 sits within ~3% of the authors' 665.

### Count integrity
- **seqtab_clean.rds: 647 field samples x 1348 ASVs, 35,716,239 integer reads.**
- samples per ASV: median 1, max 318 (gate requires NOT every ASV in <=1 sample -> PASS).
- ASVs per sample: median 12, max 69; sparsity 0.9889.

### Species assignment on the corrected table (NOT a concordance -- published 68-matrix is Mac-side)
- distinct ASSIGNED species among surviving ASVs: **33** (D3 relabels applied; was 37 ASV-level pre-fix).
- **Largest assigned taxon: Salmo_macrostigma (45.9% of clean reads).** Salmo_macrostigma (Adriatic/Med trout, cannot occur in N. Wales) is a
  DATABASE NEAREST-NEIGHBOUR ARTIFACT = local Salmo (trutta/salar). The macrostigma->trutta relabel is a
  Mac-side D3 step (per PHASE2_ISSUES_FOR_PC.md P4), NOT a PC fix -- flagged here so the count is not mistaken
  for a finding. species_table.rds (43.4M reads) is the PRE-decontam artifact; seqtab_clean.rds is post-decontam.
- sample axis keyed by ID (E##_R#_T##); ASV-level (taxonomy.tsv ships the species map).
