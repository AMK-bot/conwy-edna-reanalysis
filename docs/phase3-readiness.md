# Phase-3 Readiness Report — Conwy River 12S eDNA Reanalysis

**Date:** 2026-06-28
**Basis:** synthesis of six deep-reader audits + adversarial verification
**Scope:** Readiness of the WSL2 Phase-2 handback (external `handback/` drive) to start Phase 3
(occumb community occupancy + Layer E compositional baseline + Layer B reclassification) on the Mac.
**Status convention:** only CONFIRMED or PARTIAL findings are treated as established. Refuted claims are
demoted; corrected numbers from the verifiers are used throughout.

> **Current-status note (2026-07-18):** this is a historical hand-back audit, not
> the current release dashboard. The corrupt count table described in preserved
> red-flag sections was replaced and Phase 3 has run. Gate 1 used the documented
> fallback (drop gradient; E13 donor for supporting E14 climate metadata); Gate 2
> added structural NA and depth-aware detection, then failed the recovery claim
> gate, so rare-species direction is withheld. R7 is resolved. Use
> `README.md`, `eDNAplan.md`, and `docs/phase3-claim-status.md` for current claims.

---

## 1. Executive summary & verdict

### VERDICT: GO-WITH-CAVEATS (the two original blockers are cleared; Phase 3 may begin)

**Update 2026-06-28 — both load-bearing blockers from the original NO-GO are now closed.** The original
verdict (preserved below for the record) rested on two confirmed failures: a structurally degenerate
`seqtab_clean.rds` (tag-jump = 1.0 bug) and an unrun concordance gate. Both are resolved:

1. **The count table is fixed (Blocker 0 cleared).** The ASV/species count table was regenerated on the
   WSL2 PC with a sane tag-jump cutoff. The validation gate now runs on a healthy
   `species_x_sample_counts_rawnames.rds` lineage (647 libraries, ~35.7M reads, biologically sensible
   recurring detection incidence) — **not** the degenerate 173,471-read / one-sample-per-ASV `seqtab_clean.rds` that
   triggered the NO-GO. The singleton-detection-history pathology is gone.

2. **The concordance / validation gate is CLOSED (PASS-WITH-CAVEATS).** `R/08_validate_concordance.R` ran,
   was reproduced from scratch and re-run live by four adversarial lenses, and every headline reproduces
   exactly. See `docs/phase2-validation-gate.md`. Of the 16 source-study labels detected in the 642 shared
   libraries, **11 are recovered (68.8%)**; all 5 source-only labels occur in 1–2 libraries and remain
   trace taxonomic disagreements rather than confirmed local species; composition agrees well (mean Bray-Curtis
   similarity 0.84, cosine similarity 0.90, genus recovery 66.7%).

**The caveat carried into Phase 3** was a modelling constraint: reprocessed detection-cell recall is
**~84% of the source-study total over recovered taxa** (selected label-specific differences 23–29%), and
misses co-vary with library depth, source-study relative abundance, and the river gradient. This creates a
recoverability/confounding risk; it does not establish the direction of occupancy bias. Phase 3 tested a
depth-aware detection term and structural NA, but the known-truth recovery MCMC failed convergence and an
18-library comparison mask remained inconsistent for vectors that are zero only within the 14-label
harmonised input. Sparse-label occupancy direction and magnitude are withheld.

**Observed residual QC pattern:** primer-trim retention is low (median ~24.8% kept) and structured along
the river gradient (E01 7.8% → E14 32.6%). This documents target-read loss; without independent biomass or
an orthogonal technical experiment, it does not identify the cause as biology, inhibition, or configuration.
The release preserves the observed table and models the resulting detection risk rather than assigning a cause.

This was **GO-WITH-CAVEATS**, not unconditional GO, because the `occumb` detection model needed to test a
read-depth covariate and structural missingness before any rare ψ estimate could be considered (see §5).
The completed recovery gate did not pass, so those estimates were not promoted. **Update 2026-06-28 — Phase 2 is now fully closed:** the detection offset was confirmed to be
measured on the regenerated ASV table (not a proxy); the precision side was added (cell precision 0.99, no
over-calling); the top 20 ASVs within the 5.9% Unassigned bin were BLASTed and were non-target vertebrate
bycatch, a bounded result that does not exclude the unexamined tail; and the detection design was rebuilt to
the 647 survivors (`R/lib_design.R` `reprocessed` mask + `R/test_reprocessed_design.R`). Current caveats
also include failed recovery convergence and the all-zero comparator-mask mismatch, so no rare-species
occupancy magnitude or direction is retained.

---

<details>
<summary><strong>Superseded original verdict (2026-06-28, pre-fix): NO-GO</strong> — kept for provenance</summary>

Phase 3's single load-bearing input — an ASV x sample integer **count** table that occumb can model — was
**not usable as delivered**, and the fix lived on the WSL2 PC, not the Mac. Two independent, adversarially
**confirmed** failures made this unambiguous:

1. **`seqtab_clean.rds` was structurally degenerate.** A 375 x 829 integer matrix holding only
   **173,471 reads**, in which **every one of the 829 ASVs occurred in exactly one sample**. Every per-ASV
   detection history was a singleton → occupancy/detection parameters non-identifiable by construction.

2. **The table was internally contradicted by its sibling file.** `species_table.rds` carried
   **43,416,098 reads** with sensible occupancy — a **250x** read discrepancy. Cause: tag-jump applied at an
   impossible cutoff of **1.0000**, which zeroed **354 healthy field libraries** and ~575 undisclosed ASVs,
   retaining just **0.35% of field reads**.

The corrective action (regenerate the cleaned count table on the PC, re-run the audit) has since been
completed, which is why this verdict is superseded above. The earlier interpretation of primer-trim
attrition as proven biology is also superseded: only the structured retention pattern is established;
its biological versus technical cause remains unresolved.

</details>

---

## 2. Deliverable inventory vs plan

| Handback file | What it is | Status vs plan |
|---|---|---|
| `seqtab_clean.rds` | Intended ASV x sample integer count table (375 x 829). | **PRESENT but UNUSABLE / corrupt** — degenerate (1 sample per ASV), 173,471 reads, 0.35% of field signal. Must be regenerated on PC. |
| `species_table.rds` | Assigned-label-level (not ASV-level) count matrix, 889 libs x 36 cols (35 assigned labels + Unassigned). 43.4M reads. | **PARTIAL** — usable as a Layer-E count baseline only after stripping 134 control rows; label-collapsed, so it cannot feed ASV-level occupancy. |
| `taxonomy.tsv` | ASV → taxon map, 1,443 ASVs in the broader field+control classification universe; MitoHelper 12S **Nov-2025** build. | **PRESENT (with defects)** — CRLF+CR terminators (needs `dos2unix`); DB-version drift; not a valid 2020 baseline. |
| `repseqs.fasta` | 1,443 representative ASV sequences from the broader field+control universe (line-wrapped). | **PRESENT** — link ASVs via `asv_map.tsv$sequence`, not naive FASTA line matching. |
| `asv_map.tsv` | ASV id ↔ sequence ↔ read totals (1,443-ASV field+control universe). | **PRESENT** — authoritative ASV linkage; 829 seqtab seqs are a clean subset (829/829). |
| `controls_audit.md` | Control role accounting. | **PRESENT (defects)** — negative set stated 79; role total exclude=522 not reproducible from file (file shows 42); off-by issues. No adversarial verification pass was run on this dimension. |
| `decontam_audit.md` | decontam + tag-jump audit. | **PRESENT but MISLEADING** — frames loss as "removed 134 controls → 375" and conceals 380 field-library loss; applied tag-jump cutoff = 1.0 (bug). Negative set stated 77. |
| `dada2_track.tsv` | Per-library read tracking, 889 rows. | **PRESENT** — `reads_postchimera`/`reads_final` columns are monotone & consistent (the "chimera bug" was **REFUTED**); but neither column equals seqtab rowSums, so do not seed counts from it. |
| `trimming_summary.tsv` | Primer-trim accounting, 892 rows. | **PRESENT** — documents structured retention; biological versus technical cause remains unresolved. |
| `qc_report.html` | DADA2 `plotQualityProfile` substitute + worst-15 trim table. | **PARTIAL** — declared PV4 swap; not FastQC/MultiQC; too thin to localize per-sample read loss. |
| `validation_concordance.md` | Concordance reproduction note. | **PARTIAL / PENDING** — concordance not run (published table not staged on PC) AND ran against 2026 DB, so not a valid 2020 gate. |
| 3-lens adversarial review (statistical/code/fidelity) | Plan Phase-2 closing requirement. | **MISSING** — not delivered in handback. |
| Documented 2020-baseline concordance reproduction | Plan Phase-2 checkpoint requirement. | **MISSING** — see above. |

---

## 3. Quality red-flags (severity-ordered)

### R1 — CRITICAL: `seqtab_clean.rds` is degenerate; every ASV occurs in exactly one sample
**Verdict: CONFIRMED.** 375 x 829 integer matrix; `sum(st>0)=829`; `table(colSums(st>0))` returns only the
value 1 (count 829); rows hold 1–11 ASVs (median 2). Every detection history is a singleton → occupancy and
detection probability are non-identifiable. **This is not a depth problem; the table is structurally broken.**
*Corrected phrasing:* the degeneracy is **one sample per ASV**, not "one ASV per sample" (samples carry
median 2, range 1–11 ASVs).
**Action:** Do not run occumb on this object. Regenerate the true cleaned ASV x sample count matrix on the PC.

### R2 — CRITICAL: tag-jump cutoff = 1.0000 is a calibration bug that destroyed the data
**Verdict: CONFIRMED.** A cross-talk cutoff of 1.0 = "treat 100% of reads as index-hopping" is degenerate;
real Illumina rates are ~0.01–2%. The audit's own sweep shows sane behaviour (0.001 → 77.5% retained,
0.01 → 37.65%); 1.0 retains 905/13,013 cells = 6.95%. It zeroed **354 field libraries that had real reads**
(median 41,828, max 608,745) — only 26 of the 380 dropped were truly empty — retaining **173,471 of
50,178,127 field reads = 0.346%**. *Caveat (uncertain):* the precise root-cause mechanism (a `max()` over a
proportion including the positive control's own mackerel/mackerel = 1.0 cell) is a plausible hypothesis only;
the calibration code is on the PC and not in the handback.
**Action:** Re-run tag-jump calibration with the positive-control's own cell excluded; sanity-check against
the sweep. Expect to recover toward ~729 field libraries with >0 post-chimera reads.

### R3 — CRITICAL/HIGH: ~250x read discrepancy + undisclosed ASV collapse → handback is internally inconsistent
**Verdict: CONFIRMED.** `species_table.rds` = 43,416,098 reads vs `seqtab_clean.rds` = 173,471 reads
(**ratio 250.3x**). ASV universe collapsed 1443 → 829 (lost 614), but the audit discloses only 39 removals
(13 contaminant + 26 mackerel); **575 ASVs are silently lost** (all-zero columns dropped after the 1.0
tag-jump zeroed them). The `decontam_audit.md` "removed 134 control rows → 375 field samples" framing
conceals the 380 field-library loss entirely.
**Action:** Treat `seqtab_clean.rds` as the wrong/corrupt export. Re-export and re-audit on the PC; track
field-library survival **separately** from control removal in the regenerated audit.

### R4 — HIGH: sample attrition 375/755 (49.7%), spatially structured upstream loss
**Verdict: CONFIRMED.** Only 375/755 field libraries present; authors retained 665/755 (88.1%). Loss is
spatially structured: E01=4 reps, E02=6, E04=7, E03=11 vs E13=44, E14=43. Of 252 (site,timepoint) cells,
**70 are empty, only 68 retain a full triplicate, 184/252 have <3 reps**; E01 has detections at only 3 of 18
timepoints. **Note:** most of this 375-number is an artifact of the R2 tag-jump bug. The residual trim-level
loss is structured but its biological versus technical cause is unresolved (R6).
**Action:** Do not assume the nominal 14x18x3 balanced design. Rebuild the realised design grid from the
**corrected** table; budget for unmodellable upstream cells and a read-depth detection covariate.

### R5 — HIGH historical concern: only 37 assigned taxonomy labels recovered
**Current verdict: the original biological interpretation is REFUTED.** The 37 assignments and the
publisher's 68-label universe are not a valid recall numerator and denominator. The fish previously called
"genuinely present but missing" have zero detections in both tables over the shared comparison and therefore
cannot support local-presence or recovery claims. The validated denominator is the 16 source-study taxa
actually detected in 642 shared libraries: 11 are recovered, and five one-to-two-library labels remain trace
review cases. Exact-name counts are also distorted by database namespace differences.

### R6 — HIGH (re-scoped): primer-trim retention is structured; cause unresolved
**Verdict: the structure is CONFIRMED; a biological or technical cause is not identified.** Field median pct_kept = 24.8–26.6%;
catastrophic single-digit cases exist (E01_R1_T08 219621→7) but cluster upstream (E04 frac<100 = 0.15) and
in summer (T07/T09/T10). pct_kept rises monotonically (broadly) along the river gradient E01 7.8% → E14
32.6%; `cor(pct_kept, log10 reads_in) = 0.178`. Site/season structure alone does not distinguish biology,
inhibition, primer mismatch, or a configuration interaction. The validated local release retains this as a
detection/QC limitation. A future full raw-read rerun should include a cutadapt configuration sensitivity;
the current analysis uses structural missingness and does not convert low-retention libraries into claimed
biological absences.

### R7 — RESOLVED 2026-07-18: 79 manifest negatives; 77 decontam inputs
**Historical disagreement:** The decontam reader argued
**77** is correct (54 T1 [ENC 45 + EBlank 9] + 23 BLANK; well-blanks = 12). The controls reader argues
**79** is correct (54 T1 + 25 T2; well-blanks = 14) — but that dimension had **no adversarial verification
pass**. The entire disagreement is the well-blank count (`BLANK_[CD]#_well`: 12 vs 14). This is a
**transparency/reproducibility defect, not a detection-bias risk** (both readers agree the negative set is
clean of cross-river/other-study contamination). The controls reader also flags `exclude=522` in
`controls_audit.md` as not reproducible from the file (file shows 42; 522 = 480 never-trimmed off-target +
42 T3 LB, a manifest-scope number).
**Resolution:** both counts were correct at different stages. The manifest contains 14 well blanks; two
(`BLANK_D1_A1`, `BLANK_C3_A1`) have no post-DADA2 ASV-table row, so 12 well blanks and 77 total negatives
enter decontam. See `data/metadata/control_reconciliation.md`.

### R8 — MEDIUM: DB-version drift — *Salmo* species label unresolved
**Historical verdict: PARTIAL; current disposition: resolved by claim demotion.** The Nov-2025 build assigns
283 ASVs to the label *Salmo macrostigma*, including ASV0001 (18.5M reads). Only **2 of 283** have
pident=100, so a universal exact-tie mechanism was overstated. No deterministic
`macrostigma→trutta` relabel was staged. Layer B instead compared current resources and a direct witness
set: genus *Salmo* is stable, while the 170-bp species label is not. The marker does not establish a local
or native species identity, and no replacement species label is treated as ground truth. The unarchived
Sep-2020 MitoFish build prevents a faithful historical replay.

### R9 — MEDIUM: two high-abundance no-hit ASVs + 19 hit-but-Unassigned ASVs
**Verdict: CONFIRMED.** ASV0007 (1,016,158 reads) and ASV0008 (887,798) are Unassigned with no BLAST hit =
**4.35%** of all 43.8M reads and dominate the Unassigned pool. Separately, 19 Unassigned ASVs *do* have a hit
(some at pident=100, e.g. ASV0016/AB972099, ASV1181) with no documented exclusion reason.
**Action:** BLAST the two no-hit ASVs against full GenBank nt before deciding whether to keep an "Unassigned"
column in any model input; document/repair the 19 hit-but-Unassigned cases.

### R10 — MEDIUM: `species_table.rds` still holds 134 control rows; Scomber spike absent
**Verdict: CONFIRMED.** 889 rows = 755 field + 134 controls (23 BLANK, 9 EBlank, 45 ENC, 16 PC, ~41 LB*) →
**pre control-removal**. 95 libraries have depth 0, 164 < 100 reads. **Scomber_scombrus** (mackerel positive
spike) is absent from the 36 assigned-label columns — either already removed, unassigned, or folded into "Unassigned".
**Action:** Strip the 134 control rows before any Layer-E baseline; reconcile the missing Scomber control;
decide handling for zero/low-depth libraries.

### R11 — MEDIUM: tag-jump + decontam reduce the matrix to near-noise sparsity
**Verdict: CONFIRMED.** Even among the 375 survivors, depth is min 1 / median 202 / max 6,287; 10th
percentile = 13 reads. Rows with 1–13 reads carry no detection power. (Subsumed by R2 — fixing the tag-jump
is the prerequisite.)
**Action:** Post-fix, impose a minimum read-depth filter for detection-history construction.

### R12 — LOW: data-completeness notes
**Verdict: CONFIRMED.** One field library genuinely missing: **E02_R3_T17** (755 = 756 − 1). Three libraries
have NA trim records: **E05_R3_T18, ENC_R1_T19, PC_C4_H6**. taxonomy.tsv has **CRLF+CR** terminators (run
`dos2unix` or strip `\r` before any parse, or 557 Unassigned ASVs silently mis-bin as assigned).
Lamprey/European flounder are a shared 12S-MiFish blind spot (absent from both tables) — report occupancy as
a *partial* community.

### Demoted / refuted (do not act on)
- **"dada2_track chimera column is buggy/misaligned"** — **REFUTED.** `reads_final ≤ reads_postchimera` in
  889/889 rows (zero reversals); the column is monotone and consistent with a legitimate downstream
  control-blanking filter (52 of 57 zero-final rows are controls). *But* neither column equals seqtab
  rowSums, so still do not seed Phase-3 counts from this file.
- **"primer-trim is a global cutadapt bug requiring Phase-2 rework"** — **REFUTED** (see R6).
- **"06_occupancy_occumb.R would code lost libs as observed-zero as written"** — **PARTIAL/overstated;** 06 is
  a bare `stop()` stub mentioning `observable_mask` only in a comment. The false-zero risk is a *latent*
  hazard once implemented, not current behaviour.

---

## 4. Reconciled key numbers (verified)

| Quantity | Reconciled value |
|---|---|
| seqtab_clean.rds dims / reads | 375 x 829 integer; **173,471 reads** |
| seqtab non-zero structure | every ASV (829/829) in exactly 1 sample; samples hold median 2 (range 1–11) ASVs |
| seqtab sparsity / depth | 0.9973; per-sample min 1, median 202, max 6,287; 10th pct = 13 reads |
| species_table.rds dims / reads | 889 x 36 (35 assigned labels + Unassigned); **43,416,098 reads** |
| Read discrepancy (species_table : seqtab) | **250.3x** |
| Field libraries surviving in seqtab | **375 / 755 (49.7%)**; authors retained 665/755 (88.1%) |
| Field libraries dropped | 380; 26 truly empty, **354 had reads** (median 41,828, max 608,745) |
| Field reads surviving into seqtab | **173,471 / 50,178,127 = 0.346%** |
| Tag-jump cutoff applied | **1.0000 (bug)**; retained 905/13,013 cells = 6.95% (sweep: 0.001→77.5%, 0.01→37.65%) |
| ASV universe | 1443 → 829 (lost 614; only 39 disclosed; **575 undisclosed**) |
| decontam prevalence | P*=0.1 → 13 contaminant ASVs (reasonable) |
| Negative-set size | **RESOLVED:** 79 manifest negatives; 77 decontam inputs because 2 zero-read well blanks have no ASV-table row |
| Taxonomy labels in broader field+control universe | **37 distinct**; 886/1,443 ASVs assigned (61.4%); read-weighted 93.3% |
| Unassigned ASVs | 557 (38.6%); 538 no hit, 19 with hit |
| Top no-hit ASVs | ASV0007 1,016,158 + ASV0008 887,798 = **4.35%** of reads |
| Salmo_macrostigma | 283 ASVs (only 2 at pident=100), ASV0001=18.5M reads; **uncorrected, no relabel staged** |
| Primer-trim field median pct_kept | 24.8–26.6%; gradient E01 7.8% → E14 32.6%; cor(pct_kept, log10 reads_in)=0.178 |
| Design grid (252 cells) | 70 empty, 57 with 1 rep, 57 with 2, **68 full triplicate** |
| Missing field library | **E02_R3_T17** |
| NA trim libraries | E05_R3_T18, ENC_R1_T19, PC_C4_H6 |
| E14 covariates | river_gradient_percent / temp_monthly_average / rainfall_monthly_average **blank in 54/54 rows**; Catchment_flow is an 18-value time series (2.53…57.64), **not constant 2.53** |
| taxonomy.tsv encoding | CRLF + CR (needs dos2unix) |
| Concordance file (Mac) | `data/metadata/source_data_12S_normalised.tsv`, 665 rows x 68 source-study label columns; **665/665 IDs join** to species_table |

---

## 5. Phase-3 entry gates & blockers

**BLOCKER 0 — Count table regenerated on the WSL2 PC. ✅ CLEARED.**
The degenerate `seqtab_clean.rds` has been replaced by a regenerated count table with a sane tag-jump cutoff;
the validation gate now runs on the healthy 647-library / ~35.7M-read lineage with sensible recurring
occupancy. The singleton-history pathology (R1–R3) is resolved. occumb may now be fit on the rebuilt design.
Primer trimming was left unchanged, but neither trimming nor a biological mechanism is established as the
cause of the remaining discordance. The detection offset was re-measured on this regenerated lineage and is
already reflected in the Phase-3 validation; it is a recoverability warning, not a ground-truth error rate.

**Gate 1 — E14 covariates (RESOLVED IN PHASE 3).**
`river_gradient_percent`, `temp_monthly_average`, and `rainfall_monthly_average` are empty for E14 in the
source covariates. The documented fallback was used: drop `river_gradient_percent`, keep
`Dist_from_lake + season`, and use an E13 donor only for supporting temperature/rainfall metadata.
Catchment-wide flow was not added to the headline spatial formula.

**Gate 2 — Detection / occasions design (IMPLEMENTED; FINAL CLAIM GATE FAILED).**
`lib_design.R::observable_mask("reprocessed")` is keyed to the 647 regenerated count-table rows and equality
tested, so unavailable cells are structural missingness rather than false zeros. The model includes read
depth on detection and supports a minimum-depth filter. These safeguards address the design requirement but
do not guarantee recovery. The known-truth MCMC failed the Rhat threshold, and 18 observed libraries that
are zero within the 14-label harmonised Layer-A input exposed a comparator-mask inconsistency; those
libraries contain reads under raw taxonomy labels. Sparse-label ψ direction and magnitude are therefore
withheld rather than described as corrected.

**Gate 3 — Concordance reproduction. ✅ DONE (PASS-WITH-CAVEATS).**
Run via `R/08_validate_concordance.R` against the staged published `source_data_12S_normalised.tsv` (665 x 68)
and the regenerated count table; bridged across the two DB vintages with the consensus crosswalk
(`data/derived/crosswalk.tsv`); D3 applied to both namespaces. Reproduced exactly by four adversarial lenses.
**642 shared reps; species recovery 68.8% (11/16; all 5 misses are one-to-two-library detections and the
one-library *Leuciscus* alternative remains unresolved);
genus recovery 66.7%; Bray-Curtis similarity 0.835/0.864; cosine similarity 0.903; Spearman(abund) ~1.00 (0.996).** Full report:
`docs/phase2-validation-gate.md`. Note this is a **2026-DB concordance**, not a valid 2020 baseline (the
Sep-2020 MitoFish build is unarchived, PV5). **Due diligence DONE (2026-06-28):** the precision pass was
added (cell precision 0.99, F1 0.91 — little over-calling) and the top 20 ASVs within the 5.9% Unassigned
bin were BLASTed against NCBI nt. Those leading ASVs are non-target vertebrate bycatch; this characterises
the dominant bin content but does not prove absence from its unexamined tail. Phase-2 validation is
**closed (PASS-WITH-CAVEATS)**.

**Gate 4 — Layer B historical-baseline decision. ✅ RESOLVED.**
The Sep-2020 MitoFish build is unarchived, so a faithful historical replay is not claimed. The executed
Layer B is current-resource concordance: a same-reference classifier comparison, an independent MIDORI2
GenBank271 comparison, and a direct current public-record witness set. It reports genus-level stability and
species-level uncertainty rather than a corrected local identity or temporal novelty.

**Historical non-blocking list — final disposition:** R7 is reconciled (79 manifest / 77 decontam); D3 and
the many-to-one crosswalk are explicit and read-conserving; no deterministic *Salmo* replacement label is
used; the leading Unassigned ASVs are characterised with a bounded conclusion; Phase-3 claim status is
documented in `docs/phase3-claim-status.md`.

---

## 6. Historical Phase-3 first-task list — final disposition

1. **Count table:** regenerated after the tag-jump error; current table is 647 libraries × 1,348 ASVs / 35.7 M reads.
2. **Denominators/provenance:** ASV, control, and hand-back counts are reconciled; the 79-manifest/77-decontam boundary is explicit.
3. **Taxonomy parsing:** CRLF/Unassigned handling is checked by the validation scripts.
4. **Concordance:** executed over 642 shared libraries; status is PASS-WITH-CAVEATS.
5. **Purported missing well-detected labels:** the original interpretation was refuted; those labels are zero in both comparison tables.
6. **Detection design:** `observable_mask("reprocessed")` and structural missingness are implemented and tested.
7. **Cleaning and assignment:** source-study bookkeeping is reproduced for concordance, while Layer B keeps *Salmo* species unresolved rather than forcing a replacement label.
8. **E14:** source gaps remain explicit; climate, flow, and gradient are excluded from the executed headline formula.
9. **Layer E:** implemented with site/time aggregation and restricted-series inference.
10. **Layer A:** implemented as a stress test; failed recovery convergence and a comparator-mask mismatch trigger numerical suppression.
11. **Layer B:** implemented as current-resource concordance because the Sep-2020 MitoFish build is unarchived.

---

*Bottom line (updated 2026-07-18): the PC regenerated a healthy count table and the concordance gate is
closed PASS-WITH-CAVEATS. Across recovered taxa and 642 shared libraries, pooled detection-cell recall is
~84% and is structured by depth and source-study relative abundance. Phase 3 tested structural-NA coding,
a depth-aware detection model, and known-truth recovery; the recovery MCMC failed convergence and an
18-library harmonised-input-zero comparison-mask mismatch remained. Consequently no sparse-label correction magnitude or direction
is retained. Layer E uses site/time-aware inference, Layer B reports genus stability/species uncertainty,
and the portfolio edition is independent of production MCMC.*
