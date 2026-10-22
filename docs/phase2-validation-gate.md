# Phase-2 Validation Gate — Conwy River 12S eDNA Reanalysis

**Date:** 2026-06-28
**Basis:** This note synthesises the R/08 concordance run and four verification lenses: recompute, crosswalk audit, detection offset, and honesty.
**Question under test:** Does the reprocessed pipeline (MitoHelper 12S, Nov-2025 build) recover the
detections published by Perry et al. 2024 (68-species normalised relative-abundance table)?
**Script under test:** `R/08_validate_concordance.R`
**Status convention:** only CONFIRMED or PARTIAL verifier findings are treated as established; corrected
values from the verifiers are used throughout. Where a headline was rounded, the unrounded value is given.

---

## 1. Verdict

### VERDICT: PASS-WITH-CAVEATS

The reprocessing reproduces the dominant/common community composition with caveats, but has materially
lower pooled detection-cell agreement across recovered taxa; that shortfall — not the recovery percentage — is the load-bearing result for
Phase 3. Every headline number reproduced exactly under an independent rewrite (all 17 figures), the
crosswalk that bridges the two database vintages is arithmetically sound (reads conserved, no
double-counting, every one of the 35 reprocessed raw columns handled exactly once), and the bounded recovery result holds:
of the 16 source-study taxa actually detected in the 642 shared replicates, the comparison splits into
**11 recovered** and **5 missed, all one-to-two-library detections**. Four are marine/estuarine or
trans-Atlantic nearest-neighbour candidates; the one-library *Leuciscus* label remains unresolved. Every
taxon the authors saw in ≥3 replicates was recovered (10/10). The recovery figure of 68.8% is therefore
**conservative**, not spun. The one caveat that survives into Phase 3 is **(i) a systematic ~16% pooled
detection-cell recall deficit (23–29% for the *Salmo salar*, *Anguilla*, and *Barbatula* labels,
which have intermediate source-study relative abundance)** that the
read-weighted agreement metrics structurally mask because they condition on co-detection — a recall
(sensitivity) shortfall, not a precision problem. The three caveats the first draft listed are now
**closed** (see §6, resolved 2026-06-28): (ii) the precision/false-positive side was measured —
**cell-level precision 0.99** (14 reprocessed-only versus 382 source-only cells descriptively shows the
asymmetry; these repeated taxon×library cells are not treated as independent biological inference); (iii) the **top 20 ASVs within the 5.9% Unassigned bin
were BLASTed against NCBI nt** and were non-target vertebrate bycatch (livestock 4.8%,
wildlife/amphibian/bird 0.4%, one marine flatfish stray 0.4%). This supports dominance by non-targets but
does not establish absence from the unexamined tail; and (iv) the offset
was in fact measured on the regenerated ASV count table (`seqtab_clean.rds` → species), **not** the
pre-decontam species-table — so it already reflects the true Phase-3 input.

---

## 2. Method — why a naive name-join is invalid, and the crosswalk that replaces it

The two database vintages can assign **different nearest-neighbour names to compatible short-marker sequence evidence**. The Nov-2025 build
assigns the dominant Conwy *Salmo* sequence the label `Salmo_macrostigma` and assigns Conwy minnow sequences `Phoxinus_phoxinus` /
`Phoxinus_sp._Morella` / `Phoxinus_bigerri`, whereas the source-study table names the corresponding minnow signal
`Phoxinus_ujmonensis`. An exact-string join across these namespaces would manufacture false discordance
(every relabelled sequence signal would appear as both a source-only miss and a reprocessed-only extra) and would understate
recovery to near zero. Naive name-overlap between the raw tables is only 15 taxa and is meaningless.

The gate instead bridges namespaces with a **consensus crosswalk** (built through three independent audit passes plus
reconciliation; persisted to `data/derived/crosswalk.tsv`) and applies the authors' own D3 cleaning to both
sides:

- **Authors side (D3 cleaning, 68 → 55):** drop 10 marine/freshwater strays + the mackerel spike
  (`DROP_A` lists 11, but `Favonigobius_gymnauchen` is absent from the published table, so only 10 fire),
  then 3 relabels (`Salmo_obtusirostris`→`Salmo_trutta`, `Anguilla_rostrata`→`Anguilla_anguilla`,
  `Salvelinus_fontinalis_x_S.malma`→`Salvelinus_malma`) which **add** abundance into the existing target
  column (not overwrite). Net 68 − 10 strays − 3 merged = 55.
- **Reprocessed side (crosswalk, 35 → 14):** 21 raw columns map (many-to-one) onto 14 source-study namespace targets
  (e.g. 3 Phoxinus → `Phoxinus_ujmonensis`; `Salmo_macrostigma` + `Salmo_trutta` → `Salmo_trutta`;
  `Barbatula_aff._barbatula_lineage_I` + `Barbatula_barbatula` → `Barbatula_barbatula`); the remaining 14
  columns (marine strays + the Unassigned bin) are dropped.
- **Conservation gates (all pass):** reads are perfectly conserved through the many-to-one collapse
  (`sum(Mk) == sum(Mh)`; mapped + dropped = grand total); the collapse is a disjoint function over columns
  so no double-counting is possible; and `setequal(colnames(M), c(names(MAP), DROP_M))` enforces that every
  one of the 35 reprocessed columns is handled exactly once.
- **Comparison space:** tables are intersected to the **642 shared replicates** (E## sample-ID namespace —
  the E##-vs-C## gotcha is not in play); detection incidence is counted as values > 0 per library; read-weighted
  metrics renormalise each row to 1 before comparison (so only **compositional**, not absolute-magnitude,
  agreement is tested — appropriate, since the reprocessed side is counts and the source-study side is relative abundance).

**Robustness notes (non-blocking, from the crosswalk audit):** the 55-taxa gate passes by *coincidence* —
`DROP_A` names 11 source-study labels but only 10 exist, and the only guardrail is the hardcoded `==55`; the same latent
risk exists for `DROP_M` (matched by membership, no presence assertion). A future table/DB revision could
silently no-op a drop. Suggest adding `stopifnot(sum(DROP_A %in% colnames(A)) == 10)`. Separately, the
recompute lens read the crosswalk *from disk* (an R/08 output), so it independently verifies the downstream
arithmetic but takes the **mapping decisions themselves as given** — the biology of the crosswalk was not
re-derived in verification.

---

## 3. Results

All values below reproduced exactly under an independent rewrite (recompute lens) and a live re-run
(crosswalk-audit lens). Rounded headlines are shown with their unrounded value where it matters.

### 3.1 Recovery and set agreement

| Metric | Value | Note |
|---|---|---|
| Shared replicates | **642** | reprocessed 647 ∩ source-study 665 |
| Source-study cleaned taxa | **55** | 68 → 55 via D3 |
| Reprocessed harmonised labels | **14** | 35 raw-collapse columns (34 assigned labels + `Unassigned`) → 14 via crosswalk |
| Shared analysis labels | **14** | all 14 reprocessed harmonised labels occur in the source-study 55-label set |
| Source-study taxa **detected** in 642 libraries | **16** | the correct recall denominator |
| Species recovery, permissive | **68.8%** (11/16, =0.6875) | recall vs realised detection set |
| Species recovery, strict (≥2 reads) | **68.8%** (=0.6875) | identical → present-calls are not single-read noise |
| Jaccard (detected sets) | **0.58** (=0.5789, 11/19) | both 11, source-study-only 5, reprocessed-only 3 |
| Genus recovery | **66.7%** (=0.6667, 10/15) | the decisive 12S pass |

**Recovery framing (important):** 68.8% = 11 of the **16** authors taxa actually present in the shared
replicates — *not* 11/55 (=20%) of the cleaned universe nor 11/68 of the published table. The denominator
of 16 is the methodologically correct recall denominator (computing against the 55 would be wrong, since 39
of them are "neither"). Because all 5 misses sit at occ_authors 1–2 and every ≥3-rep taxon is recovered,
the aggregate percentage hides incidence: all source-study taxa detected in ≥3 libraries were recovered.

### 3.2 Read-weighted agreement (14 shared analysis labels, per-replicate, renormalised)

| Metric | Value | Note |
|---|---|---|
| Bray-Curtis similarity, mean / median | **0.835 / 0.864** | per-replicate, shared labels |
| Cosine similarity, mean | **0.903** | |
| Spearman of per-taxon total rel. abundance | **~1.00 (=0.9956)** | NOT exactly 1; over 14 taxa dominated by ~6 labels — near-tautological, report as ~1.00 |
| Genus-level Bray-Curtis similarity, mean | **0.855** | |

**Caveat:** these metrics are conditioned on co-detection within the 14 shared analysis labels and therefore
**cannot see** the detection-cell recall deficit in §5. They evidence that *the same shared labels dominate in the same
proportions*, not that detection sensitivity matches. They should not be the headline.

### 3.3 The recovered set — BOTH (11 taxa; source-study / reprocessed detection counts over 642 libraries)

| Taxon | Source study | Reprocessed | Reprocessed/source ratio |
|---|---|---|---|
| Salmo_trutta | 642 | 604 | 0.941 |
| Phoxinus_ujmonensis | 454 | 408 | 0.899 |
| Salmo_salar | 400 | 310 | 0.775 |
| Barbatula_barbatula | 381 | 270 | 0.709 |
| Anguilla_anguilla | 251 | 179 | 0.713 |
| Gasterosteus_aculeatus | 142 | 122 | 0.859 |
| Oncorhynchus_mykiss_x_Salmo_salar | 37 | 39 | 1.054 |
| Pungitius_pungitius | 26 | 27 | 1.038 |
| Cottus_perifretum | 6 | 3 | 0.500 |
| Salvelinus_malma | 4 | 4 | 1.000 |
| Perca_fluviatilis | 1 | 3 | 3.000 |

---

## 4. Divergence map interpreted

### 4.1 Source-study-only (5 labels — all trace, 1–2 libraries)

Every source-study-only label occurs in one or two libraries. Several are marine, estuarine, or
trans-Atlantic nearest-neighbour candidates; none establishes a Conwy population. The
source-study data contain **no** grayling, pike, gudgeon, chub, or tench detections (all are
"neither", 0 libraries in both tables). The five source-only labels are trace observations and do not
support claims about established Conwy populations.

| Taxon | Reps | Classification | Reasoning |
|---|---|---|---|
| Myoxocephalus_scorpius | 2 | **trace marine-associated label; review** | the short-marker label and two-library detection do not establish organism origin or a local population |
| Chelon_labrosus | 1 | **trace estuarine-associated label; review** | the one-library detection does not establish organism origin or a local population |
| Alosa_sapidissima | 1 | **geographic/taxonomic review candidate** | a trans-Atlantic label from one library; local identity and alternative assignments remain unresolved |
| Osmerus_mordax | 1 | **geographic/taxonomic review candidate** | a North American species label from one library; local identity and alternative assignments remain unresolved |
| Leuciscus_burdigalensis | 1 | **geographic/taxonomic review candidate** | no Leuciscus/Squalius/dace label occurs on the reprocessed side; a single short-marker observation cannot resolve the alternatives |

**The one residual the gate cannot fully foreclose:** 12S barely separates *Leuciscus* spp., native UK
*L. leuciscus* (dace) exists, and a one-library dace detection would be indistinguishable from a
reference-label mismatch at this marker. It could in principle sit in the unexamined tail of the 5.9%
Unassigned bin, so the alternatives remain unresolved.

### 4.2 Reprocessed-only (3 extras — all trace)

The gate applies the same skepticism to reprocessed-only labels. All three are short-marker review
candidates, not established local records. Notably, **Abramis brama occurs in 5 libraries**, more detection cells than
any source-only label (maximum 2), so it cannot be dismissed solely because it is reprocessed-only.

| Taxon | Reps | Classification | Reasoning |
|---|---|---|---|
| Abramis_brama | 5 | **trace label; local population not established** | highest detection count of any divergent label, but the current audit cannot establish local identity |
| Barbus_barbus | 1 | **trace label; identity unresolved** | one-library detection; no replacement identity is inferred |
| Perca_flavescens | 1 | **geographic/taxonomic review candidate** | North American species label from one library; kept as a distinct observed label, with no inferred local identity or replacement |

### 4.3 Join exclusions (verified harmless)

The 642-library intersection excludes 23 source-study libraries and 5 reprocessed libraries. Checked: the 23 excluded
authors reps contain only BOTH-set taxa (Salmo_trutta, Barbatula, Phoxinus, Salmo_salar) plus the dropped
marine stray *Scomber_scombrus*. No non-trace taxon is hidden by the exclusion, so recall is not biased — but
the script asserts the count without reporting what was dropped (a transparency gap).

---

## 5. The systematic detection-offset finding and its Phase-3 implications

**Finding (CONFIRMED):** per-taxon detection counts across shared libraries are systematically **lower**
on the reprocessed side for mid/high-incidence taxa. This is the load-bearing validation result, and it is
masked by the reassuring read-weighted metrics in §3.2.

**Quantified shortfall (use the pooled figure, not the mean ratio):**

- **Pooled detection-count ratio = Σ(reprocessed)/Σ(source study) over the 11 recovered taxa = 0.840 →
  16.0% net under-detection.**
- Concentrated, not flat: **Salmo_salar −22.5%, Anguilla −28.7%, Barbatula −29.1%**; Phoxinus −10%, even
  Salmo_trutta −6%. Trace taxa match or exceed the authors (Oncorhynchus 1.05, Pungitius 1.04, Perca 3.0).
- **382 source-positive/reprocessed-negative detection cells** across the recovered label set (Barbatula 112, Salmo_salar 90,
  Anguilla 72).
- **Do NOT report the "mean per-taxon ratio" (1.044)** — it is inflated by Perca_fluviatilis 1→3 (ratio 3.0)
  and reads as *no shortfall*. Median ratio is 0.899.

**Cause — both originally-proposed causes are REFUTED:**

- *Not* a minimum-read detection threshold: the validation input has minimum nonzero cell = **1 read**;
  every positive cell is counted as a detection.
- *Not* the superseded tag-jump = 1.0 table: the concordance uses the regenerated 35.7 M-read
  `seqtab_clean.rds` lineage after the 0.001 tag-jump calibration. The local `R/03_decontam.R` placeholder
  is not claimed as the clean-clone generator for the separate-PC hand-back.
- **Observed association — raw-name-table count total and source-study relative abundance track sensitivity loss.** An
  unclustered descriptive logistic diagnostic gives coefficients **−1.87** for log10(the 35-column
  reprocessed raw-name-table count total) and **−3.45** for log10(source-study relative abundance).
  This processed total includes `Unassigned` and is not whole-library sequencing depth. Missed cells have **7.3× lower**
  source-study relative abundance than matched cells (median 0.029 vs 0.208). These taxon×library cells
  repeat taxa and site-times, so the coefficients and raw medians diagnose structure rather than provide
  cluster-aware biological inference. Precision remains high; the discrepancy is primarily lower recall.

**Phase-3 occupancy (ψ) implication — a recoverability risk, not a known bias correction:**

The missed detections are structured by the raw-name-table count total and source-study relative
abundance, both of which co-vary with the upstream-to-downstream biological gradient. This can confound
an occupancy covariate with the detection process. A separately defined whole-library read-depth term
and structural-NA coding were therefore required tests, not
guaranteed fixes. In the completed Phase-3 stress test, the known-truth recovery MCMC failed convergence
and 18 observed libraries that are zero within the 14-label harmonised input are dropped by the
community model while the naïve comparator retains them; each contains positive reads under raw taxonomy
labels. Consequently, no sparse-label occupancy
magnitude or direction is retained as a result.

**Safeguards tested (necessary but not sufficient):**
1. **Read depth as a detection covariate** in the occumb detection sub-model.
2. **Dropped / low-depth libraries coded as structural NA**, not observed-zero.
3. A **minimum read-depth filter** for detection-history construction.

**Lineage status:** the 16% pooled and 23–29% selected-label shortfalls were measured on
`species_x_sample_counts_rawnames.rds`, produced by collapsing the regenerated 647 × 1,348
`seqtab_clean.rds` table through `taxonomy.tsv`. They therefore describe the actual Phase-3 input lineage,
not a pre-decontamination proxy. They remain detection-agreement statistics, not occupancy estimates.

---

## 6. Open items

### 6.1 Characterise the leading Unassigned ASVs — ✅ DONE (2026-06-28)
The top 20 Unassigned ASVs (5.59% of the clean table, within the 5.9% Unassigned bin, incl. ASV0007 + ASV0008) were
BLASTed against NCBI nt (`data/derived/unassigned_blast_summary.tsv`). All 20 hit at 99–100% identity to
**non-target vertebrates**, none a freshwater fish: **livestock 4.82%** (cattle ASV0007 886k + sheep ASV0008
725k + pig — a Welsh grazing catchment), **wildlife** (water vole, field vole, red fox, mallard, crow),
**amphibian** (palmate newt *Lissotriton helveticus*, 0.21%), and **one marine flatfish stray**
(*Platichthys stellatus*, which the authors themselves dropped in D3). These audited ASVs show that the
bin is dominated by non-target vertebrate bycatch; they do not prove that the smaller unexamined tail
contains no fish. The one-library *Leuciscus* alternative therefore remains explicit. The observed recall
deficit is not explained by the audited leading Unassigned ASVs.

### 6.2 Re-measure the detection-offset on the regenerated ASV table — ✅ RESOLVED (was a misread)
The offset was **already** measured on the regenerated ASV table: `species_x_sample_counts_rawnames.rds`
is `seqtab_clean.rds` (647 × 1348 corrected ASVs) collapsed via `taxonomy.tsv` — **not** the pre-decontam
`species_table.rds`. So the 16% figure reflects the true Phase-3 input. The only residual robustness item is
a **tag-jump-cutoff sensitivity sweep** (how the offset moves at floors 0 / 0.0005 / 0.001); the sweep lives
PC-side and is a nice-to-have, not a blocker.

### 6.3 Measure the precision/false-positive side — ✅ DONE (2026-06-28)
Added to `R/08` (§5b). Cell-level over the 14 shared analysis labels × 642 libraries: **precision 0.99, recall 0.84,
F1 0.91**; 1,962 both / **14 reprocessed-only** / 382 source-only. This raw asymmetry is descriptive;
McNemar's unclustered cell p-value is not used as biological inference. The reprocessed-only cells are
trace and concentrate in *Abramis_brama* (5), a taxonomic review candidate.
Conclusion: over-calling is rare; the dominant discordance is lower recall on the reprocessed side.

### 6.4 Harden the drop gates — ✅ DONE (2026-06-28)
`R/08` now asserts `sum(DROP_A %in% names(A)) == 10` and `all(names(RELAB_A) %in% names(A))`, plus the
existing `setequal` gate on the 35 reprocessed columns (34 assigned labels + `Unassigned`), so the 55-/14-label universes validate rather than pass by
coincidence.

### 6.5 Phase-3 disposition — ✅ TESTED AND CLAIM-GATED
Phase 3 included read depth in the `occumb` detection sub-model and used structural NA for unavailable
libraries. The known-truth recovery MCMC nevertheless failed convergence and the 18-library harmonised-input comparison mask
was inconsistent, so rare-species occupancy magnitudes and directions remain exploratory rather than final.

---

*Bottom line (Phase-2 close-out complete): composition agrees well (mean Bray-Curtis similarity 0.84, cosine similarity 0.90, all
≥3-library taxa recovered), over-calling is rare (precision 0.99), and the Unassigned bin is dominated by
non-target vertebrate bycatch. Reprocessed detection-cell recall is ~84% of the source-study table and is
structured along the river gradient. Phase 3 therefore tested depth-aware detection and structural NA, then
demoted rare-species occupancy when the recovery gate failed. **PASS-WITH-CAVEATS — Phase 2 is closed.***
