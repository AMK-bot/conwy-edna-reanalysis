# Accession-backed *Salmo* multi-marker benchmark

**Frozen specification:** `docs/impact-sprint-preregistration.md`, section 4  
**Machine-readable contract:** `data/metadata/impact_sprint_manifest.json`  
**NCBI retrieval date:** 2026-07-19  
**Evidence scope:** submitted-label discrimination in a small public-record panel

## Question and claim boundary

Within a provenance-audited complete-mitogenome panel, what is the finest submitted
taxonomic rank consistently discriminated by the MiFish 12S product, a 648-bp 5-prime
COI window and a 1,140-bp cyt-b window?

Public-record organism names are submitted annotations, not independently verified
biological identities. This benchmark does **not** establish local occurrence, stock,
ancestry, hybrid status, abundance or a species identity for `ASV0001`. Its marker
windows are independent of the existing Conwy query fragment.

## Network-only acquisition

`scripts/fetch_salmo_multimarker.py` queried NCBI Nucleotide E-utilities with the exact
submitted organism label, the mitochondrial filter, a 15,000--20,000-bp length screen
and `NOT srcdb_refseq[prop]`. Results were accession-sorted and capped at 50 per label.
The official ESearch responses, plus ESummary responses when a query returned records,
are cached under `data/public/salmo_multimarker/query_receipts/`. Zero-result queries
have no invented ESummary receipt. Each candidate was then inspected in
accession-version order from its raw GenBank record until five records passed.

A record passed only when it had the exact submitted organism annotation, a complete
primary mitochondrial record, a unique unambiguous MiFish-U product, one annotated
COI/COX1 CDS with at least 648 unambiguous coding bases, one annotated CYTB CDS with at
least 1,140 unambiguous coding bases, and no duplicate voucher/isolate key already in
the selected label panel. RefSeq derivatives and partial-marker substitutes were not
used. Records without a voucher/isolate/strain/haplotype/clone qualifier remain distinct
at accession-record level but are flagged as having undocumented specimen independence.

The official submitted definition for `LT617630.1` contains the typographical wording
“mitochondrionial geneome”. Its 16,677-bp circular record, source annotation and full
marker features passed; the typo is retained in the raw record rather than converted
into an exclusion.

### Coverage audit

| Submitted label | NCBI total | Returned (cap 50) | Passing | Rank-gate eligible (minimum 2) |
|---|---:|---:|---:|---|
| *Salmo trutta* | 30 | 30 | 5 | yes |
| *Salmo salar* | 21 | 21 | 5 | yes |
| *Salmo carpio* | 0 | 0 | 0 | no |
| *Salmo marmoratus* | 0 | 0 | 0 | no |
| *Salmo macedonicus* | 0 | 0 | 0 | no |
| *Salmo macrostigma* | 0 | 0 | 0 | no |
| *Salmo obtusirostris* | 1 | 1 | 1 | no |
| *Salvelinus alpinus* | 111 | 50 | 5 | yes |
| *Oncorhynchus mykiss* | 41 | 41 | 5 | yes |

The NCBI totals describe the frozen query receipts on the retrieval date, not timeless
database counts. Fourteen of the 21 selected records have a documented voucher,
isolate or strain-like source key; seven have accession-level independence only. The
species panel therefore contains two eligible *Salmo* labels, below the preregistered
minimum of four. No partial COI, cyt-b or 12S records were introduced to repair that
shortfall.

### Selected primary accessions

| Submitted label | Accession versions |
|---|---|
| *S. trutta* | `JQ390057.1`, `LT617630.1`, `LT617631.1`, `LT617632.1`, `MF621760.1` |
| *S. salar* | `AF133701.1`, `JQ390055.1`, `JQ390056.1`, `KF792729.1`, `LC012541.1` |
| *S. obtusirostris* | `LT617633.1` |
| *Salvelinus alpinus* | `AF154851.1`, `MF621740.1`, `MF621741.1`, `MF621743.1`, `MN530962.1` |
| *Oncorhynchus mykiss* | `DQ288268.1`, `DQ288269.1`, `DQ288270.1`, `DQ288271.1`, `KP013084.1` |

Every accession has an official NCBI Nucleotide URL, cached raw-record path and SHA-256
digest in `data/metadata/impact_salmo_accessions.tsv`.

### Assessed exclusions before the five-record cap

| Submitted query | Accession | Retained reason |
|---|---|---|
| *S. trutta* | `AM910409.1` | organism annotation is *S. trutta trutta*, not the exact frozen label |
| *S. trutta* | `KP171723.1`, `KT633607.1`, `KT634053.1`, `LC137015.1` | organism annotation is *S. trutta fario*, not the exact frozen label |
| *S. alpinus* | `CM071963.1` | organism annotation is *S. alpinus oquassa*, not the exact frozen label |
| *S. alpinus* | `CM098937.1`, `CM134157.1` | no uniquely recognised annotated COI/COX1 CDS |

Candidates after the fifth passing accession were left unassessed and are explicitly
marked `not_assessed_after_five_passing` in the accession manifest.

## Marker extraction

- **MiFish 12S:** the unique product between
  `GTCGGTAAAACTCGTGCCAGC` and the reverse binding site of
  `CATAGTGGGGTATCTAATCCCAGTTTG`, allowing at most three substitutions per
  primer and a 140--250-bp inner product. All 21 selected products were 170 bp,
  with one forward-primer substitution and zero reverse-primer substitutions.
- **COI:** the first 648 unambiguous coding bases of the single annotated
  mitochondrial COI/COX1 CDS, after respecting feature strand and `codon_start`.
- **cyt-b:** the first 1,140 unambiguous coding bases of the single annotated
  CYTB CDS, after respecting feature strand and `codon_start`.

All three markers come from the same 21 selected records. Sequences and extraction
coordinates are in `data/public/salmo_multimarker/marker_sequences.fasta` and
`data/metadata/impact_salmo_marker_manifest.tsv`.

## Offline metrics and gates

`R/12_salmo_multimarker_benchmark.R` makes no network call. It verifies every cached
receipt, raw GenBank file and marker sequence against its SHA-256 digest before
analysis. Global Levenshtein edit distance is divided by the longer sequence length;
exact haplotypes require identical strings.

For every eligible submitted *Salmo* label, the barcode gap is the nearest
between-label distance minus the maximum within-label distance. Leave-one-accession-out
classification retains **all** tied nearest neighbours; a species call is correct only
when every tie has the query's submitted species annotation. The genus comparison asks
whether each eligible *Salmo* query is strictly nearer to another eligible *Salmo*
record than to either eligible outgroup genus.

The species gate requires at least four eligible *Salmo* labels, 100% all-ties
classification, no cross-label exact haplotype and a positive gap for every eligible
label. The genus gate requires at least two eligible *Salmo* labels, both outgroup
genera with at least two records, species not supported, and strict *Salmo*-versus-each-
outgroup separation for every query.

## Results

| Marker | Eligible candidate *Salmo* labels | Narrowed-panel all-ties LOO | Cross-label exact haplotypes | Barcode gaps | Species gate | Genus boundary | Finest supported rank |
|---|---:|---:|---:|---|---|---|---|
| MiFish 12S | 2 of 7; minimum 4 | 10/10 | 0 | both positive | coverage-inadequate | 10/10, pass | genus |
| COI 5-prime | 2 of 7; minimum 4 | 10/10 | 0 | both positive | coverage-inadequate | 10/10, pass | genus |
| cyt-b | 2 of 7; minimum 4 | 10/10 | 0 | both positive | coverage-inadequate | 10/10, pass | genus |

The positive barcode gaps were 0.02353 and 0.00588 for MiFish 12S, 0.05864 and
0.06019 for COI, and 0.05000 and 0.05351 for cyt-b (*S. trutta* and *S. salar*,
respectively; normalised distance). These are two-label diagnostics, not evidence that
any marker consistently separates the preregistered seven-label *Salmo* candidate set.

The defensible result is therefore:

> In this accession-backed submitted-label panel, MiFish 12S, COI and cyt-b each
> passed the preregistered record-level genus boundary, but only two of seven candidate
> labels were eligible and the species gate required at least four. The longer-marker
> species comparison is therefore withheld.

## Limitations and evidence needed next

- Five of seven candidate *Salmo* labels lacked two exact-label complete primary
  mitogenomes; four had none and *S. obtusirostris* had one.
- The eligible species diagnostic is restricted to *S. trutta* versus *S. salar*.
  Its perfect all-ties result cannot be extrapolated to the uncovered labels.
- Seven selected records—including all five *S. salar* records—lack specimen-level
  source identifiers, so duplicate biological sampling cannot be excluded beyond
  distinct primary accession records. The genus result is consequently a record-level
  boundary, not proof based on independently vouched specimens.
- Mitochondrial markers cannot resolve nuclear ancestry, hybrids, introgression,
  incomplete lineage sorting or stock provenance.
- Database labels may be synonymous, outdated or misidentified; this analysis measures
  discrimination of submitted annotations, not taxonomic truth.
- Species-rank evaluation would require at least two auditable complete primary
  mitogenomes for at least two additional exact *Salmo* labels, followed by the same
  frozen all-ties, exact-sharing and barcode-gap tests.

## Reproduction and verification

Network acquisition (only when deliberately refreshing the public cache):

```bash
python3 scripts/fetch_salmo_multimarker.py --retrieval-date 2026-07-19
```

Offline benchmark and tests:

```bash
Rscript R/12_salmo_multimarker_benchmark.R
Rscript R/test_salmo_multimarker.R
Rscript R/12_salmo_multimarker_benchmark.R --validate-only
```

On the test laptop, the offline build ran in 6.95--8.69 seconds. The test suite performs
two additional complete builds and verifies byte-identical SHA-256 digests for six TSVs,
the PNG and the SVG. Runtime telemetry is printed only and is not included in canonical
outputs.

## Output map

- acquisition receipts and raw records: `data/public/salmo_multimarker/`
- acquisition provenance: `data/metadata/impact_salmo_queries.tsv`,
  `impact_salmo_accessions.tsv`, `impact_salmo_marker_manifest.tsv`
- panel coverage: `data/metadata/impact_salmo_panel_coverage.tsv`
- distances, haplotypes, gaps, LOO and gates: `data/derived/impact_salmo_*.tsv`
- figure: `figures/impact_salmo_multimarker.png` and `.svg`
- offline analysis: `R/12_salmo_multimarker_benchmark.R`
- invariant/regression tests: `R/test_salmo_multimarker.R`
