# Salmo and fish-traits validation receipt

Validation date: 2026-07-18

## Machine checks passed

- `fish_traits.tsv` has exactly 14 unique analysis labels and exactly matches the columns of the harmonised label matrix.
- Every semicolon-delimited `source_id` in the trait table resolves to `fish_traits_sources.tsv`; all 28 source rows carry the 2026-07-18 access date.
- `salmo_reference_amplicons.fasta` contains the query plus eight references; every sequence is exactly 170 nt.
- Every aligned sequence in `salmo_sequence_alignment.tsv` is identical to its FASTA sequence. Match, mismatch, identity and edit-distance fields recalculate from the bases.
- `salmo_edit_distance_matrix.tsv` is 9 × 9, symmetric, has a zero diagonal and exactly equals recomputed position-wise distances.
- Four unique haplotypes occur in `salmo_haplotype_nodes.tsv`; all three edge weights in `salmo_haplotype_edges.tsv` recalculate from their sequences.
- All 42 taxon × section rows in `fish_species_section_prevalence.tsv` reproduce from `species_x_sample_counts_clean.rds` after replicate OR-combination and restriction to the 11 timepoints present at all 14 sites.
- All 21 guild × section rows in `fish_trait_section_summary.tsv` reproduce from the deterministic mapping in `fish_trait_analysis_groups.tsv`. The balanced panel contains 154 site × timepoint units: 55 upper, 55 middle and 44 lower.

The machine-check outputs were:

```text
PASS: 14 trait rows; 28 sources; 9 aligned sequences; symmetric 9x9 distance matrix; 4 haplotypes/3 MST edges
PASS: 42 taxon-section prevalences and 21 guild-section summaries reproduce exactly from 14 sites x 11 common timepoints (154 units)
```

## Source-record check

The eight cited nucleotide records were fetched together from NCBI EFetch in GenBank format on 2026-07-18. Accession/version, record length, organism, isolate/voucher, geographic qualifier, collection date where present and LOCUS date were checked against `salmo_reference_provenance.tsv`. The current NCBI records confirm the exact sequence labels used in the alignment.

The `PX646298.1` and `PX646330.1` LOCUS dates are 2026-05-14, later than the local MIDORI2 GenBank271 release dated 2026-04-07. They are therefore marked only as contemporary exact-match witnesses.

## Visual check

`salmo_molecular_detective.svg` is valid XML and is backed by the committed alignment, haplotype-node and haplotype-edge tables. A 1600 × 1000 PNG was rendered from the SVG and visually inspected at original resolution. Labels, all seven mismatch marks, the three edge weights and the five exact species labels are visible without overlap.

## Reproducibility boundary

The standalone base-R generator is `scripts/build_biological_evidence.R`. It is deliberately not coupled to `scripts/build_website.py`; release orchestration can call it explicitly before building the website bundle.

```bash
Rscript scripts/build_biological_evidence.R
Rscript scripts/build_biological_evidence.R --validate-only
Rscript scripts/build_biological_evidence.R --no-rds --validate-only
```

The default invocation rebuilds both balanced trait-summary tables when the RDS and covariates are available, validates all Salmo sequences/distances/haplotype edges, and regenerates the SVG and 1600 × 1000 PNG from committed evidence tables. `--validate-only` is read-only. `--no-rds` exercises a clean-clone fallback that validates committed prevalence arithmetic and guild incidences while still validating the Salmo evidence; without the RDS it cannot recover guild co-occurrence unions from marginal summary rows.

The generator uses base R only. Byte-identical raster output can still depend on the platform's Cairo/font implementation; the script therefore treats sequence, distance, summary arithmetic and figure source-data consistency as the portable semantic contract.

## Limits not resolved by validation

- Validation establishes internal arithmetic and provenance consistency, not the truth of public reference organism names.
- Current NCBI exact matches do not reconstruct the original 2020 database.
- The missing historical MIDORI target accession remains unrecoverable from committed outputs.
- Guild assignments are literature-based species associations; eDNA does not observe realised habitat use, diet, size class, migration or abundance.
- Section summaries are descriptive and inherit assignment uncertainty. Existing restricted-series inferential tests did not cross the 0.05 support threshold.
