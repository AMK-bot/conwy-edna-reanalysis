# Salmo molecular detective: methods, result and limits

## Question

Can the dominant 170-bp `Salmo` 12S ASV (`ASV0001`) support a species-level call, and why did the local classifiers report different *Salmo* species-level labels?

## Evidence chain

`ASV0001` is 170 bp long and has 18,531,288 reads in the broader field+control `data/derived/asv_map.tsv` universe; its field-only count is 14,123,368 reads (39.5% of the 647-field-library table). The **local Nov-2025 MitoHelper baseline** first hit in `data/derived/taxonomy.tsv` is `Salmo macrostigma`, 100% identity, accession [PQ275758](https://www.ncbi.nlm.nih.gov/nuccore/PQ275758). This is not the source study's unarchived Sep-2020 MitoFish output. Reclassification against the same Nov-2025 MitoFish reference vintage with a different classifier reported `S. trutta` at 100%; because the database vintage was shared, this is classifier/tie-breaking evidence, not database independence. An independent MIDORI2 GenBank271 run reported `S. carpio` at 100%. The committed MIDORI provenance did not retain the selected accession, so its exact target cannot now be reconstructed from repository outputs.

The current NCBI witness comparison is more direct. The ASV is an exact 170/170 match to records labelled:

- `S. macrostigma` — [PQ275758](https://www.ncbi.nlm.nih.gov/nuccore/PQ275758)
- `S. trutta` — [LC385098](https://www.ncbi.nlm.nih.gov/nuccore/LC385098)
- `S. carpio` — [PX646298](https://www.ncbi.nlm.nih.gov/nuccore/PX646298)
- `S. marmoratus` — [OP931658](https://www.ncbi.nlm.nih.gov/nuccore/OP931658)
- `S. macedonicus` — [PX646330](https://www.ncbi.nlm.nih.gov/nuccore/PX646330)

All NCBI records were accessed 2026-07-18. The two `PX` records carry a GenBank LOCUS date of 2026-05-14, after the local MIDORI2 GenBank271 release dated 2026-04-07. They therefore demonstrate contemporary reference-label non-uniqueness; they cannot be invoked as the actual references behind the earlier MIDORI call.

Nearby comparison records are one substitution away (`S. trutta`, [OP931685](https://www.ncbi.nlm.nih.gov/nuccore/OP931685)), two substitutions away (`S. obtusirostris`, [NC_037939](https://www.ncbi.nlm.nih.gov/nuccore/NC_037939)) and four substitutions away (`S. salar`, [PQ285811](https://www.ncbi.nlm.nih.gov/nuccore/PQ285811)). The exact sequences, provenance, pairwise distances and figure inputs are committed in `data/metadata/salmo_*` files.

## Sequence procedure

1. Read `ASV0001` from `data/derived/asv_map.tsv`.
2. Retrieve each cited NCBI nucleotide record in FASTA form.
3. For each record, select the 170-nt forward-strand window with the fewest substitutions to `ASV0001`. All selected windows aligned without a gap.
4. Count position-wise substitutions. With equal-length ungapped windows, Hamming distance equals Levenshtein edit distance.
5. Collapse identical windows into haplotypes and calculate all pairwise distances. The displayed network uses the minimum-spanning edges in `salmo_haplotype_edges.tsv`; it is a distance diagram, not a phylogenetic tree.

The output was independently checked for 170 bases per sequence, matrix symmetry, zero diagonal, sequence/TSV agreement and consistency of every stated mismatch position. The checks and reproducibility boundary are recorded in `docs/salmo_fish_traits_validation_receipt.md`.

## Result

The strongest defensible molecular result is:

> `ASV0001` is robustly assignable to genus `Salmo`, but this 170-bp fragment is not species-diagnostic. A single exact haplotype occurs in public records carrying at least five species labels.

Ecological and geographic context may make a brown-trout-complex interpretation a plausible prior or trait proxy in the Conwy, but a prior cannot convert a non-diagnostic sequence into species-level molecular proof. The molecular call supported by this dominant haplotype is `Salmo spp.`; a brown-trout-complex interpretation is not a co-equal sequence identification.

This conclusion is also consistent with broader taxonomic work: Tougard et al. (2018) concluded that `S. macrostigma` is not a valid species but a lineage within the brown trout complex ([doi:10.1111/jfb.13751](https://doi.org/10.1111/jfb.13751)). MIDORI2 itself is a quality-controlled, regularly updated, GenBank-derived mitochondrial reference system, so database changes and source-label accuracy remain part of assignment uncertainty (Leray, Knowlton & Machida 2022, [doi:10.1002/edn3.303](https://doi.org/10.1002/edn3.303)).

## Limits that must travel with the result

- Public-record organism names are submitted identifications, not independent ground truth.
- The selected records are a small witness set designed to test identifiability, not an exhaustive phylogeographic analysis.
- Exact identity can reflect shared mitochondrial haplotypes, incomplete lineage sorting, introgression, synonymy/taxonomic inflation or reference misidentification; this analysis does not discriminate among those mechanisms.
- A short mitochondrial marker cannot distinguish hybrids, nuclear ancestry or population provenance.
- The current NCBI witness set is not a reconstruction of the original 2020 reference database.
- The query read total is not abundance or biomass and should not weight taxonomic certainty.
- The four-substitution `S. salar` comparison concerns `ASV0001` only. Other local ASVs labelled `S. salar` require their own reference audit.
- The local MIDORI top-hit accession was not retained. That missing provenance prevents exact historical reproduction of the `S. carpio` call.

## File map

- `data/metadata/salmo_reference_amplicons.fasta` — query and the eight cited 170-bp comparison windows
- `data/metadata/salmo_sequence_alignment.tsv` — sequence-level matches, identities and mismatch positions
- `data/metadata/salmo_edit_distance_matrix.tsv` — complete pairwise edit-distance matrix
- `data/metadata/salmo_reference_provenance.tsv` — accession, record, locality/date notes and historical role
- `data/metadata/salmo_classifier_provenance.tsv` — local Nov-2025 baseline, same-reference-vintage, independent-resource and current direct-sequence evidence
- `data/metadata/salmo_haplotype_nodes.tsv` and `salmo_haplotype_edges.tsv` — source data for the compact distance network
- `figures/salmo_molecular_detective.svg` — visual summary
- `scripts/build_biological_evidence.R` — deterministic base-R validator and table/figure generator
