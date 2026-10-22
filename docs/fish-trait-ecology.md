# Fish traits and site-aware ecology

## What the balanced panel shows

The clearest descriptive pattern is **downstream guild broadening, not taxon replacement**. Rheophilic-freshwater-phase labels occur in 53/55 upper, 55/55 middle and 44/44 lower site × timepoint units. Eurytopic labels, by contrast, occur in 0/55 upper, 7/55 middle and 42/44 lower units; limnophilic labels occur in 0/55, 1/55 and 11/44 units. Thus downstream samples retain the rheophilic component while adding labels whose conditional trait proxies include eel, stickleback and other habitat-generalist or slower-water associations. These are analysis-label detections; ecological coherence does not validate local species identities or occurrence.

The trophic view tells the same, cautious story. A broad benthic-feeding group occurs in 0/55 upper, 37/55 middle and 44/44 lower units. Small-invertebrate feeders occur in 0/55, 1/55 and 32/44 units. Upper units are instead almost entirely the `Salmo trutta`-complex signal, with one `Phoxinus` occurrence. This is a widening of detected trait space downstream, not evidence that feeding interactions or fish size were observed.

These numbers use presence/absence at the biological sampling unit. Replicate libraries were OR-combined within site × timepoint, then restricted to the 11 timepoints represented at every site: 5 upper sites × 11 = 55 units, 5 middle × 11 = 55 and 4 lower × 11 = 44. Raw reads were not treated as abundance.

The ecological vocabulary follows the standard European rheophilic/limnophilic/eurytopic and river-zonation framework at [freshwaterecology.info](https://www.freshwaterecology.info/), with definitions documented on its [rheophily](https://www.freshwaterecology.info/fwe_info.php?p=b2c9ZmlzaCNwYXJhbT0yNzg%3D) and [fish-region](https://www.freshwaterecology.info/fwe_info.php?lang=2&p=b2c9ZmlzaCNwYXJhbT0yNTQ%3D) pages (all accessed 2026-07-18). Species natural-history evidence and every caveat are itemised in `data/metadata/fish_traits_sources.tsv`.

## Which labels drive the pattern

- Upper: `Salmo_trutta` complex occurs in 53/55 balanced units; `Phoxinus` occurs once. The molecular detective analysis shows why this must remain a complex-level `Salmo` interpretation.
- Middle: the `Salmo_trutta` label is ubiquitous, with `Phoxinus` labels in 51/55, stone-loach labels in 36/55, `Salmo_salar` labels in 35/55 and `Anguilla_anguilla` labels in 6/55 units.
- Lower: `Salmo_trutta` and `Phoxinus` labels occur in every unit, accompanied by `Anguilla_anguilla` (41/44), `Salmo_salar` (40/44), stone-loach (36/44), three-spined-stickleback (32/44), the raw-rainbow-trout proxy (11/44) and `Pungitius` (9/44) labels.

The `Anguilla_anguilla`-label pattern is ecologically coherent with a catadromous fish using connected estuary–river–lake habitats ([FishBase species synthesis](https://www.fishbase.se/summary/Anguilla-anguilla.html)); three-spined-stickleback literature includes both resident and anadromous populations in shallow streams, lakes and estuaries ([FishBase species synthesis](https://www.fishbase.se/summary/Gasterosteus-aculeatus.html)). These are conditional coherence checks, not occurrence or identity validation independent of the same eDNA data.

## What is and is not supported

This interpretation is descriptive. In the repository's restricted-series tests:

- section PERMANOVA gives large descriptive `R²` values (0.879 Bray–Curtis, 0.753 Jaccard and 0.564 robust Aitchison), but each exact restricted-series `p = 0.0714`;
- no analysis label passes the documented multi-method differential-abundance rule plus the later spatial-series gate;
- no section indicator has restricted-series `p ≤ 0.05` (the smallest reported value is 0.077).

Accordingly, the pattern can be described as observed longitudinal turnover/broadening in this panel. It should not be promoted to a confirmatory habitat-effect, indicator-species or differential-abundance claim.

## Assignment sensitivity

Two labels are excluded from every guild summary: the one-ASV/four-read `Perca_flavescens` geographic/taxonomic review label and the incoherent `Salvelinus_malma` collapse. The former is excluded because support is trace and identity is unresolved; no mis-hit, local occurrence or alternative `Perca` identity is inferred. The harmonised `Oncorhynchus_mykiss_x_Salmo_salar` name is never treated as a real hybrid; its lower-river label detections use only the raw `O. mykiss` ecological proxy. Removing that proxy does not change any section-level guild-presence result above because other rheophilic and ontogenetically piscivorous labels already occur in those units.

`Phoxinus_ujmonensis` is interpreted only as `Phoxinus cf. phoxinus`; `Pungitius_pungitius` remains a `P. pungitius/laevis` complex; and the dominant `Salmo` haplotype remains `Salmo` / `S. trutta` complex rather than species proof. Official Welsh pages establish that Llyn Conwy is regulated as trout water ([Natural Resources Wales](https://naturalresources.wales/guidance-and-advice/business-sectors/fisheries/angling-byelaws/when-what-and-where-you-can-fish/non-migratory-trout-arctic-charr-coarse-fish-and-eels-open-seasons/?lang=en)), but ecological plausibility cannot repair a non-diagnostic marker.

## Management relevance without overreach

The defensible management connection is connectivity and monitoring priority: European eel is CITES Appendix II ([JNCC](https://jncc.gov.uk/resources/20810981-e500-4207-a1df-b37cd9eb47e1)); Atlantic salmon is a Habitats Directive Annex II species requiring clean gravel spawning habitat ([JNCC](https://sac.jncc.gov.uk/species/S1106/)); and Natural Resources Wales identifies salmon, sea trout and eel as current conservation priorities in its [2026–2031 fisheries strategy](https://naturalresources.wales/guidance-and-advice/business-sectors/fisheries/climate-resilient-fisheries-a-strategy-for-the-sustainable-management-of-migratory-and-freshwater-fishes-and-fisheries-of-wales-2026-2031/?lang=en) (accessed 2026-07-18). The eDNA panel can motivate targeted follow-up; it cannot estimate stock status, passage success, biomass or recruitment.

## Audit trail

- `data/metadata/fish_traits.tsv` — 14-label evidence table
- `data/metadata/fish_traits_sources.tsv` — cited source register with URLs, access dates and caveats
- `data/metadata/fish_trait_analysis_groups.tsv` — deterministic label-to-guild mapping and exclusions
- `data/metadata/fish_species_section_prevalence.tsv` — all label × section balanced-unit prevalence
- `data/metadata/fish_trait_section_summary.tsv` — guild-presence and taxon-incidence summaries
- `docs/fish_traits_evidence_notes.md` — decision rules and high-priority admissions
- `data/derived/layerE_permanova.tsv`, `layerE_da_consensus.tsv` and `layerE_indicator_species.tsv` — existing inferential receipts
- `scripts/build_biological_evidence.R` — standalone base-R rebuild/validation entry point

The site-aware summaries were recomputed from `species_x_sample_counts_clean.rds` and `conwy_covariates.tsv`; they do not modify the existing Layer E analysis.
