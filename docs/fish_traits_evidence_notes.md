# Fish traits evidence notes

## Purpose

`data/metadata/fish_traits.tsv` is a conservative interpretation layer for the 14 labels in `data/derived/species_x_sample_counts_clean.rds`. It supports an portfolio-edition comparison between plausible ecological signal and what the molecular data can actually resolve. It is not a checklist of confirmed Conwy fish populations. Every species natural-history entry is conditional on the reported/harmonised label: compatibility with a trait profile does not establish local species identity, occurrence, origin or conservation status.

The companion `data/metadata/fish_traits_sources.tsv` records the evidence source, URL or repository path, access date, scope and caveat. Web sources were restricted to an established European trait database, FishBase species syntheses, peer-reviewed database documentation, UK statutory conservation pages and Natural Resources Wales guidance/strategy. FishBase is explicitly treated as a secondary synthesis rather than local occurrence proof.

## Decision rules

1. Trait assignments follow the ecological taxon that can defensibly be interpreted, not necessarily the harmonised string.
2. `NA` is used when harmonisation created an incoherent biological entity (`Salvelinus_malma`) or when trace support and geographic/taxonomic uncertainty make guild interpretation unsafe (`Perca_flavescens`).
3. Complex/proxy labels are explicit for `Cottus_perifretum`, `Phoxinus_ujmonensis`, `Pungitius_pungitius`, `Salmo_trutta` and the harmonised rainbow-trout label.
4. Native status is separated from catchment occurrence. A fish described as native in Welsh rules is not thereby shown to occur naturally in the Conwy.
5. A detection count is the number of libraries with a non-zero harmonised count. It is not occupancy, abundance, biomass, persistence, spawning evidence or a population trend.
6. Global IUCN categories reported by FishBase are not substituted for local stock condition. Where official Welsh conservation statements exist (eel, salmon and sea trout), they are reported as Wales-wide relevance, not as a Conwy estimate.

## River-zone and flow vocabulary

The rheophilic/limnophilic/eurytopic and trout/grayling/barbel/bream-zone vocabulary follows freshwaterecology.info. These categories are tendencies. They should be used to structure hypotheses or explain community contrasts, not to infer measured hydraulic conditions from a single taxon.

## High-priority admissions

- The clean matrix contains analysis labels, not a uniformly species-validated inventory.
- `Oncorhynchus_mykiss_x_Salmo_salar` is a biologically incoherent harmonisation label with no evidence of a real hybrid here; only an `O. mykiss` ecological proxy is supplied.
- `Perca_flavescens` is a one-ASV/four-read geographic/taxonomic review label and receives no ecological traits. This does not prove a mis-hit or any replacement identity.
- `Phoxinus_ujmonensis` belongs to a same-genus label cluster. Use `Phoxinus spp.`; `P. phoxinus` literature supplies traits but is not proposed as local identity.
- `Pungitius_pungitius` requires genus-level caution in north-western Europe. Use a `Pungitius` trait envelope and propose no replacement identity.
- `Salvelinus_malma` combines incompatible raw labels and receives no ecological traits.
- The dominant `Salmo` 12S haplotype is not species-diagnostic. See `docs/salmo_molecular_detective_methods.md` and the associated sequence evidence.

## Appropriate use

Suitable uses include a qualitative traits appendix, pre-specified guild summaries with uncertainty, or sensitivity analyses that exclude low-confidence labels. Unsuitable uses include species-richness claims without a label audit, local conservation-status assertions from eDNA detections alone, or habitat reconstruction that treats trait categories as direct environmental measurements.
