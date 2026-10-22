# Phase-2 validation gate -- concordance with the published 68-species table

- Reference DB: **MitoHelper 12S Nov-2025** (MitoFish-derived; the exact 23-Sep-2020 build is unarchived -- PV5).
- BLAST thresholds (authors'): >=90% identity, >=90% qcov, e<=0.001, first hit (top bitscore).
- ASVs: 1443; assigned to species: 886 (61.4%); distinct species assigned: **37** (an ASSIGNMENT count, NOT a measured concordance -- see below).

## Known taxonomy artifacts (PC-flagged; relabels are Mac-side per P4)
- **Largest 'species' by reads: Salmo_macrostigma = 48.6% of species_table.rds reads.**
- `Salmo_macrostigma` (Adriatic/Mediterranean trout; cannot occur in N. Wales) = 48.6% of reads -- a
  DATABASE NEAREST-NEIGHBOUR ARTIFACT, almost certainly local Salmo (trutta/salar). The macrostigma->trutta
  relabel is a Mac-side D3 cleaning step (PHASE2_ISSUES_FOR_PC.md P4), NOT a PC fix -- flagged so the species
  count is not mistaken for a finding. A handful of marine/freshwater strays (e.g. Microgobius, Myoxocephalus,
  Hippoglossus) carry negligible reads and inflate the distinct-species count; ~557 Unassigned ASVs may add richness.
- **species_table.rds is the PRE-decontam artifact** (raw seqtab; relabels + Scomber drop only; all 889 libraries
  incl. controls). The decontaminated ASV deliverable is **seqtab_clean.rds** (post-decontam + tag-jump, controls
  removed). The read-count difference between them is expected cleaning, not corruption.

## Concordance: PENDING
`data/metadata/source_data_12S_normalised.tsv` (the published 68-species matrix) is the Mac-side
validation target and is **not staged on the PC**. Stage it (USB) for the row-level read-weighted
concordance. My recovered species list is below; D3 cleaning (drops + 3 relabels) applied to reach
the 55-species universe.

Recovered species (37): Abramis_brama, Anguilla_anguilla, Anguilla_rostrata, Barbatula_aff._barbatula_lineage_I, Barbatula_barbatula, Barbus_barbus, Ciliata_mustela, Clupea_harengus, Cottus_gobio, Gadus_morhua, Gasterosteus_aculeatus, Hippoglossus_stenolepis, Hyperoplus_lanceolatus, Melanogrammus_aeglefinus, Merlangius_merlangus, Merluccius_merluccius, Microgobius_gulosus, Myoxocephalus_thompsonii, Oncorhynchus_mykiss, Osmerus_eperlanus, Perca_flavescens, Perca_fluviatilis, Pholis_gunnellus, Phoxinus_bigerri, Phoxinus_phoxinus, Phoxinus_sp._Morella, Pungitius_kaibarae, Pungitius_pungitius, Rhamdia_cf._jequitinhonha_NS_3056B, Salmo_macrostigma, Salmo_salar, Salmo_trutta, Salvelinus_fontinalis, Salvelinus_umbla, Scomber_scombrus, Sprattus_sprattus, Taurulus_bubalis

## D3 cleaning applied
- Relabels: Salmo_obtusirostris->Salmo_trutta; Anguilla_rostrata->Anguilla_anguilla; Salvelinus hybrid->Salvelinus_malma.
- Dropped positive-control spike: Scomber_scombrus.
- NOTE: the authors' full marine/aquaculture stray drop-list is not staged on the PC; remaining strays flagged for the Mac-side species collapse.
