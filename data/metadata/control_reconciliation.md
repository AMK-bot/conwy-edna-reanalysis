# Phase-2 negative-control reconciliation

**Resolved:** 2026-07-18

The earlier 79-versus-77 discrepancy describes two different stages, not two
competing negative-control definitions:

- **79 manifest negatives:** 54 tier-1 controls plus 25 tier-2 process blanks,
  including 14 `BLANK_[CD]#_well` wells.
- **77 decontam inputs:** two manifest well blanks, `BLANK_D1_A1` and
  `BLANK_C3_A1`, have no row in the post-DADA2 ASV table because they had zero
  usable reads. The remaining 12 well blanks enter the table, giving 77 total
  negative libraries used by decontam.

The local evidence is `data/derived/controls_audit.md` (manifest accounting) and
`data/derived/dada2_track.tsv` (the 12 surviving `BLANK_[CD]#_well` aliases).
Those Phase-2 hand-back files are intentionally git-ignored and are not claimed
to regenerate from a clean clone. This included note preserves the resolved
interpretation and the exact two excluded aliases.

Use the numbers with their denominators: **79 negatives in the manifest; 77
negative libraries in the decontam input**. Do not collapse them into one
unqualified control count.
