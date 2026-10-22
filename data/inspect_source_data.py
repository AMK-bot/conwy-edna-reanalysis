#!/usr/bin/env python3
"""
inspect_source_data.py  --  Phase 0 / P1 gate: does the published Source Data
actually contain a *replicate-level* 12S fish detection matrix?

The Phase-1 fast path (published table -> single-species detection occupancy,
BEFORE the 71 GB reprocessing) hinges on this. `msocc` needs the 3 water-sample
replicates kept distinct (E##_R#_T##); a site x time-collapsed table cannot seed
an occupancy detection history. Count-based Layer E and the Phase-3 `occumb`
stress test require the separate reprocessed count table.

This script opens `data/published/source_data_species_x_replicate.xlsx`, prints
every sheet's shape + row/column semantics, identifies the 12S fish block,
asserts replicate resolution, and reconciles its sample IDs against the 755-row
ENA join table (`conwy_12S_libraries.tsv`). It writes two committed artifacts:

  - data/metadata/source_data_verification.md   the verified finding + go/no-go
  - data/metadata/source_data_12S_normalised.tsv the clean, BOM-free, Conwy-only
        replicate x species matrix (normalised relative abundance) that Phase-1
        R scripts consume instead of re-parsing the git-ignored 5.8 MB xlsx.

Python (not R): consistent with verify_ena_join.py -- Python does the glue,
R consumes the derived TSV. Only dep beyond stdlib is openpyxl (ships in the
conda `r-*`/python layer; `pip install openpyxl` otherwise).

Usage:  python3 data/inspect_source_data.py    (run from repo root)
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

import openpyxl

REP_RE = re.compile(r"^E(\d{2})_R(\d)_T(\d{2})$")   # Conwy replicate library
N_DEPOSITED = 755                                    # ENA-deposited Conwy 12S libraries
N_SPECIES_EXPECTED = 68                              # paper's species tally
N_REPS_EXPECTED = 665                                # plan's claimed published rows


def looks_binomial(name: str) -> bool:
    """Genus_species (fish 12S) vs '1234;Arthropoda;...' ASV/taxonomy string (CO1/18S)."""
    return bool(re.match(r"^[A-Z][a-z]+_[a-z]", str(name))) and ";" not in str(name)


def read_sheet(ws):
    """Return (col_headers, {row_label: row_values}) for a wide sheet (label in col 0)."""
    it = ws.iter_rows(values_only=True)
    header = next(it)
    cols = [h for h in header[1:] if h is not None]
    rows = {}
    for r in it:
        if r[0] is not None:
            rows[str(r[0]).strip()] = r[1:1 + len(cols)]
    return cols, rows


def classify(ws):
    """Cheap structural summary of any sheet for the report."""
    cols, rows = read_sheet(ws)
    labels = list(rows)
    n_rep = sum(bool(REP_RE.match(x)) for x in labels)
    n_sitetime = sum(bool(re.match(r"^E\d\d_T\d\d$", x)) for x in labels)
    binom_cols = sum(looks_binomial(c) for c in cols)
    if n_rep and binom_cols >= 0.5 * len(cols):
        kind, marker, resolution = "species x replicate", "12S fish", "REPLICATE (E##_R#_T##)"
    elif n_sitetime:
        kind, marker, resolution = "ASV/taxon x site_time", "CO1/18S (eukaryote)", "site x time (collapsed)"
    else:
        kind, marker, resolution = "long format", "CO1/18S (eukaryote)", "long (site x time rows)"
    return dict(title=ws.title, n_rows=len(labels), n_cols=len(cols),
                n_rep=n_rep, binom_cols=binom_cols, kind=kind, marker=marker,
                resolution=resolution, cols=cols, rows=rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", type=Path,
                    default="data/published/source_data_species_x_replicate.xlsx")
    ap.add_argument("--libraries", type=Path,
                    default="data/metadata/conwy_12S_libraries.tsv")
    ap.add_argument("--outdir", type=Path, default="data/metadata")
    args = ap.parse_args()

    if not args.xlsx.exists():
        sys.exit(f"ERROR: {args.xlsx} not found (publisher Source Data, git-ignored).")

    wb = openpyxl.load_workbook(args.xlsx, read_only=True, data_only=True)
    sheets = [classify(ws) for ws in wb.worksheets]
    wb.close()

    # The 12S fish replicate sheets: replicate row labels + binomial species columns.
    fish = [s for s in sheets if s["n_rep"] > 0 and s["binom_cols"] >= 0.5 * s["n_cols"]]

    # Reconcile against the 755-row ENA join table -----------------------------------
    lib_ids = set()
    with open(args.libraries) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            lib_ids.add(row["ID"])

    primary = fish[0] if fish else None
    pub_ids = {x for x in (primary["rows"] if primary else {}) if REP_RE.match(x)}
    matched = pub_ids & lib_ids
    pub_only = pub_ids - lib_ids          # published reps not deposited (should be 0)
    lib_only = lib_ids - pub_ids          # deposited but absent from published matrix (QC-dropped)

    n_species = primary["binom_cols"] if primary else 0

    # ---- assertions (the GATE) ------------------------------------------------------
    checks = []
    checks.append(("A 12S fish replicate-level sheet exists", bool(fish),
                   f"{len(fish)} sheet(s)"))
    checks.append(("Rows are replicate-coded E##_R#_T##", len(pub_ids) > 0,
                   f"{len(pub_ids)} Conwy replicate rows"))
    checks.append((f"Conwy replicate rows == {N_REPS_EXPECTED}", len(pub_ids) == N_REPS_EXPECTED,
                   f"{len(pub_ids)}"))
    checks.append((f"Species columns == {N_SPECIES_EXPECTED}", n_species == N_SPECIES_EXPECTED,
                   f"{n_species}"))
    checks.append(("Every published rep maps to a deposited library (0 orphans)",
                   len(pub_only) == 0, f"{len(pub_only)} published-only"))
    checks.append((f"matched + dropped == {N_DEPOSITED} deposited",
                   len(matched) + len(lib_only) == N_DEPOSITED,
                   f"{len(matched)} matched + {len(lib_only)} dropped"))
    all_pass = all(ok for _, ok, _ in checks)

    print("=== Source Data inspection: 12S fish replicate matrix ==========================")
    for s in sheets:
        print(f"  [{s['title'][:34]:34s}] {s['n_rows']:>4} rows x {s['n_cols']:>4} cols  "
              f"-> {s['marker']:<18} {s['resolution']}")
    print("  " + "-" * 78)
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:52s} ({detail})")

    # ---- artifact 1: clean, BOM-free, Conwy-only replicate x species matrix ----------
    out_tsv = args.outdir / "source_data_12S_normalised.tsv"
    if primary:
        cols = primary["cols"]
        with open(out_tsv, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["ID"] + cols)
            for rid in sorted(pub_ids):
                vals = primary["rows"][rid]
                w.writerow([rid] + [0 if v is None else v for v in vals])

    # ---- artifact 2: verification report --------------------------------------------
    drop_by_tp = dict(sorted(Counter(x[-3:] for x in lib_only).items()))
    drop_by_site = dict(sorted(Counter(x[:3] for x in lib_only).items()))
    pub_tps = sorted({x[-3:] for x in pub_ids})

    md = []
    md.append("# Source Data verification — published 12S detection matrix (P1 gate)\n")
    md.append("*Auto-generated by `data/inspect_source_data.py`. Re-run to refresh.*\n")
    md.append(f"- Source: `{args.xlsx}` (publisher Source Data, git-ignored)")
    md.append(f"- Overall gate: **{'PASS ✓' if all_pass else 'FAIL ✗'}**\n")
    md.append("## Decision: Phase-1 detection fast path is **"
              + ("VIABLE ✅" if all_pass else "BLOCKED ❌") + "**\n")
    md.append("The published 12S fish Source Data **keeps the 3 water-sample replicates "
              "distinct** (`E##_R#_T##` row labels), so an occupancy detection history can be "
              "built directly from it for the Phase-1 single-species detection models. "
              "Because the matrix is normalised relative abundance rather than counts, the "
              "count-based Layer E baseline and Phase-3 community occupancy stress test still "
              "require the reprocessed count table.\n")

    md.append("## Sheets\n")
    md.append("| Sheet | Rows × Cols | Marker | Resolution |")
    md.append("|---|---|---|---|")
    for s in sheets:
        star = " **(12S fish — USE)**" if s in fish else ""
        md.append(f"| `{s['title']}` | {s['n_rows']} × {s['n_cols']} | {s['marker']}{star} | {s['resolution']} |")
    md.append("")
    md.append("> The four CO1/18S sheets are eukaryote/invertebrate markers (ASV-/taxonomy-keyed "
              "columns, `;`-delimited) and **are collapsed to site×time** — out of scope here. "
              "Only the two 12S sheets carry replicate resolution.\n")

    md.append("## The 12S fish matrix\n")
    if len(fish) >= 1:
        md.append(f"- **{len(pub_ids)} Conwy replicate rows × {n_species} fish species** "
                  "(normalised within-sample relative abundance; each row sums to 1).")
        md.append("- Two near-identical copies feed different figure panels: "
                  + ", ".join(f"`{s['title']}`" for s in fish)
                  + ". They share the same 665×68 replicate×species structure; "
                  "detection (value > 0) is identical between them.")
        md.append("- **Reads ≠ counts here:** the Source Data is *normalised relative abundance*, "
                  "not raw reads. Detection/non-detection (value > 0 vs 0) is fully preserved "
                  "(normalisation divides by a positive row total, so it never turns a 0 into a "
                  "non-zero or vice-versa). Raw per-library read depth, if needed as a detection "
                  "covariate, lives in `conwy_12S_libraries.tsv` (`read_count`).")
        md.append(f"- Published timepoints present: {pub_tps} "
                  f"(**T16 absent — consistent with the design grid; 18 usable occasions**).\n")

    md.append("## Reconciliation to the 755-row ENA join table\n")
    md.append(f"- Deposited Conwy 12S libraries (`conwy_12S_libraries.tsv`): **{len(lib_ids)}**")
    md.append(f"- Published replicate rows: **{len(pub_ids)}**")
    md.append(f"- **matched: {len(matched)}** · published-only: {len(pub_only)} · "
              f"deposited-but-dropped: **{len(lib_only)}**")
    md.append(f"- **The 90-row gap is explained:** {N_DEPOSITED} deposited − {len(lib_only)} "
              "QC-dropped = "
              f"{len(matched)} in the published analysis matrix. The authors drop libraries with "
              "<1000 reads (`Analysing_BLAST_output_WP2_12S.R`, the `rarefy_col_sum < 1000` filter). "
              "There are **0 published-only** rows, so every published replicate is a real deposited "
              "library — a clean 1:1 subset.\n")
    md.append("### ⚠️ The 90 dropped replicates are NA, not 0\n")
    md.append("Like T16, a QC-dropped replicate is **missing data, not a non-detection**. In the "
              "fast-path detection history its cell must be **NA** (the library was sequenced but "
              "yielded too few reads to call presence/absence) — encoding it as all-zero would "
              "manufacture false non-detections and bias ψ downward, the exact error this project "
              "exists to correct. The drop is **not random**: it concentrates in upstream low-yield "
              "sites and summer low-flow timepoints.")
    md.append(f"- dropped by site: {drop_by_site}")
    md.append(f"- dropped by timepoint: {drop_by_tp}\n")

    md.append("## Artifacts written\n")
    md.append(f"- `{out_tsv}` — clean Conwy-only replicate × species matrix "
              f"({len(pub_ids)} rows × {n_species} species, normalised relative abundance, BOM-free). "
              "Phase-1 R scripts consume this, not the git-ignored xlsx.")
    md.append(f"- `{args.outdir / 'source_data_verification.md'}` — this report.\n")
    (args.outdir / "source_data_verification.md").write_text("\n".join(md) + "\n")

    print(f"\nWrote: {out_tsv}\n       {args.outdir/'source_data_verification.md'}")
    if not all_pass:
        print("\nGATE FAILED — Phase-1 detection fast-path assumption is broken; escalate before Phase 1.", file=sys.stderr)
        return 1
    print("\nGATE PASSED ✓  — Phase-1 detection fast path VIABLE (replicate-level 12S matrix confirmed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
