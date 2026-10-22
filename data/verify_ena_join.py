#!/usr/bin/env python3
"""
verify_ena_join.py  --  Phase 0 "verify-on-download" gate for the Conwy 12S reanalysis.

Reads the ENA read_run metadata for study PRJEB48362 (Perry et al. 2024, Nat Commun)
and asserts the facts the execution plan depends on BEFORE any analysis begins:

  1. The Conwy 12S subset is the set of sample_alias matching  E##_R#_T##
     (E## = River Conwy sites; C## is a DIFFERENT river -- the #1 risk to avoid).
  2. That subset is exactly 755 read-runs  (798 theoretical - 43 absent).
  3. The ID (E##_R#_T##) <-> ENA run_accession join is strictly 1:1.
  4. Every run is paired-end (exactly two FASTQ files).

It then writes three committed artifacts into data/metadata/:
  - conwy_12S_libraries.tsv   the clean join table (ID -> run -> FASTQ URLs/md5/bytes)
  - design_gaps.tsv           the missing site x timepoint x rep cells
  - ena_join_verification.md  a human-readable, reproducible verification report

Usage:
    python3 data/verify_ena_join.py
Run from the repo root, or set --tsv / --outdir explicitly.

This is intentionally Python (not R): per the plan, Python/bash handle ENA download +
glue; R consumes the resulting conwy_12S_libraries.tsv. No third-party deps beyond pandas.
"""
from __future__ import annotations
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

# --- Design constants (verified against the archive; see data/README.md) -------------
SITES = [f"E{i:02d}" for i in range(1, 15)]      # E01 .. E14  (14 Conwy sites)
TIMEPOINTS = [f"T{i:02d}" for i in range(1, 20)]  # T01 .. T19  (19 sampling occasions)
REPS = ["R1", "R2", "R3"]                          # 3 water-sample replicates per site x time
N_THEORETICAL = len(SITES) * len(TIMEPOINTS) * len(REPS)  # 14 * 19 * 3 = 798
N_EXPECTED = 755                                   # deposited 12S Conwy replicate-libraries
ALIAS_RE = re.compile(r"^E(\d{2})_R(\d)_T(\d{2})$")


def load(tsv: Path) -> pd.DataFrame:
    if not tsv.exists():
        sys.exit(f"ERROR: {tsv} not found. Run data/download_ena.sh --metadata-only first.")
    return pd.read_csv(tsv, sep="\t", dtype=str).fillna("")


def conwy_subset(df: pd.DataFrame) -> pd.DataFrame:
    mask = df["sample_alias"].str.match(ALIAS_RE)
    con = df[mask].copy()
    parts = con["sample_alias"].str.extract(ALIAS_RE)
    con["ID"] = con["sample_alias"]
    con["site"] = "E" + parts[0]
    con["rep"] = "R" + parts[1]
    con["timepoint"] = "T" + parts[2]
    con["n_fastq"] = con["fastq_ftp"].apply(lambda x: len([u for u in x.split(";") if u]))
    return con.sort_values(["site", "timepoint", "rep"]).reset_index(drop=True)


def design_gaps(con: pd.DataFrame) -> list[tuple[str, str, str]]:
    have = set(zip(con["site"], con["timepoint"], con["rep"]))
    full = {(s, t, r) for s in SITES for t in TIMEPOINTS for r in REPS}
    return sorted(full - have)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="data/ena/PRJEB48362_read_run.tsv", type=Path)
    ap.add_argument("--outdir", default="data/metadata", type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = load(args.tsv)
    con = conwy_subset(df)
    gaps = design_gaps(con)

    # ---- assertions (the verification GATE) -----------------------------------------
    checks: list[tuple[str, bool, str]] = []
    n = len(con)
    checks.append(("Conwy 12S (E##_R#_T##) run count == 755", n == N_EXPECTED, f"{n}"))
    checks.append(("ID is unique (no alias maps to >1 run)",
                   con["ID"].is_unique, f"{con['ID'].nunique()} unique IDs"))
    checks.append(("run_accession is unique (1:1)",
                   con["run_accession"].is_unique, f"{con['run_accession'].nunique()} unique runs"))
    checks.append(("All runs are paired-end (2 FASTQ each)",
                   bool((con["n_fastq"] == 2).all()), f"{sorted(con['n_fastq'].unique())}"))
    checks.append(("Exactly 43 design cells absent",
                   len(gaps) == (N_THEORETICAL - N_EXPECTED), f"{len(gaps)} of {N_THEORETICAL}"))
    checks.append(("14 distinct sites E01-E14",
                   sorted(con["site"].unique()) == SITES, f"{con['site'].nunique()} sites"))

    print("=== ENA join verification: PRJEB48362 Conwy 12S =================================")
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:48s} ({detail})")

    # ---- artifact 1: clean join table -----------------------------------------------
    cols = ["ID", "site", "timepoint", "rep", "run_accession", "experiment_accession",
            "sample_accession", "n_fastq", "read_count", "base_count",
            "fastq_ftp", "fastq_md5", "fastq_bytes"]
    cols = [c for c in cols if c in con.columns]
    join_path = args.outdir / "conwy_12S_libraries.tsv"
    con[cols].to_csv(join_path, sep="\t", index=False)

    # ---- artifact 2: design gaps ----------------------------------------------------
    gap_df = pd.DataFrame(gaps, columns=["site", "timepoint", "rep"])
    gap_path = args.outdir / "design_gaps.tsv"
    gap_df.to_csv(gap_path, sep="\t", index=False)

    # ---- artifact 3: verification report (computed, reproducible) -------------------
    by_tp = Counter(t for _, t, _ in gaps)
    by_site = Counter(s for s, _, _ in gaps)
    total_bytes = con["fastq_bytes"].apply(lambda x: sum(int(b) for b in x.split(";") if b)).sum()
    total_reads = con["read_count"].apply(lambda x: int(x) if x else 0).sum()
    all_pass = all(ok for _, ok, _ in checks)

    md = []
    md.append("# ENA join verification — PRJEB48362 (Conwy 12S)\n")
    md.append("*Auto-generated by `data/verify_ena_join.py`. Re-run to refresh.*\n")
    md.append(f"- Source: `{args.tsv}` ({len(df)} total read-runs in the study)")
    md.append(f"- Conwy 12S subset (`E##_R#_T##`): **{n}** read-runs")
    md.append(f"- Overall gate: **{'PASS ✓' if all_pass else 'FAIL ✗'}**\n")
    md.append("## Checks\n")
    md.append("| Check | Result | Detail |")
    md.append("|---|---|---|")
    for name, ok, detail in checks:
        md.append(f"| {name} | {'PASS ✓' if ok else 'FAIL ✗'} | {detail} |")
    md.append("\n## Download footprint (Phase 2)\n")
    md.append(f"- FASTQ files: **{int(con['n_fastq'].sum())}** ({n} runs × 2, paired-end)")
    md.append(f"- Total size: **{total_bytes/1e9:.1f} GB**  ·  total reads: **{total_reads:,}**")
    md.append("\n## The 43 absent design cells (NOT random QC dropout)\n")
    md.append(f"- Theoretical design: {len(SITES)} sites × {len(TIMEPOINTS)} timepoints × "
              f"{len(REPS)} reps = **{N_THEORETICAL}**; deposited = **{n}**; absent = **{len(gaps)}**.")
    md.append(f"- Absent **by timepoint**: {dict(sorted(by_tp.items()))}")
    md.append(f"- Absent **by site**: {dict(sorted(by_site.items()))}")
    md.append("- **Interpretation:** timepoint **T16 is absent for all 14 sites** (42 cells) plus "
              "`E02_T17_R3` (1 cell). T16 is effectively a whole missing sampling occasion — so there "
              "are **18 usable timepoints, not 19**. This is structured, not scattered, missingness and "
              "must be handled deliberately when building occupancy detection histories (treat T16 as "
              "unsurveyed across all sites rather than as non-detections).")
    md.append("\n## Artifacts written\n")
    md.append(f"- `{join_path}` — clean ID → run → FASTQ join table ({n} rows)")
    md.append(f"- `{gap_path}` — the {len(gaps)} absent cells\n")
    (args.outdir / "ena_join_verification.md").write_text("\n".join(md) + "\n")

    print(f"\nWrote: {join_path}\n       {gap_path}\n       {args.outdir/'ena_join_verification.md'}")
    if not all_pass:
        print("\nGATE FAILED — do not proceed; investigate the failing check above.", file=sys.stderr)
        return 1
    print("\nGATE PASSED ✓")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
