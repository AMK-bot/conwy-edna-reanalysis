#!/usr/bin/env python3
"""
build_occasions.py  --  Phase 0 / P5: emit the canonical occasion map so the
"18 usable timepoints; T16 = unsurveyed (NA, NOT 0)" rule is a machine-readable
artifact, not prose that gets re-derived from memory weeks later (and mishandled
as all-zero detections, which would bias occupancy ψ downward — the exact error
this project exists to expose).

Emits data/metadata/occasions.tsv: one row per cell of the full 14×19×3 = 798
design grid, flagged so detection-history construction can mark NA correctly at
two distinct levels of missingness:

  surveyed            : a library exists (ENA). FALSE => UNSURVEYED design cell
                        => detection-history cell is NA, never 0. (43 cells:
                        T16×42 + E02_T17_R3.)
  in_published_12S    : the library survived the authors' <1000-read QC into the
                        published 665-row matrix. surveyed & !in_published => in
                        the FAST-PATH (published-matrix) history the cell is NA
                        (sequenced but no usable call); Phase-2 reprocessing may
                        recover it. (90 cells.)
  usable_timepoint    : timepoint != T16 (T16 is unsurveyed for ALL sites, so it
                        is a whole missing occasion → 18 usable, not 19).

Stdlib only; consistent with verify_ena_join.py. Run from repo root.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

SITES = [f"E{i:02d}" for i in range(1, 15)]
TIMEPOINTS = [f"T{i:02d}" for i in range(1, 20)]
REPS = ["R1", "R2", "R3"]
UNUSABLE_TP = {"T16"}   # absent for all 14 sites → not a usable occasion


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--libraries", type=Path, default="data/metadata/conwy_12S_libraries.tsv")
    ap.add_argument("--published", type=Path, default="data/metadata/source_data_12S_normalised.tsv")
    ap.add_argument("--outdir", type=Path, default="data/metadata")
    args = ap.parse_args()

    surveyed = set()
    with open(args.libraries) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            surveyed.add(r["ID"])

    published = set()
    if args.published.exists():
        with open(args.published) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                published.add(r["ID"])

    rows = []
    for s in SITES:
        for t in TIMEPOINTS:
            for rep in REPS:
                cid = f"{s}_{rep}_{t}"
                rows.append({
                    "site": s, "timepoint": t, "rep": rep, "ID": cid,
                    "surveyed": int(cid in surveyed),
                    "in_published_12S": int(cid in published),
                    "usable_timepoint": int(t not in UNUSABLE_TP),
                })

    out = args.outdir / "occasions.tsv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    n_surv = sum(r["surveyed"] for r in rows)
    n_pub = sum(r["in_published_12S"] for r in rows)
    usable_tps = sorted({r["timepoint"] for r in rows if r["usable_timepoint"]})
    n_unsurv = len(rows) - n_surv

    print("=== Occasion map (P5) ==========================================================")
    print(f"  design grid cells (14×19×3)         : {len(rows)}")
    print(f"  surveyed (library exists)           : {n_surv}")
    print(f"  UNSURVEYED → NA, never 0            : {n_unsurv}  (T16×42 + E02_T17_R3)")
    print(f"  in published 12S matrix (665)       : {n_pub}")
    print(f"  surveyed-but-QC-dropped → NA in     : {n_surv - n_pub}  (fast-path history only)")
    print(f"    the fast-path history")
    print(f"  usable timepoints                   : {len(usable_tps)}  {usable_tps}")
    assert len(usable_tps) == 18, "expected 18 usable timepoints (T16 excluded)"
    assert n_unsurv == 43, "expected 43 unsurveyed design cells"
    print(f"\nWrote: {out}")
    print("\nGATE PASSED ✓  — 18 usable occasions; 43 cells flagged unsurveyed (NA, not 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
