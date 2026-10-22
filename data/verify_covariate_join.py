#!/usr/bin/env python3
"""
verify_covariate_join.py  --  Phase 0 / P3 gate: make the covariate join a
committed, BOM-safe, re-runnable assertion (the single most-cited integration
fact in the docs previously had no backing code).

Does four things:
  1. Reads supplementary_data_1_covariates.csv with encoding="utf-8-sig" so the
     UTF-8 BOM on the `ID` header (EF BB BF) is stripped -- a naive read.csv in R
     names the column `i..ID` and a join on `ID` silently yields ZERO matches.
  2. Asserts the 755 deposited Conwy libraries (conwy_12S_libraries.tsv) join 1:1
     to covariate rows on `ID` (755 matched / 0 left-only); the CSV is a superset
     holding the full 798-cell design grid, of which exactly 43 are unsequenced
     (cross-checks the ENA gate's 43 absent cells).
  3. Audits NA in the prereg's locked modelling covariates -- and reports the
     structure honestly. (Spoiler: 3 of 8 are NOT complete; see the report.)
  4. Writes a clean, BOM-free, Conwy-only modelling table conwy_covariates.tsv
     that R scripts consume instead of the raw publisher CSV.

Python (not R), consistent with verify_ena_join.py. Stdlib only.

Usage:  python3 data/verify_covariate_join.py    (run from repo root)
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

DESIGN_RE = re.compile(r"^E\d\d_R\d_T\d\d$")   # a Conwy design-grid replicate cell

# Prereg-locked modelling covariates (plan §2.2). Order = output column order.
PSI_COVS = ["Dist_from_lake", "river_gradient_percent", "river_section",
            "season", "Latitude", "Longitude"]
DET_COVS = ["week_before_flow_by_catchment", "Water_temp_oC"]
LOCKED = PSI_COVS + DET_COVS
# Resolutions to the locked-covariate gaps (dated 2026-06-19 in the project record;
# repository history does not independently establish pre-fit timing; see prereg §1):
#   temp_monthly_average  <- substitutes the 100%-NA Water_temp_oC
#   Catchment_flow_m3_s    -> on-hand per-timepoint catchment-flow series for E13/E14 (their
#                            *_by_catchment is NA) -- a time-series, NOT a scalar (2.53 is only T01)
EXTRA = ["temp_monthly_average", "Catchment_flow_m3_s", "daily_flow_by_catchment",
         "monthly_flow_by_catchment", "seasonal_flow_by_catchment",
         "rainfall_monthly_average", "pH", "conductivity_uS_cm"]
OUT_COLS = ["ID"] + LOCKED + EXTRA

NA_TOKENS = {"", "NA", "N/A", "na", "nan"}


def is_na(v: str) -> bool:
    return v.strip() in NA_TOKENS


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path,
                    default="data/metadata/supplementary_data_1_covariates.csv")
    ap.add_argument("--libraries", type=Path,
                    default="data/metadata/conwy_12S_libraries.tsv")
    ap.add_argument("--outdir", type=Path, default="data/metadata")
    args = ap.parse_args()
    if not args.csv.exists():
        sys.exit(f"ERROR: {args.csv} not found (publisher covariate CSV, git-ignored).")

    # ---- 1. BOM-safe read -----------------------------------------------------------
    with open(args.csv, "rb") as fb:
        bom = fb.read(3) == b"\xef\xbb\xbf"
    cov = {}
    with open(args.csv, encoding="utf-8-sig") as f:
        rd = csv.DictReader(f)
        assert rd.fieldnames and rd.fieldnames[0] == "ID", \
            f"first column is {rd.fieldnames[0]!r}, not 'ID' — BOM not stripped?"
        for r in rd:
            cov[r["ID"]] = r

    lib_ids = []
    with open(args.libraries) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            lib_ids.append(r["ID"])

    # ---- 2. join integrity ----------------------------------------------------------
    lib_set, cov_set = set(lib_ids), set(cov)
    matched = lib_set & cov_set
    left_only = lib_set - cov_set        # libraries with no covariate row (must be 0)
    # The CSV is a SUPERSET: it carries the full 798-cell Conwy design grid (plus
    # other rivers + controls). The design rows it holds that were NOT sequenced
    # must be exactly the 43 absent cells from the ENA gate (T16×42 + E02_T17_R3).
    design_rows = {i for i in cov_set if DESIGN_RE.match(i)}
    unsequenced = design_rows - lib_set

    # ---- 3. NA audit on the 755 ----------------------------------------------------
    def na_ids(col):
        return [i for i in lib_ids if col not in cov[i] or is_na(cov[i][col])]
    na_counts = {c: na_ids(c) for c in LOCKED + EXTRA}

    def status(c):
        n = len(na_counts[c])
        return "EMPTY" if n == len(lib_ids) else ("COMPLETE" if n == 0 else "PARTIAL")

    complete_locked = [c for c in LOCKED if status(c) == "COMPLETE"]
    incomplete_locked = [c for c in LOCKED if status(c) != "COMPLETE"]

    # ---- assertions (the GATE: join integrity + the BOM) ---------------------------
    checks = []
    checks.append(("BOM present on CSV and stripped (utf-8-sig)", bom, "EF BB BF" if bom else "no BOM"))
    checks.append(("All 755 libraries have a covariate row (0 left-only)",
                   len(left_only) == 0, f"{len(matched)} matched, {len(left_only)} left-only"))
    checks.append(("Join is exactly 755", len(matched) == 755, f"{len(matched)}"))
    checks.append(("CSV covers the 798 design grid; 43 unsequenced (= ENA gate)",
                   len(design_rows) == 798 and len(unsequenced) == 43,
                   f"{len(design_rows)} design rows, {len(unsequenced)} unsequenced"))
    join_ok = all(ok for _, ok, _ in checks)

    print("=== Covariate join verification: 755 Conwy 12S libraries =======================")
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name:54s} ({detail})")
    print("  " + "-" * 78)
    print("  Locked covariate completeness (NA over the 755):")
    for c in LOCKED:
        print(f"    [{status(c):8s}] {c:32s} {len(na_counts[c]):>4} NA")
    if incomplete_locked:
        print("  note: 3 locked covariates had gaps; the executed headline formula excludes the unresolved E14 fields. See report.")

    # ---- 4. clean, BOM-free modelling table ----------------------------------------
    out_tsv = args.outdir / "conwy_covariates.tsv"
    with open(out_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(OUT_COLS)
        for i in lib_ids:
            row = cov[i]
            w.writerow([i] + [("" if (c not in row or is_na(row[c])) else row[c].strip())
                              for c in OUT_COLS[1:]])

    # ---- verification report --------------------------------------------------------
    wb_by_tp = dict(sorted(Counter(x[-3:] for x in na_counts["week_before_flow_by_catchment"]).items()))
    wb_by_site = dict(sorted(Counter(x[:3] for x in na_counts["week_before_flow_by_catchment"]).items()))
    rg_by_site = dict(sorted(Counter(x[:3] for x in na_counts["river_gradient_percent"]).items()))

    # Catchment_flow_m3_s is a per-timepoint catchment series (constant across the 14 sites at
    # each timepoint), NOT a scalar -- derive its real structure from the data so the report can
    # never drift back to the documented-wrong "= 2.53" claim (eDNAplan.md D8).
    def _flow(i):
        v = cov[i].get("Catchment_flow_m3_s")
        return v.strip() if isinstance(v, str) and not is_na(v) else None
    flow_vals = sorted({float(_flow(i)) for i in lib_ids if _flow(i) is not None})
    flow_t01 = next((_flow(i) for i in lib_ids if i[-3:] == "T01" and _flow(i) is not None), None)
    flow_desc = (f"a per-timepoint catchment series ({len(flow_vals)} distinct values, "
                 f"~{flow_vals[0]:g}–{flow_vals[-1]:g} m³/s, constant across sites at each "
                 f"timepoint; {flow_t01} is only T01)") if flow_vals else \
                "the on-hand catchment-flow series (per timepoint, not a scalar)"

    md = []
    md.append("# Covariate join verification — 755 Conwy 12S libraries (P3 gate)\n")
    md.append("*Auto-generated by `data/verify_covariate_join.py`. Re-run to refresh.*\n")
    md.append(f"- Source: `{args.csv}` (publisher Supplementary Data 1, git-ignored)")
    md.append(f"- **UTF-8 BOM on the `ID` header: {'present (EF BB BF) — stripped via utf-8-sig' if bom else 'none'}.** "
              "In R, always read this with `readr::read_csv()` or "
              "`read.csv(..., fileEncoding = 'UTF-8-BOM')`; a naive `read.csv()` names the first "
              "column `ï..ID` and any join on `ID` silently returns zero rows.")
    md.append(f"- Join integrity: **{'PASS ✓' if join_ok else 'FAIL ✗'}**\n")
    md.append("## Join (`ID` = `E##_R#_T##`)\n")
    md.append("| Check | Result | Detail |")
    md.append("|---|---|---|")
    for name, ok, detail in checks:
        md.append(f"| {name} | {'PASS ✓' if ok else 'FAIL ✗'} | {detail} |")
    md.append(f"\nEvery deposited library has a covariate row; the CSV is a superset covering the "
              f"full 798-cell design grid, of which **{len(unsequenced)}** are unsequenced (= the ENA "
              "gate's 43 absent cells). The headline **755/755** integration fact is now backed by code.\n")

    md.append("## Locked-covariate completeness (NA over the 755)\n")
    md.append("| Covariate | Role | Status | NA |")
    md.append("|---|---|---|---|")
    for c in PSI_COVS:
        md.append(f"| `{c}` | ψ | {status(c)} | {len(na_counts[c])} |")
    for c in DET_COVS:
        md.append(f"| `{c}` | p | {status(c)} | {len(na_counts[c])} |")
    md.append("")
    md.append(f"- **Complete (0 NA):** {', '.join('`'+c+'`' for c in complete_locked)}")
    md.append(f"- **Incomplete:** {', '.join('`'+c+'`' for c in incomplete_locked)}\n")

    md.append("## Locked-covariate gaps — project-recorded decisions (dated 2026-06-19; timing limitation applies)\n")
    md.append("Supplementary Data 1 does not fully supply 3 historically locked covariates. The project "
              "record documents the following resolutions, but git history does not independently establish "
              "that they preceded every model fit. The executed analysis uses no silent imputation:\n")
    md.append(f"1. **`Water_temp_oC` is 100% NA ({len(na_counts['Water_temp_oC'])}/755)** → **RESOLVED: "
              "substitute `temp_monthly_average`** "
              f"({len(na_counts['temp_monthly_average'])} NA, all E14) — the column the authors actually "
              "used. Carried in `conwy_covariates.tsv`; prereg detection model updated.")
    md.append(f"2. **`week_before_flow_by_catchment`: {len(na_counts['week_before_flow_by_catchment'])} NA** "
              f"(and daily/monthly/seasonal flow) — by site {wb_by_site}: only the two lowest sites "
              f"**E13, E14**. → **RESOLVED: use the on-hand `Catchment_flow_m3_s`** for these sites "
              f"(carried in `conwy_covariates.tsv`) — {flow_desc}, *not* a single scalar.")
    md.append(f"3. **Site E14 (Trefriw)** has NA `river_gradient_percent` ({len(na_counts['river_gradient_percent'])}, "
              f"by site {rg_by_site}), `temp_monthly_average`, and `rainfall_monthly_average`, and source values "
              "are **not in any file on hand**. Final local disposition: a separately generated table uses an "
              "explicitly flagged, unvalidated E13 donor for supporting climate values; gradient remains NA. "
              "The executed headline occupancy formula excludes climate, flow, and gradient "
              "(`Dist_from_lake + season`), so these gaps do not block it.\n")

    md.append("## Artifacts written\n")
    md.append(f"- `{out_tsv}` — BOM-free Conwy modelling table: `ID` + {len(LOCKED)} locked + "
              f"{len(EXTRA)} documented-alternative covariates, 755 rows. Empty cells are blank "
              "(true NA). **R scripts read this, not the raw publisher CSV.**")
    md.append(f"- `{args.outdir / 'covariate_join_verification.md'}` — this report.\n")
    (args.outdir / "covariate_join_verification.md").write_text("\n".join(md) + "\n")

    print(f"\nWrote: {out_tsv}\n       {args.outdir/'covariate_join_verification.md'}")
    if not join_ok:
        print("\nGATE FAILED — join integrity broken; investigate before Phase 1.", file=sys.stderr)
        return 1
    print("\nGATE PASSED ✓  — 755/755 join sound (BOM-safe). "
          "E14 source gaps are explicit and excluded from the executed headline formula (see report).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
