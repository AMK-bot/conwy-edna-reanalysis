#!/usr/bin/env python3
"""Run the complete offline Phase-5 scientific-impact workflow.

Run this same command twice.  The first run records canonical output hashes; the
second refuses to pass unless every deterministic scientific artifact is
byte-identical.  Use ``--release`` after reader-facing edits to render and
verify the complete portable release as well.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "data" / "metadata"
HASHES = META / "impact_sprint_canonical_hashes.tsv"
RECEIPT = META / "impact_sprint_run_receipt.json"
EXPECTED_PREREG = "5d4ab9ce7d1f0dc695646f6e3c781b4931fad5fecc1513ddd79046c3b316e567"
EXPECTED_MANIFEST = "0f6e6464ec54492cb23b8adbdc81b3b05a9ec98e73e1b04b413c764c03194843"

CANONICAL = [
    # Shared portable input.
    "data/derived/impact_library_counts.tsv",
    "data/metadata/impact_library_export_receipt.tsv",
    # Synthetic occupancy atlas (telemetry receipt excluded).
    "data/derived/impact_occupancy_scenarios.tsv",
    "data/derived/impact_occupancy_simulation_audit.tsv",
    "data/derived/impact_occupancy_fit_results.tsv",
    "data/derived/impact_occupancy_atlas.tsv",
    "data/derived/impact_occupancy_generator_fixture.tsv",
    "data/derived/impact_occupancy_manifest.json",
    "figures/impact_occupancy_atlas.png",
    "figures/impact_occupancy_atlas.svg",
    # Cached, accession-backed Salmo benchmark.
    "data/metadata/impact_salmo_accessions.tsv",
    "data/metadata/impact_salmo_marker_manifest.tsv",
    "data/metadata/impact_salmo_panel_coverage.tsv",
    "data/metadata/impact_salmo_queries.tsv",
    "data/derived/impact_salmo_pairwise_distances.tsv",
    "data/derived/impact_salmo_exact_haplotypes.tsv",
    "data/derived/impact_salmo_barcode_gaps.tsv",
    "data/derived/impact_salmo_loo_classification.tsv",
    "data/derived/impact_salmo_marker_gates.tsv",
    "figures/impact_salmo_multimarker.png",
    "figures/impact_salmo_multimarker.svg",
    # Layer-E multiverse (runtime/determinism receipts excluded).
    "data/derived/impact_layerE_specifications.tsv",
    "data/derived/impact_layerE_library_draws.tsv",
    "data/derived/impact_layerE_draw_results.tsv",
    "data/derived/impact_layerE_results.tsv",
    "data/derived/impact_layerE_influence.tsv",
    "data/derived/impact_layerE_gate.tsv",
    "figures/impact_layerE_multiverse.png",
    "figures/impact_layerE_multiverse.svg",
    # Trait fixed-margin null.
    "data/derived/impact_trait_gower.tsv",
    "data/derived/impact_trait_unit_rao.tsv",
    "data/derived/impact_trait_null_draws.tsv",
    "data/derived/impact_trait_null_summary.tsv",
    "data/metadata/impact_trait_null_receipt.tsv",
    "figures/impact_trait_null_model.png",
    "figures/impact_trait_null_model.svg",
    # Cluster-aware pipeline disagreement.
    "data/metadata/impact_pipeline_disagreement_cells.tsv",
    "data/derived/impact_pipeline_model.tsv",
    "data/derived/impact_pipeline_bootstrap.tsv",
    "data/derived/impact_pipeline_leaveout.tsv",
    "data/derived/impact_pipeline_rates.tsv",
    "data/derived/impact_pipeline_summary.tsv",
    "data/metadata/impact_pipeline_receipt.tsv",
    "figures/impact_pipeline_disagreement.png",
    "figures/impact_pipeline_disagreement.svg",
    # Checked synthesis.
    "data/metadata/claim_survival.tsv",
    "website/claim-survival.html",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def frozen_input_gate() -> None:
    actual = {
        "docs/impact-sprint-preregistration.md": sha256(ROOT / "docs/impact-sprint-preregistration.md"),
        "data/metadata/impact_sprint_manifest.json": sha256(ROOT / "data/metadata/impact_sprint_manifest.json"),
    }
    expected = {
        "docs/impact-sprint-preregistration.md": EXPECTED_PREREG,
        "data/metadata/impact_sprint_manifest.json": EXPECTED_MANIFEST,
    }
    if actual != expected:
        raise SystemExit(f"Frozen Phase-5 specification hash mismatch: {actual}")


def resolve_executable(explicit: str | None, env_name: str, name: str,
                       sibling_of: Path | None = None) -> str:
    candidates: list[str] = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get(env_name):
        candidates.append(os.environ[env_name])
    if shutil.which(name):
        candidates.append(shutil.which(name) or "")
    if sibling_of:
        candidates.append(str(sibling_of.parent / name))
    for candidate in candidates:
        path = Path(candidate).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            return str(path.resolve())
    raise SystemExit(f"Could not find {name}; pass --{name.lower()} or set {env_name}")


def run(label: str, command: list[str], env: dict[str, str] | None = None) -> None:
    print(f"\n[Phase 5] {label}\n$ {' '.join(command)}", flush=True)
    completed = subprocess.run(command, cwd=ROOT, env=env, check=False)
    if completed.returncode:
        raise SystemExit(f"{label} failed with exit code {completed.returncode}")


def read_previous_hashes() -> dict[str, str] | None:
    if not HASHES.exists():
        return None
    with HASHES.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"artifact", "sha256", "bytes"}:
        raise SystemExit("Existing Phase-5 canonical hash manifest has an invalid schema")
    return {row["artifact"]: row["sha256"] for row in rows}


def collect_hashes() -> dict[str, str]:
    missing = [rel for rel in CANONICAL if not (ROOT / rel).is_file()]
    if missing:
        raise SystemExit("Missing canonical Phase-5 artifacts: " + ", ".join(missing))
    return {rel: sha256(ROOT / rel) for rel in CANONICAL}


def write_hashes(hashes: dict[str, str]) -> None:
    with HASHES.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["artifact", "sha256", "bytes"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for rel in CANONICAL:
            writer.writerow({"artifact": rel, "sha256": hashes[rel],
                             "bytes": (ROOT / rel).stat().st_size})


def run_science(rscript: str, validate_only: bool) -> None:
    if validate_only:
        run("validate Salmo benchmark", [rscript, "R/12_salmo_multimarker_benchmark.R", "--validate-only"])
        run("validate trait null", [rscript, "R/14_trait_null_model.R", "--validate-only"])
        run("validate pipeline disagreement", [rscript, "R/15_pipeline_disagreement.R", "--validate-only"])
    else:
        run("export portable impact input", [rscript, "scripts/export_impact_inputs.R"])
        atlas_env = dict(os.environ)
        atlas_env.setdefault("IMPACT_ATLAS_CORES", "4")
        run("build occupancy recoverability atlas",
            [rscript, "R/11_occupancy_recoverability_atlas.R", "--verify-rebuild"],
            env=atlas_env)
        run("build cached Salmo benchmark", [rscript, "R/12_salmo_multimarker_benchmark.R"])
        run("build Layer-E multiverse", [rscript, "R/13_layerE_multiverse.R"])
        run("build trait fixed-margin null", [rscript, "R/14_trait_null_model.R"])
        run("build cluster-aware pipeline disagreement", [rscript, "R/15_pipeline_disagreement.R"])

    run("build checked claim-survival synthesis", [sys.executable, "scripts/build_claim_survival.py"])
    for label, test in (
        ("test occupancy atlas", "R/test_impact_atlas.R"),
        ("test Salmo benchmark", "R/test_salmo_multimarker.R"),
        ("test Layer-E multiverse", "R/test_layerE_multiverse.R"),
        ("test trait null", "R/test_trait_null_model.R"),
        ("test pipeline disagreement", "R/test_pipeline_disagreement.R"),
    ):
        run(label, [rscript, test])


def render_release(rscript: str, quarto_arg: str | None, pandoc_arg: str | None) -> None:
    quarto = resolve_executable(quarto_arg, "EDNA_QUARTO", "quarto", Path(rscript))
    pandoc = resolve_executable(pandoc_arg, "EDNA_PANDOC", "pandoc", Path(rscript))
    quarto_env = dict(os.environ)
    qbin = Path(quarto).resolve().parent
    conda_root = qbin.parent
    # The conda-forge Quarto package stores its tools beside the wrapper rather
    # than below bin/tools. Supply the official override variables when that
    # layout is detected; system/app Quarto installations need no overrides.
    overrides = {
        "QUARTO_SHARE_PATH": conda_root / "share" / "quarto",
        "QUARTO_DENO": qbin / "deno",
        "QUARTO_DENO_DOM": conda_root / "lib" / "deno_dom.dylib",
        "QUARTO_PANDOC": Path(pandoc),
        "QUARTO_DART_SASS": qbin / "sass",
        "QUARTO_ESBUILD": qbin / "esbuild",
    }
    if all(path.exists() for path in overrides.values()):
        quarto_env.update({key: str(path) for key, path in overrides.items()})
    quarto_env["PATH"] = os.pathsep.join(
        [str(Path(rscript).resolve().parent), str(qbin), quarto_env.get("PATH", "")]
    )
    run("render technical report", [quarto, "render", "report"], env=quarto_env)
    run("render self-contained research brief", [
        pandoc, "website/research_brief.md", "--standalone", "--embed-resources",
        "--resource-path=website", "--css=website/research_brief.css",
        "--output=website/research_brief.html",
    ])
    run("build website figures and PDFs", [sys.executable, "scripts/build_website.py"])


def package_and_verify_release() -> None:
    run("build and clean-extract portable website", [sys.executable, "scripts/package_website.py"])
    run("verify complete release", [sys.executable, "scripts/verify_release.py"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", help="Path to the Rscript executable")
    parser.add_argument("--quarto", help="Path to Quarto (used with --release)")
    parser.add_argument("--pandoc", help="Path to Pandoc (used with --release)")
    parser.add_argument("--validate-only", action="store_true",
                        help="Run validators/tests without rebuilding the five analyses")
    parser.add_argument("--release", action="store_true",
                        help="Also render, package and verify the complete offline release")
    args = parser.parse_args()

    started = time.monotonic()
    frozen_input_gate()
    previous = read_previous_hashes()
    rscript = resolve_executable(args.rscript, "EDNA_RSCRIPT", "Rscript")
    run_science(rscript, args.validate_only)
    current = collect_hashes()
    mismatches: list[dict[str, str]] = []
    determinism = "NOT_EVALUATED_FIRST_RUN"
    if previous is not None:
        for rel in CANONICAL:
            if previous.get(rel) != current[rel]:
                mismatches.append({"artifact": rel, "previous": previous.get(rel, "MISSING"),
                                   "current": current[rel]})
        determinism = "PASS" if not mismatches else "FAIL"
    write_hashes(current)

    if args.release and not mismatches:
        render_release(rscript, args.quarto, args.pandoc)

    receipt = {
        "status": "PASS" if not mismatches else "FAIL",
        "scope": "offline Phase-5 scientific-impact workflow",
        "frozen_specification": "PASS",
        "analysis_mode": "validate-only" if args.validate_only else "full-rebuild",
        "release_mode": args.release,
        "rscript": rscript,
        "canonical_artifacts": len(CANONICAL),
        "deterministic_against_previous_run": determinism,
        "mismatches": mismatches,
        "runtime_seconds": round(time.monotonic() - started, 3),
        "runtime_scope": (
            "scientific rebuild and release rendering before packaging; packaging and "
            "verification excluded so this finalized receipt can be embedded in the archive"
            if args.release else "scientific rebuild and tests"
        ),
        "receipt_finalized_before_packaging": True,
        "expected_ordinary_laptop_runtime": "approximately 4-7 minutes per full scientific rebuild; release rendering adds about 1 minute",
        "claim_boundary": "Synthetic occupancy results do not estimate Conwy occupancy; public sequence labels are annotations; pipeline disagreement is not a biological false-negative rate.",
    }
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")

    if args.release and not mismatches:
        package_and_verify_release()

    print(json.dumps(receipt, indent=2))
    return 0 if not mismatches else 1


if __name__ == "__main__":
    raise SystemExit(main())
