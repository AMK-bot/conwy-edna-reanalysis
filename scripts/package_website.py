#!/usr/bin/env python3
"""Build and verify a portable, offline ZIP of the project website.

The archive preserves the repository-relative paths used by the static website,
collects every local file linked by its HTML pages, adds small unlinked website
inputs, and verifies the result after a clean extraction.  It deliberately does
not include raw FASTQs, heavyweight model checkpoints, a hosted service, or a
claim that the bundle is a complete clean-clone analysis environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from html.parser import HTMLParser
from pathlib import Path, PurePosixPath
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
BUNDLE_SLUG = "conwy-edna-project-website"
RELEASE_DATE = "2026-07-19"
DEFAULT_OUTPUT = ROOT / "output" / "share" / f"{BUNDLE_SLUG}.zip"

ENTRY_HTML = (
    "website/index.html",
    "website/editor.html",
    "website/deliverables.html",
    "website/research_brief.html",
    "website/visual-brief.html",
    "website/detection-lab.html",
    "website/claim-survival.html",
    "website/poster.html",
    "_site/conwy-reanalysis.html",
)

# Phase-5 files that define, execute, explain, and expose the six impact-sprint
# extensions.  The broader globs below also retain every machine table and run
# receipt, while this explicit core prevents a partially generated sprint from
# being packaged as complete.
PHASE5_REQUIRED_FILES = (
    "docs/impact-sprint-preregistration.md",
    "data/metadata/impact_sprint_manifest.json",
    "data/metadata/impact_sprint_preregistration_receipt.tsv",
    "data/metadata/impact_sprint_canonical_hashes.tsv",
    "data/metadata/impact_sprint_run_receipt.json",
    "scripts/run_impact_sprint.py",
    "scripts/export_impact_inputs.R",
    "scripts/fetch_salmo_multimarker.py",
    "scripts/build_claim_survival.py",
    "R/11_occupancy_recoverability_atlas.R",
    "R/12_salmo_multimarker_benchmark.R",
    "R/13_layerE_multiverse.R",
    "R/14_trait_null_model.R",
    "R/15_pipeline_disagreement.R",
    "R/test_impact_atlas.R",
    "R/test_salmo_multimarker.R",
    "R/test_layerE_multiverse.R",
    "R/test_trait_null_model.R",
    "R/test_pipeline_disagreement.R",
    "docs/impact_occupancy_atlas_methods.md",
    "docs/impact_salmo_multimarker_methods.md",
    "docs/impact_layerE_multiverse_methods.md",
    "docs/impact_trait_null_methods.md",
    "docs/impact_pipeline_disagreement_methods.md",
    "data/metadata/claim_survival.tsv",
    "website/claim-survival.html",
    "website/assets/claim-survival.png",
    "data/derived/impact_occupancy_scenarios.tsv",
    "data/derived/impact_occupancy_simulation_audit.tsv",
    "data/derived/impact_occupancy_fit_results.tsv",
    "data/derived/impact_occupancy_atlas.tsv",
    "data/derived/impact_occupancy_generator_fixture.tsv",
    "data/derived/impact_occupancy_manifest.json",
    "data/derived/impact_occupancy_receipt.json",
    "data/derived/impact_occupancy_test_receipt.tsv",
    "data/metadata/impact_salmo_queries.tsv",
    "data/metadata/impact_salmo_accessions.tsv",
    "data/metadata/impact_salmo_marker_manifest.tsv",
    "data/metadata/impact_salmo_panel_coverage.tsv",
    "data/derived/impact_salmo_exact_haplotypes.tsv",
    "data/derived/impact_salmo_pairwise_distances.tsv",
    "data/derived/impact_salmo_barcode_gaps.tsv",
    "data/derived/impact_salmo_loo_classification.tsv",
    "data/derived/impact_salmo_marker_gates.tsv",
    "data/public/salmo_multimarker/README.md",
    "data/public/salmo_multimarker/marker_sequences.fasta",
    "data/derived/impact_library_counts.tsv",
    "data/metadata/impact_library_export_receipt.tsv",
    "data/derived/impact_layerE_specifications.tsv",
    "data/derived/impact_layerE_library_draws.tsv",
    "data/derived/impact_layerE_draw_results.tsv",
    "data/derived/impact_layerE_results.tsv",
    "data/derived/impact_layerE_influence.tsv",
    "data/derived/impact_layerE_gate.tsv",
    "data/metadata/impact_layerE_run_receipt.tsv",
    "data/metadata/impact_layerE_determinism_receipt.tsv",
    "data/derived/impact_trait_gower.tsv",
    "data/derived/impact_trait_unit_rao.tsv",
    "data/derived/impact_trait_null_draws.tsv",
    "data/derived/impact_trait_null_summary.tsv",
    "data/metadata/impact_trait_null_receipt.tsv",
    "data/metadata/impact_pipeline_disagreement_cells.tsv",
    "data/derived/impact_pipeline_bootstrap.tsv",
    "data/derived/impact_pipeline_leaveout.tsv",
    "data/derived/impact_pipeline_model.tsv",
    "data/derived/impact_pipeline_rates.tsv",
    "data/derived/impact_pipeline_summary.tsv",
    "data/metadata/impact_pipeline_receipt.tsv",
    "figures/impact_occupancy_atlas.png",
    "figures/impact_occupancy_atlas.svg",
    "figures/impact_salmo_multimarker.png",
    "figures/impact_salmo_multimarker.svg",
    "figures/impact_layerE_multiverse.png",
    "figures/impact_layerE_multiverse.svg",
    "figures/impact_trait_null_model.png",
    "figures/impact_trait_null_model.svg",
    "figures/impact_pipeline_disagreement.png",
    "figures/impact_pipeline_disagreement.svg",
)

PHASE5_GLOBS = (
    "data/derived/impact_*",
    "data/metadata/impact_*",
    "docs/impact*",
    "figures/impact_*",
)

PHASE5_TREES = ("data/public/salmo_multimarker",)

# These are small presentation inputs or fallbacks used to understand/rebuild a
# linked deliverable but are not all referenced through an HTML href/src.
EXTRA_FILES = (
    "EDIT_WEBSITE.command",
    "EDIT_WEBSITE.md",
    "scripts/edit_website.py",
    "website/content-overrides.json",
    "website/research_brief.md",
    "website/research_brief.css",
    "website/graphical_abstract.svg",
    "website/assets/research-timeline.png",
    "website/assets/models-overview.png",
    "website/assets/model-community-multiverse.png",
    "website/assets/model-occupancy-atlas.png",
    "website/assets/model-pipeline-bootstrap.png",
    "website/assets/model-salmo-benchmark.png",
    "website/assets/model-trait-null.png",
    "website/data/community-gradient.csv",
    "website/data/salmo-calls.csv",
    "figures/salmo_molecular_detective.png",
    "data/metadata/fish_trait_analysis_groups.tsv",
    "data/metadata/salmo_haplotype_edges.tsv",
    "data/metadata/salmo_haplotype_nodes.tsv",
    "data/derived/layerE_alpha_site_time.tsv",
    "data/derived/reclassB_partB_provenance.txt",
    "data/derived/reclassB_partC_provenance.txt",
    "figures/phase3_occupancy_cov_corrected_vs_raw.png",
    "figures/phase3_occupancy_cov_design.png",
    "figures/phase3_occupancy_cov_dist_effect.png",
    "figures/phase3_occupancy_detection_fnr.png",
    "figures/phase3_occupancy_psi_forest.png",
    "R/09_e14_covariates.R",
    "report/_quarto.yml",
    "scripts/apply_layerA_claim_gate.R",
    "scripts/package_website.py",
    *PHASE5_REQUIRED_FILES,
)

FORBIDDEN_PARTS = {
    ".git",
    ".venv",
    ".venv-website",
    "__pycache__",
    ".DS_Store",
}
FORBIDDEN_SUFFIXES = {".fastq", ".fq", ".rds"}
NETWORK_URL_RE = re.compile(r"url\(\s*['\"]?(?:https?:)?//", re.IGNORECASE)
NETWORK_HREF_RE = re.compile(
    r"(?:href|xlink:href)\s*=\s*['\"](?:https?:)?//", re.IGNORECASE
)


class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids: set[str] = set()
        self.links: list[tuple[str, str]] = []
        self.network_dependencies: list[str] = []

    def handle_starttag(self, tag, attrs):
        attr = dict(attrs)
        if attr.get("id"):
            self.ids.add(attr["id"])
        for key in ("href", "src"):
            if attr.get(key):
                raw = attr[key]
                self.links.append((key, raw))
                parsed = urlsplit(raw)
                dependency = key == "src" or (
                    tag == "link"
                    and any(
                        token in (attr.get("rel") or "").lower().split()
                        for token in ("stylesheet", "preload", "modulepreload", "icon")
                    )
                )
                if dependency and (parsed.scheme or parsed.netloc) and parsed.scheme != "data":
                    self.network_dependencies.append(raw)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_html(path: Path) -> LinkParser:
    parser = LinkParser()
    parser.feed(path.read_text(encoding="utf-8"))
    return parser


def has_network_asset_reference(path: Path) -> bool:
    if path.suffix.lower() not in {".css", ".svg"}:
        return False
    text = path.read_text(encoding="utf-8")
    return bool(
        NETWORK_URL_RE.search(text)
        or NETWORK_HREF_RE.search(text)
        or (path.suffix.lower() == ".css" and re.search(r"@import\b", text, re.IGNORECASE))
    )


def local_target(source: Path, raw: str, root: Path) -> tuple[Path | None, str]:
    parsed = urlsplit(raw)
    if parsed.scheme or parsed.netloc or raw.startswith(("mailto:", "javascript:")):
        return None, ""
    target = source if not parsed.path else source.parent / unquote(parsed.path)
    target = target.resolve()
    try:
        target.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"link escapes bundle root: {source}: {raw}") from exc
    return target, parsed.fragment


def collect_files() -> tuple[list[Path], set[str]]:
    pending = [ROOT / rel for rel in ENTRY_HTML]
    files = {ROOT / rel for rel in EXTRA_FILES}
    for pattern in PHASE5_GLOBS:
        files.update(path for path in ROOT.glob(pattern) if path.is_file())
    for rel in PHASE5_TREES:
        tree = ROOT / rel
        if not tree.is_dir():
            raise FileNotFoundError(f"missing Phase-5 public-input cache: {rel}")
        files.update(path for path in tree.rglob("*") if path.is_file())
    external_links: set[str] = set()
    parsed_html: set[Path] = set()

    while pending:
        page = pending.pop()
        page = page.resolve()
        if page in parsed_html:
            continue
        if not page.is_file():
            raise FileNotFoundError(f"missing website page: {page.relative_to(ROOT)}")
        parsed_html.add(page)
        files.add(page)
        parser = parse_html(page)
        if parser.network_dependencies:
            raise ValueError(
                f"network-dependent website asset: {page.relative_to(ROOT)} -> "
                f"{sorted(set(parser.network_dependencies))}"
            )
        for _, raw in parser.links:
            parsed = urlsplit(raw)
            if parsed.scheme or parsed.netloc:
                if parsed.scheme not in {"data"}:
                    external_links.add(raw)
                continue
            target, _ = local_target(page, raw, ROOT)
            if target is None:
                continue
            if not target.exists():
                raise FileNotFoundError(
                    f"broken website dependency: {page.relative_to(ROOT)} -> {raw}"
                )
            if not target.is_file():
                raise ValueError(
                    f"website dependency is not a file: {page.relative_to(ROOT)} -> {raw}"
                )
            files.add(target)
            if target.suffix.lower() in {".html", ".htm"}:
                pending.append(target)

    clean_files: list[Path] = []
    casefolded_paths: dict[str, Path] = {}
    for source in sorted(files):
        rel = source.resolve().relative_to(ROOT.resolve())
        if source.is_symlink():
            raise ValueError(f"symbolic-link source is not allowed: {rel}")
        if FORBIDDEN_PARTS & set(rel.parts):
            raise ValueError(f"forbidden path requested by website: {rel}")
        suffixes = {suffix.lower() for suffix in source.suffixes}
        if FORBIDDEN_SUFFIXES & suffixes:
            raise ValueError(f"heavy checkpoint/raw file requested by website: {rel}")
        if not source.is_file():
            raise FileNotFoundError(f"missing bundle source: {rel}")
        if has_network_asset_reference(source):
            raise ValueError(f"network-dependent CSS/SVG asset: {rel}")
        folded = rel.as_posix().casefold()
        if folded in casefolded_paths and casefolded_paths[folded] != rel:
            raise ValueError(
                f"case-insensitive path collision: {casefolded_paths[folded]} and {rel}"
            )
        casefolded_paths[folded] = rel
        clean_files.append(source)
    return clean_files, external_links


def gateway_html() -> str:
    return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta http-equiv="refresh" content="0; url=website/index.html">
  <title>Open River Conwy fish-eDNA project website</title>
  <style>
    body { max-width: 48rem; margin: 12vh auto; padding: 2rem; color: #173f3a;
      background: #f5f0e6; font: 18px/1.55 system-ui, sans-serif; }
    a { color: #0b6658; font-weight: 700; }
  </style>
</head>
<body>
  <h1>River Conwy fish-eDNA reanalysis</h1>
  <p>Opening the offline project website…</p>
  <p><a href="website/index.html">Continue to the project summary</a></p>
</body>
</html>
"""


def readme_text() -> str:
    return """RIVER CONWY FISH-eDNA PROJECT WEBSITE — READ ME FIRST

1. Extract the entire ZIP file. Do not send or move an individual HTML file by itself.
2. Open OPEN_WEBSITE.html (or index.html) in a current web browser.
3. The project summary, claim-survival map, self-contained research brief, figures,
   simulator, report, PDFs, evidence, cached public-sequence inputs, and code work locally.
   External citation and reference links require an internet connection, but none is
   needed for the bundled evidence or local deliverables.
4. To personalize selected narrative text, double-click EDIT_WEBSITE.command (or run
   python3 scripts/edit_website.py), then follow EDIT_WEBSITE.md. Python is needed
   for saving edits, but not for viewing the website or using the synthetic simulator.

This is a static website/review bundle, not a hosted service, peer-reviewed publication,
or complete clean-clone raw-read/MCMC environment. The rare-label recovery result remains
withheld because its convergence gate failed.

Integrity check (optional):
    python3 VERIFY_BUNDLE.py

Reuse and citation:
- Project reuse status is stated in LICENSE (all rights reserved).
- CITATION.cff records citation metadata and the source study.
- THIRD_PARTY_NOTICES.txt records source-study attribution and bundled dependency boundaries.
- Source-study and other third-party materials retain their own terms and provenance.
"""


def third_party_notices_text() -> str:
    return """THIRD-PARTY NOTICES

Source study and supplementary information
------------------------------------------
Perry, W. B. et al. (2024). "An integrated spatio-temporal view of riverine
biodiversity using environmental DNA metabarcoding." Nature Communications 15,
4372. https://doi.org/10.1038/s41467-024-48640-3

The source article is recorded by this project as Creative Commons Attribution
4.0 (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/
The bundled docs/supplementary_information.pdf is an unmodified supporting copy
and is kept separate from project-authored conclusions.

Embedded report libraries
-------------------------
The self-contained Quarto technical report embeds browser libraries and fonts.
Their original embedded notices and headers are preserved; this packager does not
strip or minify them further.

Reference-derived sequence evidence
-----------------------------------
Reference accessions, classifier comparisons, and database-derived sequences retain
their upstream terms. See data/metadata/impact_salmo_accessions.tsv,
data/metadata/impact_salmo_queries.tsv, data/metadata/salmo_reference_provenance.tsv,
and data/metadata/salmo_classifier_provenance.tsv for the exact evidence provenance.

No notice in this file changes the repository-level reuse status in LICENSE.
"""


def verifier_wrapper() -> str:
    return """#!/usr/bin/env python3
import sys
from pathlib import Path

sys.dont_write_bytecode = True
from scripts.package_website import verify_extracted

result = verify_extracted(Path(__file__).resolve().parent, verify_manifest=True)
print(__import__("json").dumps(result, indent=2))
"""


def copy_sources(files: list[Path], destination: Path) -> None:
    for source in files:
        rel = source.resolve().relative_to(ROOT.resolve())
        target = destination / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)


def write_manifest(bundle_root: Path) -> int:
    manifest_path = bundle_root / "SHA256SUMS.txt"
    rows = []
    for path in sorted(item for item in bundle_root.rglob("*") if item.is_file()):
        if path == manifest_path:
            continue
        rel = path.relative_to(bundle_root).as_posix()
        rows.append(f"{sha256(path)}  {rel}")
    manifest_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return len(rows)


def verify_manifest_file(bundle_root: Path, failures: list[str]) -> int:
    manifest = bundle_root / "SHA256SUMS.txt"
    if not manifest.is_file():
        failures.append("missing SHA256SUMS.txt")
        return 0
    checked = 0
    listed: set[str] = set()
    for number, line in enumerate(manifest.read_text(encoding="utf-8").splitlines(), 1):
        if not line.strip():
            continue
        try:
            expected, rel = line.split("  ", 1)
        except ValueError:
            failures.append(f"malformed SHA256SUMS.txt line {number}")
            continue
        listed.add(rel)
        path = bundle_root / rel
        try:
            path.resolve().relative_to(bundle_root.resolve())
        except ValueError:
            failures.append(f"manifest path escapes bundle: {rel}")
            continue
        if not path.is_file():
            failures.append(f"manifest file missing: {rel}")
            continue
        actual = sha256(path)
        if actual != expected:
            failures.append(f"checksum mismatch: {rel}")
        checked += 1
    actual_files = {
        path.relative_to(bundle_root).as_posix()
        for path in bundle_root.rglob("*")
        if path.is_file() and path != manifest
    }
    for rel in sorted(actual_files - listed):
        failures.append(f"file absent from manifest: {rel}")
    for rel in sorted(listed - actual_files):
        failures.append(f"manifest entry absent from bundle: {rel}")
    return checked


def verify_extracted(bundle_root: Path, verify_manifest: bool = True) -> dict[str, object]:
    bundle_root = bundle_root.resolve()
    failures: list[str] = []
    local_links_checked = 0
    html_pages = 0

    for required in (
        "OPEN_WEBSITE.html",
        "index.html",
        "README_FIRST.txt",
        "THIRD_PARTY_NOTICES.txt",
        "BUNDLE_INFO.json",
        "VERIFY_BUNDLE.py",
        "EDIT_WEBSITE.command",
        "EDIT_WEBSITE.md",
        "website/index.html",
        "website/editor.html",
        "website/content-overrides.json",
        "website/assets/editor.css",
        "website/assets/editor.js",
        "website/assets/apply-content.js",
        "website/assets/content-overrides.js",
        "website/deliverables.html",
        "website/research_brief.html",
        "website/detection-lab.html",
        "website/claim-survival.html",
        "website/graphical_abstract.svg",
        "website/assets/research-timeline.png",
        "website/assets/claim-survival.png",
        "website/assets/models-overview.png",
        "website/assets/model-community-multiverse.png",
        "website/assets/model-occupancy-atlas.png",
        "website/assets/model-pipeline-bootstrap.png",
        "website/assets/model-salmo-benchmark.png",
        "website/assets/model-trait-null.png",
        "_site/conwy-reanalysis.html",
        "output/pdf/RESEARCH_BRIEF.pdf",
        "output/pdf/POSTER.pdf",
        "scripts/edit_website.py",
    ):
        if not (bundle_root / required).is_file():
            failures.append(f"missing required bundle file: {required}")

    for required in PHASE5_REQUIRED_FILES:
        if not (bundle_root / required).is_file():
            failures.append(f"missing required Phase-5 bundle file: {required}")
    public_cache = bundle_root / PHASE5_TREES[0]
    if not public_cache.is_dir():
        failures.append(f"missing Phase-5 public-input cache: {PHASE5_TREES[0]}")
    else:
        genbank_records = list((public_cache / "genbank").glob("*.gb"))
        query_receipts = list((public_cache / "query_receipts").glob("*.json"))
        if len(genbank_records) != 29:
            failures.append(
                f"Phase-5 public-input cache expected 29 GenBank records, found {len(genbank_records)}"
            )
        if len(query_receipts) != 14:
            failures.append(
                f"Phase-5 public-input cache expected 14 NCBI query receipts, found {len(query_receipts)}"
            )

    for page in sorted(bundle_root.rglob("*.html")):
        html_pages += 1
        parser = parse_html(page)
        if parser.network_dependencies:
            failures.append(
                f"network-dependent assets: {page.relative_to(bundle_root)} -> "
                f"{sorted(set(parser.network_dependencies))}"
            )
        for _, raw in parser.links:
            try:
                target, fragment = local_target(page, raw, bundle_root)
            except ValueError as exc:
                failures.append(str(exc))
                continue
            if target is None:
                continue
            local_links_checked += 1
            if not target.exists():
                failures.append(f"broken local link: {page.relative_to(bundle_root)} -> {raw}")
                continue
            if fragment and target.suffix.lower() in {".html", ".htm"}:
                if fragment not in parse_html(target).ids:
                    failures.append(
                        f"missing HTML anchor: {page.relative_to(bundle_root)} -> {raw}"
                    )

    for path in bundle_root.rglob("*"):
        if path.is_symlink():
            failures.append(f"symbolic link is not allowed in bundle: {path.relative_to(bundle_root)}")
        rel = path.relative_to(bundle_root)
        if FORBIDDEN_PARTS & set(rel.parts):
            failures.append(f"forbidden path in bundle: {rel}")
        suffixes = {suffix.lower() for suffix in path.suffixes}
        if path.is_file() and FORBIDDEN_SUFFIXES & suffixes:
            failures.append(f"heavy checkpoint/raw file in bundle: {rel}")
        if path.is_file() and has_network_asset_reference(path):
            failures.append(f"network-dependent CSS/SVG asset in bundle: {rel}")

    checksums_checked = verify_manifest_file(bundle_root, failures) if verify_manifest else 0
    result = {
        "status": "PASS" if not failures else "FAIL",
        "path_contains_spaces": " " in str(bundle_root),
        "html_pages": html_pages,
        "local_links_checked": local_links_checked,
        "checksums_checked": checksums_checked,
        "failures": failures,
    }
    if failures:
        raise RuntimeError(json.dumps(result, indent=2))
    return result


def safe_extract(archive: Path, destination: Path) -> Path:
    with zipfile.ZipFile(archive) as handle:
        for member in handle.infolist():
            member_path = PurePosixPath(member.filename)
            if member_path.is_absolute() or ".." in member_path.parts:
                raise ValueError(f"unsafe ZIP member: {member.filename}")
            if member.is_dir():
                continue
            unix_mode = member.external_attr >> 16
            if unix_mode and (unix_mode & 0o170000) == 0o120000:
                raise ValueError(f"symbolic-link ZIP member: {member.filename}")
        handle.extractall(destination)
    return destination / BUNDLE_SLUG


def write_zip(bundle_root: Path, archive: Path) -> None:
    archive.parent.mkdir(parents=True, exist_ok=True)
    temporary = archive.with_suffix(archive.suffix + ".tmp")
    if temporary.exists():
        temporary.unlink()
    fixed_timestamp = (2026, 7, 19, 0, 0, 0)
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as handle:
        for path in sorted(item for item in bundle_root.rglob("*") if item.is_file()):
            rel = Path(BUNDLE_SLUG) / path.relative_to(bundle_root)
            info = zipfile.ZipInfo(rel.as_posix(), date_time=fixed_timestamp)
            info.compress_type = zipfile.ZIP_DEFLATED
            mode = 0o100755 if path.name == "EDIT_WEBSITE.command" else 0o100644
            info.external_attr = (mode & 0xFFFF) << 16
            handle.writestr(info, path.read_bytes())
    os.replace(temporary, archive)


def build(output: Path) -> dict[str, object]:
    sources, external_links = collect_files()
    with tempfile.TemporaryDirectory(prefix="conwy-edna-site-") as temporary:
        temporary_root = Path(temporary)
        bundle_root = temporary_root / BUNDLE_SLUG
        bundle_root.mkdir()
        copy_sources(sources, bundle_root)

        gateway = gateway_html()
        (bundle_root / "OPEN_WEBSITE.html").write_text(gateway, encoding="utf-8")
        (bundle_root / "index.html").write_text(gateway, encoding="utf-8")
        (bundle_root / "README_FIRST.txt").write_text(readme_text(), encoding="utf-8")
        (bundle_root / "THIRD_PARTY_NOTICES.txt").write_text(
            third_party_notices_text(), encoding="utf-8"
        )
        (bundle_root / "VERIFY_BUNDLE.py").write_text(verifier_wrapper(), encoding="utf-8")

        info = {
            "title": "River Conwy fish-eDNA reanalysis — portable project website",
            "bundle_version": RELEASE_DATE,
            "entrypoint": "OPEN_WEBSITE.html",
            "offline": True,
            "scope": (
                "Static project/project website and every local deliverable directly linked "
                "from its HTML catalogue, plus the complete Phase-5 machine audit, cached public "
                "sequence inputs, presentation inputs, and integrity records."
            ),
            "not_in_scope": [
                "raw FASTQs or a full raw-read pipeline rerun",
                "production MCMC or heavyweight model checkpoints",
                "hosted deployment, DOI, peer review, or publication",
            ],
            "external_links": sorted(external_links),
            "source_file_count": len(sources),
        }
        (bundle_root / "BUNDLE_INFO.json").write_text(
            json.dumps(info, indent=2) + "\n", encoding="utf-8"
        )
        manifested_files = write_manifest(bundle_root)
        staged_check = verify_extracted(bundle_root, verify_manifest=True)
        write_zip(bundle_root, output)

    with tempfile.TemporaryDirectory(prefix="conwy-edna-unzip-") as extracted:
        extraction_parent = Path(extracted) / "clean extraction with spaces"
        extraction_parent.mkdir()
        extracted_root = safe_extract(output, extraction_parent)
        if " " not in str(extracted_root):
            raise RuntimeError("clean-extraction path unexpectedly contains no spaces")
        extracted_check = verify_extracted(extracted_root, verify_manifest=True)
        bundled_verifier = subprocess.run(
            [sys.executable, str(extracted_root / "VERIFY_BUNDLE.py")],
            cwd=extracted_root,
            check=False,
            capture_output=True,
            text=True,
        )
        if bundled_verifier.returncode != 0:
            raise RuntimeError(
                "bundled VERIFY_BUNDLE.py failed after clean extraction:\n"
                + bundled_verifier.stdout
                + bundled_verifier.stderr
            )
        try:
            bundled_verifier_result = json.loads(bundled_verifier.stdout)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                "bundled VERIFY_BUNDLE.py did not emit a JSON result"
            ) from exc
        if bundled_verifier_result.get("status") != "PASS":
            raise RuntimeError("bundled VERIFY_BUNDLE.py did not report PASS")

    digest = sha256(output)
    checksum_path = output.with_suffix(output.suffix + ".sha256")
    checksum_path.write_text(f"{digest}  {output.name}\n", encoding="utf-8")
    try:
        archive_label = output.relative_to(ROOT).as_posix()
    except ValueError:
        archive_label = str(output)
    report = {
        "status": "PASS",
        "archive": archive_label,
        "archive_bytes": output.stat().st_size,
        "archive_sha256": digest,
        "source_files": len(sources),
        "manifested_files": manifested_files,
        "staged_verification": staged_check,
        "clean_extraction_verification": extracted_check,
        "bundled_verifier": {
            "status": bundled_verifier_result["status"],
            "returncode": bundled_verifier.returncode,
            "path_contains_spaces": bundled_verifier_result.get("path_contains_spaces"),
        },
    }
    report_path = output.with_suffix(".bundle.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a self-contained, offline ZIP of the project website and linked deliverables."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"archive path (default: {DEFAULT_OUTPUT.relative_to(ROOT)})",
    )
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if output.suffix.lower() != ".zip":
        parser.error("--output must end in .zip")
    report = build(output.resolve())
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
