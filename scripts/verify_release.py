#!/usr/bin/env python3
"""Verify the local website release without network access or raw-data work."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

try:
    from pypdf import PdfReader
except ImportError as exc:
    raise SystemExit(
        "Missing release-check dependency. Activate the START_HERE environment "
        "and install requirements-website.txt."
    ) from exc


ROOT = Path(__file__).resolve().parents[1]

MODEL_DIAGRAM_ASSETS = (
    "website/assets/model-community-multiverse.png",
    "website/assets/model-occupancy-atlas.png",
    "website/assets/model-pipeline-bootstrap.png",
    "website/assets/model-salmo-benchmark.png",
    "website/assets/model-trait-null.png",
)

WEBSITE_BUILD_OUTPUTS = (
    "website/assets/community-gradient.png",
    "website/assets/salmo-resolution.png",
    "website/assets/detection-simulation.png",
    "website/assets/claim-survival.png",
    "website/assets/research-timeline.png",
    "website/assets/models-overview.png",
    *MODEL_DIAGRAM_ASSETS,
    "output/pdf/RESEARCH_BRIEF.pdf",
    "output/pdf/POSTER.pdf",
)

WEBSITE_PDF_PAGES = {
    "output/pdf/RESEARCH_BRIEF.pdf": 4,
    "output/pdf/POSTER.pdf": 1,
}

EDITOR_REQUIRED_ARTIFACTS = (
    "EDIT_WEBSITE.command",
    "EDIT_WEBSITE.md",
    "website/editor.html",
    "website/content-overrides.json",
    "website/assets/editor.css",
    "website/assets/editor.js",
    "website/assets/apply-content.js",
    "website/assets/content-overrides.js",
    "scripts/edit_website.py",
)

# Narrative fields are deliberately enumerated. Expanding this contract requires
# a conscious release change so evidence, results, links, and verdict labels do
# not accidentally become editable presentation copy.
ALLOWED_EDITABLE_TEXT_KEYS = frozenset({
    "claims.summary",
    "claims.title",
    "contribution.body",
    "contribution.title",
    "hero.eyebrow",
    "hero.summary",
    "hero.title",
    "methodology.diagram_caption",
    "methodology.step1.body",
    "methodology.step1.title",
    "methodology.step2.body",
    "methodology.step2.title",
    "methodology.step3.body",
    "methodology.step3.title",
    "methodology.step4.body",
    "methodology.step4.title",
    "methodology.step5.body",
    "methodology.step5.title",
    "methodology.summary",
    "methodology.title",
    "models.diagram_caption",
    "models.summary",
    "models.title",
    "next_evidence.title",
    "question.biological.body",
    "question.biological.title",
    "question.identity.body",
    "question.identity.title",
    "question.observation.body",
    "question.observation.title",
    "question.summary",
    "question.title",
    "routes.brief.body",
    "routes.brief.title",
    "routes.catalogue.body",
    "routes.catalogue.title",
    "routes.release.body",
    "routes.release.title",
    "routes.report.body",
    "routes.report.title",
    "stress_tests.diagram_caption",
    "stress_tests.summary",
    "stress_tests.title",
    "verdict.title",
})

REQUIRED_EDITABLE_TEXT_KEYS = frozenset({
    "hero.eyebrow",
    "hero.summary",
    "question.title",
    "question.summary",
    "methodology.title",
    "methodology.summary",
    "contribution.title",
    "contribution.body",
})

FROZEN_IMPACT_HASHES = {
    "docs/impact-sprint-preregistration.md":
        "5d4ab9ce7d1f0dc695646f6e3c781b4931fad5fecc1513ddd79046c3b316e567",
    "data/metadata/impact_sprint_manifest.json":
        "0f6e6464ec54492cb23b8adbdc81b3b05a9ec98e73e1b04b413c764c03194843",
}

PHASE5_REQUIRED = (
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
    "data/metadata/claim_survival.tsv",
    "website/claim-survival.html",
    "website/assets/claim-survival.png",
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


class PageParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids: list[str] = []
        self.links: list[tuple[str, str]] = []
        self.inputs: dict[str, dict[str, str | None]] = {}
        self.labels_for: list[str] = []
        self.tag_counts: dict[str, int] = {}
        self.images: list[dict[str, str | None]] = []
        self.buttons: list[dict[str, str | None]] = []
        self.html_attrs: dict[str, str | None] = {}
        self.text_chunks: list[str] = []
        self.edit_keys: list[str] = []

    def handle_starttag(self, tag, attrs):
        attr = dict(attrs)
        self.tag_counts[tag] = self.tag_counts.get(tag, 0) + 1
        if attr.get("id"):
            self.ids.append(attr["id"])
        for key in ("href", "src"):
            if attr.get(key):
                self.links.append((key, attr[key]))
        if tag == "input" and attr.get("id"):
            self.inputs[attr["id"]] = attr
        if tag == "label" and attr.get("for"):
            self.labels_for.append(attr["for"])
        if tag == "img":
            self.images.append(attr)
        if tag == "button":
            self.buttons.append(attr)
        if tag == "html":
            self.html_attrs = attr
        if attr.get("data-edit-key"):
            self.edit_keys.append(attr["data-edit-key"])

    def handle_data(self, data):
        self.text_chunks.append(data)


def local_target(source: Path, raw: str) -> tuple[Path | None, str]:
    parsed = urlsplit(raw)
    if parsed.scheme or parsed.netloc or raw.startswith(("mailto:", "javascript:")):
        return None, ""
    path = unquote(parsed.path)
    target = source if not path else (source.parent / path).resolve()
    return target, parsed.fragment


def check_html(path: Path, failures: list[str]):
    parser = PageParser()
    parser.feed(path.read_text(encoding="utf-8"))
    duplicates = sorted({item for item in parser.ids if parser.ids.count(item) > 1})
    if duplicates:
        failures.append(f"{path.relative_to(ROOT)} duplicate ids: {duplicates}")
    for key, raw in parser.links:
        target, fragment = local_target(path, raw)
        if target is None:
            continue
        if not target.exists():
            failures.append(f"{path.relative_to(ROOT)} broken {key}: {raw}")
            continue
        if fragment and target.suffix.lower() in {".html", ".htm"}:
            other = PageParser()
            other.feed(target.read_text(encoding="utf-8"))
            if fragment not in other.ids:
                failures.append(f"{path.relative_to(ROOT)} missing anchor: {raw}")
    for input_id in parser.inputs:
        if input_id not in parser.labels_for:
            failures.append(f"{path.relative_to(ROOT)} input #{input_id} has no label")
    return parser


def check_markdown(path: Path, failures: list[str]):
    text = path.read_text(encoding="utf-8")
    for raw in re.findall(r"(?<!!)\[[^\]]+\]\(([^)]+)\)", text):
        raw = raw.strip().split(maxsplit=1)[0].strip("<>")
        target, _ = local_target(path, raw)
        if target is not None and not target.exists():
            failures.append(f"{path.relative_to(ROOT)} broken link: {raw}")


def tsv_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_editor(
    failures: list[str], html_parsers: dict[str, PageParser]
) -> dict[str, int]:
    """Verify the editor's narrow content contract and local-write boundary."""

    landing_path = ROOT / "website/index.html"
    landing = html_parsers.get("website/index.html")
    actual_keys: set[str] = set()
    if landing and landing_path.is_file():
        duplicate_keys = sorted(
            key for key, count in Counter(landing.edit_keys).items() if count > 1
        )
        if duplicate_keys:
            failures.append(f"website landing duplicate editable keys: {duplicate_keys}")
        actual_keys = set(landing.edit_keys)
        unsupported = sorted(actual_keys - ALLOWED_EDITABLE_TEXT_KEYS)
        if unsupported:
            failures.append(
                "website landing exposes unapproved editable fields: "
                + ", ".join(unsupported)
            )
        missing_core = sorted(REQUIRED_EDITABLE_TEXT_KEYS - actual_keys)
        if missing_core:
            failures.append(
                "website landing lacks core personalization fields: "
                + ", ".join(missing_core)
            )

        landing_text = landing_path.read_text(encoding="utf-8")
        editable_tags = re.findall(
            r"<[^>]*\bdata-edit-key=\"[^\"]+\"[^>]*>", landing_text
        )
        protected = [
            tag for tag in editable_tags
            if re.search(r"\bdata-(?:claim-id|status|decision)=", tag)
        ]
        if protected:
            failures.append(
                "website landing makes a protected claim/status element editable"
            )

        script_sources = [raw for key, raw in landing.links if key == "src"]
        required_scripts = (
            "assets/content-overrides.js",
            "assets/apply-content.js",
        )
        for source in required_scripts:
            if source not in script_sources:
                failures.append(f"website landing lacks editor content script: {source}")
        if all(source in script_sources for source in required_scripts):
            positions = [script_sources.index(source) for source in required_scripts]
            if positions != sorted(positions):
                failures.append(
                    "website landing loads apply-content.js before content-overrides.js"
                )

    override_path = ROOT / "website/content-overrides.json"
    override_values: dict[str, str] | None = None
    if override_path.is_file():
        try:
            override_state = json.loads(override_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError, UnicodeError) as exc:
            failures.append(f"website/content-overrides.json is unreadable: {exc}")
            override_state = None
        if not isinstance(override_state, dict):
            failures.append("website/content-overrides.json must contain an object")
        else:
            expected_fields = {"version", "revision", "updated_at", "values"}
            if set(override_state) != expected_fields:
                failures.append(
                    "website/content-overrides.json fields differ: expected "
                    f"{sorted(expected_fields)}, found {sorted(override_state)}"
                )
            if type(override_state.get("version")) is not int or override_state.get("version") != 1:
                failures.append("website/content-overrides.json version must be integer 1")
            updated_at = override_state.get("updated_at")
            if updated_at is not None and not isinstance(updated_at, str):
                failures.append("website/content-overrides.json updated_at must be null or text")
            values = override_state.get("values")
            if not isinstance(values, dict):
                failures.append("website/content-overrides.json values must be an object")
            else:
                override_values = values
                unknown = sorted(set(values) - actual_keys)
                if unknown:
                    failures.append(
                        "website/content-overrides.json contains unknown/non-editable keys: "
                        + ", ".join(unknown)
                    )
                non_text = sorted(key for key, value in values.items() if not isinstance(value, str))
                if non_text:
                    failures.append(
                        "website/content-overrides.json has non-text values: "
                        + ", ".join(non_text)
                    )
                canonical = json.dumps(
                    values,
                    ensure_ascii=True,
                    sort_keys=True,
                    separators=(",", ":"),
                ).encode("utf-8")
                expected_revision = hashlib.sha256(canonical).hexdigest()
                if override_state.get("revision") != expected_revision:
                    failures.append(
                        "website/content-overrides.json revision does not match its values"
                    )

    content_js_path = ROOT / "website/assets/content-overrides.js"
    if content_js_path.is_file():
        content_js = content_js_path.read_text(encoding="utf-8")
        match = re.fullmatch(
            r"\s*window\.WEBSITE_CONTENT\s*=\s*Object\.freeze\((\{.*\})\);\s*",
            content_js,
            flags=re.DOTALL,
        )
        if not match:
            failures.append(
                "website/assets/content-overrides.js lacks the generated "
                "window.WEBSITE_CONTENT marker"
            )
        else:
            try:
                js_values = json.loads(match.group(1))
            except json.JSONDecodeError as exc:
                failures.append(f"website/assets/content-overrides.js has invalid JSON: {exc}")
            else:
                if override_values is not None and js_values != override_values:
                    failures.append(
                        "content-overrides.js does not match content-overrides.json values"
                    )

    apply_js_path = ROOT / "website/assets/apply-content.js"
    if apply_js_path.is_file():
        apply_js = apply_js_path.read_text(encoding="utf-8")
        for marker in (
            "window.WEBSITE_CONTENT",
            'querySelectorAll("[data-edit-key]")',
            "textContent",
        ):
            if marker not in apply_js:
                failures.append(f"apply-content.js lacks required marker: {marker}")

    editor_path = ROOT / "website/editor.html"
    editor = html_parsers.get("website/editor.html")
    if editor and editor_path.is_file():
        required_ids = {
            "discardButton",
            "undoButton",
            "exportButton",
            "importButton",
            "resetButton",
            "saveButton",
            "openSiteLink",
            "selectionLabel",
            "saveStatus",
            "sitePreview",
            "importFile",
        }
        missing_ids = sorted(required_ids - set(editor.ids))
        if missing_ids:
            failures.append(f"website editor lacks controls: {missing_ids}")
        editor_assets = {raw for _, raw in editor.links}
        for asset in ("assets/editor.css", "assets/editor.js"):
            if asset not in editor_assets:
                failures.append(f"website editor lacks asset: {asset}")
        editor_text = editor_path.read_text(encoding="utf-8")
        if not re.search(
            r"<iframe\b(?=[^>]*\bid=\"sitePreview\")"
            r"(?=[^>]*\bsrc=\"index\.html\?editor-preview=1\")[^>]*>",
            editor_text,
        ):
            failures.append("website editor lacks its local editable-preview iframe")

    editor_js_path = ROOT / "website/assets/editor.js"
    if editor_js_path.is_file():
        editor_js = editor_js_path.read_text(encoding="utf-8")
        for marker in (
            "/api/content",
            "X-Editor-Token",
            "base_revision",
            "contentEditable",
            "textContent",
            "beforeunload",
        ):
            if marker not in editor_js:
                failures.append(f"editor.js lacks required behavior marker: {marker}")

    server_path = ROOT / "scripts/edit_website.py"
    if server_path.is_file():
        server_text = server_path.read_text(encoding="utf-8")
        for marker in (
            "127.0.0.1",
            "ThreadingHTTPServer",
            "secrets.token_urlsafe",
            "X-Editor-Token",
            "base_revision",
            "409",
            "Cache-Control",
            "no-store",
            "Origin",
            "Content-Length",
            "application/json",
            "os.replace",
        ):
            if marker not in server_text:
                failures.append(f"edit_website.py lacks security marker: {marker}")
        if "0.0.0.0" in server_text:
            failures.append("edit_website.py must not bind to all network interfaces")

    launcher_path = ROOT / "EDIT_WEBSITE.command"
    if launcher_path.is_file():
        launcher_text = launcher_path.read_text(encoding="utf-8")
        if "scripts/edit_website.py" not in launcher_text:
            failures.append("EDIT_WEBSITE.command does not launch the local editor server")

    return {
        "editable_text_fields": len(actual_keys),
        "saved_text_overrides": len(override_values or {}),
    }


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "pass"}


def check_phase5(failures: list[str], html_parsers: dict[str, PageParser]) -> dict[str, int]:
    """Check the frozen Phase-5 contract, exact denominators, and claim boundaries."""

    def rows(rel: str) -> list[dict[str, str]]:
        path = ROOT / rel
        if not path.is_file():
            return []
        try:
            return tsv_rows(path)
        except (csv.Error, OSError, UnicodeError) as exc:
            failures.append(f"{rel}: cannot read TSV: {exc}")
            return []

    def json_file(rel: str) -> dict[str, object]:
        path = ROOT / rel
        if not path.is_file():
            return {}
        try:
            value = json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError, UnicodeError) as exc:
            failures.append(f"{rel}: cannot read JSON: {exc}")
            return {}
        if not isinstance(value, dict):
            failures.append(f"{rel}: expected a JSON object")
            return {}
        return value

    def exact(actual: object, expected: object, label: str):
        if actual != expected:
            failures.append(f"{label}: expected {expected!r}, found {actual!r}")

    def number(value: object, label: str) -> float:
        try:
            return float(str(value))
        except (TypeError, ValueError):
            failures.append(f"{label}: expected a number, found {value!r}")
            return float("nan")

    # The receipt is not merely self-consistent: its two source digests are the
    # hashes frozen before any new Phase-5 analysis ran.
    freeze_rows = rows("data/metadata/impact_sprint_preregistration_receipt.tsv")
    freeze_by_path = {row.get("artifact", ""): row for row in freeze_rows}
    exact(len(freeze_rows), 2, "impact preregistration receipt row count")
    for rel, expected_digest in FROZEN_IMPACT_HASHES.items():
        path = ROOT / rel
        if path.is_file():
            exact(sha256_file(path), expected_digest, f"frozen impact digest for {rel}")
        receipt = freeze_by_path.get(rel)
        if not receipt:
            failures.append(f"impact preregistration receipt lacks {rel}")
            continue
        exact(receipt.get("sha256"), expected_digest, f"recorded impact digest for {rel}")
        exact(receipt.get("frozen_date"), "2026-07-19", f"freeze date for {rel}")
        if "not externally registered" not in receipt.get("timing_scope", ""):
            failures.append(f"impact receipt overstates registration timing for {rel}")

    manifest = json_file("data/metadata/impact_sprint_manifest.json")
    exact(manifest.get("schema_version"), 1, "impact manifest schema")
    exact(manifest.get("specification_date"), "2026-07-19", "impact specification date")
    if "not externally registered" not in str(manifest.get("timing_scope", "")):
        failures.append("impact manifest does not disclose its repository-local timing scope")
    exact(
        manifest.get("seeds"),
        {
            "occupancy_atlas": 2026071901,
            "layerE_multiverse": 2026071902,
            "trait_null": 2026071903,
            "pipeline_bootstrap": 2026071904,
        },
        "impact fixed seeds",
    )
    atlas_spec = manifest.get("occupancy_atlas", {})
    if isinstance(atlas_spec, dict):
        exact(atlas_spec.get("expected_scenarios"), 648, "atlas scenario specification")
        exact(atlas_spec.get("simulations_per_scenario"), 100, "atlas simulation specification")
        exact(atlas_spec.get("n_units"), 120, "atlas unit specification")
        exact(
            atlas_spec.get("estimators"),
            ["naive_occurrence", "replication_aware", "depth_aware"],
            "atlas estimator specification",
        )
    layer_spec = manifest.get("layerE_multiverse", {})
    if isinstance(layer_spec, dict):
        exact(layer_spec.get("expected_specification_families"), 72, "Layer-E family specification")
        exact(layer_spec.get("one_library_draws"), 100, "Layer-E draw specification")
        exact(layer_spec.get("expected_reference_influence_rows"), 75, "Layer-E influence specification")
    trait_spec = manifest.get("trait_null", {})
    if isinstance(trait_spec, dict):
        exact(trait_spec.get("balanced_units"), 154, "trait-null balanced units")
        exact(trait_spec.get("null_draws"), 9999, "trait-null draw specification")
    pipeline_spec = manifest.get("pipeline_disagreement", {})
    if isinstance(pipeline_spec, dict):
        for field, expected in (
            ("shared_libraries", 642),
            ("source_positive_cells", 2344),
            ("recovered_cells", 1962),
            ("nonrecovered_cells", 382),
            ("bootstrap_draws", 999),
        ):
            exact(pipeline_spec.get(field), expected, f"pipeline specification {field}")

    canonical_rows = rows("data/metadata/impact_sprint_canonical_hashes.tsv")
    exact(len(canonical_rows), 47, "impact canonical-hash rows")
    for row in canonical_rows:
        rel = row.get("artifact", "")
        path = (ROOT / rel).resolve()
        try:
            path.relative_to(ROOT.resolve())
        except ValueError:
            failures.append(f"impact canonical-hash path escapes the project: {rel}")
            continue
        if not path.is_file():
            failures.append(f"impact canonical-hash artifact is missing: {rel}")
            continue
        exact(sha256_file(path), row.get("sha256"), f"impact canonical digest {rel}")
        exact(str(path.stat().st_size), row.get("bytes"), f"impact canonical bytes {rel}")
    sprint_receipt = json_file("data/metadata/impact_sprint_run_receipt.json")
    exact(sprint_receipt.get("status"), "PASS", "impact sprint run status")
    exact(sprint_receipt.get("frozen_specification"), "PASS", "impact sprint frozen-specification gate")
    exact(sprint_receipt.get("analysis_mode"), "full-rebuild", "impact sprint full-rebuild mode")
    exact(sprint_receipt.get("release_mode"), True, "impact sprint release mode")
    exact(
        sprint_receipt.get("receipt_finalized_before_packaging"),
        True,
        "impact sprint receipt finalized before packaging",
    )
    exact(sprint_receipt.get("canonical_artifacts"), 47, "impact sprint canonical artifact count")
    exact(sprint_receipt.get("mismatches"), [], "impact sprint deterministic mismatches")
    exact(
        sprint_receipt.get("deterministic_against_previous_run"),
        "PASS",
        "impact sprint second-run determinism gate",
    )

    # Known-truth occupancy recoverability atlas.
    scenarios = rows("data/derived/impact_occupancy_scenarios.tsv")
    atlas = rows("data/derived/impact_occupancy_atlas.tsv")
    fit_rows = rows("data/derived/impact_occupancy_fit_results.tsv")
    audit_rows = rows("data/derived/impact_occupancy_simulation_audit.tsv")
    fixture_rows = rows("data/derived/impact_occupancy_generator_fixture.tsv")
    exact(len(scenarios), 648, "occupancy atlas scenario rows")
    exact(len(atlas), 1944, "occupancy atlas estimator-scenario rows")
    exact(len(fit_rows), 194400, "occupancy atlas fit rows")
    exact(len(audit_rows), 64800, "occupancy atlas simulation-audit rows")
    exact(len(fixture_rows), 1080, "occupancy atlas generator-fixture rows")
    if atlas:
        exact(
            Counter(row.get("gate_class") for row in atlas),
            Counter({"recoverable": 613, "unreliable": 732, "non-identifiable": 599}),
            "occupancy atlas gate-class totals",
        )
        exact({row.get("data_scope") for row in atlas}, {"synthetic_known_truth_only"}, "occupancy atlas scope")
        exact({row.get("simulations") for row in atlas}, {"100"}, "occupancy atlas simulations per row")
        exact(
            {row.get("estimator") for row in atlas},
            {"naive_occurrence", "replication_aware", "depth_aware"},
            "occupancy atlas estimators",
        )
    if fit_rows:
        counts = Counter((row.get("scenario_id"), row.get("estimator")) for row in fit_rows)
        if len(counts) != 1944 or set(counts.values()) != {100}:
            failures.append("occupancy fit table is not exactly 100 simulations for each of 1,944 estimator-scenarios")
        exact({row.get("data_scope") for row in fit_rows}, {"synthetic_known_truth_only"}, "occupancy fit scope")
    if audit_rows:
        audit_counts = Counter(row.get("scenario_id") for row in audit_rows)
        if len(audit_counts) != 648 or set(audit_counts.values()) != {100}:
            failures.append("occupancy audit is not exactly 100 simulations for each of 648 scenarios")
    occupancy_manifest = json_file("data/derived/impact_occupancy_manifest.json")
    exact(occupancy_manifest.get("layerA_real_inputs_read"), [], "occupancy atlas real-input exclusion")
    if "Cannot estimate Conwy occupancy" not in str(occupancy_manifest.get("claim_boundary", "")):
        failures.append("occupancy atlas manifest lacks the real-Conwy claim boundary")
    canonical_paths = {
        "scenarios": "data/derived/impact_occupancy_scenarios.tsv",
        "simulation_audit": "data/derived/impact_occupancy_simulation_audit.tsv",
        "fit_results": "data/derived/impact_occupancy_fit_results.tsv",
        "atlas": "data/derived/impact_occupancy_atlas.tsv",
        "generator_fixture": "data/derived/impact_occupancy_generator_fixture.tsv",
        "figure_png": "figures/impact_occupancy_atlas.png",
        "figure_svg": "figures/impact_occupancy_atlas.svg",
    }
    canonical = occupancy_manifest.get("canonical_output_sha256", {})
    if isinstance(canonical, dict):
        for key, rel in canonical_paths.items():
            if (ROOT / rel).is_file():
                exact(canonical.get(key), sha256_file(ROOT / rel), f"occupancy canonical digest {key}")
    occupancy_receipt = json_file("data/derived/impact_occupancy_receipt.json")
    exact(occupancy_receipt.get("deterministic_rebuild_pass"), True, "occupancy deterministic rebuild")
    exact(occupancy_receipt.get("layerA_real_inputs_read"), [], "occupancy receipt real-input exclusion")
    test_receipt = rows("data/derived/impact_occupancy_test_receipt.tsv")
    if test_receipt:
        exact(test_receipt[0].get("tests_passed"), "59", "occupancy atlas tests passed")
        for field in ("unit_invariant_regression_tests", "full_dimension_denominator_provenance_tests", "deterministic_rebuild"):
            exact(test_receipt[0].get(field), "PASS", f"occupancy atlas test receipt {field}")

    # Public-record Salmo multi-marker benchmark and its complete offline cache.
    salmo_queries = rows("data/metadata/impact_salmo_queries.tsv")
    salmo_accessions = rows("data/metadata/impact_salmo_accessions.tsv")
    marker_manifest = rows("data/metadata/impact_salmo_marker_manifest.tsv")
    panel_coverage = rows("data/metadata/impact_salmo_panel_coverage.tsv")
    salmo_gates = rows("data/derived/impact_salmo_marker_gates.tsv")
    salmo_pairs = rows("data/derived/impact_salmo_pairwise_distances.tsv")
    salmo_loo = rows("data/derived/impact_salmo_loo_classification.tsv")
    salmo_gaps = rows("data/derived/impact_salmo_barcode_gaps.tsv")
    exact(len(salmo_queries), 9, "Salmo query rows")
    exact(len(salmo_accessions), 143, "Salmo accession audit rows")
    exact(sum(row.get("selected") == "yes" for row in salmo_accessions), 21, "selected Salmo/outgroup accessions")
    exact(len(marker_manifest), 63, "Salmo marker-sequence rows")
    exact(len(panel_coverage), 27, "Salmo marker-panel coverage rows")
    exact(len(salmo_gates), 3, "Salmo marker gate rows")
    exact(len(salmo_pairs), 630, "Salmo pairwise-distance rows")
    exact(len(salmo_loo), 60, "Salmo leave-one-accession-out rows")
    exact(len(salmo_gaps), 6, "Salmo barcode-gap rows")
    expected_lengths = {"mifish_12S": "170", "coi_5prime": "648", "cytb": "1140"}
    if marker_manifest:
        exact({row.get("marker") for row in marker_manifest}, set(expected_lengths), "Salmo marker IDs")
        for marker, expected_length in expected_lengths.items():
            lengths = {row.get("sequence_length") for row in marker_manifest if row.get("marker") == marker}
            exact(lengths, {expected_length}, f"{marker} sequence length")
    for row in salmo_queries:
        for receipt_field, digest_field in (
            ("esearch_receipt", "esearch_sha256"),
            ("esummary_receipt", "esummary_sha256"),
        ):
            rel = row.get(receipt_field, "")
            if not rel:
                continue
            path = ROOT / rel
            if not path.is_file():
                failures.append(f"Salmo query receipt missing: {rel}")
            else:
                exact(sha256_file(path), row.get(digest_field), f"Salmo query-receipt digest {rel}")
    public_root = (ROOT / "data/public/salmo_multimarker").resolve()
    for row in salmo_accessions:
        rel = row.get("raw_genbank_path", "")
        if not rel:
            continue
        path = (ROOT / rel).resolve()
        try:
            path.relative_to(public_root)
        except ValueError:
            failures.append(f"Salmo raw-record path escapes the public cache: {rel}")
            continue
        if not path.is_file():
            failures.append(f"Salmo cached GenBank record missing: {rel}")
        else:
            exact(sha256_file(path), row.get("raw_genbank_sha256"), f"Salmo GenBank digest {rel}")
    exact(len(list((public_root / "genbank").glob("*.gb"))), 29, "cached Salmo GenBank files")
    exact(len(list((public_root / "query_receipts").glob("*.json"))), 14, "cached Salmo query-receipt files")
    fasta_path = public_root / "marker_sequences.fasta"
    if fasta_path.is_file() and marker_manifest:
        fasta: dict[str, str] = {}
        current = ""
        for line in fasta_path.read_text(encoding="ascii").splitlines():
            if line.startswith(">"):
                current = line[1:].split()[0]
                fasta[current] = ""
            elif current:
                fasta[current] += line.strip()
        exact(len(fasta), 63, "Salmo cached FASTA records")
        for row in marker_manifest:
            sequence = fasta.get(row.get("sequence_id", ""), "")
            if not sequence:
                failures.append(f"Salmo FASTA lacks {row.get('sequence_id', '')}")
            else:
                exact(
                    hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                    row.get("sequence_sha256"),
                    f"Salmo marker digest {row.get('sequence_id', '')}",
                )
    for row in salmo_gates:
        marker = row.get("marker", "unknown marker")
        expected = {
            "passing_salmo_accessions": "11",
            "candidate_salmo_labels": "7",
            "eligible_salmo_labels": "2",
            "required_salmo_labels_for_species": "4",
            "species_panel_coverage_adequate": "no",
            "species_loo_queries": "10",
            "species_loo_all_ties_accuracy": "1",
            "cross_label_exact_haplotypes": "0",
            "species_gate_status": "coverage-inadequate",
            "eligible_outgroup_labels": "2",
            "genus_panel_coverage_adequate": "yes",
            "genus_loo_queries": "10",
            "salmo_queries_nearer_than_each_outgroup": "10",
            "genus_gate_status": "pass",
            "finest_consistently_supported_rank": "genus",
        }
        for field, wanted in expected.items():
            exact(row.get(field), wanted, f"{marker} gate {field}")
        if "2 of 7 candidate Salmo labels eligible; minimum 4 required" not in row.get("species_gate_reason", ""):
            failures.append(f"{marker} species gate obscures its candidate denominator and threshold")
        boundary = row.get("claim_boundary", "")
        for phrase in ("annotations, not truth", "No local occurrence", "ASV0001 species identity"):
            if phrase not in boundary:
                failures.append(f"{marker} gate lacks claim boundary: {phrase}")

    # Layer-E specification multiverse.
    library_counts = rows("data/derived/impact_library_counts.tsv")
    export_receipt = rows("data/metadata/impact_library_export_receipt.tsv")
    layer_specs = rows("data/derived/impact_layerE_specifications.tsv")
    layer_library_draws = rows("data/derived/impact_layerE_library_draws.tsv")
    layer_draw_results = rows("data/derived/impact_layerE_draw_results.tsv")
    layer_results = rows("data/derived/impact_layerE_results.tsv")
    layer_influence = rows("data/derived/impact_layerE_influence.tsv")
    layer_gate = rows("data/derived/impact_layerE_gate.tsv")
    exact(len(library_counts), 647, "portable library-count rows")
    if export_receipt and any(row.get("status") != "PASS" for row in export_receipt):
        failures.append("impact library export receipt contains a failed check")
    exact(len(layer_specs), 72, "Layer-E specification families")
    exact(len(layer_library_draws), 23400, "Layer-E one-library selections")
    exact(len(layer_draw_results), 3600, "Layer-E one-library result rows")
    exact(len(layer_results), 72, "Layer-E family result rows")
    exact(len(layer_influence), 75, "Layer-E reference influence rows")
    exact(len(layer_gate), 5, "Layer-E gate rows")
    layer_by_criterion = {row.get("criterion", ""): row for row in layer_gate}
    layer_expectations = {
        "valid_specification_families": ("64", "72", "FALSE"),
        "positive_valid_family_contrasts": ("64", "64", "TRUE"),
        "positive_one_library_draw_contrasts": ("3600", "3600", "TRUE"),
        "reference_leaveout_sign_reversals": ("0", "75", "TRUE"),
        "overall_gate": ("NA", "NA", "FALSE"),
    }
    for criterion, (numerator, denominator, passed) in layer_expectations.items():
        row = layer_by_criterion.get(criterion)
        if not row:
            failures.append(f"Layer-E gate lacks {criterion}")
            continue
        exact(row.get("numerator"), numerator, f"Layer-E {criterion} numerator")
        exact(row.get("denominator"), denominator, f"Layer-E {criterion} denominator")
        exact(row.get("pass"), passed, f"Layer-E {criterion} pass")
    if layer_results:
        exact(sum(truthy(row.get("family_valid")) for row in layer_results), 64, "valid Layer-E families")
        exact(sum(not truthy(row.get("family_valid")) for row in layer_results), 8, "invalid Layer-E families")
    if layer_draw_results:
        exact(sum(truthy(row.get("valid")) for row in layer_draw_results), 3200, "valid Layer-E draw results")
        exact(
            sum(
                number(
                    row.get("richness_contrast_lower_minus_upper"),
                    "Layer-E draw richness contrast",
                ) > 0
                for row in layer_draw_results
            ),
            3600,
            "positive Layer-E draw contrasts",
        )
    if layer_influence:
        exact(sum(truthy(row.get("sign_reversal")) for row in layer_influence), 0, "Layer-E sign reversals")
    layer_determinism = rows("data/metadata/impact_layerE_determinism_receipt.tsv")
    exact(len(layer_determinism), 8, "Layer-E determinism receipt rows")
    for row in layer_determinism:
        rel = row.get("artifact", "")
        if not truthy(row.get("byte_identical_to_preexisting")) or row.get("status") != "PASS":
            failures.append(f"Layer-E deterministic rebuild failed for {rel}")
        path = ROOT / rel
        if path.is_file():
            exact(sha256_file(path), row.get("rebuilt_sha256"), f"Layer-E deterministic digest {rel}")

    # Fixed-margin trait null model.
    trait_gower = rows("data/derived/impact_trait_gower.tsv")
    trait_units = rows("data/derived/impact_trait_unit_rao.tsv")
    trait_draws = rows("data/derived/impact_trait_null_draws.tsv")
    trait_summary = rows("data/derived/impact_trait_null_summary.tsv")
    exact(len(trait_gower), 169, "trait Gower rows")
    exact(len(trait_units), 308, "trait Rao-Q unit rows")
    exact(len(trait_draws), 19998, "trait-null draw rows")
    exact(len(trait_summary), 2, "trait-null summary rows")
    if trait_draws:
        exact(Counter(row.get("analysis") for row in trait_draws), Counter({"primary_12_labels": 9999, "strict_5_labels": 9999}), "trait-null draws per panel")
        if any(not truthy(row.get("margin_pass")) or not truthy(row.get("swap_target_pass")) for row in trait_draws):
            failures.append("trait-null draw table contains a margin or movable-block swap failure")
    trait_by_analysis = {row.get("analysis", ""): row for row in trait_summary}
    trait_expected = {
        "primary_12_labels": ("12", "0.3159"),
        "strict_5_labels": ("5", "0.1142"),
    }
    for analysis, (labels, p_value) in trait_expected.items():
        row = trait_by_analysis.get(analysis)
        if not row:
            failures.append(f"trait-null summary lacks {analysis}")
            continue
        exact(row.get("n_labels"), labels, f"{analysis} label count")
        exact(row.get("n_units"), "154", f"{analysis} unit count")
        exact(row.get("null_draws"), "9999", f"{analysis} null draws")
        exact(row.get("monte_carlo_p"), p_value, f"{analysis} Monte Carlo p")
        exact(row.get("gate_pass"), "FALSE", f"{analysis} gate pass")
        exact(row.get("claim_gate"), "FAILED_WITHHOLD_BEYOND_RICHNESS", f"{analysis} claim gate")

    # Cluster-aware pipeline-disagreement analysis.
    pipeline_cells = rows("data/metadata/impact_pipeline_disagreement_cells.tsv")
    pipeline_bootstrap = rows("data/derived/impact_pipeline_bootstrap.tsv")
    pipeline_leaveout = rows("data/derived/impact_pipeline_leaveout.tsv")
    pipeline_model = rows("data/derived/impact_pipeline_model.tsv")
    pipeline_rates = rows("data/derived/impact_pipeline_rates.tsv")
    pipeline_summary = rows("data/derived/impact_pipeline_summary.tsv")
    exact(len(pipeline_cells), 2344, "pipeline source-positive cells")
    if pipeline_cells:
        exact(sum(int(row.get("miss", "0")) for row in pipeline_cells), 382, "pipeline nonrecovered cells")
        if any(
            int(float(row.get("reprocessed_assigned_fish_read_total", "-1")))
            - int(float(row.get("reprocessed_count", "-1")))
            != int(float(row.get("reprocessed_assigned_fish_reads_excluding_focal_label", "-2")))
            for row in pipeline_cells
        ):
            failures.append("pipeline focal-excluded total does not exactly remove the focal count")
        exact(
            {row.get("evidence_scope") for row in pipeline_cells},
            {"pipeline disagreement relative to source table; source is not truth"},
            "pipeline evidence scope",
        )
    exact(len(pipeline_bootstrap), 7992, "pipeline bootstrap term rows")
    exact(len(pipeline_leaveout), 344, "pipeline leave-out term rows")
    exact(len(pipeline_rates), 18, "pipeline grouped-rate rows")
    exact(len(pipeline_summary), 1, "pipeline summary rows")
    if pipeline_summary:
        summary = pipeline_summary[0]
        expected = {
            "shared_libraries": "642",
            "harmonised_labels": "14",
            "source_positive_labels": "11",
            "source_positive_cells": "2344",
            "recovered_cells": "1962",
            "nonrecovered_cells": "382",
            "collinearity_gate_pass": "FALSE",
            "bootstrap_draws": "999",
            "finite_bootstrap_fits": "999",
            "leaveout_groups": "43",
            "promoted_associations": "2",
            "sensitivity_bootstrap_draws": "999",
            "sensitivity_finite_bootstrap_fits": "999",
            "sensitivity_leaveout_groups": "43",
            "depth_association_survives_focal_exclusion": "TRUE",
            "depth_association_sensitivity_status": "SURVIVES_FOCAL_EXCLUSION",
            "overall_gate": "ASSOCIATION_EVIDENCE",
        }
        for field, wanted in expected.items():
            exact(summary.get(field), wanted, f"pipeline summary {field}")
        if "neither pipeline is biological truth" not in summary.get("evidence_scope", ""):
            failures.append("pipeline summary lacks the comparator-relative claim boundary")
        for definition_marker in (
            "all 35 raw-name columns",
            "including the focal harmonised-label count and Unassigned",
        ):
            if definition_marker not in summary.get("primary_depth_definition", ""):
                failures.append(f"pipeline summary obscures primary total: {definition_marker}")
        for definition_marker in (
            "34 assigned raw-name fish columns excluding Unassigned",
            "mapped to the focal harmonised label",
        ):
            if definition_marker not in summary.get("sensitivity_depth_definition", ""):
                failures.append(f"pipeline summary obscures focal-excluded total: {definition_marker}")
    pipeline_by_term = {
        row.get("term", ""): row
        for row in pipeline_model
        if row.get("model") == "abundance_depth_primary"
    }
    pipeline_expected = {
        "z_log_source_ra": (-1.57389682597455, -2.94986805287256, -0.806717037664016),
        "z_log_depth": (-1.11238760339895, -1.97158410446502, -0.719453694460575),
    }
    for term, (estimate, lower, upper) in pipeline_expected.items():
        row = pipeline_by_term.get(term)
        if not row:
            failures.append(f"pipeline primary model lacks {term}")
            continue
        for field, wanted in (("estimate.x", estimate), ("bootstrap_q025", lower), ("bootstrap_q975", upper)):
            actual = number(row.get(field), f"pipeline {term} {field}")
            if abs(actual - wanted) > 1e-12:
                failures.append(f"pipeline {term} {field}: expected {wanted}, found {row.get(field)!r}")
        if row.get("gate_status") != "ASSOCIATION_PASSED" or not truthy(row.get("gate_pass")):
            failures.append(f"pipeline scoped association did not pass for {term}")
        exact(row.get("leaveout_sign_stability"), "1", f"pipeline {term} leave-out sign stability")
    primary_leaveouts = [
        row for row in pipeline_leaveout
        if row.get("term") in pipeline_expected
    ]
    exact(len(primary_leaveouts), 86, "pipeline primary-effect leave-out rows")
    if primary_leaveouts and any(
        not truthy(row.get("fit_ok"))
        or number(row.get("estimate"), "pipeline leave-out estimate") >= 0
        for row in primary_leaveouts
    ):
        failures.append("pipeline primary effect has a failed or non-negative leave-out fit")
    sensitivity_rows = [
        row for row in pipeline_model
        if row.get("model") == "focal_excluded_assigned_fish_sensitivity"
        and row.get("term") == "z_log_focal_excluded_assigned_fish_reads"
    ]
    exact(len(sensitivity_rows), 1, "pipeline focal-excluded sensitivity rows")
    if sensitivity_rows:
        row = sensitivity_rows[0]
        for field, wanted in (
            ("estimate.x", -1.14939008320629),
            ("bootstrap_q025", -3.08646314340788),
            ("bootstrap_q975", -0.579192692358123),
            ("leaveout_min", -1.7375914842437),
            ("leaveout_max", -1.05429213637694),
        ):
            actual = number(row.get(field), f"pipeline focal-excluded {field}")
            if abs(actual - wanted) > 1e-12:
                failures.append(
                    f"pipeline focal-excluded {field}: expected {wanted}, found {row.get(field)!r}"
                )
        for field in (
            "bootstrap_finite_fraction",
            "leaveout_sign_stability",
        ):
            exact(row.get(field), "1", f"pipeline focal-excluded {field}")
        if not truthy(row.get("gate_pass")) or not truthy(row.get("no_leaveout_reversal")):
            failures.append("pipeline focal-excluded sensitivity gate did not pass")

    figure_stems = (
        "impact_occupancy_atlas",
        "impact_salmo_multimarker",
        "impact_layerE_multiverse",
        "impact_trait_null_model",
        "impact_pipeline_disagreement",
    )
    for stem in figure_stems:
        png = ROOT / "figures" / f"{stem}.png"
        svg = ROOT / "figures" / f"{stem}.svg"
        if png.is_file():
            header = png.read_bytes()[:24]
            if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
                failures.append(f"figures/{stem}.png is not a valid PNG header")
            else:
                width = int.from_bytes(header[16:20], "big")
                height = int.from_bytes(header[20:24], "big")
                if width <= 0 or height <= 0:
                    failures.append(f"figures/{stem}.png has invalid dimensions {width}x{height}")
        if svg.is_file():
            svg_text = svg.read_text(encoding="utf-8")
            if "<svg" not in svg_text:
                failures.append(f"figures/{stem}.svg lacks an SVG root")
            if re.search(r"(?:href|xlink:href)\s*=\s*['\"](?:https?:)?//", svg_text):
                failures.append(f"figures/{stem}.svg has a network-dependent asset")

    # Canonical claim-survival register and accessible, fully local rendering.
    claims = rows("data/metadata/claim_survival.tsv")
    expected_claims = {
        "river-community": "SURVIVES_CONDITIONAL",
        "imperfect-detection": "WITHHELD",
        "pipeline-disagreement": "SURVIVES",
        "taxonomic-resolution": "WITHHELD",
        "trait-diversity": "FAILED_GATE",
    }
    exact(len(claims), 5, "claim-survival register rows")
    claims_by_id = {row.get("claim_id", ""): row for row in claims}
    exact(set(claims_by_id), set(expected_claims), "claim-survival IDs")
    for claim_id, status in expected_claims.items():
        row = claims_by_id.get(claim_id)
        if not row:
            continue
        exact(row.get("status"), status, f"claim-survival status {claim_id}")
        for field in ("machine_evidence", "methods", "figure"):
            raw = row.get(field, "")
            target = (ROOT / "website" / raw).resolve()
            try:
                target.relative_to(ROOT.resolve())
            except ValueError:
                failures.append(f"claim register {claim_id} {field} escapes the project: {raw}")
                continue
            if not target.is_file():
                failures.append(f"claim register {claim_id} has missing {field}: {raw}")
    claim_boundaries = {
        "river-community": ("does not establish a section centroid shift", "full-multiverse robustness label is withheld"),
        "imperfect-detection": ("cannot estimate Conwy occupancy", "do not repair or validate"),
        "pipeline-disagreement": (
            "Neither pipeline is biological ground truth",
            "associated with pipeline non-recovery",
            "Neither count total is total library sequencing depth",
            "focal label",
        ),
        "taxonomic-resolution": (
            "cannot establish ASV0001 species identity",
            "longer-marker species comparison is withheld",
            "all five S. salar records",
        ),
        "trait-diversity": ("says nothing about abundance", "neither was unusually large"),
    }
    for claim_id, phrases in claim_boundaries.items():
        row = claims_by_id.get(claim_id, {})
        combined = " ".join(row.values())
        for phrase in phrases:
            if phrase not in combined:
                failures.append(f"claim-survival {claim_id} lacks boundary wording: {phrase}")

    claim_page = html_parsers.get("website/claim-survival.html")
    if claim_page:
        exact(claim_page.html_attrs.get("lang"), "en", "claim map document language")
        for tag, expected in (("main", 1), ("article", 5), ("details", 5), ("summary", 5), ("button", 5), ("noscript", 1)):
            exact(claim_page.tag_counts.get(tag, 0), expected, f"claim map <{tag}> count")
        if any(button.get("type") != "button" or button.get("aria-pressed") not in {"true", "false"} for button in claim_page.buttons):
            failures.append("claim-map filters lack button type or aria-pressed state")
        external = [raw for _, raw in claim_page.links if urlsplit(raw).scheme or urlsplit(raw).netloc]
        if external:
            failures.append(f"claim map has network dependencies/links: {external}")
        claim_html = (ROOT / "website/claim-survival.html").read_text(encoding="utf-8")
        for marker in ('class="skip-link" href="#claim-map"', 'aria-live="polite"', "Failed gates are results", "<noscript>"):
            if marker not in claim_html:
                failures.append(f"claim map lacks accessibility/offline marker: {marker}")
        exact(claim_html.count('class="claim-path"'), 5, "claim map rendered claim paths")
        css_path = ROOT / "website/assets/style.css"
        if css_path.is_file():
            css = css_path.read_text(encoding="utf-8").lower()
            if "@import" in css or re.search(r"url\(\s*['\"]?(?:https?:)?//", css):
                failures.append("claim map stylesheet has a network dependency")

    for rel in (
        "README.md",
        "START_HERE.md",
        "report/conwy-reanalysis.qmd",
        "docs/preprint-draft.md",
        "website/research_brief.md",
    ):
        path = ROOT / rel
        if path.is_file():
            text = path.read_text(encoding="utf-8")
            if "claim-survival" not in text and "claim_survival" not in text:
                failures.append(f"{rel} does not reference the canonical claim-survival audit")

    return {
        "impact_required_artifacts": len(PHASE5_REQUIRED),
        "occupancy_atlas_rows": len(atlas),
        "salmo_marker_rows": len(marker_manifest),
        "layerE_families": len(layer_results),
        "trait_null_draws": len(trait_draws),
        "pipeline_cells": len(pipeline_cells),
        "claim_rows": len(claims),
    }


def main():
    failures: list[str] = []

    required = [
        "README.md",
        "eDNAplan.md",
        "START_HERE.md",
        "LICENSE",
        "CITATION.cff",
        "website/index.html",
        "website/deliverables.html",
        "website/research_brief.html",
        "website/research_brief.md",
        "website/research_brief.css",
        "website/graphical_abstract.svg",
        "website/visual-brief.html",
        "website/detection-lab.html",
        "website/poster.html",
        "website/data/community-gradient.csv",
        "website/assets/community-gradient.png",
        "website/assets/salmo-resolution.png",
        "website/assets/detection-simulation.png",
        "website/assets/research-timeline.png",
        "website/assets/models-overview.png",
        *MODEL_DIAGRAM_ASSETS,
        "figures/salmo_molecular_detective.svg",
        "figures/salmo_molecular_detective.png",
        "docs/fish-trait-ecology.md",
        "docs/salmo_molecular_detective_methods.md",
        "docs/phase3-claim-status.md",
        "scripts/build_biological_evidence.R",
        "scripts/apply_layerA_claim_gate.R",
        "scripts/package_website.py",
        "scripts/verify_release.py",
        *EDITOR_REQUIRED_ARTIFACTS,
        "report/conwy-reanalysis.qmd",
        "docs/preprint-draft.md",
        "data/metadata/fish_traits.tsv",
        "data/metadata/fish_traits_sources.tsv",
        "data/metadata/fish_trait_section_summary.tsv",
        "data/metadata/salmo_sequence_alignment.tsv",
        "data/metadata/control_reconciliation.md",
        "data/metadata/layerA_recovery_status.tsv",
        "data/metadata/layerA_figure_gate.tsv",
        "output/pdf/RESEARCH_BRIEF.pdf",
        "output/pdf/POSTER.pdf",
        "_site/conwy-reanalysis.html",
        *PHASE5_REQUIRED,
    ]
    for rel in required:
        if not (ROOT / rel).exists():
            failures.append(f"missing required artifact: {rel}")

    html_parsers = {}
    for rel in (
        "website/index.html",
        "website/editor.html",
        "website/deliverables.html",
        "website/research_brief.html",
        "website/visual-brief.html",
        "website/detection-lab.html",
        "website/claim-survival.html",
        "website/poster.html",
    ):
        path = ROOT / rel
        if path.exists():
            html_parsers[rel] = check_html(path, failures)

    markdown_paths = {
        *ROOT.glob("*.md"),
        *(ROOT / "website").glob("*.md"),
        *(ROOT / "docs").glob("*.md"),
        *(ROOT / "data" / "metadata").glob("*.md"),
    }
    for path in sorted(markdown_paths):
        check_markdown(path, failures)

    editor_counts = check_editor(failures, html_parsers)
    phase5_counts = check_phase5(failures, html_parsers)

    pdf_specs = WEBSITE_PDF_PAGES
    for rel, pages in pdf_specs.items():
        path = ROOT / rel
        if not path.exists():
            continue
        reader = PdfReader(path)
        if len(reader.pages) != pages:
            failures.append(f"{rel}: expected {pages} pages, found {len(reader.pages)}")
        pdf_text = "\n".join(page.extract_text() or "" for page in reader.pages)
        for stale in (
            "Layer E hypothesis-test p-values until clustered correction",
            "I re-examined",
            "+0.30",
            "+0.07",
        ):
            if stale in pdf_text:
                failures.append(f"{rel}: stale text present: {stale}")
    brief_path = ROOT / "output/pdf/RESEARCH_BRIEF.pdf"
    brief_pages: list[str] = []
    brief_link_targets: list[str] = []
    if brief_path.exists():
        brief_reader = PdfReader(brief_path)
        brief_pages = [page.extract_text() or "" for page in brief_reader.pages]
        for page in brief_reader.pages:
            for annotation_ref in page.get("/Annots") or []:
                annotation = annotation_ref.get_object()
                action = annotation.get("/A")
                if action and action.get("/S") == "/URI" and action.get("/URI"):
                    brief_link_targets.append(str(action.get("/URI")))

    # ReportLab emits selectable headings and captions separately from the raster
    # diagrams. Normalize typography and whitespace so the release gate checks
    # scientific content, not a particular minus glyph or text-extraction wrap.
    normalized_brief_pages = [
        re.sub(
            r"\s+",
            " ",
            re.sub(r"\s*=\s*", "=", page.replace("\u2212", "-").replace("\u2013", "-")),
        ).strip()
        for page in brief_pages
    ]
    brief_text = " ".join(normalized_brief_pages)
    if brief_text and "0.0714" not in brief_text:
        failures.append("RESEARCH_BRIEF.pdf lacks corrected Layer E p=0.0714 boundary")
    if brief_text:
        for marker in (
            "Phase-5",
            "64/72",
            "613/1,944",
            "2 of 7",
            "-1.149",
            "999/999",
            "43/43",
            "0.3159",
            "0.1142",
            "questions, tests, and verdict rules were written in a dated project plan before analysis",
        ):
            if marker not in brief_text:
                failures.append(f"RESEARCH_BRIEF.pdf lacks required evidence marker: {marker}")

    page_markers = {
        1: (
            "Question, evidence base, and research trail",
            "Our question",
            "Our methodology + key research documentation moments",
            "WHAT RE-RUNNING PROVES",
            "the same committed inputs produced the same files",
            "does not add new biological replication or independent field evidence",
        ),
        2: (
            "Five claims enter the audit",
            "Each hypothesis receives one decision label and one reason",
            "Our five initial claims",
            "Stressing those claims with more tests",
            "Decision tally | 1 accepted | 1 rejected | 3 non-conclusive",
        ),
        3: (
            "The modelling laboratory",
            "The models and simulations used to evaluate these claims",
            "MODEL 01 | NON-CONCLUSIVE | Reason: 64/72 valid (<90%); direction stable",
            "MODEL 02 | NON-CONCLUSIVE | Reason: real fit failed; atlas is synthetic",
            "MODEL 03 | ACCEPTED | Reason: 999/999 finite; 43/43 leave-outs negative",
        ),
        4: (
            "What the evidence decided - and what I added",
            "the same three-state decision rule, applied to identity and traits",
            "MODEL 04 | NON-CONCLUSIVE | Reason: 2 of 7 labels eligible; minimum 4",
            "MODEL 05 | REJECTED | Reason: p=0.3159 / 0.1142; both > 0.05",
            "Our final verdict",
            "ACCEPTED=rule met with adequate evidence",
            "REJECTED=adequate test did not support",
            "NON-CONCLUSIVE=insufficient to decide",
            "showing which claims are accepted, rejected, or non-conclusive - and why",
        ),
    }
    page_footer_markers = {
        1: "Rules written before testing | unsuccessful models and invalid analyses remain visible",
        2: "Decision tally | 1 accepted | 1 rejected | 3 non-conclusive",
        3: "M1-M3 methods, machine results, full figures and executable code are bundled offline",
        4: "Phase-5 methods, data, code, figures and clean-extraction checks bundled offline",
    }
    claim_decision_markers = (
        "River community NON-CONCLUSIVE Reason: 64/72 valid (<90%); all 64 valid directions positive",
        "Real occupancy NON-CONCLUSIVE Reason: Real fit failed; recovery atlas is synthetic only",
        "Pipeline association ACCEPTED Reason: 999/999 finite; all 43 leave-outs negative",
        "Species resolution NON-CONCLUSIVE Reason: Only 2 of 7 labels eligible; minimum was 4",
        "Trait broadening REJECTED Reason: p=0.3159 and 0.1142; both above 0.05",
    )
    if len(normalized_brief_pages) == 4:
        for page_number, markers in page_markers.items():
            page_text = normalized_brief_pages[page_number - 1]
            for marker in markers:
                if marker not in page_text:
                    failures.append(
                        f"RESEARCH_BRIEF.pdf page {page_number} lacks selectable marker: {marker}"
                    )
            footer = f"{page_number} / 4"
            if footer not in page_text:
                failures.append(
                    f"RESEARCH_BRIEF.pdf page {page_number} lacks page footer: {footer}"
                )
            footer_marker = page_footer_markers[page_number]
            if footer_marker not in page_text:
                failures.append(
                    f"RESEARCH_BRIEF.pdf page {page_number} lacks footer marker: "
                    f"{footer_marker}"
                )
        verdict_page = normalized_brief_pages[3]
        for marker in claim_decision_markers:
            if marker not in verdict_page:
                failures.append(
                    "RESEARCH_BRIEF.pdf page 4 lacks selectable claim decision: "
                    f"{marker}"
                )

    expected_model_links = {
        "../../website/index.html#model-community",
        "../../website/index.html#model-detection",
        "../../website/index.html#model-pipeline",
        "../../website/index.html#model-taxonomy",
        "../../website/index.html#model-traits",
    }
    missing_model_links = expected_model_links - set(brief_link_targets)
    if missing_model_links:
        failures.append(
            "RESEARCH_BRIEF.pdf lacks clickable model/evidence links: "
            + ", ".join(sorted(missing_model_links))
        )

    poster_path = ROOT / "output/pdf/POSTER.pdf"
    if poster_path.exists():
        poster_text = "\n".join(page.extract_text() or "" for page in PdfReader(poster_path).pages)
        if "10.1038/s41467-024-48640-3" not in poster_text:
            failures.append("POSTER.pdf lacks the source-study DOI")

    bundle_path = ROOT / "output/share/conwy-edna-project-website.zip"
    bundle_checksum_path = ROOT / "output/share/conwy-edna-project-website.zip.sha256"
    bundle_report_path = ROOT / "output/share/conwy-edna-project-website.bundle.json"
    actual_bundle_digest = None
    actual_bundle_bytes = None
    if bundle_path.exists():
        actual_bundle_digest = hashlib.sha256(bundle_path.read_bytes()).hexdigest()
        actual_bundle_bytes = bundle_path.stat().st_size
    if bundle_path.exists() and bundle_checksum_path.exists():
        checksum_parts = bundle_checksum_path.read_text(encoding="utf-8").strip().split()
        if len(checksum_parts) != 2 or checksum_parts[1] != bundle_path.name:
            failures.append("website bundle checksum sidecar is malformed")
        else:
            if actual_bundle_digest != checksum_parts[0]:
                failures.append("website bundle SHA-256 does not match its sidecar")
        try:
            with zipfile.ZipFile(bundle_path) as bundle:
                names = set(bundle.namelist())
                for required_member in (
                    "conwy-edna-project-website/OPEN_WEBSITE.html",
                    "conwy-edna-project-website/SHA256SUMS.txt",
                    "conwy-edna-project-website/VERIFY_BUNDLE.py",
                    "conwy-edna-project-website/EDIT_WEBSITE.command",
                    "conwy-edna-project-website/EDIT_WEBSITE.md",
                    "conwy-edna-project-website/website/index.html",
                    "conwy-edna-project-website/website/editor.html",
                    "conwy-edna-project-website/website/content-overrides.json",
                    "conwy-edna-project-website/website/assets/editor.css",
                    "conwy-edna-project-website/website/assets/editor.js",
                    "conwy-edna-project-website/website/assets/apply-content.js",
                    "conwy-edna-project-website/website/assets/content-overrides.js",
                    "conwy-edna-project-website/website/deliverables.html",
                    "conwy-edna-project-website/website/research_brief.html",
                    "conwy-edna-project-website/website/claim-survival.html",
                    "conwy-edna-project-website/website/assets/model-community-multiverse.png",
                    "conwy-edna-project-website/website/assets/model-occupancy-atlas.png",
                    "conwy-edna-project-website/website/assets/model-pipeline-bootstrap.png",
                    "conwy-edna-project-website/website/assets/model-salmo-benchmark.png",
                    "conwy-edna-project-website/website/assets/model-trait-null.png",
                    "conwy-edna-project-website/_site/conwy-reanalysis.html",
                    "conwy-edna-project-website/data/metadata/claim_survival.tsv",
                    "conwy-edna-project-website/data/metadata/impact_sprint_manifest.json",
                    "conwy-edna-project-website/data/public/salmo_multimarker/marker_sequences.fasta",
                    "conwy-edna-project-website/scripts/run_impact_sprint.py",
                    "conwy-edna-project-website/scripts/edit_website.py",
                ):
                    if required_member not in names:
                        failures.append(f"website bundle lacks required member: {required_member}")
                missing_phase5_members = [
                    rel for rel in PHASE5_REQUIRED
                    if f"conwy-edna-project-website/{rel}" not in names
                ]
                if missing_phase5_members:
                    failures.append(
                        "website bundle lacks required Phase-5 members: "
                        + ", ".join(missing_phase5_members)
                    )
                genbank_members = [
                    name for name in names
                    if name.startswith("conwy-edna-project-website/data/public/salmo_multimarker/genbank/")
                    and name.endswith(".gb")
                ]
                query_receipt_members = [
                    name for name in names
                    if name.startswith("conwy-edna-project-website/data/public/salmo_multimarker/query_receipts/")
                    and name.endswith(".json")
                ]
                if len(genbank_members) != 29:
                    failures.append(f"website bundle expected 29 cached GenBank records, found {len(genbank_members)}")
                if len(query_receipt_members) != 14:
                    failures.append(f"website bundle expected 14 cached query receipts, found {len(query_receipt_members)}")
                unsafe = [
                    name for name in names
                    if name.startswith("/") or ".." in name.split("/") or ".git" in name.split("/")
                ]
                if unsafe:
                    failures.append(f"website bundle has unsafe members: {sorted(unsafe)}")
        except zipfile.BadZipFile:
            failures.append("website bundle is not a valid ZIP archive")
    if bundle_report_path.exists():
        bundle_report = json.loads(bundle_report_path.read_text(encoding="utf-8"))
        if bundle_report.get("status") != "PASS":
            failures.append("website bundle report is not PASS")
        if bundle_report.get("archive") != "output/share/conwy-edna-project-website.zip":
            failures.append("website bundle report names a different archive")
        if actual_bundle_digest is not None and bundle_report.get("archive_sha256") != actual_bundle_digest:
            failures.append("website bundle report SHA-256 is not bound to the current archive")
        if actual_bundle_bytes is not None and bundle_report.get("archive_bytes") != actual_bundle_bytes:
            failures.append("website bundle report byte count is not bound to the current archive")
        extracted = bundle_report.get("clean_extraction_verification", {})
        if extracted.get("status") != "PASS" or extracted.get("failures"):
            failures.append("website bundle clean-extraction verification is not PASS")
        if extracted.get("path_contains_spaces") is not True:
            failures.append("website bundle was not verified from a clean extraction path containing spaces")
        bundled_verifier = bundle_report.get("bundled_verifier", {})
        if (
            bundled_verifier.get("status") != "PASS"
            or bundled_verifier.get("returncode") != 0
            or bundled_verifier.get("path_contains_spaces") is not True
        ):
            failures.append("bundled VERIFY_BUNDLE.py did not pass from the clean path containing spaces")
    else:
        failures.append("missing website bundle verification report")

    receipt_path = ROOT / "website/build-receipt.json"
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt.get("status") != "PASS":
            failures.append("website/build-receipt.json is not PASS")
        receipt_outputs = receipt.get("outputs")
        expected_outputs = list(WEBSITE_BUILD_OUTPUTS)
        if receipt_outputs != expected_outputs:
            failures.append(
                "website build receipt output contract differs: "
                f"expected {expected_outputs}, found {receipt_outputs}"
            )
        if isinstance(receipt_outputs, list) and len(receipt_outputs) != len(set(receipt_outputs)):
            failures.append("website build receipt contains duplicate output paths")
        if receipt.get("pdf_pages") != WEBSITE_PDF_PAGES:
            failures.append(
                "website build receipt PDF page contract differs: "
                f"expected {WEBSITE_PDF_PAGES}, found {receipt.get('pdf_pages')}"
            )
        if receipt.get("failures") != []:
            failures.append(
                f"website build receipt contains failures: {receipt.get('failures')}"
            )
        scenario = receipt.get("synthetic_scenario", {})
        params = scenario.get("parameters", {})
        expected_params = {
            "occupancy": 0.35,
            "detection": 0.25,
            "replicates": 3,
            "relative_depth": 1.0,
            "depth_gradient_coupling": 0.70,
            "seed": 20260718,
            "units": 800,
        }
        if params != expected_params:
            failures.append(f"build receipt synthetic parameters differ: {params}")
        if abs(float(scenario.get("correction_error", 99)) - 0.0710271687) > 1e-8:
            failures.append("build receipt synthetic correction does not match fixed-seed reference")
        panel = receipt.get("community_gradient_panel", {})
        if panel.get("sites") != 14 or len(panel.get("common_timepoints", [])) != 11:
            failures.append(f"build receipt community panel is not 14 sites x 11 common timepoints: {panel}")
        if panel.get("libraries_per_unit_range") != [1, 3]:
            failures.append(f"build receipt community effort range is not [1, 3]: {panel}")
    else:
        failures.append("missing website/build-receipt.json")

    gradient_path = ROOT / "website/data/community-gradient.csv"
    if gradient_path.exists():
        with gradient_path.open(newline="", encoding="utf-8") as handle:
            gradient_rows = list(csv.DictReader(handle))
        if len(gradient_rows) != 14:
            failures.append(f"community hero expected 14 sites, found {len(gradient_rows)}")
        if {row.get("n_common_timepoints") for row in gradient_rows} != {"11"}:
            failures.append("community hero is not restricted to 11 common timepoints")
        effort = {
            int(row[key])
            for row in gradient_rows
            for key in ("min_libraries", "max_libraries")
            if row.get(key)
        }
        if not effort or min(effort) != 1 or max(effort) != 3:
            failures.append(f"community hero does not disclose the verified 1-3-library effort range: {sorted(effort)}")
    else:
        failures.append("missing website/data/community-gradient.csv")

    lab = html_parsers.get("website/detection-lab.html")
    if lab:
        expected_defaults = {"occ": "0.35", "det": "0.25", "reps": "3", "depth": "1", "conf": "0.70"}
        for input_id, value in expected_defaults.items():
            actual = (lab.inputs.get(input_id) or {}).get("value")
            if actual != value:
                failures.append(f"simulator #{input_id} default {actual!r}, expected {value!r}")
        lab_text = (ROOT / "website/detection-lab.html").read_text(encoding="utf-8")
        for svg_id in ("bars", "miss"):
            if lab_text.count(f'id="{svg_id}-title"') < 2 or lab_text.count(f'id="{svg_id}-desc"') < 2:
                failures.append(f"simulator redraw does not preserve accessible {svg_id} title/description")
        for stale in ("Average true occupancy", "Depth-biological correlation", "fits detection as a function of depth"):
            if stale in lab_text:
                failures.append(f"simulator contains inaccurate parameter wording: {stale}")

    landing = html_parsers.get("website/index.html")
    if landing:
        hrefs = {raw for key, raw in landing.links if key == "href"}
        required_landing_links = (
            "research_brief.html",
            "deliverables.html",
            "../_site/conwy-reanalysis.html",
            "../START_HERE.md",
            "../scripts/build_website.py",
            "../data/metadata/data_dictionary.md",
            "claim-survival.html",
            "../docs/impact-sprint-preregistration.md",
            "../data/metadata/impact_sprint_manifest.json",
            "../data/metadata/impact_sprint_run_receipt.json",
            "../data/metadata/claim_survival.tsv",
            "../scripts/run_impact_sprint.py",
            "../scripts/verify_release.py",
            "../figures/impact_layerE_multiverse.svg",
            "../docs/impact_layerE_multiverse_methods.md",
            "../data/derived/impact_layerE_gate.tsv",
            "../R/13_layerE_multiverse.R",
            "../_site/conwy-reanalysis.html#layer-e-robustness-multiverse",
            "../figures/impact_occupancy_atlas.svg",
            "../docs/impact_occupancy_atlas_methods.md",
            "../data/derived/impact_occupancy_atlas.tsv",
            "../R/11_occupancy_recoverability_atlas.R",
            "../_site/conwy-reanalysis.html#synthetic-occupancy-recoverability-atlas",
            "detection-lab.html",
            "../figures/impact_pipeline_disagreement.svg",
            "../docs/impact_pipeline_disagreement_methods.md",
            "../data/derived/impact_pipeline_model.tsv",
            "../R/15_pipeline_disagreement.R",
            "../_site/conwy-reanalysis.html#cluster-aware-pipeline-disagreement",
            "../figures/impact_salmo_multimarker.svg",
            "../docs/impact_salmo_multimarker_methods.md",
            "../data/derived/impact_salmo_marker_gates.tsv",
            "../R/12_salmo_multimarker_benchmark.R",
            "../_site/conwy-reanalysis.html#accession-backed-salmo-multi-marker-benchmark",
            "../figures/impact_trait_null_model.svg",
            "../docs/impact_trait_null_methods.md",
            "../data/derived/impact_trait_null_summary.tsv",
            "../R/14_trait_null_model.R",
            "../_site/conwy-reanalysis.html#functional-trait-fixed-margin-null",
        )
        for required_link in required_landing_links:
            if required_link not in hrefs:
                failures.append(f"website landing lacks explicit technical link: {required_link}")

        landing_text = (ROOT / "website/index.html").read_text(encoding="utf-8")
        expected_sections = [
            "question",
            "methodology",
            "initial-claims",
            "stress-tests",
            "models",
            "verdict",
        ]
        actual_sections = re.findall(r'<section id="([^"]+)"', landing_text)
        if actual_sections != expected_sections:
            failures.append(
                f"website landing section order is {actual_sections}, expected {expected_sections}"
            )
        for tag, expected in (("main", 1), ("h1", 1), ("h2", 6)):
            actual = landing.tag_counts.get(tag, 0)
            if actual != expected:
                failures.append(f"website landing has {actual} <{tag}> elements, expected {expected}")
        for heading_marker in (
            "Our question:",
            "Our methodology:",
            "Our initial claims:",
            "Stressing those claims:",
            "The models and simulations used:",
            "Our final verdict:",
        ):
            if heading_marker not in landing_text:
                failures.append(f"website landing lacks six-part heading: {heading_marker}")

        required_visuals = {
            "graphical_abstract.svg",
            "assets/research-timeline.png",
            "assets/claim-survival.png",
            "assets/models-overview.png",
        }
        landing_image_sources = {attrs.get("src") for attrs in landing.images}
        missing_visuals = sorted(required_visuals - landing_image_sources)
        if missing_visuals:
            failures.append(f"website landing lacks required overview visuals: {missing_visuals}")
        for index, attrs in enumerate(landing.images, start=1):
            if not (attrs.get("alt") or "").strip():
                failures.append(f"website landing image {index} lacks alternative text")

        expected_claim_statuses = {
            "river-community": "SURVIVES_CONDITIONAL",
            "imperfect-detection": "WITHHELD",
            "pipeline-disagreement": "SURVIVES",
            "taxonomic-resolution": "WITHHELD",
            "trait-diversity": "FAILED_GATE",
        }
        actual_claim_statuses = dict(re.findall(
            r'data-claim-id="([^"]+)"\s+data-status="([^"]+)"',
            landing_text,
        ))
        if actual_claim_statuses != expected_claim_statuses:
            failures.append(
                f"website landing claim/status contract is {actual_claim_statuses}, "
                f"expected {expected_claim_statuses}"
            )

        expected_claim_decisions = {
            "river-community": "NON_CONCLUSIVE",
            "imperfect-detection": "NON_CONCLUSIVE",
            "pipeline-disagreement": "ACCEPTED",
            "taxonomic-resolution": "NON_CONCLUSIVE",
            "trait-diversity": "REJECTED",
        }
        actual_claim_decisions = dict(re.findall(
            r'data-claim-id="([^"]+)"[^>]*\bdata-decision="([^"]+)"',
            landing_text,
        ))
        if actual_claim_decisions != expected_claim_decisions:
            failures.append(
                f"website landing claim/decision contract is {actual_claim_decisions}, "
                f"expected {expected_claim_decisions}"
            )
        decision_tally = Counter(actual_claim_decisions.values())
        expected_decision_tally = Counter({
            "ACCEPTED": 1,
            "REJECTED": 1,
            "NON_CONCLUSIVE": 3,
        })
        if decision_tally != expected_decision_tally:
            failures.append(
                f"website landing decision tally is {dict(decision_tally)}, "
                f"expected {dict(expected_decision_tally)}"
            )

        landing_plain = " ".join(" ".join(landing.text_chunks).split())
        for marker in (
            "1 accepted",
            "1 rejected",
            "3 non-conclusive",
            "64/72",
            "613/1,944",
            "−1.149",
            "2 of 7",
            "0.3159",
            "0.1142",
            "p = 0.0714",
            "not a biological false-negative rate",
            "stable descriptive downstream broadening",
            "Layer A coefficients or survey utility remain not claimed",
        ):
            if marker.lower() not in landing_plain.lower():
                failures.append(f"website landing lacks claim-boundary marker: {marker}")

    research_brief = html_parsers.get("website/research_brief.html")
    if research_brief:
        for tag, expected in (("h1", 1), ("h2", 6), ("blockquote", 1), ("table", 4), ("img", 4)):
            actual = research_brief.tag_counts.get(tag, 0)
            if actual != expected:
                failures.append(f"research brief has {actual} <{tag}> elements, expected {expected}")
        for index, attrs in enumerate(research_brief.images, start=1):
            if not (attrs.get("alt") or "").strip():
                failures.append(f"research brief image {index} lacks alternative text")
            if not (attrs.get("src") or "").startswith("data:image/"):
                failures.append(f"research brief image {index} is not embedded")
        brief_html = (ROOT / "website/research_brief.html").read_text(encoding="utf-8")
        brief_text = " ".join(" ".join(research_brief.text_chunks).split())
        for marker in (
            "Our question",
            "Our methodology + key research documentation moments",
            "Our initial claims from the first part of the project",
            "Stressing those claims with more tests",
            "The models and simulations used to evaluate the claims",
            "Our final verdict",
            "p = 0.0714",
            "WITHHELD_UNCONVERGED",
            "170-bp",
            "1,348",
            "613/1,944",
            "64/72",
            "2 of 7",
            "−1.149",
            "0.3159",
            "not a biological false-negative rate",
            "not a publication",
        ):
            if marker not in brief_text:
                failures.append(f"research brief lacks required claim-boundary marker: {marker}")
        for marker in ("@page", "page-break-before"):
            if marker not in brief_html:
                failures.append(f"research brief lacks embedded print-layout marker: {marker}")
        hrefs = {raw for key, raw in research_brief.links if key == "href"}
        for required_link in (
            "index.html",
            "detection-lab.html",
            "deliverables.html",
            "../output/pdf/RESEARCH_BRIEF.pdf",
            "../_site/conwy-reanalysis.html",
            "../docs/phase2-validation-gate.md",
            "../docs/phase3-claim-status.md",
            "../START_HERE.md",
            "claim-survival.html",
            "../docs/impact-sprint-preregistration.md",
        ):
            if required_link not in hrefs:
                failures.append(f"research brief lacks evidence or route link: {required_link}")

    catalogue = html_parsers.get("website/deliverables.html")
    if catalogue:
        hrefs = {raw for key, raw in catalogue.links if key == "href"}
        expected_catalogue_links = {
            "index.html",
            "research_brief.html",
            "research_brief.md",
            "claim-survival.html",
            "../data/metadata/claim_survival.tsv",
            "visual-brief.html",
            "../output/pdf/RESEARCH_BRIEF.pdf",
            "poster.html",
            "../output/pdf/POSTER.pdf",
            "detection-lab.html",
            "../figures/salmo_molecular_detective.svg",
            "../docs/salmo_molecular_detective_methods.md",
            "../docs/fish-trait-ecology.md",
            "../_site/conwy-reanalysis.html",
            "../report/conwy-reanalysis.qmd",
            "../docs/preprint-draft.md",
            "../docs/phase2-validation-gate.md",
            "../docs/phase3-claim-status.md",
            "../START_HERE.md",
            "../README.md",
            "../eDNAplan.md",
            "../scripts/build_website.py",
            "../scripts/build_biological_evidence.R",
            "../scripts/package_website.py",
            "../scripts/verify_release.py",
            "../scripts/run_impact_sprint.py",
            "../docs/impact-sprint-preregistration.md",
            "../docs/impact_occupancy_atlas_methods.md",
            "../docs/impact_salmo_multimarker_methods.md",
            "../docs/impact_layerE_multiverse_methods.md",
            "../docs/impact_trait_null_methods.md",
            "../docs/impact_pipeline_disagreement_methods.md",
            "../data/derived/impact_occupancy_atlas.tsv",
            "../data/derived/impact_salmo_marker_gates.tsv",
            "../data/derived/impact_layerE_gate.tsv",
            "../data/derived/impact_trait_null_summary.tsv",
            "../data/derived/impact_pipeline_model.tsv",
            "build-receipt.json",
            "../shiny/app.R",
            "../docs/prereg.md",
            "../CITATION.cff",
            "../LICENSE",
        }
        missing_catalogue = sorted(expected_catalogue_links - hrefs)
        if missing_catalogue:
            failures.append(f"deliverables catalogue lacks required links: {missing_catalogue}")

    traits_path = ROOT / "data/metadata/fish_traits.tsv"
    if traits_path.exists() and len(tsv_rows(traits_path)) != 14:
        failures.append("fish_traits.tsv must have exactly 14 evidence rows")

    alignment_path = ROOT / "data/metadata/salmo_sequence_alignment.tsv"
    if alignment_path.exists():
        alignment = tsv_rows(alignment_path)
        if len(alignment) != 9:
            failures.append(f"Salmo alignment expected 9 rows, found {len(alignment)}")
        sequence_keys = ("aligned_amplicon_5to3", "aligned_sequence", "sequence")
        sequence_key = next(
            (key for key in sequence_keys if alignment and key in alignment[0]),
            "",
        )
        if not sequence_key:
            failures.append("Salmo alignment has no recognised sequence column")
        lengths = {len(row.get(sequence_key, "").replace("-", "")) for row in alignment}
        if lengths != {170}:
            failures.append(f"Salmo aligned sequence lengths are {sorted(lengths)}")

    current_files = [
        ROOT / "README.md",
        ROOT / "START_HERE.md",
        ROOT / "eDNAplan.md",
        ROOT / "report/conwy-reanalysis.qmd",
        ROOT / "docs/preprint-draft.md",
        ROOT / "docs/phase2-validation-gate.md",
        ROOT / "docs/phase3-readiness.md",
        ROOT / "docs/phase3-claim-status.md",
        ROOT / "docs/prereg.md",
        ROOT / "docs/salmo_molecular_detective_methods.md",
        ROOT / "R/05_occupancy_msocc.R",
        ROOT / "R/06_occupancy_occumb.R",
        ROOT / "R/07_reclass_B.R",
        ROOT / "R/09_e14_covariates.R",
        ROOT / "R/10_sim_recovery.R",
        ROOT / "data/metadata/data_dictionary.md",
        ROOT / "data/metadata/covariate_join_verification.md",
        ROOT / "data/metadata/salmo_classifier_provenance.tsv",
        ROOT / "data/derived/reclassB_drift_catalogue.tsv",
        ROOT / "data/derived/reclassB_partB_provenance.txt",
        ROOT / "data/derived/reclassB_partC_provenance.txt",
        ROOT / "shiny/app.R",
        ROOT / "scripts/build_website.py",
        *list((ROOT / "website").glob("*.html")),
    ]
    forbidden = {
        "+0.30": "unconverged Layer A magnitude",
        "+0.07": "unconverged Layer A comparator magnitude",
        "ESTABLISHED SYNTHETIC FAILURE MODE": "stale Layer A claim label",
        "Layer E hypothesis-test p-values until clustered correction": "stale Layer E status",
        "PERMANOVA by river\nsection is highly significant": "stale Layer E conclusion",
        "per shared library": "incorrect unit for pooled detection-cell deficit",
        "per-library recall": "incorrect unit for pooled detection-cell deficit",
        "per-library detection-cell": "incorrect unit for pooled detection-cell deficit",
        "per-replicate detection sensitivity": "incorrect unit for pooled detection-cell deficit",
        "depth-biomass": "unsupported biomass confound",
        "biomass gradient": "unsupported biomass gradient",
        "HadUK-Grid": "unsupported E14 climate provenance",
        "biogeographically impossible": "unsupported geographic certainty",
        "geographically implausible": "unsupported local-identity verdict",
        "genuine trace stray": "unsupported trace-label verdict",
        "DB-naming artifact": "unsupported deterministic taxonomy verdict",
        "no hidden native fish": "overclaim from bounded Unassigned audit",
        "no hidden freshwater fish": "overclaim from bounded Unassigned audit",
        "large upward corrections": "unconverged Layer A direction",
        "all-fish-zero": "incorrect Layer A harmonised-zero denominator",
        "all-zero comparison mask": "ambiguous Layer A harmonised-zero denominator",
        "Original pipeline output": "mislabelled local Nov-2025 classifier provenance",
        "well-detected species": "analysis labels presented as independently validated species",
        "into their true species": "source-study relabels presented as biological truth",
        "fixed before any model fit": "unsupported pre-fit timing claim",
        "pre-fit substitutions": "unsupported pre-fit timing claim",
        "predeclared multi-method": "unsupported decision-timing claim",
        "depth_biology_correlation": "inaccurate simulator parameter name",
        "not fully final": "stale Phase-3 status",
        "Fukaya 2025": "incorrect occumb citation year",
    }
    for path in current_files:
        text = path.read_text(encoding="utf-8")
        plain_text = re.sub(r"[*_`]", "", text)
        for needle, meaning in forbidden.items():
            if needle in text or needle in plain_text:
                failures.append(f"{path.relative_to(ROOT)} contains {meaning}: {needle}")

    r10_text = (ROOT / "R/10_sim_recovery.R").read_text(encoding="utf-8")
    self_check_tail = r10_text.split("# ---- self-check", 1)[-1]
    if "if (convergence_pass) {" not in self_check_tail or "CLAIM WITHHELD" not in self_check_tail:
        failures.append("R/10 console diagnostics are not guarded by the convergence claim gate")
    r06_text = (ROOT / "R/06_occupancy_occumb.R").read_text(encoding="utf-8")
    if "write_layerA_withheld_figures(layerA_status)" not in r06_text:
        failures.append("R/06 rerun path can resurrect Layer-A numerical figures")

    drift_path = ROOT / "data/derived/reclassB_drift_catalogue.tsv"
    if drift_path.exists():
        drift_rows = tsv_rows(drift_path)
        required_audit_fields = {"audit_category", "audit_hypothesis", "evidence_flag", "audit_note"}
        fields = set(drift_rows[0]) if drift_rows else set()
        if not required_audit_fields <= fields:
            failures.append("Layer B drift catalogue lacks neutral audit schema")
        if {"corrected_label", "flag", "note"} & fields:
            failures.append("Layer B drift catalogue still exposes verdict/corrected-identity schema")
    for rel in ("data/derived/reclassB_partB_provenance.txt", "data/derived/reclassB_partC_provenance.txt"):
        path = ROOT / rel
        if path.exists() and "cached-results" not in path.read_text(encoding="utf-8"):
            failures.append(f"{rel} does not disclose cached-results status")

    status_path = ROOT / "data/metadata/layerA_recovery_status.tsv"
    if status_path.exists():
        status_rows = tsv_rows(status_path)
        if len(status_rows) != 1 or status_rows[0].get("claim_status") != "WITHHELD_UNCONVERGED":
            failures.append("Layer A recovery receipt must be one WITHHELD_UNCONVERGED row")

    figure_gate_path = ROOT / "data/metadata/layerA_figure_gate.tsv"
    if figure_gate_path.exists():
        gate_rows = tsv_rows(figure_gate_path)
        if len(gate_rows) != 5:
            failures.append(f"Layer A figure gate expected 5 rows, found {len(gate_rows)}")
        for row in gate_rows:
            rel = row.get("figure", "")
            path = ROOT / rel
            if row.get("claim_status") != "WITHHELD_UNCONVERGED":
                failures.append(f"Layer A figure gate has non-withheld row: {rel}")
            if not path.exists():
                failures.append(f"Layer A gated figure is missing: {rel}")
                continue
            digest = hashlib.md5(path.read_bytes()).hexdigest()
            if digest != row.get("md5"):
                failures.append(f"Layer A gated figure changed after suppression: {rel}")

    site_path = ROOT / "_site/conwy-reanalysis.html"
    report_source = ROOT / "report/conwy-reanalysis.qmd"
    if site_path.exists():
        site_text = site_path.read_text(encoding="utf-8")
        for required_text in (
            "p = 0.0714",
            "WITHHELD_UNCONVERGED",
            "1,944 estimator",
            "613 as recoverable",
            "64/72",
            "0.3159",
            "claim-survival.html",
        ):
            if required_text not in site_text:
                failures.append(f"rendered report lacks current evidence marker: {required_text}")
        stale_site_phrases = (
            "highly significant",
            "Adriatic endemic",
            "Across three databases",
            "Phase 3's multi-chain `occumb` firms this up",
            "Phase 3’s multi-chain occumb firms this up",
            "the WAIC <em>ranking</em> is safe",
        )
        for stale in stale_site_phrases:
            if stale in site_text:
                failures.append(f"rendered report contains stale claim: {stale}")
        if report_source.exists() and site_path.stat().st_mtime < report_source.stat().st_mtime:
            failures.append("rendered report is older than report/conwy-reanalysis.qmd")

    result = {
        "status": "PASS" if not failures else "FAIL",
        "checks": {
            "required_artifacts": len(required),
            "html_pages": len(html_parsers),
            "pdfs": len(pdf_specs),
            "trait_rows": len(tsv_rows(traits_path)) if traits_path.exists() else 0,
            "salmo_alignment_rows": len(tsv_rows(alignment_path)) if alignment_path.exists() else 0,
            **editor_counts,
            **phase5_counts,
        },
        "failures": failures,
    }
    print(json.dumps(result, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
