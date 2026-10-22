#!/usr/bin/env python3
"""Build the deterministic Phase-5 claim register and offline reader view.

The five records are not a second analysis layer.  They are a checked synthesis
of the machine-readable gate outputs created by the Phase-5 analyses.
"""

from __future__ import annotations

import csv
import html
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DERIVED = ROOT / "data" / "derived"
META = ROOT / "data" / "metadata"
OUT_TSV = META / "claim_survival.tsv"
OUT_HTML = ROOT / "website" / "claim-survival.html"

FIELDS = [
    "claim_id",
    "initial_claim",
    "credibility_audit",
    "gate_result",
    "status",
    "surviving_wording",
    "limitation",
    "next_evidence",
    "machine_evidence",
    "methods",
    "figure",
]
ALLOWED_STATUSES = {
    "SURVIVES",
    "SURVIVES_CONDITIONAL",
    "WITHHELD",
    "FAILED_GATE",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"Claim-survival build gate failed: {message}")


def checked_records() -> list[dict[str, str]]:
    atlas = read_tsv(DERIVED / "impact_occupancy_atlas.tsv")
    require(len(atlas) == 648 * 3, f"occupancy atlas has {len(atlas)} rows, expected 1,944")
    atlas_classes = Counter(row["gate_class"] for row in atlas)
    require(set(atlas_classes) == {"recoverable", "unreliable", "non-identifiable"},
            f"unexpected atlas classes: {sorted(atlas_classes)}")

    layer_gate = {row["criterion"]: row for row in read_tsv(DERIVED / "impact_layerE_gate.tsv")}
    require(layer_gate["overall_gate"]["pass"] == "FALSE", "Layer-E overall gate must remain failed")
    require(layer_gate["valid_specification_families"]["numerator"] == "64", "Layer-E valid numerator changed")
    require(layer_gate["valid_specification_families"]["denominator"] == "72", "Layer-E valid denominator changed")
    require(layer_gate["positive_valid_family_contrasts"]["numerator"] == "64", "Layer-E sign numerator changed")
    require(layer_gate["positive_one_library_draw_contrasts"]["numerator"] == "3600", "Layer-E draw numerator changed")

    salmo = read_tsv(DERIVED / "impact_salmo_marker_gates.tsv")
    require(len(salmo) == 3, f"Salmo benchmark has {len(salmo)} markers, expected three")
    require(all(row["species_gate_status"] == "coverage-inadequate" for row in salmo),
            "a species gate was unexpectedly promoted")
    require(all(row["candidate_salmo_labels"] == "7" for row in salmo),
            "Salmo candidate-label denominator changed")
    require(all(row["eligible_salmo_labels"] == "2" and
                row["required_salmo_labels_for_species"] == "4" for row in salmo),
            "Salmo species-coverage numerator or threshold changed")
    require(all(row["genus_gate_status"] == "pass" for row in salmo),
            "a genus gate did not pass")
    require(all(row["genus_loo_queries"] == "10" and row["salmo_queries_nearer_than_each_outgroup"] == "10" for row in salmo),
            "Salmo genus-boundary denominator changed")

    trait = {row["analysis"]: row for row in read_tsv(DERIVED / "impact_trait_null_summary.tsv")}
    require(set(trait) == {"primary_12_labels", "strict_5_labels"}, "trait sensitivity panels changed")
    require(all(row["gate_pass"] == "FALSE" for row in trait.values()), "trait claim was unexpectedly promoted")
    require(all(row["null_draws"] == "9999" for row in trait.values()), "trait null denominator changed")

    pipeline_summary = read_tsv(DERIVED / "impact_pipeline_summary.tsv")
    require(len(pipeline_summary) == 1, "pipeline summary must contain one analysis row")
    pipeline_summary = pipeline_summary[0]
    require(pipeline_summary["source_positive_cells"] == "2344", "pipeline cell denominator changed")
    require(pipeline_summary["nonrecovered_cells"] == "382", "pipeline disagreement numerator changed")
    require(pipeline_summary["finite_bootstrap_fits"] == "999", "pipeline bootstrap denominator changed")
    pipeline_terms = {
        row["term"]: row
        for row in read_tsv(DERIVED / "impact_pipeline_model.tsv")
        if row["model"] == "abundance_depth_primary" and row["term"] in {"z_log_source_ra", "z_log_depth"}
    }
    require(set(pipeline_terms) == {"z_log_source_ra", "z_log_depth"}, "pipeline primary terms changed")
    require(all(row["gate_pass"] == "TRUE" for row in pipeline_terms.values()), "pipeline association gate failed")
    pipeline_sensitivity = [
        row for row in read_tsv(DERIVED / "impact_pipeline_model.tsv")
        if row["model"] == "focal_excluded_assigned_fish_sensitivity"
        and row["term"] == "z_log_focal_excluded_assigned_fish_reads"
    ]
    require(len(pipeline_sensitivity) == 1, "pipeline focal-excluded sensitivity row changed")
    pipeline_sensitivity = pipeline_sensitivity[0]
    require(pipeline_sensitivity["gate_pass"] == "TRUE" and
            pipeline_sensitivity["no_leaveout_reversal"] == "TRUE" and
            pipeline_sensitivity["bootstrap_finite_fraction"] == "1",
            "pipeline focal-excluded assigned-fish sensitivity gate failed")
    require(pipeline_summary["depth_association_sensitivity_status"] ==
            "SURVIVES_FOCAL_EXCLUSION",
            "pipeline summary did not retain the focal-exclusion sensitivity")

    records = [
        {
            "claim_id": "river-community",
            "initial_claim": "Detected fish analysis-label richness broadens downstream in the balanced River Conwy series.",
            "credibility_audit": "A 72-family multiverse varied panel, unit, taxonomy, threshold, effort handling and distance; it also retained every declared one-library draw and reference leave-out.",
            "gate_result": "Full robustness gate failed: 64/72 families were valid (88.9%, below 90%). Direction stayed positive in 64/64 valid families and 3,600/3,600 draws; 0/75 reference leave-outs reversed it.",
            "status": "SURVIVES_CONDITIONAL",
            "surviving_wording": "The downstream descriptive direction is highly stable across valid specifications, but the preregistered full-multiverse robustness label is withheld; the restricted reference PERMANOVA remains p = 2/28 = 0.0714.",
            "limitation": "Eight rCLR families were undefined because some real 154-unit communities were empty. Detected-label richness is not true richness, and the analysis does not establish a section centroid shift or causation.",
            "next_evidence": "Predefine an empty-community treatment or obtain balanced biological observations that separate ecological structure from library effort.",
            "machine_evidence": "../data/derived/impact_layerE_gate.tsv",
            "methods": "../docs/impact_layerE_multiverse_methods.md",
            "figure": "../figures/impact_layerE_multiverse.svg",
        },
        {
            "claim_id": "imperfect-detection",
            "initial_claim": "A correction model can recover the true occupancy of rare Conwy labels.",
            "credibility_audit": "A 648-scenario, 100-simulation known-truth atlas compared naive occurrence, replication-aware and depth-aware estimators under depth coupling, missingness and mask error.",
            "gate_result": f"Synthetic atlas: {atlas_classes['recoverable']:,}/1,944 estimator-scenarios recoverable, {atlas_classes['unreliable']:,} unreliable and {atlas_classes['non-identifiable']:,} non-identifiable. The real Conwy recovery fit still fails its convergence/mask gate.",
            "status": "WITHHELD",
            "surviving_wording": "Known-truth simulations identify design regions where fast occupancy estimators can or cannot recover their target; they do not repair or validate the unconverged real Conwy correction.",
            "limitation": "The atlas is synthetic and conditional on its generator. It cannot estimate Conwy occupancy, fill unsurveyed cells, or convert a failed real fit into biological evidence.",
            "next_evidence": "Collect a design with independent biological truth, balanced replication and depth separated from the river gradient, then validate recovery before interpreting rare-label coefficients.",
            "machine_evidence": "../data/derived/impact_occupancy_atlas.tsv",
            "methods": "../docs/impact_occupancy_atlas_methods.md",
            "figure": "../figures/impact_occupancy_atlas.svg",
        },
        {
            "claim_id": "pipeline-disagreement",
            "initial_claim": "Source-positive cells are more likely to be absent from the reprocessed table when source signal or the reprocessed raw-name-table count total is low.",
            "credibility_audit": "The analysis used exactly 2,344 source-positive cells, paired 999-draw site/time/taxon block bootstraps and 43 leave-outs. Because the preregistered 35-column total includes Unassigned and the focal count, a focal-excluded 34-column assigned-fish sensitivity removed that direct part-whole coupling.",
            "gate_result": "Source signal passed in the preregistered model: coefficient -1.574 (bootstrap 95% interval -2.950 to -0.807). The focal-excluded assigned-fish total also passed: -1.149 (-3.086 to -0.579), 999/999 finite fits and 43/43 leave-outs negative.",
            "status": "SURVIVES",
            "surviving_wording": "Within the verified shared-library comparison, lower source-study relative abundance and a lower reprocessed assigned-fish count total excluding the focal label are associated with pipeline non-recovery.",
            "limitation": "Neither count total is total library sequencing depth, and both are processed, potentially endogenous quantities. The outcome is disagreement relative to the source table. Neither pipeline is biological ground truth, so this is not a causal effect, abundance effect or biological false-negative rate.",
            "next_evidence": "Use mock communities, technical replicates and independently quantified templates to partition extraction, sequencing, filtering and assignment losses.",
            "machine_evidence": "../data/derived/impact_pipeline_model.tsv",
            "methods": "../docs/impact_pipeline_disagreement_methods.md",
            "figure": "../figures/impact_pipeline_disagreement.svg",
        },
        {
            "claim_id": "taxonomic-resolution",
            "initial_claim": "Longer mitochondrial markers consistently resolve the relevant public Salmo species labels more finely than the MiFish fragment.",
            "credibility_audit": "Accession-backed MiFish 12S, COI and cyt-b panels measured exact sharing, barcode gaps and all-ties leave-one-accession-out classification, with two outgroup genera and explicit coverage gates.",
            "gate_result": "Species comparison was coverage-inadequate for every marker: 2 of 7 candidate Salmo labels were eligible; the preregistered minimum was 4. In that narrowed two-label diagnostic, all 10/10 queries per marker retained their submitted label; the separate genus boundary passed 10/10 against each outgroup for all three markers.",
            "status": "WITHHELD",
            "surviving_wording": "The longer-marker species comparison is withheld. A separate, narrower record-level result is that MiFish 12S, COI and cyt-b each passed the audited Salmo-versus-outgroup genus boundary.",
            "limitation": "Seven of 21 selected records, including all five S. salar records, lack voucher/isolate/strain source keys, so their specimen independence is unverified. Submitted database labels are annotations, not truth; the benchmark cannot establish ASV0001 species identity, local occurrence, stock, ancestry, hybrid status or abundance.",
            "next_evidence": "Add independently curated, voucher-backed records for at least four relevant Salmo labels and re-run the frozen all-ties species gate.",
            "machine_evidence": "../data/derived/impact_salmo_marker_gates.tsv",
            "methods": "../docs/impact_salmo_multimarker_methods.md",
            "figure": "../figures/impact_salmo_multimarker.svg",
        },
        {
            "claim_id": "trait-diversity",
            "initial_claim": "Downstream functional-trait diversity broadens more than expected from detected-label richness and prevalence alone.",
            "credibility_audit": "Rao's Q used a documented equal-axis Gower matrix and 9,999 fixed-margin checkerboard null draws, preserving unit richness and label prevalence within each timepoint; a strict five-label sensitivity repeated the test.",
            "gate_result": f"Failed in both panels: primary one-sided p = {float(trait['primary_12_labels']['monte_carlo_p']):.4f}; strict sensitivity p = {float(trait['strict_5_labels']['monte_carlo_p']):.4f}.",
            "status": "FAILED_GATE",
            "surviving_wording": "Both observed downstream trait-diversity contrasts were slightly above their null means, but neither was unusually large under the fixed-margin null; the proposed broadening beyond richness and prevalence is unsupported.",
            "limitation": "The result is conditional on analysis-label assignments and cited trait proxies. It says nothing about abundance, diet, realised habitat use, causal filtering or independently verified occurrence.",
            "next_evidence": "Use verified species occurrences, population-specific traits and independent abundance or habitat measurements before proposing a causal functional mechanism.",
            "machine_evidence": "../data/derived/impact_trait_null_summary.tsv",
            "methods": "../docs/impact_trait_null_methods.md",
            "figure": "../figures/impact_trait_null_model.svg",
        },
    ]
    require([row["claim_id"] for row in records] == [
        "river-community", "imperfect-detection", "pipeline-disagreement",
        "taxonomic-resolution", "trait-diversity"
    ], "required claim IDs changed")
    require(all(row["status"] in ALLOWED_STATUSES for row in records), "invalid status")
    return records


def write_tsv(records: list[dict[str, str]]) -> None:
    with OUT_TSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)


def esc(value: str) -> str:
    return html.escape(value, quote=True)


def claim_article(row: dict[str, str], number: int) -> str:
    status_label = row["status"].replace("_", " ").title()
    return f'''      <article class="claim-path" data-status="{esc(row['status'])}">
        <header class="claim-path-heading">
          <p class="number">{number:02d} / {esc(row['claim_id'].replace('-', ' ').upper())}</p>
          <span class="claim-status status-{esc(row['status'].lower().replace('_', '-'))}">{esc(status_label)}</span>
        </header>
        <div class="claim-flow" aria-label="Claim evidence path">
          <section><p class="claim-step">Initial claim</p><p>{esc(row['initial_claim'])}</p></section>
          <span class="claim-arrow" aria-hidden="true">&#8594;</span>
          <section><p class="claim-step">Credibility audit</p><p>{esc(row['credibility_audit'])}</p></section>
          <span class="claim-arrow" aria-hidden="true">&#8594;</span>
          <section><p class="claim-step">Surviving wording</p><p><strong>{esc(row['surviving_wording'])}</strong></p></section>
        </div>
        <details>
          <summary>Inspect gate, limitation and evidence</summary>
          <dl class="claim-detail">
            <div><dt>Gate result</dt><dd>{esc(row['gate_result'])}</dd></div>
            <div><dt>Limitation</dt><dd>{esc(row['limitation'])}</dd></div>
            <div><dt>Evidence still needed</dt><dd>{esc(row['next_evidence'])}</dd></div>
          </dl>
          <p class="link-row"><a href="{esc(row['machine_evidence'])}">Machine-readable evidence</a><a href="{esc(row['methods'])}">Methods</a><a href="{esc(row['figure'])}">Figure</a></p>
        </details>
      </article>'''


def write_html(records: list[dict[str, str]]) -> None:
    articles = "\n".join(claim_article(row, i) for i, row in enumerate(records, 1))
    document = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="Offline evidence map tracing five River Conwy eDNA claims through credibility audits, gates, limitations and surviving wording.">
  <title>Claim-survival map | River Conwy eDNA</title>
  <link rel="stylesheet" href="assets/style.css">
</head>
<body>
  <a class="skip-link" href="#claim-map">Skip to the claim paths</a>
  <header>
    <nav class="site-nav" aria-label="Project website">
      <a class="brand" href="index.html">River Conwy eDNA</a>
      <div><a href="index.html">Summary</a><a href="claim-survival.html" aria-current="page">Claim map</a><a href="deliverables.html">Deliverables</a><a href="../START_HERE.md">Reproduce</a></div>
    </nav>
    <div class="catalog-hero">
      <div><p class="eyebrow">Phase 5 · credibility audit</p><h1>What survived?</h1><p class="lede">Five plausible claims enter. Every one leaves narrower: one passes its scoped evidence gate, one survives conditionally, two remain withheld, and one fails its preregistered test.</p></div>
      <aside class="hero-note" aria-label="How to read this page"><strong>Read left to right</strong><p>Initial wording &#8594; audit &#8594; surviving wording. Open any path for its exact denominator, limitation and next evidence.</p><p>Failed gates are results, not missing content.</p></aside>
    </div>
  </header>
  <main id="claim-map">
    <section aria-labelledby="map-heading">
      <p class="eyebrow">Claim-survival map</p>
      <h2 id="map-heading">Evidence changes the sentence</h2>
      <div class="claim-controls" role="group" aria-label="Filter claim paths">
        <button type="button" class="claim-filter" data-filter="ALL" aria-pressed="true">All 5</button>
        <button type="button" class="claim-filter" data-filter="SURVIVES" aria-pressed="false">Survives</button>
        <button type="button" class="claim-filter" data-filter="SURVIVES_CONDITIONAL" aria-pressed="false">Conditional</button>
        <button type="button" class="claim-filter" data-filter="WITHHELD" aria-pressed="false">Withheld</button>
        <button type="button" class="claim-filter" data-filter="FAILED_GATE" aria-pressed="false">Failed gate</button>
      </div>
      <p id="claim-count" class="filter-result" aria-live="polite">Showing all 5 claim paths.</p>
      <div class="claim-map-list">
{articles}
      </div>
      <noscript><p class="callout">Filtering needs JavaScript; all five claim paths and their evidence remain visible without it.</p></noscript>
    </section>
    <section aria-labelledby="register-heading">
      <p class="eyebrow">Portable audit trail</p><h2 id="register-heading">The page is backed by one checked register</h2>
      <p class="section-intro">The register is rebuilt from the five machine-readable gate outputs. It stops the website, report and manuscript from silently using different claim wording.</p>
      <p class="section-actions"><a class="button" href="../data/metadata/claim_survival.tsv">Open claim register</a><a class="button secondary" href="../docs/impact-sprint-preregistration.md">Read frozen specification</a></p>
    </section>
  </main>
  <footer><span>Local static release. No server or network connection is required.</span><span><a href="index.html">Project home</a> · <a href="deliverables.html">All deliverables</a> · <a href="../docs/phase3-claim-status.md">Earlier claim register</a></span></footer>
  <script>
    (() => {{
      const buttons = Array.from(document.querySelectorAll('.claim-filter'));
      const paths = Array.from(document.querySelectorAll('.claim-path'));
      const count = document.getElementById('claim-count');
      buttons.forEach((button) => button.addEventListener('click', () => {{
        const selected = button.dataset.filter;
        buttons.forEach((item) => item.setAttribute('aria-pressed', String(item === button)));
        let visible = 0;
        paths.forEach((path) => {{
          const show = selected === 'ALL' || path.dataset.status === selected;
          path.hidden = !show;
          if (show) visible += 1;
        }});
        count.textContent = selected === 'ALL' ? 'Showing all 5 claim paths.' : `Showing ${{visible}} claim path${{visible === 1 ? '' : 's'}}.`;
      }}));
    }})();
  </script>
</body>
</html>
'''
    OUT_HTML.write_text(document, encoding="utf-8")


def main() -> None:
    records = checked_records()
    write_tsv(records)
    write_html(records)
    print(f"PASS: wrote {OUT_TSV.relative_to(ROOT)} and {OUT_HTML.relative_to(ROOT)} ({len(records)} claims)")


if __name__ == "__main__":
    main()
