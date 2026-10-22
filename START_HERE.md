# Start here - five-minute quick-start capsule

This is the smallest reproducible route through the project. It rebuilds three
website figures and two PDFs from included small local tabular inputs. It does **not** read
FASTQs, download a reference database, start a server, or run Bayesian MCMC.

## Prerequisites

- Python 3.10 or newer
- approximately 100 MB of free space
- the packages pinned in `requirements-website.txt`

One-time setup:

```bash
python3 -m venv .venv-website
source .venv-website/bin/activate
python3 -m pip install -r requirements-website.txt
```

On Windows PowerShell, activate with
`.venv-website\Scripts\Activate.ps1` instead.

## One build command

From the repository root:

```bash
python3 scripts/build_website.py
```

Expected runtime is under five minutes on an ordinary laptop; the verified
reference run in this workspace took under one second.

## Full Phase-5 credibility rebuild

The larger, still offline credibility sprint has one root-level orchestrator:

```bash
python3 scripts/run_impact_sprint.py
```

It checks the frozen local specification, rebuilds and tests the synthetic occupancy
atlas, accession-backed *Salmo* benchmark, Layer-E multiverse, trait null, and
pipeline-disagreement analysis, then regenerates the
[claim-survival map](website/claim-survival.html). It does not download FASTQs,
estimate real River Conwy rare-label occupancy, or turn submitted database labels
into verified species identities.

Run the same command twice after changing an impact analysis. The first run records
SHA-256 values for 47 canonical scientific artifacts; the second run fails unless
every one is byte-identical. Add `--release` to the second run to render the report
and research brief, rebuild the PDFs, create the portable ZIP, execute its verifier
from a clean extraction path containing spaces, and run the full release gate.

The gate record is intentionally mixed. The synthetic-only atlas reports 613/1,944
estimator–scenario rows recoverable, 732 unreliable, and 599 non-identifiable; none
at baseline per-replicate detection `p = 0.20` is recoverable, and both correction
estimators are structurally non-identifiable at one replicate. The three-marker
*Salmo* species gate is coverage-inadequate (2 of 7 candidate labels eligible;
minimum 4 required), while the record-level genus boundary passes 10/10 queries
for every marker. Layer E fails its full
90% validity gate at 64/72 families even though all 64 valid families and all
3,600 one-library draws retain the positive direction and 0/75 leave-outs reverse
it. Fixed-margin trait nulls fail to support broadening beyond richness/prevalence
(`p = 0.3159` primary; `p = 0.1142` strict). The scoped pipeline-disagreement
associations pass for lower source signal (coefficient -1.574, bootstrap 95%
interval -2.950 to -0.807). The count-total direction survives removal of the focal
label: assigned-fish total coefficient -1.149 (-3.086 to -0.579), with 43/43
leave-outs negative. Neither processed total is whole-library sequencing depth;
these are comparator-relative, non-causal associations, not biological false-negative rates.

## Inputs

- `data/derived/layerE_alpha_site_time.tsv` - small site-time summary used to
  regenerate `website/data/community-gradient.csv` over the 11 timepoints
  shared by all 14 sites; the hero reports detected-label richness and retains
  the residual 1-3-library effort range, not true richness;
- `website/data/salmo-calls.csv` - three recorded ASV0001 classification
  outcomes: one same-reference-set method comparison and one independent-reference comparison;
- fixed synthetic parameters embedded in `scripts/build_website.py` for the
  detection illustration.

## Outputs

The command recreates:

1. `website/assets/community-gradient.png`
2. `website/assets/salmo-resolution.png`
3. `website/assets/detection-simulation.png`
4. `website/assets/claim-survival.png`
5. `output/pdf/RESEARCH_BRIEF.pdf` - exactly four pages, including five model architecture diagrams
6. `output/pdf/POSTER.pdf` - exactly one page
7. `website/build-receipt.json`

A successful run prints JSON containing `"status": "PASS"`, records the runtime
and balanced community-panel denominator, checks that every output is non-trivial,
and verifies the PDF page counts.

## Verify the complete local release

After the build, run the offline release audit from the same environment:

```bash
python3 scripts/verify_release.py
```

It checks required website and evidence artifacts, local HTML/Markdown links,
simulator defaults and fixed-seed values, PDF page counts, trait/alignment dimensions,
the Layer-A and Phase-5 claim gates, the claim-survival links, and the current rendered report. Success prints
`"status": "PASS"` with an empty failures list.

The cited trait tables and annotated *Salmo* molecular-detective figure have their
own base-R generator. They are already included; to regenerate or validate them in
the full `edna` R environment:

```bash
Rscript scripts/build_biological_evidence.R
Rscript scripts/build_biological_evidence.R --validate-only
```

These evidence commands are separate from the minimal Python capsule so the one
website build command remains lightweight and independent of the R/MCMC stack.

## Open the local edition

Open `website/index.html` in a browser for the project summary, then use
`website/deliverables.html` for the explained catalogue of every release output.
The canonical five-minute read is the self-contained paper-style page at
`website/research_brief.html`; the simulator is a standalone file at
`website/detection-lab.html`. None of these pages needs a server or network connection.

The research brief is already committed and does not need rebuilding for normal use. To
regenerate it after editing `website/research_brief.md`, use the full `edna` environment:

```bash
micromamba run -n edna pandoc website/research_brief.md \
  --standalone --embed-resources --resource-path=website \
  --css=website/research_brief.css \
  --output=website/research_brief.html
```

This optional render step is separate from the minimal Python capsule; reading and sending
the completed HTML require neither Pandoc nor micromamba.

## Personalize the website text

To change the website's narrative wording without editing HTML, double-click
`EDIT_WEBSITE.command`. Alternatively, run:

```bash
python3 scripts/edit_website.py
```

The command opens a private local editor. Click an outlined heading or paragraph in the
preview, type normally, and choose **Save website text**. Undo, restore-original, and
JSON import/export controls are included. Saved wording lives in a separate override file,
so the canonical page remains recoverable. Evidence numbers, decision labels, links,
figures, model results, and the generated PDFs are intentionally not editable there.
See `EDIT_WEBSITE.md` for the short guide, and press `Ctrl+C` in the launcher window when
you are finished.

## Create the sendable website file

From the repository root, create the portable ZIP with:

```bash
python3 scripts/package_website.py
```

This writes `output/share/conwy-edna-project-website.zip`, its SHA-256 checksum,
and a machine-readable verification report. The archive contains the static website,
every local deliverable directly linked from its catalogue, small presentation inputs,
an obvious `OPEN_WEBSITE.html` entry point, and a standard-library integrity checker.
The packager performs a clean extraction and checks all local HTML links and bundled
checksums before reporting `"status": "PASS"`.

The ZIP is the file to send. A recipient extracts the complete folder and opens
`OPEN_WEBSITE.html`; no Python, R, Quarto, server, or internet connection is needed
to read the local materials or use the synthetic simulator. The external source-study
DOI remains an online link. The bundle is a review/presentation package, not a claim
that the full raw-read or production-MCMC workflow is portable.

## Scope of reproducibility

This capsule reproduces the **presentation layer** from small,
auditable inputs. The full raw-read pipeline ran on a separate 64 GB PC and is
not implied by this command. The technical report records that limitation. The
headline conclusions also do not depend on the unconverged smoke-scale
community-occupancy/recovery fits; production-scale MCMC was not run.
