#!/usr/bin/env python3
"""Build the laptop-sized quick-start capsule (eleven figures + two PDFs).

The script deliberately uses only small included local tabular files. It never reads raw
FASTQs, starts a server, downloads data, or runs Bayesian MCMC.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
import sys
import textwrap
import time
from collections import Counter
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import A4, landscape
    from reportlab.lib.utils import ImageReader
    from reportlab.pdfbase.pdfmetrics import stringWidth
    from reportlab.pdfgen import canvas
except ImportError as exc:  # pragma: no cover - human-facing setup path
    raise SystemExit(
        "Missing website build dependencies. Run: "
        "python3 -m pip install -r requirements-website.txt"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "website" / "data"
ASSETS = ROOT / "website" / "assets"
PDFS = ROOT / "output" / "pdf"

INK = "#17212b"
MUTED = "#52606d"
PAPER = "#f7f6f1"
WHITE = "#ffffff"
RIVER = "#0b7285"
RIVER_DARK = "#075563"
MOSS = "#46734f"
AMBER = "#b56b16"
ROSE = "#a33d46"
LINE = "#d6d9d8"

# Reader-facing decision vocabulary. Canonical machine statuses remain in
# data/metadata/claim_survival.tsv; this translation makes the reader-facing
# brief use one consistent three-state system and always state the reason.
PRESENTATION_DECISIONS = {
    "river-community": {
        "decision": "NON-CONCLUSIVE",
        "reason": "64/72 valid (<90%); direction stable",
        "colour": AMBER,
    },
    "imperfect-detection": {
        "decision": "NON-CONCLUSIVE",
        "reason": "Real fit failed; atlas is synthetic only",
        "colour": AMBER,
    },
    "pipeline-disagreement": {
        "decision": "ACCEPTED",
        "reason": "999/999 finite; 43/43 leave-outs negative",
        "colour": MOSS,
    },
    "taxonomic-resolution": {
        "decision": "NON-CONCLUSIVE",
        "reason": "2 of 7 labels eligible; minimum 4",
        "colour": AMBER,
    },
    "trait-diversity": {
        "decision": "REJECTED",
        "reason": "Both null p-values were > 0.05",
        "colour": ROSE,
    },
}

SYNTHETIC_DEFAULTS = {
    "occupancy": 0.35,
    "detection": 0.25,
    "replicates": 3,
    "relative_depth": 1.0,
    "depth_gradient_coupling": 0.70,
    "seed": 20260718,
    "units": 800,
}


def font(size: int, bold: bool = False):
    candidates = [
        Path("/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf"),
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf" if bold else "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size)
    return ImageFont.load_default(size=size)


def read_csv(path: Path, delimiter=","):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def refresh_community_gradient():
    """Derive the hero table from the balanced 14-site x 11-timepoint panel."""
    source = ROOT / "data" / "derived" / "layerE_alpha_site_time.tsv"
    rows = read_csv(source, delimiter="\t")
    sites = sorted({row["site"] for row in rows})
    timepoint_counts = Counter(row["timepoint"] for row in rows)
    common_timepoints = sorted(
        timepoint for timepoint, count in timepoint_counts.items() if count == len(sites)
    )
    if len(sites) != 14 or len(common_timepoints) != 11:
        raise SystemExit(
            "Community hero gate failed: expected 14 sites x 11 common timepoints, "
            f"found {len(sites)} x {len(common_timepoints)}"
        )

    output_rows = []
    for site in sites:
        panel = [
            row for row in rows
            if row["site"] == site and row["timepoint"] in common_timepoints
        ]
        if len(panel) != len(common_timepoints):
            raise SystemExit(f"Community hero gate failed for {site}: {len(panel)} rows")
        output_rows.append({
            "site": site,
            "section": panel[0]["section"],
            "median_site_time_richness": f"{statistics.median(float(row['richness_q0']) for row in panel):.1f}",
            "n_common_timepoints": len(panel),
            "min_libraries": min(int(row["n_libraries"]) for row in panel),
            "max_libraries": max(int(row["n_libraries"]) for row in panel),
        })

    out = DATA / "community-gradient.csv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_rows[0].keys())
        writer.writeheader()
        writer.writerows(output_rows)
    return {
        "source": str(source.relative_to(ROOT)),
        "sites": len(sites),
        "common_timepoints": common_timepoints,
        "site_time_units": len(sites) * len(common_timepoints),
        "libraries_per_unit_range": [
            min(row["min_libraries"] for row in output_rows),
            max(row["max_libraries"] for row in output_rows),
        ],
        "interpretation": "descriptive detected-label richness; not true richness",
    }


def text(draw, xy, value, size=28, fill=INK, bold=False, anchor=None):
    draw.text(xy, value, font=font(size, bold), fill=fill, anchor=anchor)


def wrapped(draw, xy, value, width_chars, size=26, fill=INK, bold=False, spacing=8):
    value = "\n".join(textwrap.wrap(value, width=width_chars))
    draw.multiline_text(xy, value, font=font(size, bold), fill=fill, spacing=spacing)


def deterministic_scenario(params=None):
    """Python port of website/detection-lab.html's fixed-seed model."""
    p = dict(SYNTHETIC_DEFAULTS if params is None else params)
    state = int(p["seed"]) & 0xFFFFFFFF

    def random():
        nonlocal state
        state = (1664525 * state + 1013904223) & 0xFFFFFFFF
        return state / 4294967296

    def logistic(x):
        return 1 / (1 + math.exp(-x))

    def logit(x):
        return math.log(x / (1 - x))

    def clamp(x, lo, hi):
        return max(lo, min(hi, x))

    units = []
    for _ in range(int(p["units"])):
        gradient = 2 * random() - 1
        noise = (random() + random() + random() - 1.5) * 0.45
        psi = logistic(logit(p["occupancy"]) + 1.10 * gradient)
        depth = clamp(
            p["relative_depth"] * math.exp(p["depth_gradient_coupling"] * gradient + noise),
            0.08,
            5,
        )
        occupied = random() < psi
        p_true = logistic(logit(p["detection"]) + 0.62 * math.log(depth) + 0.85 * gradient)
        hits = 0
        if occupied:
            for _ in range(int(p["replicates"])):
                if random() < p_true:
                    hits += 1
        units.append((psi, depth, hits > 0))

    true_prob = sum(row[0] for row in units) / len(units)
    naive = sum(row[2] for row in units) / len(units)
    best_psi, best_ll = 0.0, -math.inf
    for n in range(1, 999):
        psi = n / 1000
        ll = 0.0
        for _, depth, seen in units:
            p_model = logistic(logit(p["detection"]) + 0.62 * math.log(depth))
            q = clamp(psi * (1 - (1 - p_model) ** int(p["replicates"])), 1e-8, 1 - 1e-8)
            ll += math.log(q if seen else 1 - q)
        if ll > best_ll:
            best_psi, best_ll = psi, ll

    return {
        "parameters": p,
        "known_mean_occupancy": true_prob,
        "naive_detected_fraction": naive,
        "depth_only_corrected": best_psi,
        "correction_error": best_psi - true_prob,
        "all_replicates_miss_risk": (1 - p["detection"]) ** int(p["replicates"]),
    }


def community_figure() -> Path:
    rows = read_csv(DATA / "community-gradient.csv")
    image = Image.new("RGB", (1400, 760), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (70, 42), "Detected fish labels broaden downstream", 48, INK, True)
    text(d, (70, 104), "Median detected-label richness across 11 common timepoints; 1-3 libraries per unit", 25, MUTED)

    left, top, right, bottom = 105, 180, 1330, 640
    max_y = 8
    section_ranges = [(0, 4, "upper", "#e3eff1"), (5, 9, "middle", "#d8e8e2"), (10, 13, "lower", "#f3e6d4")]
    step = (right - left) / (len(rows) - 1)
    for lo, hi, label, colour in section_ranges:
        x0 = max(left - 25, left + lo * step - step * 0.45)
        x1 = min(right + 25, left + hi * step + step * 0.45)
        d.rectangle((x0, top, x1, bottom), fill=colour)
        text(d, ((x0 + x1) / 2, top + 20), label.upper(), 19, MUTED, True, "mm")

    for y in range(0, max_y + 1, 2):
        py = bottom - y / max_y * (bottom - top)
        d.line((left, py, right, py), fill=LINE, width=2)
        text(d, (left - 18, py), str(y), 19, MUTED, anchor="rm")
    d.line((left, top, left, bottom), fill=INK, width=3)
    d.line((left, bottom, right, bottom), fill=INK, width=3)

    points = []
    for i, row in enumerate(rows):
        x = left + i * step
        value = float(row["median_site_time_richness"])
        y = bottom - value / max_y * (bottom - top)
        points.append((x, y))
        d.ellipse((x - 9, y - 9, x + 9, y + 9), fill=RIVER_DARK, outline=WHITE, width=2)
        text(d, (x, bottom + 32), row["site"], 17, INK, True, "mm")
    d.line(points, fill=RIVER_DARK, width=5, joint="curve")
    text(d, (700, 710), "Observed labels, not true richness; residual 1-3-library effort imbalance remains.", 20, MUTED, anchor="mm")
    out = ASSETS / "community-gradient.png"
    image.save(out, quality=95)
    return out


def salmo_figure() -> Path:
    rows = read_csv(DATA / "salmo-calls.csv")
    image = Image.new("RGB", (1400, 760), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (70, 42), "One 170 bp sequence, three species-level labels", 48, INK, True)
    text(d, (70, 104), "ASV0001 receives a 100% identity hit in each comparison - but not the same species label", 25, MUTED)

    cx, cy = 700, 345
    d.ellipse((cx - 126, cy - 72, cx + 126, cy + 72), fill=RIVER_DARK)
    text(d, (cx, cy - 11), "ASV0001", 30, WHITE, True, "mm")
    text(d, (cx, cy + 27), "170 bp MiFish 12S", 20, WHITE, anchor="mm")

    positions = [(250, 255), (1150, 255), (700, 620)]
    card_colours = ["#eaf3f4", "#edf2e9", "#fff0dc"]
    for row, (x, y), colour in zip(rows, positions, card_colours):
        if y < cy:
            d.line((cx - 90 if x < cx else cx + 90, cy - 42, x + 120 if x < cx else x - 120, y + 75), fill=RIVER, width=4)
        else:
            d.line((cx, cy + 72, x, y - 75), fill=RIVER, width=4)
        d.rounded_rectangle((x - 205, y - 86, x + 205, y + 86), radius=16, fill=colour, outline=LINE, width=2)
        text(d, (x, y - 52), row["method_or_database"], 19, MUTED, True, "mm")
        text(d, (x, y - 4), row["species_call"], 30, INK, True, "mm")
        text(d, (x, y + 40), f'{row["identity_percent"]}% identity; genus {row["genus_call"]}', 18, MUTED, anchor="mm")

    text(d, (70, 704), "Defensible inference: genus Salmo is stable; species identity is not identifiable from this marker/reference set.", 24, RIVER_DARK, True)
    out = ASSETS / "salmo-resolution.png"
    image.save(out, quality=95)
    return out


def detection_figure() -> Path:
    scenario = deterministic_scenario()
    image = Image.new("RGB", (1400, 760), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (70, 42), "Replication helps; incomplete adjustment can mislead", 48, INK, True)
    text(d, (70, 104), "Fixed-seed synthetic teaching scenario - independent of the Layer A MCMC", 25, MUTED)

    # Panel A: all-miss probability at p=0.25, matching the standalone lab default.
    left, top, right, bottom = 95, 205, 690, 625
    text(d, (left, 155), "All-miss chance at reference conditions", 25, INK, True)
    for yv in [0, .25, .5, .75, 1]:
        py = bottom - yv * (bottom - top)
        d.line((left, py, right, py), fill=LINE, width=2)
        text(d, (left - 15, py), f"{int(yv * 100)}%", 17, MUTED, anchor="rm")
    points = []
    for k in range(1, 9):
        miss = (1 - 0.25) ** k
        x = left + (k - 1) / 7 * (right - left)
        y = bottom - miss * (bottom - top)
        points.append((x, y))
        d.ellipse((x - 8, y - 8, x + 8, y + 8), fill=RIVER_DARK)
        text(d, (x, bottom + 29), str(k), 17, MUTED, anchor="mm")
    d.line(points, fill=RIVER_DARK, width=5, joint="curve")
    text(d, ((left + right) / 2, bottom + 62), "water-sample replicates", 19, MUTED, anchor="mm")

    # Panel B: deterministic standalone-lab default. These are known-truth teaching
    # values, not estimates from the unconverged Layer A recovery MCMC.
    x0, y0 = 820, 210
    text(d, (x0, 155), "Known truth versus observation and adjustment", 25, INK, True)
    labels = [
        ("Known truth", scenario["known_mean_occupancy"], MOSS),
        ("Naive detected", scenario["naive_detected_fraction"], "#7b8792"),
        ("Depth-only", scenario["depth_only_corrected"], AMBER),
    ]
    max_value = 0.50
    for i, (label, value, colour) in enumerate(labels):
        y = y0 + i * 100
        d.rounded_rectangle((x0, y, x0 + 460, y + 68), radius=10, fill=WHITE, outline=LINE, width=2)
        d.rectangle((x0, y, x0 + 460 * value / max_value, y + 68), fill=colour)
        text(d, (x0 + 18, y + 15), label, 21, INK, True)
        text(d, (x0 + 440, y + 34), f"{100 * value:.1f}%", 26, INK, True, "rm")
    wrapped(d, (x0, 525), f'At this selected scenario the depth-only calculation overshoots known occupancy by {100 * scenario["correction_error"]:.1f} percentage points.', 44, 22, MUTED)
    text(d, (70, 710), "Teaching counterexample only: the real rare-species MCMC direction is withheld after failed convergence.", 22, RIVER_DARK, True)
    out = ASSETS / "detection-simulation.png"
    image.save(out, quality=95)
    return out


def claim_survival_figure() -> Path:
    """Render the five checked Phase-5 outcomes as a compact offline summary."""
    rows = read_csv(ROOT / "data" / "metadata" / "claim_survival.tsv", delimiter="\t")
    expected = [
        "river-community", "imperfect-detection", "pipeline-disagreement",
        "taxonomic-resolution", "trait-diversity",
    ]
    if [row["claim_id"] for row in rows] != expected:
        raise SystemExit("Claim-survival figure gate failed: expected five ordered claims")

    concise = {
        "river-community": "River community",
        "imperfect-detection": "Imperfect detection",
        "pipeline-disagreement": "Pipeline disagreement",
        "taxonomic-resolution": "Taxonomic resolution",
        "trait-diversity": "Functional traits",
    }
    image = Image.new("RGB", (1400, 650), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (65, 38), "Five claims, one decision vocabulary", 45, INK, True)
    text(d, (65, 88), "ACCEPTED = adequate evidence met the decision rule | REJECTED = an adequate test did not support the claim", 18, MUTED)
    text(d, (65, 115), "NON-CONCLUSIVE = fit, validity, or coverage was insufficient to decide", 18, MUTED)
    y = 150
    for row in rows:
        claim = concise[row["claim_id"]]
        presentation = PRESENTATION_DECISIONS[row["claim_id"]]
        result = presentation["reason"]
        colour = presentation["colour"]
        d.rounded_rectangle((65, y, 1335, y + 82), radius=9, fill=WHITE, outline=LINE, width=2)
        d.rectangle((65, y, 81, y + 82), fill=colour)
        text(d, (105, y + 16), claim, 26, INK, True)
        text(d, (420, y + 18), result, 23, INK)
        status = presentation["decision"]
        text(d, (1310, y + 58), status, 18, colour, True, "rm")
        y += 94
    text(d, (65, 625), "Decision tally: 1 accepted | 1 rejected | 3 non-conclusive. Detailed decision records remain linked offline.", 19, MUTED)
    out = ASSETS / "claim-survival.png"
    image.save(out, quality=95)
    return out


def research_timeline_figure() -> Path:
    """Render the research decision trail as a minimalist documentation timeline."""
    image = Image.new("RGB", (1400, 600), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (65, 36), "How the evidence was checked and documented", 44, INK, True)
    text(d, (65, 90), "Five explained steps moved the source data toward defensible claim decisions", 24, MUTED)

    nodes = [
        (
            "CHECK SOURCE DATA",
            "Matched Conwy records,\nlinked repeated samples,\nand marked missing data\nbefore analysis.",
        ),
        (
            "REPROCESS DNA READS",
            "Re-ran raw sequencing;\n647 field libraries and\n1,348 DNA variants passed\nvalidation.",
        ),
        (
            "RECORD MODEL LIMITS",
            "The rare-label model did\nnot converge; no biological\noccupancy estimate was\nreported.",
        ),
        (
            "SET TESTS IN ADVANCE",
            "Before stress testing, wrote\nfive questions, exact settings,\nwhich runs would count, and\nverdict rules.",
        ),
        (
            "CONFIRM REPRODUCIBILITY",
            "Ran the full analysis twice;\nthen checked a clean offline\ncopy for matching, working\noutputs.",
        ),
    ]
    xs = [150, 425, 700, 975, 1250]
    line_y = 200
    d.line((xs[0], line_y, xs[-1], line_y), fill=LINE, width=8)
    for index, ((heading, body), x) in enumerate(zip(nodes, xs), start=1):
        if index < len(nodes):
            d.polygon(
                [(x + 103, line_y - 10), (x + 123, line_y), (x + 103, line_y + 10)],
                fill=RIVER_DARK,
            )
        d.ellipse((x - 30, line_y - 30, x + 30, line_y + 30), fill=RIVER_DARK)
        text(d, (x, line_y), str(index), 26, WHITE, True, "mm")
        text(d, (x, 252), heading, 18, RIVER_DARK, True, "mm")
        body_lines = body.split("\n")
        for offset, line in enumerate(body_lines):
            text(d, (x, 290 + offset * 25), line, 18, INK, anchor="mm")

    d.rounded_rectangle((65, 430, 1335, 570), radius=14, fill="#eaf3f4")
    text(d, (95, 451), "WHY THIS MATTERS", 19, RIVER_DARK, True)
    text(
        d,
        (95, 488),
        "The five questions and verdict rules were written before the new stress tests, so favourable results could not define the plan.",
        18,
        INK,
    )
    text(
        d,
        (95, 520),
        "Unsuccessful models stayed visible, and running everything twice checked that the same inputs produced the same outputs.",
        18,
        MUTED,
    )
    out = ASSETS / "research-timeline.png"
    image.save(out, quality=95)
    return out


def models_overview_figure() -> Path:
    """Render the five major stress-test engines and their bounded outcomes."""
    image = Image.new("RGB", (1400, 670), PAPER)
    d = ImageDraw.Draw(image)
    text(d, (65, 34), "Five claims, five stress-test engines", 44, INK, True)
    text(d, (65, 88), "Every engine ends with the same decision label and a stated reason", 24, MUTED)

    rows = [
        ("1", "river-community", "River community", "Restricted PERMANOVA + 72-family multiverse"),
        ("2", "imperfect-detection", "Imperfect detection", "648 scenarios x 100 simulations x 3 estimators"),
        ("3", "pipeline-disagreement", "Pipeline disagreement", "Logistic model + 999 block bootstraps + 43 leave-outs"),
        ("4", "taxonomic-resolution", "Taxonomic resolution", "MiFish / COI / cyt-b distances + all-ties LOO"),
        ("5", "trait-diversity", "Functional traits", "Gower + Rao's Q + 9,999 fixed-margin nulls"),
    ]
    y = 132
    for number, claim_id, claim, method in rows:
        presentation = PRESENTATION_DECISIONS[claim_id]
        status = presentation["decision"]
        colour = presentation["colour"]
        d.rounded_rectangle((65, y, 1335, y + 86), radius=10, fill=WHITE, outline=LINE, width=2)
        d.ellipse((88, y + 17, 140, y + 69), fill=RIVER_DARK)
        text(d, (114, y + 43), number, 23, WHITE, True, "mm")
        text(d, (170, y + 22), claim, 24, INK, True)
        text(d, (465, y + 15), method, 18, INK)
        text(d, (465, y + 50), f"Reason: {presentation['reason']}", 16, MUTED)
        d.rounded_rectangle((1110, y + 20, 1310, y + 66), radius=20, fill=colour)
        text(d, (1210, y + 43), status, 16, WHITE, True, "mm")
        y += 96

    text(
        d,
        (65, 640),
        "Decision tally: 1 accepted | 1 rejected | 3 non-conclusive.",
        20,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "models-overview.png"
    image.save(out, quality=95)
    return out


def diagram_arrow(draw, start, end, colour=RIVER_DARK, width=5, head=14):
    """Draw a clean directional connector for the model architecture figures."""
    x0, y0 = start
    x1, y1 = end
    dx, dy = x1 - x0, y1 - y0
    length = max(math.hypot(dx, dy), 1)
    ux, uy = dx / length, dy / length
    px, py = -uy, ux
    base_x, base_y = x1 - head * ux, y1 - head * uy
    draw.line((x0, y0, base_x, base_y), fill=colour, width=width)
    draw.polygon(
        [
            (x1, y1),
            (base_x + head * 0.55 * px, base_y + head * 0.55 * py),
            (base_x - head * 0.55 * px, base_y - head * 0.55 * py),
        ],
        fill=colour,
    )


def diagram_box(draw, box, heading, lines, accent=RIVER, fill=WHITE):
    """Draw one labelled analytical stage with a shared visual grammar."""
    x0, y0, x1, y1 = box
    draw.rounded_rectangle(box, radius=12, fill=fill, outline=LINE, width=2)
    draw.rectangle((x0, y0, x0 + 9, y1), fill=accent)
    text(draw, (x0 + 24, y0 + 16), heading, 19, INK, True)
    for index, line in enumerate(lines):
        text(draw, (x0 + 24, y0 + 55 + index * 27), line, 17, MUTED)


def diagram_status(draw, label, colour, x1=1360):
    width = max(170, 14 * len(label))
    x0 = x1 - width
    draw.rounded_rectangle((x0, 34, x1, 74), radius=19, fill=colour)
    text(draw, ((x0 + x1) / 2, 54), label, 15, WHITE, True, "mm")


def model_canvas(number, title_value, subtitle, status, status_colour):
    image = Image.new("RGB", (1400, 430), PAPER)
    draw = ImageDraw.Draw(image)
    text(draw, (38, 23), f"MODEL {number:02d}", 17, RIVER_DARK, True)
    text(draw, (38, 49), title_value, 30, INK, True)
    text(draw, (38, 88), subtitle, 18, MUTED)
    diagram_status(draw, status, status_colour)
    draw.line((38, 112, 1362, 112), fill=LINE, width=2)
    return image, draw


def model_community_multiverse_figure() -> Path:
    """Show how the community claim branches into 72 declared specifications."""
    image, d = model_canvas(
        1,
        "Community robustness multiverse",
        "INPUT -> MODEL FAMILY -> RESTRICTED STRESS TEST -> DECISION",
        "NON-CONCLUSIVE",
        AMBER,
    )
    diagram_box(
        d,
        (35, 138, 305, 342),
        "SPECIFICATION SPACE",
        [
            "3 panels x 2 taxonomies",
            "x 2 thresholds x 2 effort rules",
            "x 3 distances = 72 families",
        ],
        RIVER,
    )
    diagram_arrow(d, (310, 240), (355, 240))

    text(d, (365, 137), "COMPLETE FAMILY GRID", 18, INK, True)
    tile_w, tile_h, gap = 18, 18, 4
    for index in range(72):
        row, column = divmod(index, 12)
        x = 366 + column * (tile_w + gap)
        y = 172 + row * (tile_h + gap)
        colour = RIVER if index < 64 else ROSE
        d.rounded_rectangle((x, y, x + tile_w, y + tile_h), radius=3, fill=colour)
    text(d, (365, 315), "64 valid", 17, RIVER_DARK, True)
    text(d, (455, 315), "8 undefined rCLR", 17, ROSE, True)
    text(d, (365, 340), "Invalid families remain in the reported total", 15, MUTED)

    diagram_arrow(d, (645, 240), (700, 240))
    diagram_box(
        d,
        (705, 138, 1020, 342),
        "RESTRICTED INFERENCE",
        [
            "Site x timepoint units stay clustered",
            "28 river arrangements: p = 2/28",
            "3,600 draws + 75 leave-outs",
        ],
        RIVER_DARK,
    )
    diagram_arrow(d, (1025, 240), (1075, 240))
    diagram_box(
        d,
        (1080, 138, 1365, 342),
        "DECISION RULE",
        [
            "64/72 valid = 88.9% (<90%)",
            "64/64 positive; 3,600/3,600",
            "0/75 sign reversals",
        ],
        AMBER,
        "#fff7e7",
    )
    text(
        d,
        (38, 397),
        "NON-CONCLUSIVE: 64/72 valid missed the 90% rule, although every valid direction was positive.",
        16,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "model-community-multiverse.png"
    image.save(out, quality=95)
    return out


def model_occupancy_atlas_figure() -> Path:
    """Show the known-truth generator, estimators and recoverability phase gate."""
    image, d = model_canvas(
        2,
        "Known-truth occupancy recoverability atlas",
        "SYNTHETIC GENERATOR -> DETECTION HISTORIES -> THREE ESTIMATORS -> DECISION MAP",
        "NON-CONCLUSIVE",
        AMBER,
    )
    diagram_box(
        d,
        (35, 138, 310, 342),
        "648 SCENARIOS",
        [
            "psi, prevalence and baseline p",
            "1-3 replicates + depth structure",
            "missingness + missing-data coding",
        ],
        RIVER,
    )
    diagram_arrow(d, (315, 240), (365, 240))
    diagram_box(
        d,
        (370, 138, 650, 342),
        "KNOWN TRUTH",
        [
            "120 ordered site x time units",
            "latent occupancy z_i",
            "replicate detections + missing-data rules",
        ],
        RIVER_DARK,
    )
    # Small latent-state diagram.
    d.ellipse((545, 270, 580, 305), fill=MOSS)
    d.ellipse((590, 278, 615, 303), fill=RIVER)
    d.ellipse((620, 278, 645, 303), fill=RIVER)
    d.ellipse((605, 310, 630, 335), fill=RIVER)
    diagram_arrow(d, (562, 306), (601, 290), MOSS, 3, 9)
    diagram_arrow(d, (562, 306), (630, 290), MOSS, 3, 9)
    diagram_arrow(d, (562, 306), (617, 318), MOSS, 3, 9)

    diagram_arrow(d, (655, 240), (700, 240))
    d.rounded_rectangle((705, 138, 1010, 342), radius=12, fill=WHITE, outline=LINE, width=2)
    d.rectangle((705, 138, 714, 342), fill=RIVER)
    text(d, (729, 154), "ESTIMATOR LANES", 19, INK, True)
    estimators = [
        ("NAIVE OCCURRENCE", "#e8ecee"),
        ("REPLICATION-AWARE", "#dfeff1"),
        ("DEPTH-AWARE", "#dfe9e1"),
    ]
    for index, (label, fill) in enumerate(estimators):
        y = 195 + index * 46
        d.rounded_rectangle((730, y, 980, y + 33), radius=7, fill=fill, outline=LINE, width=1)
        text(d, (855, y + 17), label, 15, INK, True, "mm")
    diagram_arrow(d, (1015, 240), (1060, 240))

    d.rounded_rectangle((1065, 138, 1365, 342), radius=12, fill="#fffafa", outline=LINE, width=2)
    d.rectangle((1065, 138, 1074, 342), fill=ROSE)
    text(d, (1088, 154), "1,944 SCENARIO CELLS", 18, INK, True)
    total = 1944
    values = [(613, MOSS), (732, AMBER), (599, ROSE)]
    bar_left, bar_top, bar_width = 1088, 202, 250
    cursor = bar_left
    for value, colour in values:
        width = bar_width * value / total
        d.rectangle((cursor, bar_top, cursor + width, bar_top + 30), fill=colour)
        cursor += width
    text(d, (1088, 249), "613 recoverable", 15, MOSS, True)
    text(d, (1088, 275), "732 unreliable", 15, AMBER, True)
    text(d, (1088, 301), "599 non-identifiable", 15, ROSE, True)
    text(d, (1088, 321), "No baseline p = 0.20 cell passed", 13, MUTED)
    text(
        d,
        (38, 397),
        "NON-CONCLUSIVE: the real fit failed; the 64,800 known-truth datasets provide synthetic design evidence only.",
        16,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "model-occupancy-atlas.png"
    image.save(out, quality=95)
    return out


def model_pipeline_bootstrap_figure() -> Path:
    """Show the comparator model, cluster bootstrap and leave-out gate."""
    image, d = model_canvas(
        3,
        "Cluster-aware pipeline-disagreement model",
        "VERIFIED CELLS -> LOGISTIC ASSOCIATION -> PAIRED CLUSTER BOOTSTRAP -> DECISION",
        "ACCEPTED",
        MOSS,
    )
    diagram_box(
        d,
        (35, 138, 275, 342),
        "COMPARATOR CELLS",
        [
            "2,344 source-positive",
            "1,962 recovered",
            "382 non-recovered",
        ],
        RIVER,
    )
    diagram_arrow(d, (280, 240), (330, 240))
    diagram_box(
        d,
        (335, 138, 625, 342),
        "LOGISTIC MODEL",
        [
            "logit P(non-recovery)",
            "source signal + count context",
            "section + season",
        ],
        RIVER_DARK,
    )
    text(d, (360, 308), "Sensitivity excludes focal label + Unassigned", 14, MUTED)
    diagram_arrow(d, (630, 240), (680, 240))

    d.rounded_rectangle((685, 138, 1015, 342), radius=12, fill=WHITE, outline=LINE, width=2)
    d.rectangle((685, 138, 694, 342), fill=RIVER)
    text(d, (709, 154), "999 PAIRED BOOTSTRAPS", 19, INK, True)
    lanes = [("SITES / SECTION", RIVER), ("TIME / BLOCK LENGTH 3", AMBER), ("TAXA / WHOLE CLUSTERS", MOSS)]
    for index, (label, colour) in enumerate(lanes):
        y = 205 + index * 44
        d.line((720, y, 950, y), fill=colour, width=7)
        for offset in range(5):
            x = 735 + offset * 50
            d.ellipse((x - 7, y - 7, x + 7, y + 7), fill=colour)
        text(d, (720, y - 24), label, 14, MUTED, True)
    text(d, (709, 326), "43 leave-outs: 11 taxa + 14 sites + 18 times", 14, MUTED)
    diagram_arrow(d, (1020, 240), (1070, 240))

    d.rounded_rectangle((1075, 138, 1365, 342), radius=12, fill="#edf5ee", outline=LINE, width=2)
    d.rectangle((1075, 138, 1084, 342), fill=MOSS)
    text(d, (1099, 154), "FOCAL-EXCLUDED RESULT", 18, INK, True)
    # Mini forest interval.
    axis_left, axis_right, axis_y = 1115, 1335, 236
    d.line((axis_left, axis_y, axis_right, axis_y), fill=MUTED, width=2)
    zero_x = axis_left + (0 - (-3.5)) / 4.0 * (axis_right - axis_left)
    d.line((zero_x, axis_y - 34, zero_x, axis_y + 34), fill=ROSE, width=2)
    lo_x = axis_left + (-3.086 + 3.5) / 4.0 * (axis_right - axis_left)
    hi_x = axis_left + (-0.579 + 3.5) / 4.0 * (axis_right - axis_left)
    beta_x = axis_left + (-1.149 + 3.5) / 4.0 * (axis_right - axis_left)
    d.line((lo_x, axis_y, hi_x, axis_y), fill=MOSS, width=8)
    d.ellipse((beta_x - 9, axis_y - 9, beta_x + 9, axis_y + 9), fill=MOSS)
    text(d, (1099, 276), "beta = -1.149", 17, INK, True)
    text(d, (1099, 303), "95%: -3.086 to -0.579", 15, MUTED)
    text(d, (1099, 328), "999/999 finite; 43/43 negative", 14, MOSS, True)
    text(
        d,
        (38, 397),
        "ACCEPTED: 999/999 finite and 43/43 negative; comparator-relative and noncausal, not a biological false-negative rate.",
        15,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "model-pipeline-bootstrap.png"
    image.save(out, quality=95)
    return out


def model_salmo_benchmark_figure() -> Path:
    """Show the sequence-marker rulers and conservative submitted-label gates."""
    image, d = model_canvas(
        4,
        "Multi-marker Salmo resolution benchmark",
        "ACCESSIONS -> MARKER WINDOWS -> DISTANCE / BARCODE / ALL-TIES LOO -> DECISION",
        "NON-CONCLUSIVE",
        AMBER,
    )
    diagram_box(
        d,
        (35, 138, 255, 342),
        "PUBLIC PANEL",
        [
            "21 mitogenome records",
            "7 candidate Salmo labels",
            "2 outgroup genera",
        ],
        RIVER,
    )
    diagram_arrow(d, (260, 240), (305, 240))

    d.rounded_rectangle((310, 138, 650, 342), radius=12, fill=WHITE, outline=LINE, width=2)
    d.rectangle((310, 138, 319, 342), fill=RIVER)
    text(d, (334, 154), "MARKER RULERS", 19, INK, True)
    marker_rows = [
        ("MiFish 12S", 170, RIVER),
        ("COI", 648, MOSS),
        ("cyt-b", 1140, AMBER),
    ]
    for index, (label, length_bp, colour) in enumerate(marker_rows):
        y = 205 + index * 49
        text(d, (334, y - 18), f"{label}  {length_bp:,} bp", 15, MUTED, True)
        width = 255 * length_bp / 1140
        d.rounded_rectangle((334, y, 334 + width, y + 17), radius=7, fill=colour)
    diagram_arrow(d, (655, 240), (700, 240))

    diagram_box(
        d,
        (705, 138, 1035, 342),
        "CONSERVATIVE AUDIT",
        [
            "normalised edit distance",
            "exact-haplotype sharing + gaps",
            "all tied nearest neighbours",
        ],
        RIVER_DARK,
    )
    diagram_arrow(d, (1040, 240), (1080, 240))
    d.rounded_rectangle((1085, 138, 1365, 342), radius=12, fill="#fffafa", outline=LINE, width=2)
    d.rectangle((1085, 138, 1094, 342), fill=ROSE)
    text(d, (1109, 154), "TWO DECISION CHECKS", 18, INK, True)
    d.rounded_rectangle((1110, 201, 1340, 252), radius=8, fill="#faecee")
    text(d, (1125, 211), "SPECIES: 2 of 7 eligible", 16, ROSE, True)
    text(d, (1125, 233), "minimum 4 -> NON-CONCLUSIVE", 13, MUTED)
    d.rounded_rectangle((1110, 270, 1340, 321), radius=8, fill="#edf5ee")
    text(d, (1125, 280), "GENUS: 10/10 x 3", 16, MOSS, True)
    text(d, (1125, 302), "Salmo nearer than each outgroup", 14, MUTED)
    text(
        d,
        (38, 397),
        "NON-CONCLUSIVE: only 2 of 7 labels were eligible; submitted labels remain annotations, not local species identity.",
        15,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "model-salmo-benchmark.png"
    image.save(out, quality=95)
    return out


def model_trait_null_figure() -> Path:
    """Show the Gower-to-Rao pipeline and fixed-margin null randomisation."""
    image, d = model_canvas(
        5,
        "Fixed-margin functional-trait null model",
        "OCCURRENCE + TRAITS -> GOWER -> RAO Q -> MARGIN-PRESERVING NULL -> DECISION",
        "REJECTED",
        ROSE,
    )
    # Incidence matrix.
    d.rounded_rectangle((35, 138, 270, 342), radius=12, fill=WHITE, outline=LINE, width=2)
    d.rectangle((35, 138, 44, 342), fill=RIVER)
    text(d, (59, 154), "154 UNIT MATRICES", 19, INK, True)
    for row in range(6):
        for column in range(9):
            x = 65 + column * 19
            y = 207 + row * 19
            fill = RIVER if (row * 3 + column * 5) % 7 < 3 else "#e4e7e8"
            d.rectangle((x, y, x + 14, y + 14), fill=fill)
    text(d, (59, 326), "11 timepoint blocks", 14, MUTED)
    diagram_arrow(d, (275, 240), (320, 240))
    diagram_box(
        d,
        (325, 138, 530, 342),
        "GOWER",
        [
            "equal-weight axes",
            "flow + trophic",
            "+ movement traits",
        ],
        RIVER_DARK,
    )
    diagram_arrow(d, (535, 240), (580, 240))
    diagram_box(
        d,
        (585, 138, 770, 342),
        "RAO Q",
        [
            "unit diversity",
            "site means",
            "lower - upper",
        ],
        RIVER,
    )
    diagram_arrow(d, (775, 240), (820, 240))
    diagram_box(
        d,
        (825, 138, 1065, 342),
        "FIXED MARGINS",
        [
            "2 x 2 checkerboard swaps",
            "richness + prevalence fixed",
            "9,999 draws x 2 panels",
        ],
        AMBER,
        "#fff7e7",
    )
    diagram_arrow(d, (1070, 240), (1110, 240))

    d.rounded_rectangle((1115, 138, 1365, 342), radius=12, fill="#fffafa", outline=LINE, width=2)
    d.rectangle((1115, 138, 1124, 342), fill=ROSE)
    text(d, (1139, 154), "EMPIRICAL NULLS", 18, INK, True)

    # Plot the actual frozen fixed-margin draws rather than a decorative density.
    draw_rows = read_csv(ROOT / "data" / "derived" / "impact_trait_null_draws.tsv", delimiter="\t")
    summary_rows = {
        row["analysis"]: row
        for row in read_csv(ROOT / "data" / "derived" / "impact_trait_null_summary.tsv", delimiter="\t")
    }
    panel_specs = [
        ("primary_12_labels", "PRIMARY"),
        ("strict_5_labels", "STRICT"),
    ]
    for index, (analysis, label) in enumerate(panel_specs):
        values = [
            float(row["contrast_lower_minus_upper"])
            for row in draw_rows
            if row["analysis"] == analysis
        ]
        summary = summary_rows[analysis]
        observed = float(summary["observed_contrast"])
        p_value = summary["monte_carlo_p"]
        low, high = min(values), max(values)
        bin_count = 30
        counts = [0] * bin_count
        for value in values:
            slot = min(bin_count - 1, int((value - low) / (high - low) * bin_count))
            counts[slot] += 1
        base_y = 245 + index * 75
        left, right = 1145, 1335
        peak = max(counts)
        points = [
            (
                left + (bin_index + 0.5) / bin_count * (right - left),
                base_y - count / peak * 40,
            )
            for bin_index, count in enumerate(counts)
        ]
        d.polygon([(left, base_y), *points, (right, base_y)], fill="#dfeff1")
        d.line(points, fill=RIVER_DARK, width=4)
        observed_x = left + (observed - low) / (high - low) * (right - left)
        d.line((observed_x, base_y - 44, observed_x, base_y + 3), fill=ROSE, width=3)
        text(d, (1140, base_y - 60), f"{label}  p = {p_value}", 13, INK, True)
    text(
        d,
        (38, 397),
        "REJECTED: both p-values exceeded 0.05; no evidence of trait broadening beyond richness and prevalence.",
        16,
        RIVER_DARK,
        True,
    )
    out = ASSETS / "model-trait-null.png"
    image.save(out, quality=95)
    return out


def pdf_header(c, title, subtitle, page_width, page_height):
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.rect(0, page_height - 88, page_width, 88, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 21)
    c.drawString(36, page_height - 45, title)
    c.setFont("Helvetica", 9.5)
    c.drawString(36, page_height - 65, subtitle)


def pdf_section(c, number, title, y, page_width):
    """Draw a compact numbered section label and return the next text baseline."""
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.circle(45, y + 1, 9, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 8.5)
    c.drawCentredString(45, y - 2, str(number))
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.setFont("Helvetica-Bold", 12)
    c.drawString(62, y - 4, title)
    c.setStrokeColor(colors.HexColor(LINE))
    c.setLineWidth(0.8)
    line_start = 62 + stringWidth(title, "Helvetica-Bold", 12) + 18
    if line_start < page_width - 44:
        c.line(line_start, y - 1, page_width - 36, y - 1)
    return y - 24


def pdf_paragraph(c, value, x, y, width, size=9.2, leading=12, colour=INK, bold=False):
    c.setFillColor(colors.HexColor(colour))
    c.setFont("Helvetica-Bold" if bold else "Helvetica", size)
    words = value.split()
    lines, current = [], ""
    for word in words:
        trial = f"{current} {word}".strip()
        if stringWidth(trial, "Helvetica-Bold" if bold else "Helvetica", size) <= width:
            current = trial
        else:
            lines.append(current)
            current = word
    if current:
        lines.append(current)
    for line in lines:
        c.drawString(x, y, line)
        y -= leading
    return y


def _research_brief_pdf_two_page_legacy(figures) -> Path:
    out = PDFS / "RESEARCH_BRIEF.pdf"
    pw, ph = A4
    c = canvas.Canvas(str(out), pagesize=A4, invariant=1)
    c.setTitle("River Conwy eDNA - investigation brief")

    pdf_header(
        c,
        "From fish-eDNA data to a defensible claim",
        "Six-part claim audit | legacy compact River Conwy portfolio edition",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 1, "Our question", y, pw)
    y = pdf_paragraph(
        c,
        "What can a replicated river-eDNA dataset establish after we correct the unit of inference, test missed-detection recovery, audit pipeline disagreement and challenge taxonomic and trait interpretations?",
        36,
        y,
        pw - 72,
        9.4,
        12.2,
        INK,
        True,
    )
    y -= 5
    c.setFillColor(colors.HexColor("#eaf3f4"))
    c.rect(36, y - 27, pw - 72, 31, fill=1, stroke=0)
    pdf_paragraph(
        c,
        "14 sites | 18 usable timepoints | 647 validated field libraries | 1,348 ASVs | replicated 12S survey",
        47,
        y - 8,
        pw - 94,
        8.3,
        10,
        RIVER_DARK,
        True,
    )
    y -= 48

    y = pdf_section(c, 2, "Our methodology + key documentation moments", y, pw)
    c.drawImage(
        ImageReader(figures[4]),
        36,
        y - 150,
        width=pw - 72,
        height=150,
        preserveAspectRatio=True,
        anchor="c",
    )
    y -= 162
    y = pdf_paragraph(
        c,
        "The Phase-5 questions, seeds, denominators and gates were frozen in a repository-local decision record before the new analyses ran. Failed fits and invalid specifications stayed visible in their original denominators.",
        36,
        y,
        pw - 72,
        8.1,
        10.2,
        MUTED,
    )
    y -= 6

    y = pdf_section(c, 3, "Our five initial claims", y, pw)
    claims = [
        "Detected analysis-label richness broadens downstream robustly.",
        "A correction model can recover true occupancy of rare Conwy labels.",
        "Lower source signal and processed fish counts predict pipeline non-recovery.",
        "Longer mitochondrial markers resolve Salmo species better than MiFish.",
        "Trait diversity broadens beyond what richness and prevalence alone predict.",
    ]
    for index, claim in enumerate(claims, start=1):
        c.setFillColor(colors.HexColor(WHITE))
        c.setStrokeColor(colors.HexColor(LINE))
        c.roundRect(36, y - 36, pw - 72, 36, 5, fill=1, stroke=1)
        c.setFillColor(colors.HexColor(RIVER_DARK))
        c.circle(54, y - 18, 9, fill=1, stroke=0)
        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 8.2)
        c.drawCentredString(54, y - 21, str(index))
        pdf_paragraph(c, claim, 72, y - 13, pw - 125, 8.3, 10, INK, True)
        y -= 43
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica-Oblique", 7.2)
    c.drawString(36, 39, "These were hypotheses formulated by this reanalysis, not claims attributed to the source-study authors.")
    c.setFont("Helvetica", 7.2)
    c.drawRightString(pw - 36, 24, "1 / 2")
    c.showPage()

    pdf_header(
        c,
        "Stress tests, models, verdict",
        "One survives | one conditional | two withheld | one failed",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 4, "Stressing our claims with more tests", y, pw)
    c.drawImage(
        ImageReader(figures[3]),
        36,
        y - 220,
        width=pw - 72,
        height=220,
        preserveAspectRatio=True,
        anchor="c",
    )
    y -= 228
    y = pdf_paragraph(
        c,
        "Every result is reported against its declared denominator: 64/72 community families (reference p=0.0714), 613/1,944 recoverable estimator-scenarios, pipeline coefficient -1.149 with 999/999 finite bootstraps, 2 of 7 eligible Salmo labels, and trait-null p=0.3159 / 0.1142.",
        36,
        y,
        pw - 72,
        7.9,
        9.8,
        MUTED,
    )
    y -= 6

    y = pdf_section(c, 5, "The models and simulations used", y, pw)
    c.drawImage(
        ImageReader(figures[5]),
        36,
        y - 205,
        width=pw - 72,
        height=205,
        preserveAspectRatio=True,
        anchor="c",
    )
    y -= 215

    y = pdf_section(c, 6, "Our final verdict", y, pw)
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.roundRect(36, 61, pw - 72, y - 55, 7, fill=1, stroke=0)
    pdf_paragraph(
        c,
        "The Conwy data support a stable descriptive downstream broadening of detected fish labels and a bounded, comparator-relative pipeline association. They do not support rare-label occupancy correction, a trait mechanism beyond richness/prevalence, or species-level identity within Salmo.",
        51,
        y - 21,
        pw - 102,
        9.1,
        11.5,
        WHITE,
        True,
    )
    pdf_paragraph(
        c,
        "Scientific addition: a reproducible credibility-audit framework, pre-specified in a frozen repository-local decision record, that shows which eDNA claims are identifiable, conditional, or must be withheld.",
        51,
        101,
        pw - 102,
        8.0,
        10,
        "#dbecea",
    )
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica", 7.2)
    c.drawString(36, 39, "Phase-5 decision record, claim register, figures and clean-extraction verification are bundled offline.")
    c.drawRightString(pw - 36, 24, "2 / 2")
    c.save()
    return out


def pdf_page_base(c, page_width, page_height):
    c.setFillColor(colors.white)
    c.rect(0, 0, page_width, page_height, fill=1, stroke=0)


def pdf_footer(c, page_number, total_pages, page_width, note):
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica", 7.2)
    c.drawString(36, 24, note)
    c.drawRightString(page_width - 36, 24, f"{page_number} / {total_pages}")


def pdf_model_caption(c, y, caption, target, page_width):
    """Draw a selectable model result plus a portable relative evidence link."""
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica", 6.7)
    c.drawString(39, y, caption)
    link_label = "OPEN MODEL + EVIDENCE >"
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.setFont("Helvetica-Bold", 6.4)
    link_width = stringWidth(link_label, "Helvetica-Bold", 6.4)
    link_x = page_width - 39 - link_width
    c.drawString(link_x, y, link_label)
    c.linkURL(
        target,
        (link_x - 2, y - 2, link_x + link_width + 2, y + 8),
        relative=0,
        thickness=0,
    )


def pdf_info_card(c, x, y, width, height, heading, body, accent=RIVER):
    c.setFillColor(colors.HexColor(WHITE))
    c.setStrokeColor(colors.HexColor(LINE))
    c.roundRect(x, y, width, height, 5, fill=1, stroke=1)
    c.setFillColor(colors.HexColor(accent))
    c.rect(x, y, 5, height, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.setFont("Helvetica-Bold", 8.2)
    c.drawString(x + 14, y + height - 19, heading)
    pdf_paragraph(c, body, x + 14, y + height - 36, width - 27, 7.3, 9.2, MUTED)


def pdf_claim_row(c, index, claim, y, page_width):
    claim_ids = [
        "river-community",
        "imperfect-detection",
        "pipeline-disagreement",
        "taxonomic-resolution",
        "trait-diversity",
    ]
    c.setFillColor(colors.HexColor(WHITE))
    c.setStrokeColor(colors.HexColor(LINE))
    c.roundRect(36, y - 47, page_width - 72, 47, 5, fill=1, stroke=1)
    c.setFillColor(colors.HexColor(PRESENTATION_DECISIONS[claim_ids[index - 1]]["colour"]))
    c.rect(36, y - 47, 5, 47, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.circle(57, y - 23.5, 9, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 8.2)
    c.drawCentredString(57, y - 26.5, str(index))
    pdf_paragraph(c, claim, 75, y - 18, page_width - 155, 8.2, 10, INK, True)
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica-Bold", 6.8)
    c.drawRightString(page_width - 50, y - 27, f"HYPOTHESIS {index:02d}")
    return y - 54


def pdf_decision_row(c, top, claim, decision, reason, colour, page_width):
    """Draw one compact verdict row using the shared three-state vocabulary."""
    height = 22
    bottom = top - height
    c.setFillColor(colors.HexColor(WHITE))
    c.setStrokeColor(colors.HexColor(LINE))
    c.roundRect(36, bottom, page_width - 72, height, 4, fill=1, stroke=1)
    c.setFillColor(colors.HexColor(colour))
    c.rect(36, bottom, 4, height, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(INK))
    c.setFont("Helvetica-Bold", 6.4)
    c.drawString(48, bottom + 7.5, claim)

    pill_x, pill_width = 145, 86
    c.setFillColor(colors.HexColor(colour))
    c.roundRect(pill_x, bottom + 4, pill_width, 14, 7, fill=1, stroke=0)
    c.setFillColor(colors.white)
    c.setFont("Helvetica-Bold", 5.4)
    c.drawCentredString(pill_x + pill_width / 2, bottom + 8.5, decision)

    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica", 6.4)
    c.drawString(242, bottom + 7.5, f"Reason: {reason}")


def research_brief_pdf(figures) -> Path:
    """Build the four-page visual brief with five model diagrams."""
    out = PDFS / "RESEARCH_BRIEF.pdf"
    pw, ph = A4
    c = canvas.Canvas(str(out), pagesize=A4, invariant=1)
    c.setTitle("River Conwy eDNA - four-page investigation brief")
    c.setSubject("Six-part scientific claim audit with model and simulation diagrams")
    total_pages = 4

    # Page 1: question and documented research trail.
    pdf_page_base(c, pw, ph)
    c.bookmarkPage("question-and-method")
    c.addOutlineEntry("1-2 Question and methodology", "question-and-method", level=0)
    pdf_header(
        c,
        "Question, evidence base, and research trail",
        "River Conwy eDNA | four-page visual brief",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 1, "Our question", y, pw)
    y = pdf_paragraph(
        c,
        "What can a replicated river-eDNA dataset establish after correcting the unit of inference, testing missed-detection recovery, auditing pipeline disagreement, and challenging taxonomic and trait interpretations?",
        36,
        y,
        pw - 72,
        9.5,
        12.2,
        INK,
        True,
    )
    y -= 5
    c.setFillColor(colors.HexColor("#eaf3f4"))
    c.rect(36, y - 29, pw - 72, 33, fill=1, stroke=0)
    pdf_paragraph(
        c,
        "14 sites | 18 usable timepoints | 647 validated field libraries | 1,348 ASVs | replicated 12S survey",
        47,
        y - 8,
        pw - 94,
        8.3,
        10,
        RIVER_DARK,
        True,
    )
    y -= 49
    y = pdf_paragraph(
        c,
        "Decision target: separate a stable descriptive pattern from mechanisms or corrections that the current evidence cannot identify.",
        36,
        y,
        pw - 72,
        8.2,
        10.5,
        MUTED,
        True,
    )
    y -= 8

    y = pdf_section(c, 2, "Our methodology + key research documentation moments", y, pw)
    timeline_height = 224
    c.drawImage(
        ImageReader(figures[4]),
        36,
        y - timeline_height,
        width=pw - 72,
        height=timeline_height,
        preserveAspectRatio=True,
        anchor="c",
    )
    y -= timeline_height + 12
    card_height = 76
    pdf_info_card(
        c,
        36,
        y - card_height,
        pw - 72,
        card_height,
        "WHAT RE-RUNNING PROVES",
        "Running the workflow twice showed that the same committed inputs produced the same files and that the offline links worked. This confirms computational reproducibility; it does not add new biological replication or independent field evidence.",
        MOSS,
    )
    pdf_footer(
        c,
        1,
        total_pages,
        pw,
        "Rules written before testing | unsuccessful models and invalid analyses remain visible",
    )
    c.showPage()

    # Page 2: the five starting claims and their falsification tests.
    pdf_page_base(c, pw, ph)
    c.bookmarkPage("claims-and-tests")
    c.addOutlineEntry("3-4 Initial claims and stress tests", "claims-and-tests", level=0)
    pdf_header(
        c,
        "Five claims enter the audit",
        "Each hypothesis receives one decision label and one reason",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 3, "Our five initial claims", y, pw)
    claims = [
        "Detected analysis-label richness broadens downstream robustly.",
        "A correction model can recover true occupancy of rare Conwy labels.",
        "Lower source signal and processed fish counts predict pipeline non-recovery.",
        "Longer mitochondrial markers resolve Salmo species better than MiFish.",
        "Trait diversity broadens beyond what richness and prevalence alone predict.",
    ]
    for index, claim in enumerate(claims, start=1):
        y = pdf_claim_row(c, index, claim, y, pw)
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica-Oblique", 7.2)
    c.drawString(36, y - 3, "These were hypotheses formulated by this reanalysis, not claims attributed to the source-study authors.")
    y -= 22

    y = pdf_section(c, 4, "Stressing those claims with more tests", y, pw)
    c.drawImage(
        ImageReader(figures[3]),
        36,
        y - 232,
        width=pw - 72,
        height=232,
        preserveAspectRatio=True,
        anchor="c",
    )
    y -= 240
    pdf_paragraph(
        c,
        "Declared evidence: 64/72 community families (reference p = 0.0714); 613/1,944 recoverable occupancy cells; pipeline beta -1.149 with 999/999 finite bootstraps; 2 of 7 eligible Salmo labels; trait-null p = 0.3159 / 0.1142.",
        36,
        y,
        pw - 72,
        7.8,
        9.6,
        MUTED,
        True,
    )
    pdf_footer(
        c,
        2,
        total_pages,
        pw,
        "Decision tally | 1 accepted | 1 rejected | 3 non-conclusive",
    )
    c.showPage()

    # Page 3: detailed architecture of models 1-3.
    pdf_page_base(c, pw, ph)
    c.bookmarkPage("models-one-to-three")
    c.addOutlineEntry("5 Models 01-03", "models-one-to-three", level=0)
    pdf_header(
        c,
        "The modelling laboratory",
        "MODEL 01-03 | multiverse, known-truth simulation, and cluster-aware association",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 5, "The models and simulations used to evaluate these claims", y, pw)
    pdf_paragraph(
        c,
        "A single model could not answer five different scientific questions. Each engine below matches its own estimand, dependence structure, stress test, and decision rule.",
        36,
        y,
        pw - 72,
        8.0,
        10,
        MUTED,
    )

    model_positions = [
        (figures[6], 503, "MODEL 01 | NON-CONCLUSIVE | Reason: 64/72 valid (<90%); direction stable", "../../website/index.html#model-community"),
        (figures[7], 321, "MODEL 02 | NON-CONCLUSIVE | Reason: real fit failed; atlas is synthetic", "../../website/index.html#model-detection"),
        (figures[8], 139, "MODEL 03 | ACCEPTED | Reason: 999/999 finite; 43/43 leave-outs negative", "../../website/index.html#model-pipeline"),
    ]
    for figure, bottom, caption, target in model_positions:
        c.drawImage(
            ImageReader(figure),
            36,
            bottom,
            width=pw - 72,
            height=160,
            preserveAspectRatio=True,
            anchor="c",
        )
        pdf_model_caption(c, bottom - 10, caption, target, pw)
    pdf_footer(
        c,
        3,
        total_pages,
        pw,
        "M1-M3 methods, machine results, full figures and executable code are bundled offline",
    )
    c.showPage()

    # Page 4: models 4-5 and the final verdict.
    pdf_page_base(c, pw, ph)
    c.bookmarkPage("models-four-five-verdict")
    c.addOutlineEntry("5-6 Models 04-05 and verdict", "models-four-five-verdict", level=0)
    pdf_header(
        c,
        "What the evidence decided - and what I added",
        "MODEL 04-05 | the same three-state decision rule, applied to identity and traits",
        pw,
        ph,
    )
    y = ph - 116
    y = pdf_section(c, 5, "Models continued: identity and mechanism", y, pw)
    for figure, bottom, caption, target in [
        (figures[9], 516, "MODEL 04 | NON-CONCLUSIVE | Reason: 2 of 7 labels eligible; minimum 4", "../../website/index.html#model-taxonomy"),
        (figures[10], 334, "MODEL 05 | REJECTED | Reason: p = 0.3159 / 0.1142; both > 0.05", "../../website/index.html#model-traits"),
    ]:
        c.drawImage(
            ImageReader(figure),
            36,
            bottom,
            width=pw - 72,
            height=160,
            preserveAspectRatio=True,
            anchor="c",
        )
        pdf_model_caption(c, bottom - 10, caption, target, pw)

    y = 299
    y = pdf_section(c, 6, "Our final verdict", y, pw)
    c.setFillColor(colors.HexColor(MUTED))
    c.setFont("Helvetica", 6.1)
    c.drawString(
        36,
        267,
        "ACCEPTED = rule met with adequate evidence | REJECTED = adequate test did not support | NON-CONCLUSIVE = insufficient to decide",
    )
    decisions = [
        ("River community", "river-community", "64/72 valid (<90%); all 64 valid directions positive"),
        ("Real occupancy", "imperfect-detection", "Real fit failed; recovery atlas is synthetic only"),
        ("Pipeline association", "pipeline-disagreement", "999/999 finite; all 43 leave-outs negative"),
        ("Species resolution", "taxonomic-resolution", "Only 2 of 7 labels eligible; minimum was 4"),
        ("Trait broadening", "trait-diversity", "p=0.3159 and 0.1142; both above 0.05"),
    ]
    row_top = 259
    for claim, claim_id, reason in decisions:
        presentation = PRESENTATION_DECISIONS[claim_id]
        pdf_decision_row(
            c,
            row_top,
            claim,
            presentation["decision"],
            reason,
            presentation["colour"],
            pw,
        )
        row_top -= 25
    c.setFillColor(colors.HexColor(RIVER_DARK))
    c.roundRect(36, 61, pw - 72, 70, 6, fill=1, stroke=0)
    pdf_paragraph(
        c,
        "Scientific addition: a reproducible credibility-audit framework whose questions, tests, and verdict rules were written in a dated project plan before analysis, showing which claims are accepted, rejected, or non-conclusive - and why.",
        50,
        105,
        pw - 100,
        8.8,
        11,
        WHITE,
        True,
    )
    pdf_footer(
        c,
        4,
        total_pages,
        pw,
        "Phase-5 methods, data, code, figures and clean-extraction checks bundled offline",
    )
    c.save()
    return out


def poster_pdf(figures) -> Path:
    out = PDFS / "POSTER.pdf"
    pw, ph = landscape(A4)
    c = canvas.Canvas(str(out), pagesize=(pw, ph), invariant=1)
    c.setTitle("River Conwy eDNA - poster")
    # Paint the whole page: an unpainted PDF background is transparent and some
    # viewers display it as black, obscuring the footer.
    c.setFillColor(colors.white); c.rect(0, 0, pw, ph, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(RIVER_DARK)); c.rect(0, ph - 82, pw, 82, fill=1, stroke=0)
    c.setFillColor(colors.white); c.setFont("Helvetica-Bold", 24)
    c.drawString(30, ph - 40, "More data, less certainty? A River Conwy eDNA reanalysis")
    c.setFont("Helvetica", 10); c.drawString(30, ph - 60, "Community ecology | imperfect detection | reference-database currency")

    margin, gap = 28, 14
    col = (pw - 2 * margin - 2 * gap) / 3
    y_top = ph - 105
    titles = ["1  River pattern", "2  Molecular detective", "3  Claim stress test"]
    captions = [
        "Across 11 common timepoints, detected-label richness broadens downstream with residual 1-3-library effort; the restricted test is p=0.0714.",
        "The same 170 bp ASV receives three 100%-identity species labels; only genus is stable.",
        "One claim is accepted, one is rejected, and three are non-conclusive; every decision includes its reason.",
    ]
    for i in range(3):
        x = margin + i * (col + gap)
        c.setFillColor(colors.HexColor("#eaf3f4")); c.rect(x, 82, col, y_top - 82, fill=1, stroke=0)
        c.setFillColor(colors.HexColor(RIVER_DARK)); c.setFont("Helvetica-Bold", 13); c.drawString(x + 12, y_top - 18, titles[i])
        figure_index = (0, 1, 3)[i]
        c.drawImage(ImageReader(figures[figure_index]), x + 10, 208, width=col - 20, height=210, preserveAspectRatio=True, anchor="c")
        pdf_paragraph(c, captions[i], x + 12, 187, col - 24, 8.7, 11, INK, True)
        c.setFillColor(colors.HexColor(MUTED)); c.setFont("Helvetica", 7.4)
        c.drawString(x + 12, 94, ["NON-CONCLUSIVE - 64/72 BELOW 90%", "NON-CONCLUSIVE - SPECIES COVERAGE", "1 ACCEPTED | 1 REJECTED | 3 NON-CONCLUSIVE"][i])
    c.setFillColor(colors.HexColor(INK)); c.setFont("Helvetica-Bold", 9)
    c.drawString(margin, 54, "Central lesson: rejected and non-conclusive outcomes remain scientific results when their reasons stay visible.")
    c.setFont("Helvetica", 7.3); c.setFillColor(colors.HexColor(MUTED))
    c.drawString(margin, 34, "Local portfolio edition. Rare-species occupancy magnitudes remain exploratory; no production MCMC, publication, or hosted service is required.")
    c.drawString(margin, 19, "Source dataset: Perry et al. (2024), Nature Communications 15:4372. doi:10.1038/s41467-024-48640-3")
    c.save()
    return out


def verify(outputs, started, community_panel):
    from pypdf import PdfReader

    required = {
        str(PDFS / "RESEARCH_BRIEF.pdf"): 4,
        str(PDFS / "POSTER.pdf"): 1,
    }
    pdf_pages = {path: len(PdfReader(path).pages) for path in required}
    failures = [path for path, pages in required.items() if pdf_pages.get(path) != pages]
    for output in outputs:
        if not Path(output).exists() or Path(output).stat().st_size < 1000:
            failures.append(str(output))
    receipt = {
        "status": "PASS" if not failures else "FAIL",
        "runtime_seconds": round(time.time() - started, 3),
        "outputs": [str(Path(p).relative_to(ROOT)) for p in outputs],
        "pdf_pages": {str(Path(k).relative_to(ROOT)): v for k, v in pdf_pages.items()},
        "synthetic_scenario": deterministic_scenario(),
        "community_gradient_panel": community_panel,
        "failures": failures,
    }
    (ROOT / "website" / "build-receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    if failures:
        raise SystemExit("Website build verification failed: " + ", ".join(failures))
    print(json.dumps(receipt, indent=2))


def main():
    started = time.time()
    ASSETS.mkdir(parents=True, exist_ok=True)
    PDFS.mkdir(parents=True, exist_ok=True)
    community_panel = refresh_community_gradient()
    figures = [
        community_figure(),
        salmo_figure(),
        detection_figure(),
        claim_survival_figure(),
        research_timeline_figure(),
        models_overview_figure(),
        model_community_multiverse_figure(),
        model_occupancy_atlas_figure(),
        model_pipeline_bootstrap_figure(),
        model_salmo_benchmark_figure(),
        model_trait_null_figure(),
    ]
    outputs = figures + [research_brief_pdf(figures), poster_pdf(figures)]
    verify(outputs, started, community_panel)


if __name__ == "__main__":
    main()
