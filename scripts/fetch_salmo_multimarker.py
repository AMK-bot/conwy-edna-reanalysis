#!/usr/bin/env python3
"""Fetch the frozen Salmo multi-marker accession panel from NCBI Nucleotide.

This is the only networked step in the benchmark.  It queries accession-sorted
complete mitochondrial primary records, caches the official E-utilities
receipts and every GenBank record actually assessed, and extracts the three
preregistered marker windows for the offline R analysis.

The script deliberately does not calculate distances or evaluate claim gates.
Those operations belong to ``R/12_salmo_multimarker_benchmark.R`` and run from
the committed cache without a network connection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import ssl
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL_NAME = "conwy_edna_impact_sprint"
USER_AGENT = "conwy-edna-impact-sprint/1.0 (NCBI E-utilities acquisition)"
REQUEST_INTERVAL_SECONDS = 0.38  # safely below NCBI's unauthenticated limit

LABELS: Tuple[Tuple[str, str], ...] = (
    ("Salmo trutta", "salmo_candidate"),
    ("Salmo salar", "salmo_candidate"),
    ("Salmo carpio", "salmo_candidate"),
    ("Salmo marmoratus", "salmo_candidate"),
    ("Salmo macedonicus", "salmo_candidate"),
    ("Salmo macrostigma", "salmo_candidate"),
    ("Salmo obtusirostris", "salmo_candidate"),
    ("Salvelinus alpinus", "outgroup"),
    ("Oncorhynchus mykiss", "outgroup"),
)

MAX_QUERY_RECORDS = 50
MAX_SELECTED_RECORDS = 5
MIN_MITOGENOME_BP = 15000
MAX_MITOGENOME_BP = 20000

FORWARD_PRIMER = "GTCGGTAAAACTCGTGCCAGC"
REVERSE_PRIMER = "CATAGTGGGGTATCTAATCCCAGTTTG"
MAX_PRIMER_MISMATCHES = 3
MIN_MIFISH_BP = 140
MAX_MIFISH_BP = 250
COI_BASES = 648
CYTB_BASES = 1140


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def slug(label: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def hamming(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires equal lengths")
    return sum(a != b for a, b in zip(left, right))


def tsv_write(path: Path, fieldnames: Sequence[str], rows: Iterable[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


class NcbiClient:
    def __init__(self) -> None:
        # The framework Python on macOS can lack a populated private CA file.
        # Prefer the system CA bundle when present; TLS verification is never
        # disabled.
        system_ca = Path("/etc/ssl/cert.pem")
        self.context = (
            ssl.create_default_context(cafile=str(system_ca))
            if system_ca.exists()
            else ssl.create_default_context()
        )
        self._last_request = 0.0

    def get(self, endpoint: str, params: Dict[str, object]) -> Tuple[bytes, str]:
        elapsed = time.monotonic() - self._last_request
        if elapsed < REQUEST_INTERVAL_SECONDS:
            time.sleep(REQUEST_INTERVAL_SECONDS - elapsed)
        query = dict(params)
        query.setdefault("tool", TOOL_NAME)
        url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(query)}"
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        with urllib.request.urlopen(request, timeout=60, context=self.context) as response:
            payload = response.read()
        self._last_request = time.monotonic()
        return payload, url


@dataclass
class Feature:
    key: str
    location: str
    qualifiers: Dict[str, List[str]] = field(default_factory=dict)


@dataclass
class GenBankRecord:
    accession_version: str
    definition: str
    organism: str
    sequence: str
    features: List[Feature]


def parse_genbank(text: str) -> GenBankRecord:
    lines = text.splitlines()
    definition_parts: List[str] = []
    in_definition = False
    accession = ""
    version = ""
    features: List[Feature] = []
    current: Optional[Feature] = None
    current_qualifier: Optional[str] = None
    in_features = False
    in_origin = False
    sequence_parts: List[str] = []

    for line in lines:
        tag = line[:12].strip() if len(line) >= 12 else line.strip()
        if line.startswith("DEFINITION"):
            in_definition = True
            definition_parts = [line[12:].strip()]
            continue
        if in_definition and line.startswith("            ") and tag == "":
            definition_parts.append(line[12:].strip())
            continue
        if in_definition:
            in_definition = False
        if line.startswith("ACCESSION"):
            accession = line[12:].strip().split()[0]
            continue
        if line.startswith("VERSION"):
            version = line[12:].strip().split()[0]
            continue
        if line.startswith("FEATURES"):
            in_features = True
            in_origin = False
            continue
        if line.startswith("ORIGIN"):
            in_features = False
            in_origin = True
            if current is not None:
                features.append(current)
                current = None
            continue
        if line.startswith("//"):
            in_origin = False
            continue
        if in_origin:
            sequence_parts.append(re.sub(r"[^A-Za-z]", "", line))
            continue
        if not in_features:
            continue

        key = line[5:21].strip() if len(line) >= 21 else ""
        body = line[21:].strip() if len(line) >= 21 else ""
        if key:
            if current is not None:
                features.append(current)
            current = Feature(key=key, location=body)
            current_qualifier = None
            continue
        if current is None:
            continue
        if body.startswith("/"):
            match = re.match(r"/([^=]+)(?:=(.*))?$", body)
            if not match:
                continue
            qkey, value = match.group(1), match.group(2) or ""
            value = value.strip()
            if value.startswith('"'):
                value = value[1:]
            closed = value.endswith('"')
            if closed:
                value = value[:-1]
            current.qualifiers.setdefault(qkey, []).append(value)
            current_qualifier = None if closed or not match.group(2) else qkey
            continue
        if current_qualifier:
            continuation = body
            closed = continuation.endswith('"')
            if closed:
                continuation = continuation[:-1]
            current.qualifiers[current_qualifier][-1] += continuation
            if closed:
                current_qualifier = None
        elif body:
            current.location += body

    if current is not None:
        features.append(current)
    record_sequence = "".join(sequence_parts).upper()
    source = next((item for item in features if item.key == "source"), None)
    organism = source.qualifiers.get("organism", [""])[0] if source else ""
    accession_version = version or accession
    if not accession_version or not record_sequence:
        raise ValueError("GenBank record lacks VERSION or ORIGIN sequence")
    return GenBankRecord(
        accession_version=accession_version,
        definition=" ".join(definition_parts),
        organism=organism,
        sequence=record_sequence,
        features=features,
    )


def feature_value(feature: Feature, key: str) -> str:
    values = feature.qualifiers.get(key, [])
    return values[0] if values else ""


def extract_location(sequence: str, location: str) -> Tuple[str, str, int, int]:
    compact = re.sub(r"\s+", "", location)
    complement = compact.startswith("complement(")
    intervals: List[Tuple[int, int]] = []
    for match in re.finditer(r"<?(\d+)\.\.>?(\d+)|<?(\d+)>?", compact):
        if match.group(1):
            start, end = int(match.group(1)), int(match.group(2))
        else:
            start = end = int(match.group(3))
        intervals.append((start, end))
    if not intervals:
        raise ValueError(f"unsupported feature location: {location}")
    chunks: List[str] = []
    for start, end in intervals:
        if start <= end:
            chunks.append(sequence[start - 1 : end])
        else:  # circular wrap
            chunks.append(sequence[start - 1 :] + sequence[:end])
    extracted = "".join(chunks)
    if complement:
        extracted = reverse_complement(extracted)
    return extracted, "-" if complement else "+", intervals[0][0], intervals[-1][1]


def normalized_gene_name(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def find_cds(record: GenBankRecord, marker: str) -> Feature:
    if marker == "coi_5prime":
        gene_names = {"COI", "CO1", "COX1", "MTCO1"}
        product_pattern = re.compile(r"cytochrome c oxidase subunit (?:i|1)\b", re.I)
    elif marker == "cytb":
        gene_names = {"CYTB", "COB", "MTCYB"}
        product_pattern = re.compile(r"cytochrome b\b", re.I)
    else:
        raise ValueError(f"unknown CDS marker: {marker}")
    matches: List[Feature] = []
    for feature in record.features:
        if feature.key != "CDS":
            continue
        gene = normalized_gene_name(feature_value(feature, "gene"))
        product = feature_value(feature, "product")
        if gene in gene_names or product_pattern.search(product):
            matches.append(feature)
    if len(matches) != 1:
        raise ValueError(f"{marker}: expected one annotated CDS, found {len(matches)}")
    return matches[0]


def extract_cds_window(record: GenBankRecord, marker: str, bases: int) -> Dict[str, object]:
    feature = find_cds(record, marker)
    coding, strand, start, end = extract_location(record.sequence, feature.location)
    codon_start_text = feature_value(feature, "codon_start") or "1"
    try:
        codon_start = int(codon_start_text)
    except ValueError as exc:
        raise ValueError(f"{marker}: invalid codon_start={codon_start_text}") from exc
    coding = coding[codon_start - 1 :]
    if len(coding) < bases:
        raise ValueError(f"{marker}: annotated CDS has {len(coding)} coding bases, needs {bases}")
    window = coding[:bases].upper()
    if re.search(r"[^ACGT]", window):
        raise ValueError(f"{marker}: first {bases} coding bases contain ambiguity")
    return {
        "sequence": window,
        "source_feature": f"CDS:{feature_value(feature, 'gene') or feature_value(feature, 'product')}",
        "source_location": feature.location,
        "strand": strand,
        "start_1based": start,
        "end_1based": end,
        "forward_primer_mismatches": "",
        "reverse_primer_mismatches": "",
        "forward_primer_site_1based": "",
        "reverse_primer_site_1based": "",
        "orientation": strand,
    }


def circular_slice(sequence: str, start: int, length: int) -> str:
    n = len(sequence)
    return "".join(sequence[(start + offset) % n] for offset in range(length))


def mifish_candidates(sequence: str, orientation: str) -> List[Dict[str, object]]:
    reverse_site = reverse_complement(REVERSE_PRIMER)
    n = len(sequence)
    forward_hits: List[Tuple[int, int]] = []
    reverse_hits: List[Tuple[int, int]] = []
    for pos in range(n):
        f_window = circular_slice(sequence, pos, len(FORWARD_PRIMER))
        f_mm = hamming(f_window, FORWARD_PRIMER)
        if f_mm <= MAX_PRIMER_MISMATCHES:
            forward_hits.append((pos, f_mm))
        r_window = circular_slice(sequence, pos, len(reverse_site))
        r_mm = hamming(r_window, reverse_site)
        if r_mm <= MAX_PRIMER_MISMATCHES:
            reverse_hits.append((pos, r_mm))
    products: List[Dict[str, object]] = []
    for f_pos, f_mm in forward_hits:
        inner_start = (f_pos + len(FORWARD_PRIMER)) % n
        for r_pos, r_mm in reverse_hits:
            inner_length = (r_pos - inner_start) % n
            if MIN_MIFISH_BP <= inner_length <= MAX_MIFISH_BP:
                product = circular_slice(sequence, inner_start, inner_length)
                products.append(
                    {
                        "sequence": product,
                        "source_feature": "in_silico_MiFish-U_primer_product",
                        "source_location": f"{inner_start + 1}..{((r_pos - 1) % n) + 1}",
                        "strand": orientation,
                        "start_1based": inner_start + 1,
                        "end_1based": ((r_pos - 1) % n) + 1,
                        "forward_primer_mismatches": f_mm,
                        "reverse_primer_mismatches": r_mm,
                        "forward_primer_site_1based": f_pos + 1,
                        "reverse_primer_site_1based": r_pos + 1,
                        "orientation": orientation,
                    }
                )
    return products


def extract_mifish(record: GenBankRecord) -> Dict[str, object]:
    products = mifish_candidates(record.sequence, "+")
    reverse_products = mifish_candidates(reverse_complement(record.sequence), "-")
    products.extend(reverse_products)
    unique: Dict[Tuple[str, str, int, int], Dict[str, object]] = {}
    for product in products:
        key = (
            str(product["orientation"]),
            str(product["sequence"]),
            int(product["forward_primer_site_1based"]),
            int(product["reverse_primer_site_1based"]),
        )
        unique[key] = product
    products = list(unique.values())
    if len(products) != 1:
        raise ValueError(f"mifish_12S: expected one primer product, found {len(products)}")
    product = products[0]
    if re.search(r"[^ACGT]", str(product["sequence"])):
        raise ValueError("mifish_12S: product contains ambiguity")
    return product


def source_identity(record: GenBankRecord) -> Tuple[str, str, str]:
    source = next((feature for feature in record.features if feature.key == "source"), None)
    if source is None:
        return f"accession:{record.accession_version}", "no", "source feature missing"
    for key in ("specimen_voucher", "isolate", "strain", "haplotype", "clone"):
        value = feature_value(source, key).strip()
        if value:
            return f"{key}:{value}", "yes", key
    return (
        f"accession:{record.accession_version}",
        "no",
        "no voucher/isolate/strain/haplotype/clone qualifier; record-level independence only",
    )


def extract_all_markers(record: GenBankRecord) -> Dict[str, Dict[str, object]]:
    return {
        "mifish_12S": extract_mifish(record),
        "coi_5prime": extract_cds_window(record, "coi_5prime", COI_BASES),
        "cytb": extract_cds_window(record, "cytb", CYTB_BASES),
    }


def complete_primary_record(record: GenBankRecord, expected_label: str) -> Optional[str]:
    if record.organism != expected_label:
        return f"organism_annotation_not_exact:{record.organism or 'missing'}"
    if not (MIN_MITOGENOME_BP <= len(record.sequence) <= MAX_MITOGENOME_BP):
        return f"length_outside_{MIN_MITOGENOME_BP}_{MAX_MITOGENOME_BP}"
    definition = record.definition.lower()
    complete_wording = (
        "complete genome" in definition
        or "complete sequence" in definition
        or "mitochondrial genome" in definition
        or "mitochondrion, complete" in definition
        # LT617630 contains the submitter typos "mitochondrionial geneome";
        # its circular 16,677-bp primary record is otherwise complete and fully
        # annotated.  Preserve the official text but do not turn that typo into
        # a biological exclusion.
        or ("mitochond" in definition and "geneome" in definition)
    )
    if "mitochond" not in definition or not complete_wording:
        return "definition_not_complete_mitochondrial_record"
    if record.accession_version.startswith(("NC_", "NG_", "NM_", "NR_", "XM_", "XR_")):
        return "refseq_record_excluded_primary_records_only"
    return None


def query_term(label: str) -> str:
    # Length is a permissive candidate screen.  Completeness wording, exact
    # organism annotation and marker features are checked from each raw record.
    return (
        f'"{label}"[Organism] AND mitochondrion[filter] '
        f'AND {MIN_MITOGENOME_BP}:{MAX_MITOGENOME_BP}[SLEN] '
        "NOT srcdb_refseq[prop]"
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Repository root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--retrieval-date",
        default=date.today().isoformat(),
        help="ISO retrieval date written to deterministic provenance tables",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace this script's existing cache and metadata files",
    )
    args = parser.parse_args(argv)
    try:
        date.fromisoformat(args.retrieval_date)
    except ValueError as exc:
        parser.error(f"--retrieval-date must be ISO YYYY-MM-DD: {exc}")

    root = args.project_root.resolve()
    public_dir = root / "data" / "public" / "salmo_multimarker"
    receipt_dir = public_dir / "query_receipts"
    genbank_dir = public_dir / "genbank"
    marker_fasta = public_dir / "marker_sequences.fasta"
    query_tsv = root / "data" / "metadata" / "impact_salmo_queries.tsv"
    accession_tsv = root / "data" / "metadata" / "impact_salmo_accessions.tsv"
    marker_tsv = root / "data" / "metadata" / "impact_salmo_marker_manifest.tsv"
    generated_targets = (marker_fasta, query_tsv, accession_tsv, marker_tsv)
    existing = [path for path in generated_targets if path.exists()]
    if existing and not args.overwrite:
        names = ", ".join(str(path.relative_to(root)) for path in existing)
        parser.error(f"outputs already exist ({names}); pass --overwrite to refresh deliberately")
    if args.overwrite:
        # Remove only files owned by this acquisition script so a refreshed
        # accession-sorted selection cannot leave stale, unmanifested records.
        for path in (*receipt_dir.glob("*.json"), *genbank_dir.glob("*.gb"), *generated_targets):
            if path.exists() and path.is_file():
                path.unlink()

    receipt_dir.mkdir(parents=True, exist_ok=True)
    genbank_dir.mkdir(parents=True, exist_ok=True)
    client = NcbiClient()
    query_rows: List[Dict[str, object]] = []
    accession_rows: List[Dict[str, object]] = []
    marker_rows: List[Dict[str, object]] = []
    fasta_entries: List[Tuple[str, str, str]] = []

    for label, role in LABELS:
        label_slug = slug(label)
        term = query_term(label)
        search_params: Dict[str, object] = {
            "db": "nuccore",
            "term": term,
            "retmax": MAX_QUERY_RECORDS,
            "retmode": "json",
            "sort": "accession",
        }
        search_payload, search_url = client.get("esearch.fcgi", search_params)
        search_path = receipt_dir / f"{label_slug}_esearch.json"
        search_path.write_bytes(search_payload)
        search = json.loads(search_payload.decode("utf-8"))["esearchresult"]
        ids: List[str] = list(search.get("idlist", []))

        summary_payload = b""
        summary_url = ""
        summary_receipt = ""
        summary_sha256 = ""
        summaries: List[Dict[str, object]] = []
        if ids:
            summary_payload, summary_url = client.get(
                "esummary.fcgi",
                {"db": "nuccore", "id": ",".join(ids), "retmode": "json"},
            )
            summary_json = json.loads(summary_payload.decode("utf-8"))["result"]
            summaries = [summary_json[uid] for uid in summary_json.get("uids", [])]
            summaries.sort(key=lambda item: str(item.get("accessionversion", "")))
            summary_path = receipt_dir / f"{label_slug}_esummary.json"
            summary_path.write_bytes(summary_payload)
            summary_receipt = summary_path.relative_to(root).as_posix()
            summary_sha256 = sha256_bytes(summary_payload)
        query_rows.append(
            {
                "submitted_label": label,
                "panel_role": role,
                "query_term": term,
                "sort_order": "accession",
                "max_query_records": MAX_QUERY_RECORDS,
                "ncbi_total_count": int(search.get("count", 0)),
                "returned_count": len(ids),
                "retrieval_date": args.retrieval_date,
                "esearch_url": search_url,
                "esearch_receipt": search_path.relative_to(root).as_posix(),
                "esearch_sha256": sha256_bytes(search_payload),
                "esummary_url": summary_url,
                "esummary_receipt": summary_receipt,
                "esummary_sha256": summary_sha256,
            }
        )

        selected = 0
        selected_identities: set[str] = set()
        for rank, summary in enumerate(summaries, start=1):
            accession_version = str(summary.get("accessionversion", ""))
            uid = str(summary.get("uid", ""))
            base_row: Dict[str, object] = {
                "submitted_label": label,
                "submitted_genus": label.split()[0],
                "panel_role": role,
                "candidate_rank": rank,
                "accession_version": accession_version,
                "nuccore_uid": uid,
                "record_title": str(summary.get("title", "")),
                "record_length": int(summary.get("slen", 0) or 0),
                "assessed": "no",
                "selected": "no",
                "status": "not_assessed",
                "exclusion_reason": "not_assessed_after_five_passing",
                "organism_annotation": "",
                "independence_key": "",
                "independence_documented": "",
                "independence_note": "",
                "raw_genbank_path": "",
                "raw_genbank_sha256": "",
                "retrieval_date": args.retrieval_date,
                "ncbi_record_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession_version}",
            }
            if selected >= MAX_SELECTED_RECORDS:
                accession_rows.append(base_row)
                continue
            efetch_payload, _ = client.get(
                "efetch.fcgi",
                {"db": "nuccore", "id": uid, "rettype": "gbwithparts", "retmode": "text"},
            )
            raw_path = genbank_dir / f"{accession_version}.gb"
            raw_path.write_bytes(efetch_payload)
            base_row.update(
                {
                    "assessed": "yes",
                    "raw_genbank_path": raw_path.relative_to(root).as_posix(),
                    "raw_genbank_sha256": sha256_bytes(efetch_payload),
                }
            )
            try:
                record = parse_genbank(efetch_payload.decode("utf-8"))
                base_row["organism_annotation"] = record.organism
                if record.accession_version != accession_version:
                    raise ValueError(
                        f"summary accession {accession_version} != GenBank VERSION {record.accession_version}"
                    )
                exclusion = complete_primary_record(record, label)
                if exclusion:
                    raise ValueError(exclusion)
                identity_key, identity_documented, identity_note = source_identity(record)
                base_row.update(
                    {
                        "independence_key": identity_key,
                        "independence_documented": identity_documented,
                        "independence_note": identity_note,
                    }
                )
                if identity_key in selected_identities:
                    raise ValueError(f"duplicate_voucher_or_isolate:{identity_key}")
                marker_data = extract_all_markers(record)
            except (UnicodeDecodeError, ValueError) as exc:
                base_row.update(
                    {
                        "status": "excluded",
                        "exclusion_reason": str(exc).replace("\t", " ").replace("\n", " "),
                    }
                )
                accession_rows.append(base_row)
                continue

            selected += 1
            selected_identities.add(identity_key)
            base_row.update({"selected": "yes", "status": "selected", "exclusion_reason": ""})
            accession_rows.append(base_row)
            for marker in ("mifish_12S", "coi_5prime", "cytb"):
                item = marker_data[marker]
                marker_sequence = str(item["sequence"])
                sequence_id = f"{accession_version}__{marker}"
                fasta_entries.append((sequence_id, label, marker_sequence))
                marker_rows.append(
                    {
                        "sequence_id": sequence_id,
                        "submitted_label": label,
                        "submitted_genus": label.split()[0],
                        "panel_role": role,
                        "accession_version": accession_version,
                        "independence_key": identity_key,
                        "independence_documented": identity_documented,
                        "marker": marker,
                        "sequence_length": len(marker_sequence),
                        "sequence_sha256": sha256_text(marker_sequence),
                        "source_feature": item["source_feature"],
                        "source_location": item["source_location"],
                        "strand": item["strand"],
                        "start_1based": item["start_1based"],
                        "end_1based": item["end_1based"],
                        "forward_primer_mismatches": item["forward_primer_mismatches"],
                        "reverse_primer_mismatches": item["reverse_primer_mismatches"],
                        "forward_primer_site_1based": item["forward_primer_site_1based"],
                        "reverse_primer_site_1based": item["reverse_primer_site_1based"],
                        "retrieval_date": args.retrieval_date,
                        "raw_genbank_path": raw_path.relative_to(root).as_posix(),
                        "ncbi_record_url": f"https://www.ncbi.nlm.nih.gov/nuccore/{accession_version}",
                    }
                )

    query_fields = [
        "submitted_label",
        "panel_role",
        "query_term",
        "sort_order",
        "max_query_records",
        "ncbi_total_count",
        "returned_count",
        "retrieval_date",
        "esearch_url",
        "esearch_receipt",
        "esearch_sha256",
        "esummary_url",
        "esummary_receipt",
        "esummary_sha256",
    ]
    accession_fields = [
        "submitted_label",
        "submitted_genus",
        "panel_role",
        "candidate_rank",
        "accession_version",
        "nuccore_uid",
        "record_title",
        "record_length",
        "assessed",
        "selected",
        "status",
        "exclusion_reason",
        "organism_annotation",
        "independence_key",
        "independence_documented",
        "independence_note",
        "raw_genbank_path",
        "raw_genbank_sha256",
        "retrieval_date",
        "ncbi_record_url",
    ]
    marker_fields = [
        "sequence_id",
        "submitted_label",
        "submitted_genus",
        "panel_role",
        "accession_version",
        "independence_key",
        "independence_documented",
        "marker",
        "sequence_length",
        "sequence_sha256",
        "source_feature",
        "source_location",
        "strand",
        "start_1based",
        "end_1based",
        "forward_primer_mismatches",
        "reverse_primer_mismatches",
        "forward_primer_site_1based",
        "reverse_primer_site_1based",
        "retrieval_date",
        "raw_genbank_path",
        "ncbi_record_url",
    ]
    tsv_write(query_tsv, query_fields, query_rows)
    tsv_write(accession_tsv, accession_fields, accession_rows)
    tsv_write(marker_tsv, marker_fields, marker_rows)

    marker_fasta.parent.mkdir(parents=True, exist_ok=True)
    with marker_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for sequence_id, label, sequence in fasta_entries:
            handle.write(f">{sequence_id} submitted_label={label.replace(' ', '_')}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")

    readme = public_dir / "README.md"
    readme.write_text(
        "# Cached public input: Salmo multi-marker benchmark\n\n"
        f"Retrieved from NCBI Nucleotide on {args.retrieval_date} by "
        "`scripts/fetch_salmo_multimarker.py`. Queries were accession-sorted, "
        "limited to 50 non-RefSeq mitochondrial records of 15,000--20,000 bp per "
        "frozen submitted label, and then audited from raw GenBank records. Only "
        "records with exact organism annotations, complete-mitogenome wording, "
        "a unique MiFish-U product, and complete unambiguous COI/cyt-b windows "
        "could enter the common-accession panel.\n\n"
        "The query receipts retain every returned UID and official summary. Raw "
        "GenBank files are cached only for records assessed before five passing "
        "accessions were selected; later query candidates remain visible in the "
        "accession manifest as not assessed. Public organism names are submitted "
        "annotations, not independent biological truth. RefSeq derivatives and "
        "partial-marker substitutions are deliberately excluded.\n",
        encoding="utf-8",
    )

    # Acquisition invariants: every selected accession contributes exactly one
    # sequence per marker and no label exceeds the frozen selection cap.
    selected_rows = [row for row in accession_rows if row["selected"] == "yes"]
    selected_keys = {(row["submitted_label"], row["accession_version"]) for row in selected_rows}
    marker_keys = {
        (row["submitted_label"], row["accession_version"], row["marker"]) for row in marker_rows
    }
    expected_marker_keys = {
        (label, accession, marker)
        for label, accession in selected_keys
        for marker in ("mifish_12S", "coi_5prime", "cytb")
    }
    if marker_keys != expected_marker_keys:
        raise RuntimeError("selected-accession/marker Cartesian invariant failed")
    for label, _ in LABELS:
        count = sum(row["submitted_label"] == label for row in selected_rows)
        if count > MAX_SELECTED_RECORDS:
            raise RuntimeError(f"selection cap exceeded for {label}")
    print(
        f"Fetched {len(query_rows)} label queries; selected {len(selected_rows)} "
        f"complete accessions and wrote {len(marker_rows)} marker sequences."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
