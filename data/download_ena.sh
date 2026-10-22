#!/usr/bin/env bash
#
# download_ena.sh -- pull the River Conwy 12S read-runs from ENA study PRJEB48362.
#
# *** SAFETY: this script HARD-FILTERS to sample_alias matching E##_R#_T## ***
#     E## = River Conwy (what we want).  C## = a DIFFERENT river in the same
#     archive (the project's #1 risk). The awk filter below can never match C##.
#
# Modes:
#   ./download_ena.sh --metadata-only   refresh data/ena/PRJEB48362_read_run.tsv only (small; committed)
#   ./download_ena.sh --dry-run         list what would be downloaded + total size, fetch nothing
#   ./download_ena.sh                    download the 755 paired FASTQs (~71 GB) into $RAW
#
# Download target ($RAW) is re-targetable for the two-machine workflow: this ~71 GB
# step runs on the 64 GB Windows/WSL2 PC, NOT the 8 GB Mac. Point it at the PC's disk:
#   export EDNA_RAW=/mnt/d/eDNA_data/ena/raw   # WSL2 path to the PC's data disk
#   bash download_ena.sh
# Default (no EDNA_RAW) keeps the Mac-local data/ena/raw/ behaviour.
#
# Raw FASTQs are git-ignored (see root .gitignore). Re-runs are resumable:
# files whose md5 already matches are skipped. Needs: curl, awk, and md5sum or md5.
set -euo pipefail

STUDY="PRJEB48362"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENA_TSV="${EDNA_ENA_TSV:-$HERE/ena/${STUDY}_read_run.tsv}"
RAW="${EDNA_RAW:-$HERE/ena/raw}"
PORTAL="https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS="run_accession,experiment_accession,sample_accession,secondary_sample_accession,sample_alias,experiment_alias,run_alias,library_name,library_strategy,library_source,library_selection,experiment_title,description,tax_id,scientific_name,fastq_ftp,fastq_md5,fastq_bytes,read_count,base_count,instrument_model,first_public"
CONWY_RE='^E[0-9][0-9]_R[0-9]_T[0-9][0-9]$'   # the E## hard-filter

md5_of() { if command -v md5sum >/dev/null; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

refresh_metadata() {
  echo "Refreshing ENA metadata for $STUDY -> $ENA_TSV"
  curl -fsSL --retry 3 "${PORTAL}?accession=${STUDY}&result=read_run&fields=${FIELDS}&format=tsv" -o "$ENA_TSV"
  echo "  $(($(wc -l < "$ENA_TSV") - 1)) total read-runs in study."
}

if [[ "${1:-}" == "--metadata-only" ]]; then
  refresh_metadata
  echo "Next: python3 $HERE/verify_ena_join.py"
  exit 0
fi
DRYRUN=false; [[ "${1:-}" == "--dry-run" ]] && DRYRUN=true

[[ -f "$ENA_TSV" ]] || refresh_metadata

# Column indices for sample_alias / fastq_ftp / fastq_md5 (robust to field reordering).
read -r C_ALIAS C_FTP C_MD5 < <(head -1 "$ENA_TSV" | tr '\t' '\n' | nl -ba \
  | awk '/\tsample_alias$/{a=$1} /\tfastq_ftp$/{f=$1} /\tfastq_md5$/{m=$1} END{print a, f, m}')

# Emit "ftp_url<TAB>md5" for each FASTQ of each Conwy E## run (bash 3.2-safe; no mapfile).
emit_pairs() {
  awk -F'\t' -v A="$C_ALIAS" -v F="$C_FTP" -v M="$C_MD5" -v RE="$CONWY_RE" '
    NR>1 && $A ~ RE { n=split($F,u,";"); split($M,h,";");
      for(i=1;i<=n;i++) if(u[i]!="") print u[i]"\t"h[i] }' "$ENA_TSV"
}

n_runs=$(awk -F'\t' -v A="$C_ALIAS" -v RE="$CONWY_RE" 'NR>1 && $A ~ RE' "$ENA_TSV" | wc -l | tr -d ' ')
n_files=$(emit_pairs | wc -l | tr -d ' ')
echo "Conwy 12S (E##_R#_T##): ${n_runs} runs, ${n_files} FASTQ files."
if [ "$n_runs" -ne 755 ]; then echo "WARNING: expected 755 Conwy runs, got ${n_runs}. Investigate before downloading." >&2; fi

if $DRYRUN; then
  tot=$(awk -F'\t' -v A="$C_ALIAS" -v RE="$CONWY_RE" 'NR==1{n=split($0,H,"\t"); for(i=1;i<=n;i++) if(H[i]=="fastq_bytes") B=i}
    NR>1 && $A ~ RE {split($B,b,";"); for(i in b) s+=b[i]} END{printf "%.1f", s/1e9}' "$ENA_TSV")
  echo "DRY RUN: would download ~${tot} GB into $RAW. Nothing fetched."; exit 0
fi

mkdir -p "$RAW"
echo "Downloading into $RAW (resumable; ~71 GB total)..."
i=0
# process substitution (not a pipe) keeps the loop in the main shell -> md5 abort works.
while IFS=$'\t' read -r ftp md5; do
  url="https://${ftp}"; fn="$RAW/$(basename "$url")"; i=$((i+1))
  if [ -f "$fn" ] && [ "$(md5_of "$fn")" = "$md5" ]; then printf "[%d/%d] skip %s (md5 ok)\n" "$i" "$n_files" "$(basename "$fn")"; continue; fi
  printf "[%d/%d] get  %s\n" "$i" "$n_files" "$(basename "$fn")"
  curl -fSL --retry 3 -C - "$url" -o "$fn"
  [ "$(md5_of "$fn")" = "$md5" ] || { echo "  md5 MISMATCH for $fn" >&2; exit 1; }
done < <(emit_pairs)
echo "Done. ${n_files} files verified by md5 in $RAW"
