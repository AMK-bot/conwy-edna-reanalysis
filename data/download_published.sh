#!/usr/bin/env bash
#
# download_published.sh -- fetch the published artifacts of the source paper:
#
#   Perry, W.B. et al. (2024) "An integrated spatio-temporal view of riverine
#   biodiversity using environmental DNA metabarcoding." Nature Communications
#   15:4372.  doi:10.1038/s41467-024-48640-3  (open access, CC BY 4.0)
#
# These are third-party files (not committed; see ../.gitignore). They are the
# inputs to the fast-path Phase-1 baseline (Source Data) and the occupancy
# covariates (Supplementary Data 1), plus the authors' own R scripts.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE="https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-48640-3/MediaObjects"

get() { echo "  $2"; curl -fsSL --retry 3 "$BASE/$1" -o "$2"; }

echo "Fetching Perry et al. 2024 (Nat Commun) supplementary artifacts..."
get "41467_2024_48640_MOESM7_ESM.xlsx" "$HERE/published/source_data_species_x_replicate.xlsx"   # Source Data
get "41467_2024_48640_MOESM4_ESM.csv"  "$HERE/metadata/supplementary_data_1_covariates.csv"      # Supp Data 1 (covariates)
get "41467_2024_48640_MOESM5_ESM.zip"  "$HERE/published/supplementary_code_1.zip"                 # Supp Code 1 (authors' R)
get "41467_2024_48640_MOESM1_ESM.pdf"  "$HERE/../docs/supplementary_information.pdf"               # Supp Info

echo "Unpacking authors' R code..."
rm -rf "$HERE/published/supplementary_code_1"
unzip -o -q "$HERE/published/supplementary_code_1.zip" -d "$HERE/published/supplementary_code_1"
echo "Done. (Includes Analysing_BLAST_output_WP2_12S.R — the 12S taxonomy script.)"
