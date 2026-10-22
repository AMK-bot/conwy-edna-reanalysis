#!/usr/bin/env bash
# setup_env.sh -- one-command setup of the isolated native-arm64 `edna` R environment.
# Safe to re-run. Does NOT touch any other conda env.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MM="${MICROMAMBA:-$HOME/.local/bin/micromamba}"

command -v "$MM" >/dev/null || { echo "micromamba not found at $MM (set \$MICROMAMBA)"; exit 1; }

echo "==> Creating/locating conda env 'edna' (~/micromamba/envs/edna)"
if "$MM" env list | grep -qE '/envs/edna$|^\s*edna\b'; then
  echo "    env 'edna' already exists — skipping create (delete with: $MM env remove -n edna)"
else
  "$MM" create -y -n edna -f "$HERE/environment.yml"
fi

echo "==> R layer: CRAN top-ups + dada2/iNEXT source rebuilds + cmdstanr smoke test"
"$MM" run -n edna Rscript "$HERE/R/00_setup.R"

echo "==> Done. Use the env with:  $MM run -n edna <command>"
