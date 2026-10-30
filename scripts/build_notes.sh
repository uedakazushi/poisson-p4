#!/usr/bin/env bash
# Rebuild the existing Japanese assessment; no changes to immutable originals.
set -euo pipefail
repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
notes_source="$repo_root/research/notes/p4_poisson_status_strategy_ja_20261030_v1.tex"
notes_output="$repo_root/build/notes"
if ! command -v lualatex >/dev/null 2>&1; then
  printf '%s\n' 'LuaLaTeX is unavailable. Install a TeX distribution with luatexja and HaranoAji; the checked-in PDF remains available.' >&2
  exit 1
fi
mkdir -p "$notes_output"
cd "$repo_root/research/notes"
# The initial table of contents occupies additional space. Three passes are
# needed from a clean build to stabilize its page references in this document.
for build_pass in 1 2 3; do
  lualatex -interaction=nonstopmode -halt-on-error -output-directory="$notes_output" "$notes_source" > "$notes_output/pass-${build_pass}.stdout" 2>&1 || {
    printf 'LuaLaTeX failed; inspect %s/pass-%s.stdout\n' "$notes_output" "$build_pass" >&2
    exit 1
  }
done
printf 'Built %s/p4_poisson_status_strategy_ja_20261030_v1.pdf\n' "$notes_output"
