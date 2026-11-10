#!/usr/bin/env bash
# Compile current publication manuscripts via LuaLaTeX/latexmk (no tracked PDFs).
set -euo pipefail
cd "$(dirname "$0")/.."
if (( $# )); then
  manuscripts=("$@")
else
  manuscripts=(poisson_divisor_symmetries_ja poisson_p4_classification_ja)
fi
for manuscript in "${manuscripts[@]}"; do
  build_dir="work/tex-build/$manuscript"
  mkdir -p "$build_dir"
  latexmk -r .latexmkrc -outdir="$build_dir" "papers/$manuscript.tex" > "$build_dir/latexmk.stdout" 2>&1 || {
    tail -70 "$build_dir/latexmk.stdout" >&2; exit 1;
  }
  if command -v rg >/dev/null 2>&1 && rg -n 'Overfull|Underfull|undefined references|undefined on input|multiply defined' "$build_dir/$manuscript.log"; then
    exit 1
  fi
  echo "$build_dir/$manuscript.pdf"
done
