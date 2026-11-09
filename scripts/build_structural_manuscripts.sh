#!/usr/bin/env bash
# Build both current manuscripts twice, reusing optional local TeX packages.
set -euo pipefail
cd "$(dirname "$0")/.."
local_tex=/workspace/.local/poisson-tex/usr/share/texlive/texmf-dist
if ! kpsewhich ltjsarticle.cls >/dev/null 2>&1 && [[ -d "$local_tex" ]]; then
  export TEXMFHOME="$local_tex"
  export TEXMFVAR="$PWD/work/tex-cache/var"
  export TEXMFCONFIG="$PWD/work/tex-cache/config"
fi
if (( $# )); then
  manuscripts=("$@")
else
  manuscripts=(poisson_divisor_symmetries_ja poisson_p4_classification_ja)
fi
for manuscript in "${manuscripts[@]}"; do
  build_dir="work/tex-build/$manuscript"
  mkdir -p "$build_dir"
  for pass in 1 2; do
    lualatex -interaction=nonstopmode -halt-on-error -output-directory="$build_dir" \
      "papers/$manuscript.tex" > "$build_dir/pass$pass.stdout"
  done
  if rg -n 'Overfull|Underfull|undefined references|undefined on input|multiply defined' "$build_dir/$manuscript.log"; then
    exit 1
  fi
  printf '%s\n' "$build_dir/$manuscript.pdf"
done
