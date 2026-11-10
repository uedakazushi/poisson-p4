#!/usr/bin/env bash
# Compile all current standalone LaTeX manuscripts using the root .latexmkrc.
# Outputs remain in gitignored work/latexmk; this script does not run CI/Actions.
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ $# -gt 0 ]]; then
  manuscripts=("$@")
else
  manuscripts=(
    papers/poisson_p4_classification_ja.tex
    papers/poisson_divisor_symmetries_ja.tex
    papers/normal_actual_divisor_recognition_ja.tex
    papers/s61_first_jet_reduction_ja.tex
    papers/ca_cb_local_normalization_ja.tex
    papers/positive_weight6_component_ja.tex
    papers/poisson_subvariety_incidence_ja.tex
    research/notes/p4_poisson_status_strategy_ja_20261030_v1.tex
  )
fi
for tex in "${manuscripts[@]}"; do
  test -f "$tex" || { echo "Missing source: $tex" >&2; exit 1; }
  latexmk -r .latexmkrc "$tex"
done
echo "Generated PDFs (not tracked): work/latexmk/"
