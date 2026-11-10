#!/usr/bin/env bash
# Manual, source-only verification; no network, Library access or GitHub Actions.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUT="${1:-work/isotropic-tangent-cone-check}"
mkdir -p "$OUT/latex"
python3 code/subvariety_incidence/check_isotropic_quadratic_cone.py \
  --output "$OUT/isotropic_quadratic_cone.json" \
  --check research/subvariety_incidence/certificates/isotropic_quadratic_cone.json
for pass in 1 2; do
  lualatex -interaction=nonstopmode -halt-on-error -output-directory="$OUT/latex" \
    papers/isotropic_log_tangent_cone_ja.tex > "$OUT/latex/pass${pass}.log" 2>&1 || {
      tail -60 "$OUT/latex/pass${pass}.log"; exit 1;
    }
done
python3 - "$OUT" <<'PY'
import pathlib,re,subprocess,sys
out=pathlib.Path(sys.argv[1]);log=(out/'latex/pass2.log').read_text(errors='replace')
assert not re.search(r'LaTeX Warning: (Reference|Citation).*undefined',log)
assert 'There were undefined references' not in log
assert 'Overfull \\hbox' not in log and 'Overfull \\vbox' not in log
pdf=out/'latex/isotropic_log_tangent_cone_ja.pdf'
text=subprocess.check_output(['pdftotext',str(pdf),'-'])
assert len(text)>15000
print('PASS: exact certificate and two-pass LuaLaTeX manuscript; PDF:',pdf)
PY
