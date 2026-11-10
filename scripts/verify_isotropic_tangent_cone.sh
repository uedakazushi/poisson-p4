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
python3 code/subvariety_incidence/check_isotropic_orbit20.py --output "$OUT/isotropic_orbit20.json"
python3 - "$OUT/isotropic_orbit20.json" <<'PY'
import json,sys
from pathlib import Path
got=json.loads(Path(sys.argv[1]).read_text())
stored=json.loads(Path('research/subvariety_incidence/certificates/isotropic_orbit20.json').read_text())
assert got==stored, 'isotropic orbit certificate mismatch'
assert (got['GL5_affine_orbit_rank'],got['with_P0_scale_rank'],got['PGL5_projective_orbit_dimension'])==(21,21,20)
print('PASS: exact orbit rank and tracked certificate')
PY
python3 code/subvariety_incidence/check_isotropic_intersection_dim23.py \
  --output "$OUT/isotropic_intersection_dim23.json"
python3 - "$OUT/isotropic_intersection_dim23.json" <<'PY'
import json,sys
from pathlib import Path
got=json.loads(Path(sys.argv[1]).read_text())
expected=json.loads(Path('research/subvariety_incidence/certificates/isotropic_intersection_dim23.json').read_text())
assert got==expected,'isotropic intersection certificate mismatch'
assert got['intersection_dimension_bounds']==[23,24]
assert got['24x24_minor_determinant']=='10871635968'
assert got['gauge_plus_family_rank']==48
print('PASS: 23-dimensional explicit Poisson intersection family')
PY
python3 code/subvariety_incidence/check_universal_nilpotent_pencil.py \
  --output "$OUT/universal_nilpotent_pencil.json" \
  --check research/subvariety_incidence/certificates/universal_nilpotent_pencil.json
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
