#!/usr/bin/env bash
# Rebuild the current Poisson subvariety incidence paper and exact certificates.
# Read-only input: Git checkout. All scratch work is confined to work/.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUT="${1:-work/subvariety-incidence-check}"
mkdir -p "$OUT/certificates" "$OUT/latex"
python3 - <<'PY'
import sympy
assert sympy.__version__ == '1.14.0', f'expected SymPy 1.14.0; found {sympy.__version__}'
PY
python3 code/subvariety_incidence/verify_del_pezzo_rigidity.py \
  --output "$OUT/certificates/del_pezzo_generic_rigidity.json"
python3 code/subvariety_incidence/check_incidence_examples.py \
  --output "$OUT/certificates/incidence_examples.json" > "$OUT/incidence_examples_stdout.log"
python3 code/subvariety_incidence/check_quadric_symmetry.py \
  "$OUT/certificates/quadric_symmetry.json" > "$OUT/quadric_symmetry_stdout.log"
python3 code/subvariety_incidence/check_log_quadric.py \
  --output "$OUT/certificates/log_quadric.json" > "$OUT/log_quadric_stdout.log"
python3 code/subvariety_incidence/check_log_wall.py \
  --output "$OUT/certificates/log_wall.json" > "$OUT/log_wall_stdout.log"
python3 - "$OUT" <<'PY'
import json, pathlib, sys
out=pathlib.Path(sys.argv[1]);root=pathlib.Path.cwd()
for name in ('del_pezzo_generic_rigidity.json','incidence_examples.json','quadric_symmetry.json','log_quadric.json','log_wall.json'):
    a=json.loads((root/'research/subvariety_incidence/certificates'/name).read_text())
    b=json.loads((out/'certificates'/name).read_text())
    assert a==b, f'certificate mismatch: {name}'
assert json.loads((out/'certificates/incidence_examples.json').read_text())['check_count']==101
rig=json.loads((out/'certificates/del_pezzo_generic_rigidity.json').read_text())
assert rig['generic_quadratic_rank']==60 and rig['linear_syzygy_count']==90
quad=json.loads((out/'certificates/quadric_symmetry.json').read_text())
assert quad['rational_strict_so5_stabilizer_rank']==10
assert quad['rational_nonzero_minor']=='27980515917168640'
log=json.loads((out/'certificates/log_quadric.json').read_text())
assert log['strict_so5_sample']['rank']==9 and log['special_log_boundary']['strict_stabilizer_dimension']==2
wall=json.loads((out/'certificates/log_wall.json').read_text())
assert wall['strict_stabilizer_minor']=='13824*c/((c - 1)**4*(c + 1)**4)'
assert wall['wall_scheme']=='Spec C[t,v]/(v(t+2))'
print('EXACT CERTIFICATES: PASS (5 exact producers and five certificates)')
PY
for pass in 1 2; do
  lualatex -interaction=nonstopmode -halt-on-error \
    -output-directory="$OUT/latex" \
    papers/poisson_subvariety_incidence_ja.tex > "$OUT/latex/pass${pass}.log" 2>&1 || {
       tail -70 "$OUT/latex/pass${pass}.log"; exit 1;
    }
done
python3 - "$OUT" <<'PY'
import pathlib,subprocess,sys,re
out=pathlib.Path(sys.argv[1]); pdf=out/'latex/poisson_subvariety_incidence_ja.pdf'
info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
m=re.search(r'^Pages:\s*(\d+)',info,re.M)
assert m and int(m.group(1))>=15, 'Unexpected page count'
log=(out/'latex/pass2.log').read_text(errors='replace')
assert not re.search(r'LaTeX Warning: (Reference|Citation).*undefined',log)
assert 'There were undefined references' not in log
assert 'Overfull \\hbox' not in log and 'Overfull \\vbox' not in log
text=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'])
import hashlib
text_hash=hashlib.sha256(text).hexdigest()
assert len(text)>10000, 'PDF text extraction unexpectedly empty'
print('LuaLaTeX: PASS (%s pages; 2 passes; extracted PDF text SHA256 %s)'%(m.group(1),text_hash))
PY
printf 'PASS: incidence source + exact certificates + generated PDF at %s\n' "$OUT"
