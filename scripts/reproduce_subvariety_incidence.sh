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
python3 code/subvariety_incidence/check_log_orthogonality.py \
  --output "$OUT/certificates/log_orthogonality.json" > "$OUT/log_orthogonality_stdout.log"
python3 code/subvariety_incidence/check_log_rankone.py \
  --output "$OUT/certificates/log_rankone.json" > "$OUT/log_rankone_stdout.log"
python3 code/subvariety_incidence/check_isotropic_log_boundaries.py \
  --output "$OUT/certificates/isotropic_log_boundaries.json" > "$OUT/isotropic_log_boundaries_stdout.log"
python3 code/subvariety_incidence/check_ambient_quadratic_obstruction.py \
  --output "$OUT/certificates/ambient_quadratic_obstruction.json" > "$OUT/ambient_obstruction_stdout.log"
python3 code/subvariety_incidence/check_local_dimension28.py \
  --output "$OUT/certificates/local_dimension28.json" > "$OUT/local_dimension28_stdout.log"
python3 code/subvariety_incidence/check_rank4_long_nonconstant_membership.py \
  --output "$OUT/certificates/rank4_long_nonconstant_membership.json" > "$OUT/rank4_long_membership_stdout.log"
python3 code/subvariety_incidence/check_dp_pfaffian_conormal.py \
  --output "$OUT/certificates/dp_pfaffian_conormal.json" \
  --check research/subvariety_incidence/certificates/dp_pfaffian_conormal.json \
  > "$OUT/dp_pfaffian_conormal_stdout.log"
python3 code/subvariety_incidence/check_dp_conormal_lie_pfaffian.py \
  --output "$OUT/certificates/dp_conormal_lie_pfaffian.json" \
  --check research/subvariety_incidence/certificates/dp_conormal_lie_pfaffian.json \
  > "$OUT/dp_conormal_lie_stdout.log"
python3 - "$OUT" <<'PY'
import json, pathlib, sys
out=pathlib.Path(sys.argv[1]);root=pathlib.Path.cwd()
for name in ('del_pezzo_generic_rigidity.json','incidence_examples.json','quadric_symmetry.json','log_quadric.json','log_wall.json','log_orthogonality.json','log_rankone.json','isotropic_log_boundaries.json','ambient_quadratic_obstruction.json','local_dimension28.json','rank4_long_nonconstant_membership.json','dp_pfaffian_conormal.json','dp_conormal_lie_pfaffian.json'):
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
ortho=json.loads((out/'certificates/log_orthogonality.json').read_text())
assert ortho['status']=='PASS' and len(ortho['samples'])==5
assert [r['strict_stabilizer_dimension'] for r in ortho['samples']]==[2,1,1,1,1]
rankone=json.loads((out/'certificates/log_rankone.json').read_text())
assert rankone['status']=='PASS'
assert [r['strict_stabilizer_dimension'] for r in rankone['cases']]==[2,1]
assert rankone['cases'][0]['two_parameter_Jacobi_ideal']=='(v*(t-1))'
isotropic=json.loads((out/'certificates/isotropic_log_boundaries.json').read_text())
assert [x['exact_Jacobi_Groebner_basis'] for x in isotropic['models']] == [[], ['a*b - 4*a'], ['a*c','c**2']]
assert [x['rank_q_restriction_to_W'] for x in isotropic['models']] == [1,1,0]
amb=json.loads((out/'certificates/ambient_quadratic_obstruction.json').read_text())
assert (amb['linearization_rank'],amb['augmented_rank'])==(95,96)
assert amb['projective_tangent_dimension']==30
assert amb['obstruction_on_negative_Jacobi']=='9'
loc=json.loads((out/'certificates/local_dimension28.json').read_text())
assert loc['exact_local_dimension']==28
assert loc['logarithmic_derivative_rank']==28
assert loc['logarithmic_projective_augmented_rank']==29
assert loc['logarithmic_minor_det']=='-57177414000000'
assert loc['tangent_quadrics']==['9*t21*t29','-4*t25*t27']
assert loc['ambient_tangent_frame_rank_with_gauge_and_scale']==55
long=json.loads((out/'certificates/rank4_long_nonconstant_membership.json').read_text())
assert long['status']=='PASS'
assert long['component']=='R4_MODULAR_RANK1_LONG_NONCONSTANT_25'
assert long['source_bracket_match']=='all six exact over QQ(a)'
assert long['quadric_Hessian_det']=='-1/(32*a**4)'
assert long['affine_change_det']=='-1/a**3'
assert long['modular_Jordan_a_nonzero']==[3,1,1] and long['modular_Jordan_a_zero']==[2,2,1]
pf=json.loads((out/'certificates/dp_pfaffian_conormal.json').read_text())
assert (pf['degree5_ideal_S_dimension'],pf['degree5_ideal_S_times_ideal_C_dimension'])==(65,41)
assert pf['degree5_extra_linear_conditions_on_quintics_containing_S']==24
assert pf['smooth_del_pezzo_S_with_fixed_C_affine_fibre_dimension']==10
print('EXACT CERTIFICATES: PASS (13 exact producers and thirteen certificates)')
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
