#!/usr/bin/env python3
"""Replay the finite-Z_(32003) raw-module certificate, without GB/CRT/QQ solving.

Extract the report archive at repository root first and source Cloud activation.
All replay outputs go to a separate work directory; frozen inputs stay untouched.
"""
from pathlib import Path
import hashlib,json,subprocess,sys,shlex,time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'work/nonsquare-v62/rank2-finite-model-v59'
OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT/'work/nonsquare-v62/finite-gcd-replay-v59'
OUT.mkdir(parents=True,exist_ok=True)
started=time.monotonic();P,e=s.symbols('P e');p=32003
loadpoly=lambda terms:sum(s.Rational(c)*P**ex[0]*e**ex[1] for ex,c in terms)
model=json.loads((DATA/'primitive-input.json').read_text())
N=s.Poly(s.sympify(model['N'],locals={'P':P,'e':e}),P,e,domain=s.QQ)
M=s.Poly(s.sympify(model['primitive_minor'],locals={'P':P,'e':e}),P,e,domain=s.QQ)
cert=json.loads((DATA/'resultant-membership-certificate.json').read_text())
R,A,B=[s.Poly(loadpoly(cert[k]),P,e,domain=s.QQ) for k in ['R','A','B']]
assert A*N+B*M==R
assert s.Poly(N.as_expr(),e).degree()==4
assert s.Poly(N.as_expr(),e).LC()==72442188
assert s.Poly(R.as_expr(),P).degree()==28
for poly in [N,M,R,A,B]:
    assert all(int(s.denom(c))%p for c in poly.coeffs())
assert int(s.numer(s.Poly(R.as_expr(),P).LC()))%p
assert int(s.denom(s.Poly(R.as_expr(),P).LC()))%p
# Bind the primitive M9 to the actual rank-two exception minor.
old=json.loads((ROOT/'work/nonsquare-v62/gcd-gate-exception-v54/certificate.json').read_text())
H=s.Matrix([[loadpoly(q) for q in row] for row in old['global_H3_entries']])
assert s.Poly(loadpoly(old['N']),P,e,domain=s.QQ)==N
unit=P**6*(20*P+3*e+20)**2*(5*P-6*e-4)
Hsmall=H.applyfunc(lambda v:s.cancel(v/unit))
for val in Hsmall:s.Poly(val,P,e,domain=s.QQ)
minor=s.Poly(Hsmall[:3,:3].det(method='domain-ge'),P,e,domain=s.QQ)
content=s.Rational(model['minor_QQ_content'])
assert minor==s.Poly(P*(20*P+3*e+20)*content*M.as_expr(),P,e,domain=s.QQ)
subprocess.run([sys.executable,str(ROOT/'code/s61_nonsquare_finite_gcd_exception_algebra_20261103_v63.py'),str(OUT)],cwd=ROOT,check=True,stdout=(OUT/'algebra-verification.log').open('w'))
flags=shlex.split(subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True))
binary=OUT/'raw-block-verifier-bin'
subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'code/s61_nonsquare_finite_gcd_exception_raw_block_20261103_v63.cpp'),'-o',str(binary),*flags],cwd=ROOT,check=True)
with (OUT/'raw-block-rank.json').open('w') as f:
    subprocess.run([str(binary),str(OUT/'multiplication-matrices.txt'),str(ROOT/'work/nonsquare-v62/raw-fourquad-GF-v36/M-QQ-native.txt'),str(DATA/'finite-model-polys-native2.txt'),str(OUT/'selected-square-inverse.txt')],cwd=ROOT,check=True,stdout=f,timeout=60)
rank=json.loads((OUT/'raw-block-rank.json').read_text());alg=json.loads((OUT/'algebra-verification.json').read_text())
assert rank['raw_block_rank']==672 and rank['selected_square_inverse_identity_verified']
assert alg['all_exact_GF_checks_passed'] and alg['FULLU_exact_Bezout_identity_verified']
result=dict(all_replay_gates_passed=True,exact_QQ_resultant_membership_verified=True,monic_e4_and_punit_monic_P28_finite_model=True,
            primitive_minor_source_identity_verified=True,FULLU_unit_on_finite_V_model_by_Nakayama=True,
            full_GF_raw_block_rank=672,raw_QQ_surjective_on_entire_closed_finite_exception_cover=True,
            no_flatness_assumed=True,no_global_GF_saturation_promotion=True,no_global_QQ_rank24_or_product41_claim=True,
            pair_and_five_quadric_gates_still_open=True,seconds=time.monotonic()-started)
(OUT/'replay-report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
