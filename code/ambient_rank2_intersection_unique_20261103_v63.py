#!/usr/bin/env python3
"""Exact symplectic reduction of the critical rank2 base-flag orbit."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from ambient_weight1_image import common_plane
from weight1_structure import J,form_matrices,flatten,coefficient_map

ROOT=Path(__file__).resolve().parents[1]

def transform(vector,g):
    inverse=g.inv();mats=form_matrices(vector)
    return flatten([(inverse.T*sum((inverse[i,k]*mats[i] for i in range(4)),s.zeros(4))*inverse).applyfunc(s.cancel) for k in range(4)])

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    params,b=common_plane();t,A,B,C,D,F=params
    value=b.subs({t:0,D:0,F:0});c=A/(2*C)
    shear=s.Matrix([[1,0,0,0],[0,1,0,-c],[c,0,1,0],[0,0,0,1]])
    assert (shear.T*J*shear-J).applyfunc(s.cancel)==s.zeros(4)
    assert (transform(value,shear)-value.subs(A,0)).applyfunc(s.cancel)==s.zeros(24,1)
    a,bs=s.symbols('a b',nonzero=True);torus=s.diag(a,1/a,bs,1/bs)
    transformed=transform(value.subs(A,0),torus)
    expected=value.subs({A:0,C:C/(a*bs**2),B:B*bs/a**2},simultaneous=True)
    assert transformed==expected
    assert s.cancel((C/(a*bs**2)).subs(a,C/bs**2))==1
    assert s.cancel((B*bs/a**2).subs(a,C/bs**2))==B*bs**5/C**2
    # H rank1 implies Delta=0. In the diagonalized chart C!=0, rankK2
    # is precisely B!=0; all2x2 nonzero minors contain B orC*B.
    K=coefficient_map(value);assert K.extract([1,2],[0,2]).det()==-2*B*C
    assert coefficient_map(value.subs(B,0)).rank()==1
    canonical=b.subs(dict(zip(params,[0,0,1,1,0,0])))
    r,q,w=s.symbols('r q w');v=s.Matrix([r,1,q,w]);ell=(v.T*J).T
    T=s.Matrix.hstack(*[(v.T*m).T for m in form_matrices(canonical)]);u=list(T);tt=list(ell*ell.T)
    equations=list(coefficient_map(canonical)*v)+[s.expand(u[i]*tt[j]-u[j]*tt[i]) for i in range(16) for j in range(i)]
    gb=s.groebner(equations,r,q,w,domain=s.QQ);assert list(gb)==[r,q,w*w]
    # e_w is the only possible offchart line and fails contraction.
    off=s.Matrix([0,0,0,1]);elloff=(off.T*J).T
    co=list(s.Matrix.hstack(*[(off.T*m).T for m in form_matrices(canonical)]));tar=list(elloff*elloff.T)
    assert any(co[i]*tar[j]!=co[j]*tar[i] for i in range(16) for j in range(i))
    paths=[Path(__file__),ROOT/'code/ambient_weight1_image.py',ROOT/'code/weight1_structure.py',
        ROOT/'reports/AMBIENT_BASE_FLAG_CRITICAL_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json']
    body=dict(id='AMBIENT_RANK2_INTERSECTION_UNIQUE_20261103_v63',all_passed=True,
        symplectic_shear=[[str(x) for x in row] for row in shear.tolist()],torus=[[str(x) for x in row] for row in torus.tolist()],
        coefficient_transforms={'C':'C/(a*b^2)','B':'B*b/a^2'},normalization_choice={'a':'C/b^2','b^5':'C^2/B'},
        critical_rank2_normal_form='t=0,A=0,B=C=1,D=F=0',canonical_full_projective_flag_fiber=list(map(str,gb)),
        whole_rank_at_least2_flag_support='one geometric point; may be a double scheme point on critical rank2 orbit',
        conclusion='On the open rank(K_beta)>=2, (CA intersection CB)_red equals N_red. In particular N is the unique reduced intersection component through distinguished Q(a=1). No scheme-theoretic intersection equality.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('RANK2 REDUCED INTERSECTION PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
