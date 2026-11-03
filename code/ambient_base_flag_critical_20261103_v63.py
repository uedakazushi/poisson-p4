#!/usr/bin/env python3
"""Exact full base-flag tangent critical support; not merely rank K_beta.

Proof that this controls normality of reduced CA/CB/N is in the paired note.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_weight1_image import common_plane
from weight1_structure import J,form_matrices,coefficient_map

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    p,b=common_plane();t,A,B,C,D,F=p;r,q,w=s.symbols('r q w');v=s.Matrix([r,1,q,w])
    matrices=form_matrices(b);ell=(v.T*J).T
    T=s.Matrix.hstack(*[(v.T*M).T for M in matrices]);left=list(T);right=list(ell*ell.T)
    equations=list(coefficient_map(b)*v)+[s.expand(left[i]*right[j]-left[j]*right[i]) for i in range(16) for j in range(i)]
    derivative=s.Matrix(equations).jacobian([r,q,w]).subs({r:0,q:0,w:0})
    rows=[]
    for row in derivative.tolist():
        if row==[0,0,0] or row in rows or [-x for x in row] in rows:continue
        rows.append(row)
    matrix=s.Matrix(rows);minors=[]
    for indices in itertools.combinations(range(matrix.rows),3):
        f=s.factor(matrix[list(indices),:].det())
        if f!=0 and f not in minors and -f not in minors:minors.append(f)
    gb=s.groebner(minors,*p,domain=s.QQ);Delta=4*C*F-D*D
    assert gb.reduce(t**3)[1]==0
    central=[s.expand(f.subs(t,0)) for f in minors]
    gcentral=s.groebner(central,A,B,C,D,F,domain=s.QQ)
    gtarget=s.groebner([c*Delta for c in [A,B,C,D,F]],A,B,C,D,F,domain=s.QQ)
    assert list(gcentral)==list(gtarget)
    # These two exact ideal computations imply radical(t,Delta). Conversely
    # at t=Delta=0 the displayed 11x3 matrix has rank at most2.
    assert all(s.rem(s.Poly(f.subs(t,0),A,B,C,D,F),s.Poly(Delta,A,B,C,D,F)).as_expr()==0 for f in minors)
    # The potentially misleading coefficient rank3 condition is narrower.
    K=coefficient_map(b);coefficient_minor=s.factor(K.extract([0,1,2],[0,2,3]).det())
    assert s.expand(coefficient_minor-t*(Delta+t*t/4))==0
    witness_t=s.Matrix(rows).subs(dict(zip(p,[1,0,0,0,0,0])))
    witness_Delta=s.Matrix(rows).subs(dict(zip(p,[0,3,5,7,11,13])))
    assert witness_t.rank()==3 and witness_Delta.rank()==3
    # At t=0,Delta!=0 the full projective uniqueness proof only needs the
    # invertibility of this2x2 matrix, not a coefficient rank3 claim.
    M=s.Matrix([[2*C,D],[D,2*F]]);assert M.det()==Delta
    paths=[Path(__file__),ROOT/'code/ambient_weight1_image.py',ROOT/'code/weight1_structure.py',
        ROOT/'research/notes/AMBIENT_FULL_FIRSTJET_FIBERS_20261103_v63.md',
        ROOT/'reports/AMBIENT_RANKDROP_FLAG_FIBERS_V2_20261103_v63/report.json']
    body=dict(id='AMBIENT_BASE_FLAG_CRITICAL_20261103_v63',all_passed=True,field='QQ',
        base_parameters=list(map(str,p)),full_flag_chart_equations=list(map(str,equations)),
        flag_tangent_matrix=[[str(x) for x in row] for row in rows],
        maximal_minors=list(map(str,minors)),maximal_minors_groebner_basis=list(map(str,gb)),
        central_minors_groebner_basis=list(map(str,gcentral)),critical_support_ideal=['t','4*C*F-D**2'],
        coefficient_rank3_minor=str(coefficient_minor),
        conclusion='Full base flag is one reduced point whenever t!=0 OR 4CF-D²!=0. Exceptional support in each six-plane has codimension2. Proper smooth component resolutions are isomorphisms outside its proper image.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('BASE FLAG CRITICAL SUPPORT PASS',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
