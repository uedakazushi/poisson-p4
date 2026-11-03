#!/usr/bin/env python3
"""All-kappa constant Jacobian minors plus the missing kappa0 obstruction.

The paired-flag/properness and hypersurface-valuation proof is in the note.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector
from ambient_rank_two_flag_intersection_table_20261103_v63 import lie

ROOT=Path(__file__).resolve().parents[1]

def constant_minor(M,kappa):
    atzero=M.subs(kappa,0);cols=list(atzero.rref()[1]);rows=list(atzero[:,cols].T.rref()[1]);det=s.factor(M[rows,cols].det(method='domain-ge'))
    assert det!=0 and not det.free_symbols
    return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(det),symbolically_constant=True)

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;kappa=s.Symbol('kappa')
    geo_path=ROOT/'reports/AMBIENT_R19_GENERIC_GEOMETRY_20261103_v63/report.json';geo=json.loads(geo_path.read_text())
    dbl_path=ROOT/'reports/AMBIENT_R19_GENERIC_DOUBLE_20261103_v63/report.json';dbl=json.loads(dbl_path.read_text())
    A=s.Matrix([list(map(s.sympify,row)) for row in geo['full_CA_source_differential']]);B=s.Matrix([list(map(s.sympify,row)) for row in geo['full_CB_source_differential']]);W=A.row_join(B)
    point=s.Matrix(list(map(s.sympify,geo['point_canonical150'])));orbit=s.Matrix.hstack(*(lie(point,i,j) for i,j in itertools.product(range(5),repeat=2)))
    dirs=[vector(canonical(raw(-4*x*x,0,2*x*z*w)))]+[vector(canonical(raw(0,0,H))) for H in [w*x*x,x*x*z,x*z*z/2,2*h*x*z,z**3/3,h*z*z,3*h*h*z]]
    R=s.Matrix.hstack(*dirs).row_join(orbit)
    matrices={'CA':A,'CB':B,'tangent_union':W,'paired_flag_R_source':R};minors={name:constant_minor(M,kappa) for name,M in matrices.items()}
    assert [v['rank'] for v in minors.values()]==[29,29,37,20]
    assert [M.rank() for M in matrices.values()]==[29,29,37,20]
    difference=s.Matrix(list(map(s.sympify,dbl['actual_second_order_difference_canonical150'])));extra=s.Matrix(list(map(s.sympify,dbl['common_first_order_canonical150'])))
    zero_covector=s.zeros(150,1);zero_covector[68]=s.Rational(5,22);zero_covector[72]=1
    assert zero_covector.T*W.subs(kappa,0)==s.zeros(1,W.cols)
    assert (zero_covector.T*difference.subs(kappa,0))[0]==-18
    assert W.subs(kappa,0).row_join(difference).rank()==38
    assert R.subs(kappa,0).row_join(extra.subs(kappa,0)).rank()==21
    paths=[Path(__file__),geo_path,dbl_path,ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',
        ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py',
        ROOT/'research/notes/AMBIENT_R19_GENERIC_DOUBLE_20261103_v63.md']
    report=dict(id='AMBIENT_R19_OPEN_DOUBLE_20261103_v63',all_passed=True,field='QQ[kappa]',
        all_kappa_constant_minors=minors,kappa_zero_tangent_sum_annihilator_canonical150=list(map(str,zero_covector)),
        kappa_zero_actual_second_order_obstruction='-18',kappa_zero_augmented_tangent_union_rank=38,
        kappa_zero_extra_direction_augmented_R_rank=21,
        obstruction_cover=['kappa!=0: kappa*e72-e63 pairs-18*kappa','kappa=0: (5/22)e68+e72 pairs-18'],
        whole_normalform_family_consequence='At EVERY P_kappa and itsPGLtranslates, CA/CB/R are smooth and actualintersection completion=C[[u1,...,u19,eps]]/(eps²).',
        source_open='t*q5*q7!=0 inR8, plusitsPGLsaturation',
        intersection_depth=19,intersection_local_associated_primes=['Rprime only'],
        scope='No claim for allAF0 or t0 boundary points and no globalprimarydecomposition.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print('WHOLE AF OPEN ORDINARY DOUBLE PASS',report['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
