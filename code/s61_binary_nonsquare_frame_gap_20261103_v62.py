#!/usr/bin/env python3
"""Extend the common-frame-gap proper-image certificate to all eight nuisance coordinates."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def verify():
    started=time.time()
    path=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json'
    previous=json.loads(path.read_text())
    assert previous['all_passed']
    ids=[4,5,12,13,14,25,26,27,34,35,36,37,48,49,50,51]
    nuisance=[5,25,26,27,34,37,48,51]
    u=s.symbols('u1:62')
    e=s.Symbol('e')
    nu=s.symbols('n0:8')
    S=3*e-2
    aa=-3*S*(e+1)/50
    dd=-20/S
    beta=dict.fromkeys(ids,s.S.Zero)
    beta.update({4:aa,14:(e+1)/5,35:1,36:dd,49:e,50:-5})
    beta.update(dict(zip(nuisance,nu)))
    vector=s.Matrix([beta[i] for i in ids])
    actions=s.Matrix(previous['primitive_action']).applyfunc(s.sympify)
    actions=actions.subs(dict(zip([u[i-1] for i in ids],vector)),simultaneous=True)
    source=vector.jacobian([e,*nu])
    matrix=actions.row_join(source)
    nrows=[ids.index(i) for i in nuisance]
    assert matrix[:,11:]==s.eye(16)[:,nrows]
    keep=[i for i in range(16) if i not in nrows]
    scale=50*S**2
    reduced=(matrix[keep,:11]*scale).applyfunc(lambda x:s.cancel(x).expand())
    kernel=reduced.T.to_DM().nullspace().to_Matrix()
    assert kernel.rows==1 and kernel.cols==8
    assert (kernel*reduced).applyfunc(s.expand)==s.zeros(1,11)
    sample={e:2,**dict(zip(nu,range(1,9)))}
    evaluated=matrix.subs(sample)
    columns=list(evaluated.to_DM().rref()[1])
    rows=list(evaluated[:,columns].T.to_DM().rref()[1])
    assert len(rows)==len(columns)==15
    determinant=evaluated[rows,columns].det(method='domain-ge')
    assert determinant
    assert s.expand(50*aa+9*e**2+3*e-6)==0
    B=(50*aa+18*e**2-9*e-2)*dd+60*e-40
    assert s.cancel(B)==0
    files=[Path(__file__),path,ROOT/'research/notes/S61_BINARY_HYPERPLANE_FRAME_COVER_20261102.md',
           ROOT/'code/s61_binary_hyperplane_frame_cover_20261102.py']
    return dict(id='S61_BINARY_NONSQUARE_FRAME_GAP_20261103_v62',all_passed=True,field='QQ',
                curve=dict(a=str(aa),e='e',d=str(dd)),nuisance_indices_one_based=nuisance,
                source_dimension=9,full_source_differential_shape=[16,19],
                nuisance_identity_rows=nrows,reduced_shape=[8,11],denominator_clearing=str(scale),
                polynomial_reduced_matrix=[[str(v) for v in row] for row in reduced.tolist()],
                exact_left_kernel=[[str(v) for v in row] for row in kernel.tolist()],
                beta_image_dimension=15,
                rank15_minor=dict(rows=rows,columns=columns,point={str(k):str(v) for k,v in sample.items()},value=str(determinant)),
                consequence='For the full H11 source, the common old/new frame gap B=D0=0 on aFC!=0 has proper beta image.',
                limitations=['No other rank-defect curve is excluded.',
                             'No local t-nilpotence is asserted at the special gap points.'],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                seconds=time.time()-started)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=verify()
    with args.output.open('x') as output:
        json.dump(report,output,indent=2)
        output.write('\n')
    print(json.dumps({'all_passed':True,'all_eight_nuisance':True,'beta_image_dimension':15,'seconds':report['seconds']}))
