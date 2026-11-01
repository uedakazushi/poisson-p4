#!/usr/bin/env python3
"""Exact full-rank local certificate at an additional M0 infinity point."""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
POINT=['1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '2', '2', '3', '4/3', '0', '0', '0', '0', '0', '5', '3', '7', '11', '-19/2', '0', '10', '16/3', '0', '0', '0', '0', '-20/3', '0', '13', '3', '-12', '5', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0']

def verify():
    model=S61();point=list(map(s.Rational,POINT))
    beta={i:point[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1}
    two,polynomials=specialize(model,beta,True)
    q0=s.Matrix([point[i] for i in two]);fixed=next(i for i,v in enumerate(q0) if v)
    columns=[i for i in range(45) if i!=fixed]
    derivative=s.zeros(135,45);tcolumn=s.zeros(135,1)
    for row,f in enumerate(polynomials):
        value=0
        for mon,c in f.items():
            if len(mon)==2:
                i,j=mon;value+=c*q0[i]*q0[j]
                derivative[row,i]+=c*q0[j];derivative[row,j]+=c*q0[i]
            elif len(mon)==1:tcolumn[row,0]+=c*q0[mon[0]]
        assert value==0
    jacobian=derivative[:,columns].row_join(tcolumn)
    _,rows=jacobian[:134,:].T.to_DM().rref();rows=list(rows)
    assert len(rows)==45
    minor=jacobian[rows,:].to_DM();determinant=minor.domain.to_sympy(minor.det())
    assert determinant!=0
    _,leading_pivots=derivative.to_DM().rref();assert len(leading_pivots)==44
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']
    return dict(id='S61_MINIMAL_AXIS_LOCAL_POINT_20261101',all_passed=True,field='QQ',point61=POINT,
        full135_leading_equations_zero=True,leading_jacobian_rank=44,
        homogenized_jacobian_rank=45,minor_rows_zero_based=rows,minor_determinant=str(determinant),
        chart_fixed_u=two[fixed]+1,chart_fixed_value=str(q0[fixed]),
        conclusion='For both 134 and 135 homogenizations the local ring is C, hence t=0 and saturation removes this point.',
        proof_note='research/notes/S61_MINIMAL_AXIS_LOCAL_POINT_20261101.md',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,rank=45,minor_determinant=r['minor_determinant'])))
if __name__=='__main__':main()
