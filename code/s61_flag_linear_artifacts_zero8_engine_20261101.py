#!/usr/bin/env python3
"""Exact nonzero minors on each h=0 six-parameter source."""
import json
import sys
from pathlib import Path
import sympy as s
wd=Path(sys.argv[1]);data=json.loads((wd/'flag-eight-initial4-field.json').read_text())
values=[0,1,2,0,3,5,7,11];M=s.zeros(217,213)
for group,offset in [('rows',0),('corrected_t4',216)]:
    for key,(numerator,denominator) in data[group].items():
        assert all(not any(ex) for ex,co in denominator)
        assert all(ex[0]==0 for ex,co in numerator)
        i,j=map(int,key.split(','))
        def evaluate(poly):
            return sum(s.Rational(co)*s.prod(s.Integer(value)**e for value,e in zip(values,ex)) for ex,co in poly)
        M[i+offset,j]=evaluate(numerator)/evaluate(denominator)
checks=[]
for rows,expected in [(216,12),(217,13)]:
    A=M[:rows,:];_,columns=A.to_DM().rref();_,indices=A.T.to_DM().rref()
    assert len(columns)==len(indices)==expected
    determinant=A[list(indices),list(columns)].det()
    assert determinant
    checks.append(dict(rank=expected,rows=list(indices),columns=list(columns),determinant=str(determinant)))
result=dict(all_passed=True,source_parameters=['b4','b7','b9','b10','b11','b13'],values=[1,2,3,5,7,11],all_candidate_denominators_constant=True,checks=checks,source_dimension_with_Sp4=16,proper_closed_failure_beta_image_bound=15)
(wd/'flag-linear-zero8-result.json').write_text(json.dumps(result,indent=2)+'\n')
print('H=0 EXACT QQ MINOR WITNESS PASS',flush=True)
