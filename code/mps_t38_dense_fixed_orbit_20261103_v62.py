#!/usr/bin/env python3
"""Exact dense fixed projective orbit in the accepted T38 component."""
import hashlib,itertools,json
from pathlib import Path
import sympy as s
from s61_model import X,coefficients,vertical_rref
from mps_s61_dominance_20261103_v62 import rank_minor
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
    source=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json';P=s.zeros(5)
    for pair,expression in json.loads(source.read_text())['homogeneous_bivector'].items():
        i,j=map(int,pair.split(','));P[i,j]=s.sympify(expression);P[j,i]=-P[i,j]
    variables=sorted(set().union(*(v.free_symbols for v in P))-set(X),key=str);sample={v:i+2 for i,v in enumerate(variables)};P=P.subs(sample)
    columns=[]
    for i,j in itertools.product(range(5),repeat=2):
        g=s.zeros(5);g[i,j]=1;d=g*P+P*g.T-sum(((g*s.Matrix(X))[k]*P.diff(X[k]) for k in range(5)),s.zeros(5));columns.append(s.Matrix(coefficients(d)))
    columns.append(s.Matrix(coefficients(P)));matrix=s.Matrix.hstack(*columns);rref,pivots=vertical_rref();free=[i for i in range(150) if i not in pivots]
    quotient=matrix[free,:]-rref[:,free].T*matrix[list(pivots),:];minor=rank_minor(quotient);assert minor['rank']==25
    body=dict(id='MPS_T38_DENSE_FIXED_ORBIT_20261103_v62',all_passed=True,projective_orbit_dimension=24,affine_cone_dimension=25,
        parameter_values={str(k):str(v) for k,v in sample.items()},quotient_matrix_QQ=[[str(v) for v in row] for row in quotient.tolist()],
        nonzero_minor=minor,component_identification_dependency='Accepted general T38 source family; known projective component dimension24 gives matching upper.',
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),source,ROOT/'code/s61_model.py',ROOT/'code/mps_s61_dominance_20261103_v62.py']})
    output=ROOT/'reports/MPS_T38_DENSE_FIXED_ORBIT_20261103_v62.json';output.write_text(json.dumps(body,indent=2)+'\n');print('T38 FIXED DENSE ORBIT PASS',minor['determinant_mod_32003'])
