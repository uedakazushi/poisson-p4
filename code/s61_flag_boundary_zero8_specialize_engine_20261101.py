#!/usr/bin/env python3
"""Specialize the original polynomial equations and frame data to b8=0."""
import pickle
import sys
from pathlib import Path
import sympy as s
wd=Path(sys.argv[1]);path=wd/'flag-b1-universal-quadratic-cache.pkl'
data=pickle.loads(path.read_bytes());index=data['variables'].index('b8')
for key in ['J','T','U','A2','W2']:
    for row in data[key]:
        for j,poly in enumerate(row):
            row[j]=[[ex,coef] for ex,coef in poly if ex[index]==0]
for row in data['F']:
    for pair in row:
        pair[1]=str(s.sympify(pair[1]).subs(s.Symbol('b8'),0))
for key in ['ud','td']:
    assert not s.sympify(data[key]).has(s.Symbol('b8'))
path.write_bytes(pickle.dumps(data))
print('ORIGINAL POLYNOMIAL FRAME AND EQUATIONS SPECIALIZED AT B8=0',flush=True)
