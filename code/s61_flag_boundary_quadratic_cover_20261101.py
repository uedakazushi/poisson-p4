#!/usr/bin/env python3
"""Cover quadratic rank9 on the eight-variable flag linear frame open."""
import argparse
import hashlib
from itertools import combinations,combinations_with_replacement
import json
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'reports/S61_FLAG_BOUNDARY_QUADRATIC_COVER_20261101.json')
    args=p.parse_args();start=time.time()
    engines=[f'code/s61_flag_boundary_eight_{stage}_engine_20261101.py' for stage in ['frame','inverse','quadratic']]
    with tempfile.TemporaryDirectory(prefix='flag-quadratic-cover-') as directory:
        wd=Path(directory)
        for engine in engines:subprocess.run([sys.executable,str(ROOT/engine),str(wd)],cwd=ROOT,check=True)
        c=pickle.loads((wd/'flag-b1-universal-quadratic-cache.pkl').read_bytes())
        q=json.loads((wd/'flag-b1-universal-quadratic.json').read_text())
    b=s.symbols('b2 b4 b7 b8 b9 b10 b11 b13');x=b[0]
    def dec(poly):return sum(s.Rational(v)*s.prod(z**e for z,e in zip(b,ex)) for ex,v in poly)
    active=[j for j,m in enumerate(q['monomials']) if all(i in {4,5,7,8} for i in m)]
    A=s.Matrix([[dec(c['A2'][i][j]) for j in active] for i in q['pivot_rows']])
    assert set().union(*(v.free_symbols for v in A))=={x}
    D=(x+6)*(5*x+12)*(20*x+119)*(25*x+4)
    contents=[]
    for i in range(9):
        g=s.factor(s.gcd_list(list(A[i,:])));contents.append(g)
        assert s.cancel(g/((x+6)**2*D)).is_Rational
        for j in range(10):A[i,j]=s.cancel(A[i,j]/g)
    minors=[]
    for j in range(10):
        dm=A[:,[k for k in range(10) if k!=j]].to_DM()
        minors.append(s.factor(dm.domain.to_sympy(dm.det())))
    g=s.factor(s.gcd_list(minors))
    assert s.cancel(g/(20*x+119)**13).is_Rational
    mu=[s.cancel((-1)**j*v/g) for j,v in enumerate(minors)]
    assert (A*s.Matrix(mu)).applyfunc(s.expand)==s.zeros(9,1)
    assert s.gcd_list(mu).is_Rational
    G=s.zeros(4)
    for v,(i,j) in zip(mu,combinations_with_replacement(range(4),2)):G[i,j]=G[j,i]=v
    two=[s.factor(G.extract(i,j).det()) for i in combinations(range(4),2) for j in combinations(range(4),2)]
    gg=s.factor(s.gcd_list(two));assert s.cancel(gg/(400*x+91)).is_Rational
    exceptional=-s.Rational(91,400);assert G.subs(x,exceptional).rank()==1
    inputs=engines+['code/s61_flag_boundary_quadratic_cover_20261101.py','code/s61_model.py',
        'code/s61_flag_transverse_obstruction_20261101.py','code/primitive_jet_elliptic_curve_20261101.py',
        'code/weight1_structure.py','reports/S61_FLAG_UNIFORM_OPEN_20261101.json',
        'reports/S61_FLAG_BOUNDARY_NORMAL_FORMS_20261101.json']
    report=dict(all_passed=True,seconds=time.time()-start,linear_open=str(s.expand(D)),
        universal_quadratic_rank_on_linear_open=9,selected_quadratic_rows=q['pivot_rows'],
        active_variables=[4,5,7,8],row_contents=list(map(str,contents)),
        primitive_maximal_minors=list(map(str,minors)),maximal_minor_gcd=str(g),
        socle_multiplication_coefficients=list(map(str,mu)),
        multiplication_two_minor_gcd=str(gg),rank_one_core=str(exceptional),
        rank_at_exceptional_core=1,
        limitations=['The four linear-frame core values are excluded here.',
            'The rank-one multiplication core requires separate cubic and higher analysis.',
            'No cubic all-point rank, formal ring or saturation assertion is made.'],
        input_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print('QUADRATIC COVER AND RANK-ONE EXCEPTION PASS',report['seconds'],flush=True)

if __name__=='__main__':main()
