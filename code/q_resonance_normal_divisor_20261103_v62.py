#!/usr/bin/env python3
"""Certify normality of the actual divisor of the nonzero resonance family."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    source=ROOT/'reports/Q_RESONANCE_INTEGRABLE_DIRECTION_20261103_v62.json'
    data=json.loads(source.read_text());assert data['all_passed']
    X=s.symbols('x0:5');eps=s.Symbol('eps');lam=s.Symbol('lam',nonzero=True)
    F=s.Poly(s.sympify(data['degeneracy_polynomial']).subs(eps,1),*X)
    denominator=s.ilcm(*(c.q for c in F.coeffs()));integerF=F.mul_ground(denominator)
    # Every displayed component has exactly the required torus character.
    weights=(-2,-1,0,1,2)
    for key,text in data['integrable_R'].items():
        i,j=map(int,key.split(','));polynomial=s.Poly(s.sympify(text),*X)
        assert all(weights[i]+weights[j]-sum(a*b for a,b in zip(mon,weights))==-1
                   for mon,c in polynomial.terms())
    assert weights[0]+weights[1]-2*weights[4]==-7
    program='R=QQ[x0,x1,x2,x3,x4];\nF='+str(integerF.as_expr()).replace('**','^')+';\n'
    program+='''J=ideal(diff(vars R,F));G=gb J;
print("JAC_DIM="|toString dim(R/J));
print("JAC_DEG="|toString degree(R/J));
print("HS="|toString hilbertSeries(R/J,Reduce=>true));
print("GB_GENERATORS="|toString numColumns gens G);
assert(dim(R/J)==2);assert(degree(R/J)==40);
print "RESONANCE ACTUAL DIVISOR NORMAL PASS";exit 0;
'''
    script=args.output/'normal.m2';script.write_text(program)
    log=args.output/'normal.log';started=time.monotonic()
    with log.open('w')as stream:
        result=subprocess.run([shutil.which('M2'),'--script',str(script.resolve())],
                              stdout=stream,stderr=subprocess.STDOUT,timeout=90)
    assert result.returncode==0 and 'RESONANCE ACTUAL DIVISOR NORMAL PASS' in log.read_text()
    files=[Path(__file__),source,script,log]
    report=dict(id='Q_RESONANCE_NORMAL_DIVISOR_20261103_V62',all_passed=True,
                field='QQ, dimension preserved by extension to C',Jacobian_affine_dimension=2,
                Jacobian_projective_dimension=1,Jacobian_degree=40,
                actual_divisor_normal_for_every_nonzero_parameter=True,
                projective_equivalence='g_lam P_eps=lam^5 P_(eps*lam^-6); choose lam^6=eps',
                normality_proof='The hypersurface is S2. Jacobian affine dimension2 excludes repeated factors and intersections of distinct geometric factors, and gives singular codimension>=2. Serre R1+S2 therefore proves actual geometrically integral normal D.',
                zero_parameter_actual_divisor='5H',
                q51_identification_dependency='Would follow from bundled normal-actual-divisor theorem; no component-identification claim is added by this computation.',
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in files[:2]},
                evidence_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()for p in files[2:]},
                seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS actual normal divisor for eps nonzero',report['seconds'])


if __name__=='__main__':
    main()
