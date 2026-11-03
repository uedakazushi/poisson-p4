#!/usr/bin/env python3
"""Full135 exact fixed-beta curve through the old double junction."""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_marked_sigma_point_lift_20261103_v62 import load_model,coerce

ROOT=Path(__file__).resolve().parents[1]


def verify(folder,output):
    start=time.monotonic();L=load_model(folder.resolve())
    K=s.QQ.algebraic_field(s.sqrt(-2),alias='h');h=K.unit
    z=s.Symbol('z');F=K.frac_field(z);hf=F.convert(h,K);zf=F.from_sympy(z)
    normal=[F.zero]*13;normal[9]=hf*hf;normal[10]=hf;normal[12]=zf
    p=L.lift(hf,normal,domain=F)
    assert p['all135_zero'] and not p['t']
    assert all(not F.to_sympy(value).has(z) for value in p['beta'])
    def expression(value):return str(L.expression(value,F,s.Symbol('h')))
    expressions=[expression(x) for x in p['beta']+p['q']]
    # The formulas are polynomial over the common p-local quadratic integer ring.
    for text in expressions:
        value=s.cancel(s.sympify(text));_,den=s.fraction(value)
        assert not den.has(z,s.Symbol('h')) and int(den)%32003
    special=L.lift(h,[K.zero]*9+[h*h,h,K.zero,coerce(K,s.Rational(2,3))*h],domain=K)
    assert all(F.to_sympy(a).subs(z,s.Rational(2,3)*s.sqrt(-2))==K.to_sympy(b)
               for a,b in zip(p['beta']+p['q'],special['beta']+special['q']))
    report=dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_JUNCTION_FAMILY_20261103_v62',all_passed=True,
                family='u=h,v=z,rho=all_lights=0,h^2=-2',
                constant_beta=True,all135_identically_zero=True,t_identically_zero=True,
                discriminant='2*h^4+9*z^2=8+9*z^2',
                double_specialization='z=2*h/3',
                coefficient_denominators_prime_units=True,
                beta=expressions[:16],q=expressions[16:],
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in [Path(__file__),ROOT/'code/s61_marked_sigma_point_lift_20261103_v62.py',
                                        L.input_path,ROOT/'legacy/v02/data/symplectic_slice_reduced_equations.json',
                                        ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json']},
                seconds=time.monotonic()-start)
    output.write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model-folder',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    a=parser.parse_args();r=verify(a.model_folder,a.output);print('PASS',r['id'],r['seconds'])
