#!/usr/bin/env python3
"""Independent rational-chart leaf tests; no old leaf normalizer is imported.

Verifies the three full magnetic normal forms on the simple cubic components,
and derives degrees nine/six by intersecting leaves with two general linear
equations. Proper Hilbert extension to whole T37 is a separate written proof.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as s

t,z,w,I,J,lam = s.symbols('t z w I J lambda')


def bracket(f,g,P):
    variables = [t,z,w]
    return s.cancel(sum(s.diff(f,x)*P[i,j]*s.diff(g,y)
                        for i,x in enumerate(variables) for j,y in enumerate(variables)))


def main(output):
    f = t*(t*t-1)
    H6 = t**6/2+3*t**5/5+3*t**4/4+t**3+3*t*t/2+3*t
    models = [
        ('node', f,-t*f,-(t*t-1)*w+I+J*f,z+(w+I)/t,J/t+I/(t*t-1)),
        ('cusp_flex',f,t*f,-(t*t-1)*w+I+J*t**3,w/t-z+I/t,(J*t+I)/(t*t-1)),
        ('cusp_secant',1-t,-3*t*(1-t),3*(1-t)*z+3*lam+3*t**6,
         w+3*t*z-H6,3*(lam+1)/(1-t))]
    rows = []
    for name,a,b,c,u,coefficient in models:
        P = s.Matrix([[0,a,b],[-a,0,c],[-b,-c,0]])
        for g in [t,z,w]:
            assert s.cancel(bracket(u,g,P)+coefficient*bracket(t,g,P)) == 0, (name,g)
        poles = [0,1,-1] if name == 'node' else ([1,-1] if name == 'cusp_flex' else [1])
        residues = [s.residue(coefficient,t,p) for p in poles]
        rows.append({'source':name,'kernel_form':'d('+str(u)+') + ('+str(coefficient)+') dt',
                     'residues':dict(zip(map(str,poles),map(str,residues))),
                     'all_three_Hamiltonian_contractions_zero':True})
    # General plane equations solve x,z linearly in 1,y,w on the quartic.
    y,x,k = s.symbols('y x k')
    a0,a1,a2,b0,b1,b2 = s.symbols('a0:3 b0:3')
    wy = 2*y*y-y**4/4
    xx = a0+a1*y+a2*wy
    zz = b0+b1*y+b2*wy
    ga9 = s.Poly(s.expand(y*xx+(4-y*y)*zz-16*y**5/5+4*y**7/7-y**9/36-k),y)
    assert ga9.degree() == 9 and ga9.LC() == -s.Rational(1,36)
    # For TRIPLE, r=-x^3-3xy+3y/2 and two plane equations give
    # y=-(x^3+b1*x+b0)/(3*x+b2-3/2). The denominator is generically coprime.
    yy = -(x**3+b1*x+b0)/(3*x+b2-s.Rational(3,2))
    tt = a0+a1*x+a2*yy
    triple = s.Poly(s.cancel(tt+x*yy/3-x**4/12+yy*yy/2-k).as_numer_denom()[0],x)
    assert triple.degree() == 6
    # Scalar normalization of numerator is immaterial; a simple exact plane
    # proves the generic resultant with its denominator is not identically zero.
    denom = 3*x+b2-s.Rational(3,2)
    resultant = s.factor(s.resultant(triple.as_expr(),denom,x))
    assert resultant != 0
    result = {'scope':__doc__, 'cubic_models':rows,
        'Ga9_general_plane_leaf_degree':ga9.degree(),
        'TRIPLE_general_plane_leaf_degree':triple.degree(),
        'TRIPLE_plane_denominator_resultant':str(resultant),
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
