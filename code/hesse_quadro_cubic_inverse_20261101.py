#!/usr/bin/env python3
"""Universal cubic inverse of the Hesse elliptic-quintic quadratic map.

All identities are checked in QQ(lambda)[x0,...,x4]. This certifies the
birational map itself, not any primitive-first-jet equivariance assertion.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]


def hesse_quadrics(variables,parameter):
    return s.Matrix([variables[i]**2+parameter*variables[(i-1)%5]*variables[(i+1)%5]
                     -variables[(i-2)%5]*variables[(i+2)%5]/parameter for i in range(5)])


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    x=s.symbols('x0:5');y=s.symbols('y0:5');lam=s.Symbol('lam',nonzero=True)
    q=hesse_quadrics(x,lam);t=[0,-lam,1,-1,lam]
    phi=s.Matrix(5,5,lambda i,j:t[(j-i)%5]*x[3*(i+j)%5])
    assert phi+phi.T==s.zeros(5)
    assert (phi*q).applyfunc(s.expand)==s.zeros(5,1)
    a=(phi*s.Matrix(y)).jacobian(x)
    assert a*s.Matrix(x)==phi*s.Matrix(y)
    assert (s.Matrix(y).T*a).applyfunc(s.expand)==s.zeros(1,5)
    # Negative adj(A) column0, divided by its known factor y0.
    inverse=s.Matrix([s.cancel(-(-1)**i*a.minor_submatrix(0,i).det()/y[0]) for i in range(5)])
    assert all(s.Poly(f,*y).total_degree()==3 for f in inverse)
    assert (a*inverse).applyfunc(s.expand)==s.zeros(5,1)
    e0={y[i]:int(i==0) for i in range(5)}
    witness=s.factor(a.minor_submatrix(0,0).subs(e0).det())
    assert witness==-lam**2
    # Coordinate0 divisibility gives the degree5 projective factors.
    composed=inverse.subs(dict(zip(y,q)),simultaneous=True).applyfunc(s.expand)
    factor_x=s.cancel(composed[0]/x[0]);assert s.Poly(factor_x,*x).total_degree()==5
    assert (composed-factor_x*s.Matrix(x)).applyfunc(s.expand)==s.zeros(5,1)
    reverse=hesse_quadrics(inverse,lam).applyfunc(s.expand)
    factor_y=s.cancel(reverse[0]/y[0]);assert s.Poly(factor_y,*y).total_degree()==5
    assert (reverse-factor_y*s.Matrix(y)).applyfunc(s.expand)==s.zeros(5,1)
    assert s.cancel(factor_x.subs({x[i]:int(i==0) for i in range(5)}))==lam**2
    assert s.cancel(factor_y.subs(e0))==lam**4
    report={'id':'HESSE_QUADRO_CUBIC_INVERSE_20261101','field':'QQ(lam)','all_passed':True,
        'quadratic_map':list(map(str,q)),'syzygy_matrix':[[str(v) for v in row] for row in phi.tolist()],
        'inverse_kernel_matrix':[[str(v) for v in row] for row in a.tolist()],
        'cubic_inverse':list(map(str,inverse)),
        'composition_factor_x':str(factor_x),'composition_factor_y':str(factor_y),
        'rank4_minor_at_e0':str(witness),
        'identities':['Phi(x)*Q(x)=0','Phi(x)*y=A(y)*x','y^T*A(y)=0','A(y)*R(y)=0','R(Q(x))=factor_x*x','Q(R(y))=factor_y*y'],
        'open_sets':'lambda !=0; source factor_x !=0 and target factor_y !=0; both nonempty since factors at e0 are lambda²,lambda⁴.',
        'scope':'The quadratic map is birational with the displayed cubic inverse. No assertion that it preserves normalized primitive jets is included.',
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as f:json.dump(report,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'rank4_minor':str(witness),'output':str(args.output)}))

if __name__=='__main__':main()
