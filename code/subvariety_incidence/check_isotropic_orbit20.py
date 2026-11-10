#!/usr/bin/env python3
"""Exact PGL5 orbit dimension at the totally isotropic Log Poisson point."""
import argparse,itertools,json
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path('work/isotropic_orbit20.json'));args=ap.parse_args()
assert s.__version__=='1.14.0'
x=s.symbols('u0 u1 z w0 w1');u0,u1,z,w0,w1=x;E=s.Matrix(x)
P=s.zeros(5)
for (i,j),f in {(0,3):2*u0*z,(0,4):-2*u0*u1,(1,3):2*u1*z,(1,4):-2*u1**2,(2,3):-4*z**2,(2,4):4*u1*z,(3,4):2*(u0**2+u1*w0-2*z*w1)}.items():
 P[i,j]=f;P[j,i]=-f
pairs=list(itertools.combinations(range(5),2));mons=list(itertools.combinations_with_replacement(range(5),2))
def coeffs(B):return s.Matrix([s.Poly(B[i,j],*x).coeff_monomial(x[a]*x[b])for i,j in pairs for a,b in mons])
def action(M):
 v=M*E
 return s.Matrix([[s.expand(sum(v[k]*s.diff(P[i,j],x[k])-M[i,k]*P[k,j]-M[j,k]*P[i,k] for k in range(5)))for j in range(5)]for i in range(5)])
ms=[]
for i in range(5):
 for j in range(5):
  M=s.zeros(5);M[i,j]=1;ms.append(M)
D=s.Matrix.hstack(*[coeffs(action(M))for M in ms])
rank=DomainMatrix.from_Matrix(D).convert_to(s.QQ).rank()
extended=DomainMatrix.from_Matrix(D.row_join(coeffs(P))).convert_to(s.QQ).rank()
assert rank==extended==21
data={'date':'2026-11-10','version':'v13','status':'PASS','sympy':s.__version__,'homogeneous_bivector_coefficients':150,'GL5_generators':25,'GL5_affine_orbit_rank':21,'with_P0_scale_rank':21,'PGL5_projective_orbit_dimension':20,'intersection_dimension_bounds':[20,24],'claim':'The projective orbit is contained in the intersection of rank-2 logarithmic dimension-28 and rank-4 long dimension-25 components. The largest-dimensional logarithmic component is regular at P0 by e=1 and Nagata; its smoothness is not an output of this finite calculation.'}
args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
print('PASS: projective orbit dimension 20')
