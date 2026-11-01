"""Exact residual initial forms in a two-dimensional tangent space."""
from itertools import combinations_with_replacement
import sympy as s
from s61_flag_transverse_obstruction_20261101 import specialize

def initial_data(model,local,quartic=False):
 point=list(map(s.Rational,local['point61']));beta={i:point[i] for i,w in enumerate(model.reduced['parameter_weights']) if w==1};_,F=specialize(model,beta,True);co=local['chart_columns'];fixed=next(i for i in range(45) if i not in co);J=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);assert T.shape==(45,2);rows=list(J.T.to_DM().rref()[1]);cols=list(J.to_DM().rref()[1]);U=s.zeros(45,43);iv=J[rows,cols].inv(method='DM')
 for j,i in enumerate(cols):U[i,:]=iv[j,:]
 JU=J*U;terms=[]
 for ff in F:
  ts=[]
  for mm,c in ff.items():
   if fixed in mm:continue
   i,j=[co.index(k) for k in mm]+[44]*(2-len(mm));ts.append((i,j,c))
  terms.append(ts)
 mons={d:list(combinations_with_replacement(range(2),d)) for d in [1,2,3,4]};ix={d:{v:i for i,v in enumerate(mons[d])} for d in mons};series={1:T};residuals={};v=s.symbols('v0:2')
 for degree in range(2,5 if quartic else 4):
  raw=s.zeros(135,len(mons[degree]))
  for row,ts in enumerate(terms):
   for i,j,c in ts:
    for a in range(1,degree):
     for k,ma in enumerate(mons[a]):
      for l,mb in enumerate(mons[degree-a]):raw[row,ix[degree][tuple(sorted(ma+mb))]]+=c*series[a][i,k]*series[degree-a][j,l]
  residual=raw-JU*raw[rows,:];series[degree]=-U*raw[rows,:];residuals[degree]=residual
  if degree<3 or quartic and degree==3:assert residual==s.zeros(135,len(mons[degree]))
 final=residuals[4 if quartic else 3];dm,piv=final.to_DM().rref();basis=dm.to_Matrix()[:len(piv),:];degree=4 if quartic else 3;forms=[s.expand(sum(c*s.prod(v[j] for j in mon) for c,mon in zip(row,mons[degree]))) for row in basis.tolist()]
 if quartic:
  assert forms==[v[0]*v[1]**3,v[1]**4];g0=s.groebner([f.subs(v[0],0) for f in forms],v[1]);g1=s.groebner([f.subs(v[0],1) for f in forms],v[1]);assert list(g0)==[v[1]**4] and list(g1)==[v[1]**3]
 else:assert forms==[v[1]**3]
 first=next((i,j) for i in range(final.rows) for j in range(final.cols) if final[i,j])
 return dict(initial_degree=degree,jacobian_pivot_rows=rows,jacobian_pivot_columns=cols,initial_forms=list(map(str,forms)),first_nonzero=dict(row=first[0],column=first[1],value=str(final[first])),residuals={str(d):[[str(c) for c in row] for row in A.tolist()] for d,A in residuals.items()},transverse_length_upper=4 if quartic else 3,generic_multiplicity_upper=3,flatness_asserted=not quartic)
