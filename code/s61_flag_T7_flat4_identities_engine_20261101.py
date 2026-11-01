#!/usr/bin/env python3
"""T7 flag flat4 certificate: identities stage. Rebuild in supplied directory."""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]);WD.mkdir(parents=True,exist_ok=True)
import json,time,sys
from pathlib import Path
import sympy as s
from sympy.polys.rings import ring
sys.path.insert(0,str(ROOT/'code'))
from s61_model import S61,X,J
from s61_flag_transverse_obstruction_20261101 import specialize
st=time.time();d=json.load(open(WD/'flag-T7-universal-transport.json'));b=list(map(s.Symbol,d['variables']));e=s.Symbol('e');coeff=[s.Matrix(d['q']).applyfunc(s.sympify)]+[s.Matrix(d['coefficients'][str(i)]).applyfunc(s.sympify)for i in [1,2,3]]
R,*_=ring(b+[e],s.QQ)
den={i:s.lcm([s.denom(v)for v in coeff[i]])for i in range(4)};P={i:[R.from_expr(s.cancel(v*den[i]))for v in coeff[i]]for i in range(4)}
for i in range(4):
 p=s.Poly(den[i],*b,e);assert len(p.terms())==1
 assert all(k==0 for j,k in enumerate(p.monoms()[0])if j not in [8,9,11])
print('Laurent coefficient denominators',den,time.time()-st,flush=True)
H=s.Matrix(json.load(open(ROOT/'reports/S61_FLAG_T7_BOUNDARY_20261101.json'))['fixed_beta_kernel']).applyfunc(s.Rational);m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];_,F=specialize(m,dict(zip(one,H*s.Matrix(b))),True)
for degree in range(4):
 D=s.lcm([den[k]*den[degree-k]for k in range(degree+1)])
 for f in F:
  v=R.zero
  for mon,c in f.items():
   if len(mon)!=2:continue
   i,j=mon;pc=R.from_expr(c)
   for k in range(degree+1):
    mult=R.from_expr(s.cancel(D/(den[k]*den[degree-k])))
    v+=pc*P[k][i]*P[degree-k][j]*mult
  assert not v
 print('all135 coefficient degree',degree,'PASS',time.time()-st,flush=True)
# The second reduced parameter is supplied by the H11-preserving shear z -> z+u*w.
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
z,us,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(z[k])for(i,j),k in COORDS]);VB=vec.jacobian(us);rows=list(VB.T.rref()[1]);VBi=VB[rows,:].inv();D=VBi*infinitesimal_action(vec,sp_basis()[9])[rows,:];DB=D.jacobian(us);L=s.Matrix.hstack(*H.T.nullspace()).T;assert L*DB*H==s.zeros(5,11);assert DB**4==s.zeros(16)
assert sp_basis()[9]==s.Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]])
u=s.Symbol('u');x,y,z,w=X[1:];G=s.eye(4)+u*sp_basis()[9];Q=s.zeros(4);Q[1,3]=-x*x;Q-=Q.T;sp=s.Matrix([x,y,z,w]);transformed=G.inv()*Q.subs(dict(zip(sp,G*sp)),simultaneous=True)*G.inv().T
Pfull=X[0]**2*J;Pfull[1:,1:]+=transformed
qaction=s.Matrix(m.normalize(Pfull)['parameters'])[[i for i,wt in enumerate(m.reduced['parameter_weights'])if wt==2],:]
assert qaction[19]==u and qaction[25]==-1 and qaction[24]==0
assert all(coeff[d][i]==0 for d in [1,2,3] for i in [19,24])
# A separate regular lattice will define the endpoint algebra; these pure jets
# alone have epsilon poles and are not a closed flat extension.
out=dict(all_passed=True,all135_four_coefficients_zero=True,coefficient_denominators={str(k):str(v)for k,v in den.items()},beta_shear_matrix=[[str(v)for v in row]for row in DB.tolist()],shear_preserves_H11=True,homogenizer_zero=True,seconds=time.time()-st)
(WD/'flag-T7-universal-transport-identities.json').write_text(json.dumps(out,indent=2)+'\n');print('UNIVERSAL LAURENT JET IDENTITIES PASS',time.time()-st,flush=True)
