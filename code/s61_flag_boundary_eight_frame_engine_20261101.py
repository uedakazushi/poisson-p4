#!/usr/bin/env python3
"""Exact rational eight-variable flag calculation: frame stage.

Internal engine of s61_flag_boundary_eight_cubic_20261101.py.
All temporary products are rebuilt in the supplied working directory.
"""
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
WD=Path(sys.argv[1]); WD.mkdir(parents=True,exist_ok=True)
import sys,json,time
from pathlib import Path
sys.path.insert(0,str(ROOT/'code'))
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
from s61_flag_transverse_obstruction_20261101 import specialize
start=time.time();m=S61();case=json.load(open(ROOT / 'reports' / 'S61_FLAG_UNIFORM_OPEN_20261101.json'))['flag'];H=s.Matrix(case['beta_kernel']).applyfunc(s.Rational)
b=s.symbols('b:15');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
z,bu,beta=primitive_beta(m);vv=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS]);bs=vv.jacobian(bu);br=list(bs.T.to_DM().rref()[1]);left=bs[br,:].inv();act=s.Matrix.hstack(*[(left*infinitesimal_action(vv,A)[br,:]).applyfunc(s.expand) for A in sp_basis()]);HC=H[:,[j for j in range(15) if j!=1]];NL=s.Matrix.hstack(*HC.T.nullspace()).T;N=(NL*act).subs(dict(zip(bu,H*s.Matrix(b))),simultaneous=True).subs(b[1],0).applyfunc(s.expand)
print('boundary normal action',N,flush=True);print('b12zero rank',len(N.subs(b[12],0).to_DM().rref()[1]),flush=True)
sub={b[0]:1,b[1]:0,b[3]:0,b[5]:0,b[6]:0,b[12]:1,b[14]:0};bet=(H*s.Matrix(b)).subs(sub)
_,F=specialize(m,dict(zip(one,bet)),True);q=s.Matrix([s.Rational(case['point61'][i]) for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];J=s.zeros(135,45)
for row,f in enumerate(F):
 for mon,c in f.items():
  if len(mon)==2:
   a,bb=mon
   if a!=fixed:J[row,cols.index(a)]+=c*q[bb]
   if bb!=fixed:J[row,cols.index(bb)]+=c*q[a]
  elif len(mon)==1:J[row,44]+=c*q[mon[0]]
rr,C=J.to_DM().rref();rr=rr.to_Matrix();free=[j for j in range(45) if j not in C];T=s.zeros(45,len(free))
for j,k in enumerate(free):
 T[k,j]=1
 for i,c in enumerate(C):T[c,j]=-rr[i,k]
R=list(J[:,C].T.to_DM().rref()[1]);dm=J[R,C].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));den=s.factor(s.lcm([s.denom(v) for v in T]));print('rank',len(C),'det',det,'Tden',den,'seconds',time.time()-start,flush=True)
enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
Path(WD / 'flag-b1-normalized-frame.json').write_text(json.dumps(dict(substitutions={str(k):str(v) for k,v in sub.items()},variables=[str(v) for v in b if v not in sub],jacobian=enc(J),tangent=enc(T),pivot_rows=R,pivot_columns=list(C),free=free,minor=str(det),tangent_denominator=str(den),boundary_normal_action=enc(N)),indent=2)+'\n')
