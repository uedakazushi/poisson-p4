#!/usr/bin/env python3
"""Two exact minimal tails at the R0 double-root endpoint."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def run():
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];x,y,z,w=X[1:];tau=s.Symbol('tau');beta=s.Matrix([0,0,0,0,-s.Rational(2,3),2,2,3,5,3,7,11,13,-s.Rational(1,3),-15,10]);_,F=specialize(m,dict(zip(one,beta)),True);out={};q0=None
 for kind in ['M0','M1']:
  G=z*z*w+tau*((-18*x*x*z-3*x*x*w+2*x*z*w) if kind=='M0' else (12*x*x*z+2*x*x*w-x*z*w));Q=s.zeros(4);Q[0,1]=tau*x*x if kind=='M0' else tau*x*x/2;Q[2,3]=0 if kind=='M0' else tau*x*x/2;Q[1,2]=-s.diff(G,w)-(tau*x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);Q-=Q.T;P=X[0]**2*J;P[1:,1:]+=Q;raw=m.normalize(P)['parameters'];q=s.Matrix([raw[i] for i in two]);assert all(s.expand(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
  if q0 is None:q0=q.subs(tau,0)
  assert q.subs(tau,0)==q0;fixed=next(i for i,v in enumerate(q0) if v);cols=[i for i in range(45) if i!=fixed];dq=q.diff(tau);dq-=q0*dq[fixed]/q0[fixed];out[kind]=dict(G=str(G),Q=[[str(v) for v in row] for row in Q.tolist()],q_curve=list(map(str,q)),q_tangent=list(map(str,dq)),chart_tangent=list(map(str,dq[cols,:].col_join(s.zeros(1,1)))),generic_length_lower=1)
 p=s.zeros(61,1)
 for i,v in zip(one,beta):p[i]=v
 for i,v in zip(two,q0):p[i]=v
 Jc=s.zeros(135,45);tc=s.zeros(135,1)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if len(mon)==2:
    i,j=mon;Jc[row,i]+=c*q0[j];Jc[row,j]+=c*q0[i]
   elif len(mon)==1:tc[row]+=c*q0[mon[0]]
 full=Jc[:,cols].row_join(tc);T=s.Matrix.hstack(*[s.Matrix(v['chart_tangent']).applyfunc(s.Rational) for v in out.values()]);assert T.rank()==2 and full*T==s.zeros(135,2);kernel=s.Matrix.hstack(*full.nullspace());rows=list(kernel.T.to_DM().rref()[1]);kt=kernel[rows,:].inv()*T[rows,:];assert kernel*kt==T
 for j,kind in enumerate(out):out[kind]['local_tangent_coordinates']=list(map(str,kt[:,j]))
 paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json']
 return dict(id='S61_MINIMAL_DOUBLE_ROOT_ENDPOINT_BRANCHES_20261101',all_passed=True,field='QQ',point61=list(map(str,p)),chart_fixed_q_index=fixed,chart_columns=cols,full_jacobian_rank=len(full.to_DM().rref()[1]),tangent_basis=[[str(v) for v in row] for row in kernel.tolist()],cases=out,tail_tangent_rank=2,scope='Two actual distinct minimal fixed-beta branches through one double-root R0 endpoint; no local multiplicity upper bound or saturation conclusion.',proof_note='research/notes/S61_MINIMAL_DOUBLE_ROOT_ENDPOINT_BRANCHES_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'],tangent_rank=r['tail_tangent_rank'])))
if __name__=='__main__':main()
