#!/usr/bin/env python3
"""Universal A1 length lower bound and saturation at the rank-one M0 junction."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def tensor(a,b,d,h):
 x,y,z,w=X[1:];Q=s.zeros(4);Q[0,1]=x*x;Q[1,2]=-b*x*x-2*d*x*w-3*h*w*w;Q[1,3]=a*x*x;return Q-Q.T
def coords(m,Q):
 P=X[0]**2*J;P[1:,1:]+=Q;return m.normalize(P)['parameters']
def universal(m):
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];Q=tensor(0,0,1,1);BK=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert BK.shape==(16,9);b=s.symbols('b0:9');two,F=specialize(m,dict(zip(one,BK*s.Matrix(b))),True);p=coords(m,Q);q=s.Matrix([p[i] for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];D=s.zeros(135,45)
 for row,f in enumerate(F):
  val=0
  for mon,c in f.items():
   if len(mon)==2:
    i,j=mon;val+=c*q[i]*q[j];D[row,i]+=c*q[j];D[row,j]+=c*q[i]
  assert s.expand(val)==0
 D=D[:,cols];D0=D[:70,:];C=s.Matrix.hstack(*D0.nullspace());R=D*C;rr,piv=R.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(R.cols) if i not in piv];N=s.zeros(R.cols,len(free))
 for j,i in enumerate(free):
  N[i,j]=1
  for r,pc in enumerate(piv):N[pc,j]=-rr[r,i]
 assert (R.to_DM()*N.to_DM()).is_zero_matrix;assert len(D0.to_DM().rref()[1])+len(piv)==42 and len(free)==2
 return dict(beta_kernel=enc(BK),q_point=list(map(str,q)),chart_columns=cols,leading_jacobian=enc(D),first_kernel=enc(C),residual_matrix=enc(R),residual_kernel=enc(N),leading_rank=42)
def normalization():
 a,b,d,h=s.symbols('a b d h');V=s.Matrix(X[1:]);A=-a;B=b-2*d*a+3*h*a*a;S=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]]);Q=tensor(a,b,d,h);assert (S*J[1:,1:]*S.T-J[1:,1:]).applyfunc(s.expand)==s.zeros(4);assert (S.inv()*Q.subs(dict(zip(V,S*V)),simultaneous=True)*S.inv().T-tensor(0,0,d-3*h*a,h)).applyfunc(s.expand)==s.zeros(4)
 la,ro=s.symbols('la ro',nonzero=True);D=s.diag(ro,1/ro,1/la,la);Q=tensor(0,0,d,h);trans=(D.inv()*Q.subs(dict(zip(V,D*V)),simultaneous=True)*D.inv().T/ro**2).applyfunc(s.cancel);assert (trans-tensor(0,0,d*la**2,h*la**3/ro)).applyfunc(s.cancel)==s.zeros(4)
 return dict(shear=enc(S),transformed_d='d-3*h*a',torus=enc(D),projective_rescale='1/ro^2',transformed_d_torus='d*la^2',transformed_h_torus='h*la^3/ro',open_condition='h*(d-3*h*a) != 0',canonical_choice='la^2=1/(d-3*h*a), ro=h*la^3')
def cubic_bound(m,local,branches):
 p=list(map(s.Rational,local['point61']));beta={i:p[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};_,F=specialize(m,beta,True);cols=local['chart_columns'];fixed=next(i for i in range(45) if i not in cols);Jc=s.Matrix(local['jacobian']).applyfunc(s.Rational);T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);n=T.cols;rows=list(Jc.T.to_DM().rref()[1]);piv=list(Jc.to_DM().rref()[1]);U=s.zeros(45,len(piv));iv=Jc.extract(rows,piv).inv(method='DM')
 for i,c in enumerate(piv):U[c,:]=iv[i,:]
 JU=Jc*U;mon2=list(combinations_with_replacement(range(n),2));mon3=list(combinations_with_replacement(range(n),3));ix={a:i for i,a in enumerate(mon3)};Q=s.zeros(135,len(mon2));terms=[]
 for row,f in enumerate(F):
  ts=[]
  for mon,c in f.items():
   if fixed in mon:continue
   i,j=[cols.index(k) for k in mon]+[44]*(2-len(mon));ts.append((i,j,c))
   for h,(a,b) in enumerate(mon2):Q[row,h]+=c*(T[i,a]*T[j,b]+(T[i,b]*T[j,a] if a!=b else 0))
  terms.append(ts)
 w2=-U*Q[rows,:];A2=Q-JU*Q[rows,:];Q3=s.zeros(135,len(mon3))
 for row,ts in enumerate(terms):
  for i,j,c in ts:
   for a in range(n):
    for h,mon in enumerate(mon2):Q3[row,ix[tuple(sorted((a,)+mon))]]+=c*(T[i,a]*w2[j,h]+T[j,a]*w2[i,h])
 A3=Q3-JU*Q3[rows,:];ar=list(A2.T.to_DM().rref()[1]);ac=list(A2.to_DM().rref()[1]);C=A2[:,ac]*A2.extract(ar,ac).inv(method='DM');B=A3-C*A3[ar,:];br=list(B.T.to_DM().rref()[1]);v=s.symbols(f'v0:{n}');forms=[s.expand(sum(c*s.prod(v[i] for i in mon) for c,mon in zip(row,mons))) for M,mons in [(A2[ar,:],mon2),(B[br,:],mon3)] for row in M.tolist()];sec=s.groebner([f.subs(v[0],0) for f in forms],*v[1:],order='grevlex');assert sec.is_zero_dimensional;gb=s.groebner([f.subs(v[0],1) for f in forms],*v[1:],order='grevlex');assert gb.is_zero_dimensional;leads=[p.LM(order=gb.order).exponents for p in gb.polys];counts=[]
 for deg in range(20):
  count=0
  for inds in combinations_with_replacement(range(n-1),deg):
   e=tuple(inds.count(i) for i in range(n-1))
   if not any(all(a<=b for a,b in zip(l,e)) for l in leads):count+=1
  counts.append(count)
  if count==0:break
 assert sum(counts)==5;assert all(s.Rational(b['tangent'][0])!=0 for b in branches)
 return dict(jacobian_pivot_rows=rows,jacobian_pivot_columns=piv,quadratic_rank=len(ar),pure_cubic_rank=len(br),initial_forms=list(map(str,forms)),section_groebner=list(map(str,sec)),localized_groebner=list(map(str,gb)),standard_degrees=counts,parameter='v0',multiplicity_upper_bound=5)
def run():
 start=time.time();m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';source=json.loads(path.read_text())['cases']['M0'];local=source['local_data'];u=universal(m);normal=normalization();t=s.Symbol('t');Q=tensor(-2,s.Rational(2,9)+20*t,s.Rational(1,3)-s.Rational(15,2)*t,t);qr=s.Matrix(coords(m,Q))[two,:];p=s.Matrix(local['point61']).applyfunc(s.Rational);q=p[two,:];assert qr.subs(t,0)==q;_,F=specialize(m,dict(zip(one,list(p[one,:]))),True);assert all(s.expand(sum(c*s.prod(qr[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F);dcan=s.Rational(1,3)-s.Rational(3,2)*t;assert s.expand((s.Rational(1,3)-s.Rational(15,2)*t)-3*t*(-2)-dcan)==0
 T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);cols=local['chart_columns'];fixed=two.index(local['chart_fixed_u']-1);dq=qr.diff(t);dq-=q*dq[fixed]/q[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));rr=list(T.T.to_DM().rref()[1]);vt=T[rr,:].inv()*tv[rr,:];assert T*vt==tv;branches=[dict(kind='main',lower=3,tangent=source['main_branch_tangent']),dict(kind='rankone_A1',lower=2,tangent=list(map(str,vt)))];assert s.Matrix.hstack(*[s.Matrix(b['tangent']).applyfunc(s.Rational) for b in branches]).rank()==2;bound=cubic_bound(m,local,branches);assert sum(b['lower'] for b in branches)==bound['multiplicity_upper_bound'];assert source['rankone_beta_minor']['rank']==16
 paths=[Path(__file__),path,ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
 return dict(id='S61_MINIMAL_RANKONE_JUNCTION_20261101',all_passed=True,field='QQ and QQ(b0,...,b8)',A1_universal=u,A1_normalization=normal,A1_q_curve=list(map(str,qr)),A1_canonical_d=str(dcan),branches=branches,local_data=local,cone_bound=bound,generic_rank=5,t_nilpotent=True,generic_beta_rankone_M0_junction_stratum_removed=True,scope='The rank-one M0 main/A1 junction and its16D stratum over general beta; does not remove the entire17D A1 family or claim a free completed ring.',proof_note='research/notes/S61_MINIMAL_RANKONE_JUNCTION_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
