#!/usr/bin/env python3
"""Remove the finite M0 rank-one A1 family over sufficiently general beta."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_minimal_rankone_junction_20261101 import tensor,coords,enc,normalization
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def uniform(m):
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];Q=tensor(0,0,1,1);BK=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());bs=s.symbols('b0:9');two,F=specialize(m,dict(zip(one,BK*s.Matrix(bs))),True);p=coords(m,Q);q=s.Matrix([p[i] for i in two]);fi=next(i for i,v in enumerate(q) if v);co=[i for i in range(45) if i!=fi];D=s.zeros(135,45);tc=s.zeros(135,1)
 for row,f in enumerate(F):
  for mon,c in f.items():
   if len(mon)==2:
    i,j=mon;D[row,i]+=c*q[j];D[row,j]+=c*q[i]
   elif len(mon)==1:tc[row]+=c*q[mon[0]]
 D=D[:,co];D0=D[:70,:];pc=list(D0.rref()[1]);pr=list(D0[:,pc].T.rref()[1]);C=s.Matrix.hstack(*D0.nullspace());U=s.zeros(44,len(pc));iv=D0.extract(pr,pc).inv()
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 R=(D*C).applyfunc(s.factor);ct=(tc-D*U*tc[pr,:]).applyfunc(s.factor);E=R.row_join(ct);rr,piv=E.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(E.cols) if i not in piv];N=s.zeros(E.cols,len(free))
 for j,i in enumerate(free):
  N[i,j]=1
  for k,c in enumerate(piv):N[c,j]=-rr[k,i]
 assert (E.to_DM()*N.to_DM()).is_zero_matrix;assert len(pc)+len(piv)==43
 T=(C*N[:-1,:]-U*tc[pr,:]*N[-1,:]).col_join(N[-1,:]);full=D.row_join(tc);assert (full.to_DM()*T.to_DM()).is_zero_matrix;assert N[-1,:]==s.zeros(1,2);H=s.zeros(135,3);tf={i:T[co.index(i),:] for i in co};tf[fi]=s.zeros(1,2)
 for row,f in enumerate(F):
  for mon,c in f.items():
   fs=[tf[i] for i in mon]+[T[-1,:]]*(2-len(mon));a,b=fs;H[row,0]+=c*a[0]*b[0];H[row,1]+=c*(a[0]*b[1]+a[1]*b[0]);H[row,2]+=c*a[1]*b[1]
 H=H.applyfunc(s.cancel);Y=(H-D*U*H[pr,:]).applyfunc(s.cancel);rows=list(E[:,list(piv)].subs(dict(zip(bs,range(1,10)))).T.to_DM().rref()[1]);Crel=E[:,list(piv)]*E.extract(rows,list(piv)).inv(method='DM');res=(Y-Crel*Y[rows,:]).applyfunc(s.factor);rrh,hp=res.T.to_DM().rref();nz=list(hp);assert len(nz)==1 and res[nz,:]==s.Matrix([[0,0,-bs[7]/12]])
 # A polynomial rank minor from first pivots + residual pivots; top block constant.
 dm=E.extract(rows,list(piv)).to_DM();minor=s.factor(dm.domain.to_sympy(dm.det()));assert minor==s.Rational(1162261467,800000000)*bs[3]*bs[6]**10*bs[7]**2*(2*bs[6]-bs[7])

 assert all(s.expand(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
 return dict(beta_kernel=enc(BK),kernel=enc(T),residual=enc(E),pivot_rows=rows,pivot_columns=list(piv),jacobian_minor=str(minor),quadratic_rows=enc(res[nz,:]),full_residual=enc(res),first_constant_rank=len(pc),leading_rank=42,full_rank=43),BK
def beta_actions(m):
 z,u,B=primitive_beta(m);v=s.Matrix([B[i,j].diff(z[h]) for (i,j),h in COORDS]);M=v.jacobian(u);rr=list(M.T.rref()[1]);iv=M[rr,:].inv();return u,s.Matrix.hstack(*[iv*infinitesimal_action(v,H)[rr,:] for H in sp_basis()])
def bad_images(m,BK):
 b=s.symbols('b0:9');u,act=beta_actions(m);out=[]
 for label,sub,drop in [('b3', {b[3]:0},3),('b6',{b[6]:0},6),('b7',{b[7]:0},7),('2b6-b7',{b[7]:2*b[6]},7)]:
  params=[a for i,a in enumerate(b) if i!=drop];beta=(BK*s.Matrix(b)).subs(sub);D=beta.jacobian(params).row_join(act.subs(dict(zip(u,beta))));rr,piv=D.T.to_DM().rref();rr=rr.to_Matrix();free=[i for i in range(16) if i not in piv];N=s.zeros(16,len(free))
  for j,i in enumerate(free):
   N[i,j]=1
   for k,c in enumerate(piv):N[c,j]=-rr[k,i]
  assert len(piv)<=15 and (D.T.to_DM()*N.to_DM()).is_zero_matrix;out.append(dict(factor=label,beta_differential=enc(D),left_kernel=enc(N),image_dimension=len(piv)))
 return out
def same_curve(m,BK,base_d):
 b=s.symbols('b0:9');h=s.Symbol('h');d=base_d+b[7]*(h-1)/(2*b[6]);a=0;c=2*b[7]**2*(h-1)/(9*b[6]**2);Q=tensor(a,c,d,h);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,dict(zip(one,BK*s.Matrix(b))),True);q=s.Matrix(coords(m,Q))[two,:];assert all(s.cancel(sum(v*s.prod(q[i] for i in mon) for mon,v in f.items() if len(mon)==2))==0 for f in F);assert d.subs(h,1)==base_d
 return dict(q=list(map(str,q)),d=str(d),h='h',a='0',shape_b=str(c),all135_identities=True),q

def pure_orbit(m):
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];Q=tensor(0,0,0,1);K=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert K.shape==(16,9);curve,qr=same_curve(m,K,0);b=s.symbols('b0:9');h=s.Symbol('h');sample=dict(zip(b,range(1,10)));point=s.Matrix(coords(m,Q));beta=K*s.Matrix(range(1,10))
 for i,v in zip(one,beta):point[i]=v
 local=calculate(m,True,point);assert local['leading_jacobian_rank']==42 and local['homogenized_jacobian_rank']==43 and local['quadratic_span_rank']==1 and local['common_hessian_rank']==1;T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rad=s.Matrix(local['quadratic_common_radical']).applyfunc(s.Rational);cols=local['chart_columns'];fixed=two.index(local['chart_fixed_u']-1);q0=point[two,:];dq=qr.diff(h).subs(sample).subs(h,1);dq-=q0*dq[fixed]/q0[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));rr=list(T.T.to_DM().rref()[1]);vt=T[rr,:].inv()*tv[rr,:];assert T*vt==tv and s.Matrix.hstack(vt,rad).rank()==1
 u,act=beta_actions(m);BD=K.row_join(act.subs(dict(zip(u,beta))));dom=minor_witness(BD);assert dom['rank']==16;vv=s.Matrix(X[1:]);qd=[]
 for H in sp_basis():
  dQ=sum(((H*vv)[i]*Q.diff(vv[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q0*dq[fixed]/q0[fixed];qd.append(dq[cols,:])
 orbit=s.Matrix.hstack(*qd);assert orbit.rank()==7
 return dict(beta_kernel=enc(K),curve=curve,local_data=local,actual_curve_tangent=list(map(str,vt)),projective_orbit_differential=enc(orbit),projective_orbit_dimension=7,beta_differential=enc(BD),beta_minor=dom,closure_dimension=16,completed_local_ring='C[[u,e]]/(e^2)',t_zero=True,generic_beta_whole_orbit_stratum_removed=True)
def run():
 st=time.time();m=S61();univ,K=uniform(m);print('uniform verified',flush=True);bad=bad_images(m,K);curve,qr=same_curve(m,K,1);pure=pure_orbit(m);norm=normalization();paths=[Path(__file__),ROOT/'code/s61_minimal_rankone_junction_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_RANKONE_JUNCTION_20261101.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json']
 return dict(id='S61_MINIMAL_RANKONE_A1_REMOVAL_20261101',all_passed=True,field='QQ and QQ(b0,...,b8)',uniform=univ,bad_beta_images=bad,canonical_actual_curve=curve,normalization=norm,pure_cubic_orbit=pure,conclusion='The entire finite M0 A1 family at rank-one vertex quadrics is removed by t-saturation over sufficiently general beta.',limitations=['The projective parameter-infinity endpoint is an R0 binary locus and is not removed by this certificate.','Uses the separately proved main/A1 junction removal for h=0.','No full affine-fiber exhaustion conclusion.'],proof_note='research/notes/S61_MINIMAL_RANKONE_A1_REMOVAL_20261101.md',seconds=time.time()-st,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
