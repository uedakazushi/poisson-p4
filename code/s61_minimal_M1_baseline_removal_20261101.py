#!/usr/bin/env python3
"""Uniform open and exceptional-stratum removal for both M1 free baselines."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_rankone_junction_20261101 import enc
from s61_minimal_baseline_thickenings_20261101 import tensor,coords,normalization
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_local_initial_bound_20261101 import cubic_bound
from weight1_structure import sp_basis,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def uniform(m,source_path):
 st=time.time();kind='M1';bs=s.symbols('b0:10');src=json.loads(source_path.read_text());r=src['cases'][kind]['universal'];one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];K=s.Matrix(r['beta_kernel']).applyfunc(s.sympify);two,F=specialize(m,dict(zip(one,K*s.Matrix(bs))),True);A=s.Matrix(r['leading_jacobian']).applyfunc(s.sympify);tc=s.Matrix(r['t_column']).applyfunc(s.sympify);C=s.Matrix(r['first_kernel']).applyfunc(s.sympify);pr=r['first_pivot_rows'];pc=r['first_pivot_columns'];U=s.zeros(44,len(pc));iv=A.extract(pr,pc).inv()
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 N=s.Matrix(r['residual_kernel']).applyfunc(s.sympify);T=(C*N[:-1,:]-U*tc[pr,:]*N[-1,:]).col_join(N[-1,:]);full_frame=T.copy();assert (A.row_join(tc).to_DM()*full_frame.to_DM()).is_zero_matrix;assert full_frame[-1,:]==s.zeros(1,3);T=T[:,1:].applyfunc(s.cancel);den=[]
 for j in range(T.cols):
  d=s.factor(s.lcm([s.denom(T[i,j]) for i in range(T.rows)]));den.append(d)
  for i in range(T.rows):T[i,j]=s.cancel(T[i,j]*d)
 print('polynomial normal frame ready',flush=True);E=s.Matrix(r['residual_matrix']).applyfunc(s.sympify).row_join(s.Matrix(r['residual_t_column']).applyfunc(s.sympify));sample=dict(zip(bs,range(1,11)));num=E.subs(sample);piv=list(num.to_DM().rref()[1]);rows=list(num[:,piv].T.to_DM().rref()[1]);cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);norm=list(combinations_with_replacement(range(T.cols),2));H=s.zeros(135,len(norm));tf={i:T[cols.index(i),:] for i in cols};tf[fixed]=s.zeros(1,T.cols)
 for row,f in enumerate(F):
  for mon,c in f.items():
   forms=[tf[i] for i in mon]+[T[-1,:]]*(2-len(mon));a,b=forms
   for j,(u,v) in enumerate(norm):H[row,j]+=c*(a[u]*b[v]+(a[v]*b[u] if u!=v else 0))
 H=H.applyfunc(s.expand);Y=(H-A*U*H[pr,:]).applyfunc(s.expand);print('polynomial quadratic coefficients ready',flush=True);Em=E.extract(rows,piv);numY=Y.subs(sample);numres=numY-num[:,piv]*num.extract(rows,piv).inv()*numY[rows,:];chosen=list(numres.T.to_DM().rref()[1]);assert len(chosen)==3;assert chosen==[67,68,134];Cr=(E[chosen,piv]*Em.inv(method='DM')).applyfunc(s.cancel)
 def cancel_with_factors(poly,denominator):
  pn=s.Poly(poly,*bs);co,fac=s.factor_list(denominator);left=s.sympify(co)
  for f,e in fac:
   fp=s.Poly(f,*bs);unused=0
   for j in range(e):
    q,r=pn.div(fp)
    if r.is_zero:pn=q
    else:unused=e-j;break
   left*=f**unused
  return pn.as_expr()/left
 out=s.zeros(3,3)
 for i,row in enumerate(chosen):
  dd=s.factor(s.lcm([s.denom(v) for v in Cr[i,:]]));coeff=[s.cancel(v*dd) for v in Cr[i,:]]
  for j,(a,b) in enumerate(norm):
   numerator=s.expand(dd*Y[row,j]-sum(c*Y[r,j] for c,r in zip(coeff,rows)));out[i,j]=cancel_with_factors(numerator,dd*den[a]*den[b])

 det=s.factor(out.det(method="domain-ge"));P2=192*bs[3]*bs[9]+9*bs[5]*bs[7]-36*bs[6]*bs[7]-8*bs[8]*bs[9];P3=684*bs[3]**2*bs[9]+837*bs[3]*bs[5]*bs[7]-3348*bs[3]*bs[6]*bs[7]-176*bs[3]*bs[8]*bs[9]+63*bs[5]*bs[7]*bs[8]-252*bs[6]*bs[7]*bs[8]-28*bs[8]**2*bs[9];W2=52*bs[3]*bs[9]+9*bs[5]*bs[7]-36*bs[6]*bs[7]-4*bs[8]*bs[9];W3=89528*bs[3]**2*bs[9]+17415*bs[3]*bs[5]*bs[7]-69660*bs[3]*bs[6]*bs[7]+9424*bs[3]*bs[8]*bs[9]+1935*bs[5]*bs[7]*bs[8]-7740*bs[6]*bs[7]*bs[8]+248*bs[8]**2*bs[9];assert s.cancel(det+20*bs[7]*bs[9]**2*(9*bs[3]+bs[8])*W2**3*W3/(3*P2*P3**3))==0;dm=E.extract(rows,piv).to_DM();jd=s.factor(dm.domain.to_sympy(dm.det()));assert jd==-s.Rational(788937956456832,244140625)*bs[7]**11*(9*bs[3]+bs[8])*(bs[5]-4*bs[6])*P2*P3
 return dict(normal_denominators=list(map(str,den)),normal_frame_polynomial=enc(T),quadratic_rows=chosen,quadratic_minor_matrix=enc(out),quadratic_minor=str(det),jacobian_minor_rows=rows,jacobian_minor_columns=piv,jacobian_minor=str(jd),full_kernel_frame=enc(full_frame)),K,bs,dict(main=9*bs[3]+bs[8],P1=bs[5]-4*bs[6],P2=P2,P3=P3,W2=W2,W3=W3,b7=bs[7],b9=bs[9])
def geometry(m,K,b,factors):
 two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];u,act=beta_actions(m);Q=tensor('M1',0,0,1);q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q) if v);V=s.Matrix(X[1:]);qd=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fi]/q[fi];qd.append(dq)
 orbit=s.Matrix.hstack(*qd);stab=orbit.nullspace();assert orbit.rank()==7 and len(stab)==3;rows=list(K.T.rref()[1]);iv=K[rows,:].inv();actions=[]
 for c in stab:
  db=(act*c).subs(dict(zip(u,K*s.Matrix(b))));v=iv*db[rows,:];assert (K*v-db).applyfunc(s.expand)==s.zeros(16,1);actions.append(v)
 preservation={}
 for label,f in factors.items():
  df=[s.expand(sum(s.diff(f,bb)*a[j] for j,bb in enumerate(b))) for a in actions];assert all(s.rem(s.Poly(v,*b),s.Poly(f,*b)).is_zero for v in df);preservation[label]=list(map(str,df))
 return dict(orbit_dimension=7,stabilizer_dimension=3,stabilizer_coefficients=[list(map(str,c)) for c in stab],stabilizer_actions=[list(map(str,c)) for c in actions],factor_preservation=preservation),u,act

def restrictions(m,K,b,factors,u,act,curve):
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];a=s.Symbol('a');out=[]
 for label in ['b7','b9','W2','P1','P2','P3','W3']:
  f=factors[label];drop=b[7] if label=='b7' else b[9] if label=='b9' else b[5];assert s.degree(f,drop)==1;coef=s.diff(f,drop);ct=f.subs(drop,0);assert s.gcd(coef,ct)==1 or coef in (1,-1);sol=s.solve(f,drop)[0];pars=[bb for bb in b if bb!=drop];beta=(K*s.Matrix(b)).subs(drop,sol);D=beta.jacobian(pars).row_join(act.subs(dict(zip(u,beta))));rr,piv=D.T.to_DM().rref();rank=len(piv);entry=dict(label=label,factor=str(f),solved_variable=str(drop),solution=str(sol),primitive_linear_irreducible=True,generic_beta_rank=rank)
  if label in ['b7','b9','W2']:
   assert rank<=15;rr=rr.to_Matrix();free=[i for i in range(16) if i not in piv];N=s.zeros(16,len(free))
   for j,i in enumerate(free):
    N[i,j]=1
    for k,c in enumerate(piv):N[c,j]=-rr[k,i]
   assert (D.T.to_DM()*N.to_DM()).is_zero_matrix;entry.update(beta_differential=enc(D),left_kernel=enc(N),beta_image_dimension=rank)
  else:
   assert rank==16;sample=dict(zip(b,range(1,11)));sample[drop]=sol.subs({bb:v for bb,v in sample.items() if bb!=drop});bv=K*s.Matrix([sample[bb] for bb in b]);pt=s.Matrix(coords(m,tensor('M1',0,0,1)))
   for i,v in zip(one,bv):pt[i]=v
   local=calculate(m,True,pt);assert local['homogenized_jacobian_rank']==42 and local['leading_jacobian_rank']==41;q0=pt[two,:];cols=local['chart_columns'];fixed=two.index(local['chart_fixed_u']-1);dq=curve.diff(a);dq-=q0*dq[fixed]/q0[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);tr=list(T.T.to_DM().rref()[1]);vt=T[tr,:].inv()*tv[tr,:];assert T*vt==tv and vt[0]!=0;branches=[dict(kind='free_baseline',lower=3,tangent=list(map(str,vt)))];dom=minor_witness(D.subs(sample));assert dom['rank']==16
   if label in ['P1','W3']:upper=cubic_bound(m,local,branches,3)
   else:
    assert local['quadratic_span_rank']==3 and local['common_hessian_rank']==2;rad=s.Matrix(local['quadratic_common_radical']).applyfunc(s.Rational);assert s.Matrix.hstack(vt,rad).rank()==1;upper=dict(transverse_quadratic_length=3)
   entry.update(sample={str(k):str(v) for k,v in sample.items()},local_data=local,actual_curve_tangent=list(map(str,vt)),beta_minor=dom,initial_bound=upper,generic_rank=3,t_nilpotent=True,whole_generic_beta_stratum_removed=True)
  out.append(entry);print(label,'checked',rank,flush=True)
 return out

def run():
 st=time.time();m=S61();source_path=ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json';unif,K,b,factors=uniform(m,source_path);geo,u,act=geometry(m,K,b,factors);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];a=s.Symbol('a');curve=s.Matrix(coords(m,tensor('M1',a,0,1)))[two,:];_,F=specialize(m,dict(zip(one,K*s.Matrix(b))),True);assert all(s.expand(sum(c*s.prod(curve[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F);q0=curve.subs(a,0);cols=json.loads(source_path.read_text())['cases']['M1']['universal']['chart_columns'];fixed=next(i for i in range(45) if i not in cols);dq=curve.diff(a);dq-=q0*dq[fixed]/q0[fixed];T=s.Matrix(unif['full_kernel_frame']).applyfunc(s.sympify);assert (T[:,0]-dq[cols,:].col_join(s.zeros(1,1))).applyfunc(s.cancel)==s.zeros(45,1);restricted=restrictions(m,K,b,factors,u,act,curve)
 main_path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';main=json.loads(main_path.read_text())['cases']['M1'];A,B,P,R,U,V,KK,L=s.symbols('A B P R U V K L');bv=K*s.Matrix(b);sub={P:bv[1],R:bv[3],A:bv[4],B:bv[13],U:bv[6],V:bv[9],KK:bv[14]+5*bv[9],L:5*bv[6]-bv[15]};ma=s.factor(s.sympify(main['solution_a']).subs(sub,simultaneous=True));mb=s.factor(s.sympify(main['solution_b']).subs(sub,simultaneous=True));assert ma==(9*b[3]+b[8])/(27*b[7]) and mb==0;paired=normalization('M1')
 paths=[Path(__file__),source_path,main_path,ROOT/'code/s61_minimal_local_initial_bound_20261101.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json']
 return dict(id='S61_MINIMAL_M1_BASELINE_REMOVAL_20261101',all_passed=True,field='QQ and QQ(b0,...,b9)',uniform=unif,geometry=geo,actual_curve_q=list(map(str,curve)),factor_polynomials={k:str(v) for k,v in factors.items()},factor_cases=restricted,main_intersection=dict(a=str(ma),b=str(mb),known_k='-4/9',dependency='S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101'),paired_normalization=paired,conclusion='Both finite M1 free baseline families are removed by t-saturation over sufficiently general beta.',limitations=['Their projective parameter-infinity R0 flag endpoints are not removed here.','The M0 free baseline families are separate.','No full affine finiteness or five-point exhaustion.'],proof_note='research/notes/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.md',seconds=time.time()-st,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
