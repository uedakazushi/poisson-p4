#!/usr/bin/env python3
"""Remove both finite M0 free baselines using universal and number-field bounds."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61,X
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_rankone_junction_20261101 import enc
from s61_minimal_baseline_thickenings_20261101 import tensor,coords,normalization
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_local_quartic_bound_20261101 import quartic_bound
from weight1_structure import sp_basis,minor_witness
ROOT=Path(__file__).resolve().parents[1]
def uniform(m,source_path):
 st=time.time();kind='M0';bs=s.symbols('b0:10');src=json.loads(source_path.read_text());r=src['cases'][kind]['universal'];one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];K=s.Matrix(r['beta_kernel']).applyfunc(s.sympify);two,F=specialize(m,dict(zip(one,K*s.Matrix(bs))),True);A=s.Matrix(r['leading_jacobian']).applyfunc(s.sympify);tc=s.Matrix(r['t_column']).applyfunc(s.sympify);C=s.Matrix(r['first_kernel']).applyfunc(s.sympify);pr=r['first_pivot_rows'];pc=r['first_pivot_columns'];U=s.zeros(44,len(pc));iv=A.extract(pr,pc).inv()
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 N=s.Matrix(r['residual_kernel']).applyfunc(s.sympify);T=(C*N[:-1,:]-U*tc[pr,:]*N[-1,:]).col_join(N[-1,:]);full_frame=T.copy();assert (A.row_join(tc).to_DM()*full_frame.to_DM()).is_zero_matrix;T=T[:,1:].applyfunc(s.cancel);den=[]
 for j in range(T.cols):
  d=s.factor(s.lcm([s.denom(T[i,j]) for i in range(T.rows)]));den.append(d)
  for i in range(T.rows):T[i,j]=s.cancel(T[i,j]*d)
 print('polynomial frame',time.time()-st,flush=True);E=s.Matrix(r['residual_matrix']).applyfunc(s.sympify).row_join(s.Matrix(r['residual_t_column']).applyfunc(s.sympify));sample=dict(zip(bs,range(1,11)));num=E.subs(sample);piv=list(num.to_DM().rref()[1]);rows=list(num[:,piv].T.to_DM().rref()[1]);cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);norm=list(combinations_with_replacement(range(T.cols),2));from sympy.polys.rings import ring
 RR,*gens=ring(','.join(map(str,bs)),s.QQ)
 TP=[[RR.from_expr(T[i,j]) for j in range(T.cols)] for i in range(T.rows)]
 HP=[[RR.zero for j in norm] for i in range(135)]
 forms={i:TP[cols.index(i)] for i in cols};forms[fixed]=[RR.zero]*T.cols
 for row,f in enumerate(F):
  for mon,c in f.items():
   fs=[forms[i] for i in mon]+[TP[-1]]*(2-len(mon));aa,bb=fs;cp=RR.from_expr(c)
   for j,(uu,vv) in enumerate(norm):HP[row][j]+=cp*(aa[uu]*bb[vv]+(aa[vv]*bb[uu] if uu!=vv else RR.zero))
 print('ring H',time.time()-st,flush=True)
 AU=A*U;YP=[]
 for i in range(135):
  coeff=[RR.from_expr(s.expand(AU[i,j])) for j in range(len(pr))]
  YP.append([HP[i][k]-sum((c*HP[r][k] for c,r in zip(coeff,pr)),RR.zero) for k in range(len(norm))])
 print('ring Y',time.time()-st,flush=True)
 values=list(range(1,11))
 def ev(p):return sum((c*s.prod(values[i]**e for i,e in enumerate(mon) if e) for mon,c in p.items()),s.S.Zero)
 numY=s.Matrix([[ev(p) for p in row] for row in YP]);Em=E.extract(rows,piv);numres=numY-num[:,piv]*num.extract(rows,piv).inv()*numY[rows,:];chosen=list(numres.T.to_DM().rref()[1]);assert len(chosen)==len(norm);print('rows',chosen,flush=True);Cr=(E[chosen,piv]*Em.inv(method='DM')).applyfunc(s.cancel)
 def cancel_ring(pn,denominator):
  co,fac=s.factor_list(denominator);left=s.sympify(co)
  for ff,ee in fac:
   fp=RR.from_expr(ff);unused=0
   for j in range(ee):
    qq,rem=pn.div(fp)
    if not rem:pn=qq
    else:unused=ee-j;break
   left*=ff**unused
  return pn.as_expr()/left
 out=s.zeros(len(norm),len(norm))
 for i,row in enumerate(chosen):
  dd=s.factor(s.lcm([s.denom(v) for v in Cr[i,:]]));dp=RR.from_expr(dd);coeff=[RR.from_expr(s.cancel(v*dd)) for v in Cr[i,:]]
  for j,(aa,bb) in enumerate(norm):
   numerator=dp*YP[row][j]-sum((c*YP[r][j] for c,r in zip(coeff,rows)),RR.zero);out[i,j]=cancel_ring(numerator,dd*den[aa]*den[bb])

 rowden=[];M=out.copy()
 for i in range(M.rows):
  dd=s.factor(s.lcm([s.denom(v) for v in M[i,:]]));rowden.append(dd)
  for j in range(M.cols):M[i,j]=s.cancel(M[i,j]*dd)
 dm=M.to_DM();numdet=dm.domain.to_sympy(dm.det());det=s.factor(s.factor(numdet)/s.prod(rowden));expected=s.Rational(59049,250)*bs[7]**4*bs[8]**6*(3*bs[1]-bs[9])**2*(8*bs[1]-bs[9])**2*(44*bs[1]-3*bs[9])*(4938*bs[1]**2-727*bs[1]*bs[9]+27*bs[9]**2)/((11*bs[1]-bs[9])**4*(194*bs[1]-3*bs[9])**3);assert det==expected;dm=E.extract(rows,piv).to_DM();jd=s.factor(dm.domain.to_sympy(dm.det()));assert jd==-s.Rational(5272018014312,9765625)*bs[1]*bs[7]**12*bs[8]*(8*bs[1]-bs[9])*(11*bs[1]-bs[9])*(194*bs[1]-3*bs[9])
 return dict(full_kernel_frame=enc(full_frame),normal_denominators=list(map(str,den)),normal_frame_polynomial=enc(T),quadratic_rows=chosen,quadratic_minor_matrix=enc(out),quadratic_minor=str(det),quadratic_row_denominators=list(map(str,rowden)),polynomial_determinant_factored=str(s.factor(numdet)),jacobian_minor_rows=rows,jacobian_minor_columns=piv,jacobian_minor=str(jd)),K,bs
def geometry(m,K,b):
 two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];u,act=beta_actions(m);Q=tensor('M0',0,0,2);q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q) if v);V=s.Matrix(X[1:]);qd=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fi]/q[fi];qd.append(dq)
 orbit=s.Matrix.hstack(*qd);stab=orbit.nullspace();assert orbit.rank()==7 and len(stab)==3;rows=list(K.T.rref()[1]);iv=K[rows,:].inv();actions=[]
 for c in stab:
  db=(act*c).subs(dict(zip(u,K*s.Matrix(b))));v=iv*db[rows,:];assert (K*v-db).applyfunc(s.expand)==s.zeros(16,1);actions.append(v)
 aa=s.Symbol('a');factors={'b1':b[1],'b7':b[7],'b8':b[8],'all_ratio_hyperplanes':aa*b[1]-b[9]};preserve={}
 for label,f in factors.items():
  df=[s.expand(sum(s.diff(f,bb)*v[j] for j,bb in enumerate(b))) for v in actions];assert all(s.rem(s.Poly(v,*b),s.Poly(f,*b)).is_zero for v in df);preserve[label]=list(map(str,df))
 proper=[]
 for idx in [1,7,8]:
  params=[bb for i,bb in enumerate(b) if i!=idx];beta=(K*s.Matrix(b)).subs(b[idx],0);D=beta.jacobian(params).row_join(act.subs(dict(zip(u,beta))));rr,piv=D.T.to_DM().rref();rank=len(piv);assert rank<=15;rr=rr.to_Matrix();free=[i for i in range(16) if i not in piv];N=s.zeros(16,len(free))
  for j,i in enumerate(free):
   N[i,j]=1
   for h,c in enumerate(piv):N[c,j]=-rr[h,i]
  assert (D.T.to_DM()*N.to_DM()).is_zero_matrix;proper.append(dict(factor=str(b[idx]),image_dimension=rank,beta_differential=enc(D),left_kernel=enc(N)))
 return dict(orbit_dimension=7,stabilizer_dimension=3,stabilizer_coefficients=[list(map(str,v)) for v in stab],stabilizer_actions=[list(map(str,v)) for v in actions],preservation=preserve,proper_beta_images=proper),u,act

def rational_case(m,K,b,u,act,curve,ratio):
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];gamma=s.Symbol('gamma');sample=dict(zip(b,range(1,11)));sample[b[9]]=ratio*sample[b[1]];bv=(K*s.Matrix(b)).subs(b[9],ratio*b[1]);HK=bv.jacobian(list(b[:9]));beta=bv.subs(sample);D=HK.row_join(act.subs(dict(zip(u,beta))));dom=minor_witness(D);assert dom['rank']==16;point=s.Matrix(coords(m,tensor('M0',0,0,2)))
 for i,v in zip(one,beta):point[i]=v
 local=calculate(m,True,point);assert local['leading_jacobian_rank']<=41 and local['homogenized_jacobian_rank']==41;cols=local['chart_columns'];fixed=two.index(local['chart_fixed_u']-1);q0=point[two,:];dq=curve.diff(gamma);dq-=q0*dq[fixed]/q0[fixed];tv=dq[cols,:].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.to_DM().rref()[1]);vt=T[rr,:].inv()*tv[rr,:];assert T*vt==tv and vt[0]!=0;branches=[dict(kind='M0_free_baseline',lower=4,lower_applies_to='full135, not its t quotient',tangent=list(map(str,vt)))]
 if ratio==3:
  assert local['quadratic_span_rank']==5;upper=quartic_bound(m,local,branches,4)
 else:
  assert local['quadratic_span_rank']==6 and local['common_hessian_rank']==3;rad=s.Matrix(local['quadratic_common_radical']).applyfunc(s.Rational);assert s.Matrix.hstack(vt,rad).rank()==1;upper=dict(transverse_quadratic_length=4)
 return dict(ratio=str(ratio),sample={str(k):str(v) for k,v in sample.items()},beta_minor=dom,beta_differential=enc(D),local_data=local,actual_curve_tangent=list(map(str,vt)),initial_bound=upper,full_generic_rank=4,generic_t_squared_zero=True,t_nilpotent=True,whole_generic_beta_stratum_removed=True)

def number_field_case(m,BK,b,u,act,curve,source_path):
 a=s.Symbol('a');f=27*a*a-727*a+4938;assert s.Poly(f,a).is_irreducible and s.gcd(f,s.diff(f,a))==1;root=s.CRootOf(f,0);K=s.QQ.algebraic_field(root);alpha=K.unit;zero=K.zero;oneK=K.one
 def cv(expr):return sum((K.from_sympy(c)*alpha**mon[0] for mon,c in s.Poly(expr,a).terms()),zero)
 def dm(M):return DM([[cv(v) for v in row] for row in M.tolist()],M.shape,K)
 def null(D):
  R,pv=D.rref();R=R.to_list();free=[i for i in range(D.shape[1]) if i not in pv];N=[[zero]*len(free) for i in range(D.shape[1])]
  for j,i in enumerate(free):
   N[i][j]=oneK
   for r,c in enumerate(pv):N[c][j]=-R[r][i]
  return DM(N,(D.shape[1],len(free)),K)
 def scalar(v):
  c=v.to_list();return str(sum(s.QQ.to_sympy(x)*a**(len(c)-i-1) for i,x in enumerate(c)))
 def en(D):return [[scalar(v) for v in row] for row in D.to_list()]
 def minor(D):
  cp=list(D.rref()[1]);rp=list(D.extract(range(D.shape[0]),cp).transpose().rref()[1]);det=D.extract(rp,cp).det();assert det;return dict(rows=rp,columns=cp,determinant=scalar(det))
 one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];sample=dict(zip(b,range(1,11)));sample[b[9]]=2*a;bv=(BK*s.Matrix(b)).subs(b[9],a*b[1]);HK=bv.jacobian(list(b[:9]));beta=bv.subs(sample);_,F=specialize(m,dict(zip(one,beta)),True);r=json.loads(source_path.read_text())['cases']['M0']['universal'];A=s.Matrix(r['leading_jacobian']).applyfunc(s.sympify).subs(sample);tc=s.Matrix(r['t_column']).applyfunc(s.sympify).subs(sample);D=dm(A.row_join(tc));T=null(D);L=null(D.transpose()).transpose();assert D.rank()==41 and dm(A).rank()==41 and T.shape==(45,4);assert (D*T).is_zero_matrix;cols=r['chart_columns'];fixed=next(i for i in range(45) if i not in cols);q0=s.Matrix(r['q_point']).applyfunc(s.Rational);assert all(cv(s.expand(sum(c*s.prod(q0[i] for i in mm) for mm,c in ff.items() if len(mm)==2)))==zero for ff in F)
 ET=[[zero]*4 for i in range(46)];tt=T.to_list()
 for j,i in enumerate(cols):ET[i]=tt[j]
 ET[45]=tt[44];mon=list(combinations_with_replacement(range(4),2));raw=[[zero]*10 for i in range(135)]
 for row,ff in enumerate(F):
  for mm,expr in ff.items():
   if fixed in mm:continue
   i,j=(*mm,*([45]*(2-len(mm))));c=cv(expr)
   for h,(v,w) in enumerate(mon):raw[row][h]+=c*(ET[i][v]*ET[j][w]+(ET[i][w]*ET[j][v] if v!=w else zero))
 Q=L*DM(raw,(135,10),K);assert Q.rank()==6;hess=[]
 for row in Q.to_list():
  if not any(row):continue
  H=[[zero]*4 for i in range(4)]
  for c,(i,j) in zip(row,mon):
   if i==j:H[i][i]+=2*c
   else:H[i][j]+=c;H[j][i]+=c
  hess+=H
 HD=DM(hess,(len(hess),4),K);assert HD.rank()==3;gamma=s.Symbol('gamma');dq=curve.diff(gamma);dq-=q0*dq[fixed]/q0[fixed];tv=dm(dq[cols,:].col_join(s.zeros(1,1)));rr=list(T.transpose().rref()[1]);vt=T.extract(rr,range(4)).inv()*tv.extract(rr,[0]);assert T*vt==tv and (HD*vt).is_zero_matrix and vt.rank()==1;BD=dm(HK.row_join(act.subs(dict(zip(u,beta)))));assert BD.rank()==16
 return dict(field='QQ[a]/('+str(f)+')',covers_both_complex_embeddings=True,beta_sample=list(map(str,beta)),full_rank=41,leading_rank=41,full_minor=minor(D),tangent_basis=en(T),quadratic_rank=6,quadratic_minor=minor(Q),quadratic_nonzero_rows={str(i):list(map(scalar,row)) for i,row in enumerate(Q.to_list()) if any(row)},hessian_rank=3,hessian_minor=minor(HD),actual_curve_tangent=en(vt),restricted_beta_rank=16,restricted_beta_minor=minor(BD),transverse_length=4,completed_ring='C[[u,e1,e2,e3]]/(e1,e2,e3)^2',t_squared_zero=True,whole_generic_beta_both_strata_removed=True)

def run():
 st=time.time();m=S61();source_path=ROOT/'reports/S61_MINIMAL_BASELINE_THICKENINGS_20261101.json';unif,K,b=uniform(m,source_path);geo,u,act=geometry(m,K,b);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];gamma=s.Symbol('gamma');curve=s.Matrix(coords(m,tensor('M0',0,gamma,2)))[two,:];_,F=specialize(m,dict(zip(one,K*s.Matrix(b))),True);assert all(s.expand(sum(c*s.prod(curve[i] for i in mm) for mm,c in ff.items() if len(mm)==2))==0 for ff in F);q0=curve.subs(gamma,0);cols=json.loads(source_path.read_text())['cases']['M0']['universal']['chart_columns'];fixed=next(i for i in range(45) if i not in cols);dq=curve.diff(gamma);dq-=q0*dq[fixed]/q0[fixed];T=s.Matrix(unif['full_kernel_frame']).applyfunc(s.sympify);assert (T[:,0]+dq[cols,:].col_join(s.zeros(1,1))).applyfunc(s.cancel)==s.zeros(45,1);assert T[-1,:]!=s.zeros(1,4);cases=[]
 for ratio in [s.Integer(11),s.Rational(194,3),s.Rational(44,3),s.Integer(3)]:
  cases.append(rational_case(m,K,b,u,act,curve,ratio));print('rational stratum',ratio,'verified',flush=True)
 field=number_field_case(m,K,b,u,act,curve,source_path);main_path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';main=json.loads(main_path.read_text())['cases']['M0'];A,B,P,R,U,V,KK,L=s.symbols('A B P R U V K L');bv=K*s.Matrix(b);sub={P:bv[1],R:bv[3],A:bv[4],B:bv[13],U:bv[6],V:bv[9],KK:bv[14]+5*bv[9],L:5*bv[6]-bv[15]};ma=s.factor(s.sympify(main['solution_a']).subs(sub,simultaneous=True));mb=s.factor(s.sympify(main['solution_b']).subs(sub,simultaneous=True));assert ma==0 and mb==(8*b[1]-b[9])/(6*b[7]);paired=normalization('M0');paths=[Path(__file__),source_path,main_path,ROOT/'code/s61_minimal_local_quartic_bound_20261101.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json']
 return dict(id='S61_MINIMAL_M0_BASELINE_REMOVAL_20261101',all_passed=True,field='QQ(b0,...,b9) and QQ[a]/(27a^2-727a+4938)',uniform=unif,geometry=geo,actual_curve_q=list(map(str,curve)),rational_exceptional_strata=cases,algebraic_exceptional_strata=field,main_intersection=dict(a=str(ma),b=str(mb),known_k='-2/3'),paired_normalization=paired,conclusion='Both finite M0 free baseline families are removed by t-saturation over sufficiently general beta.',limitations=['The full generic factor length4 is not the length of its t quotient, which is generally3.','At the quartic exception no flatness or full t^2=0 assertion is used; t^2 is generically zero and finite-length torsion proves nilpotence.','Projective parameter-infinity R0 flag endpoints remain separate.','No full affine-fiber exhaustion.'],proof_note='research/notes/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.md',seconds=time.time()-st,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
