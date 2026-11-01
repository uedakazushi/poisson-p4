#!/usr/bin/env python3
"""Remove the full finite M1 triple-tail pair, including its zero orbit."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_cubic_tail_20261101 import tensor,coords,enc,normalization
from s61_minimal_semisimple_tails_removal_20261101 import proper_kernel,minor_witness
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_two_variable_initial_20261101 import initial_data
from s61_nil22_infinity_20261101 import beta_matrix
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
def uniform(m,source):
 st=time.time();r=json.loads(source.read_text())['universal'];b=s.symbols('b0:10');K=s.Matrix(r['beta_kernel']).applyfunc(s.sympify);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,dict(zip(one,K*s.Matrix(b))),True);q=s.Matrix(r['q_point']).applyfunc(s.Rational);co=r['chart_columns'];fi=next(i for i in range(45) if i not in co);A=s.Matrix(r['leading_jacobian']).applyfunc(s.sympify);C=s.Matrix(r['first_kernel']).applyfunc(s.sympify);T=s.Matrix(r['tangent']).applyfunc(s.sympify);L=s.Matrix(r['quadratic_lift']).applyfunc(s.sympify);A0=A[:70,:];pc=list(A0.rref()[1]);pr=list(A0[:,pc].T.rref()[1]);U=s.zeros(44,len(pc));iv=A0[pr,pc].inv()
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 tc=s.Matrix([s.expand(sum(c*q[mon[0]] for mon,c in ff.items() if len(mon)==1)) for ff in F]);R=A*C;ct=(tc-A*U*tc[pr,:]).applyfunc(s.cancel);E=R.row_join(ct);rf,ep=E.to_DM().rref();rf=rf.to_Matrix();free=[i for i in range(20) if i not in ep];N=s.zeros(20,len(free))
 for j,i in enumerate(free):
  N[i,j]=1
  for h,c in enumerate(ep):N[c,j]=-rf[h,i]
 assert (E.to_DM()*N.to_DM()).is_zero_matrix;assert len(ep)==18 and N[-1,:]==s.zeros(1,2);sample=dict(zip(b,range(1,11)));num=E.subs(sample);cols=list(num.to_DM().rref()[1]);rows=list(num[:,cols].T.to_DM().rref()[1]);dm=E[rows,cols].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));print('full Jacobian verified',flush=True)
 results=[]
 for which,idx in [(0,0),(1,2)]:
  v=T[:,which];w=-L[:,idx];dv=s.factor(s.lcm([s.denom(s.cancel(x)) for x in v]));dw=s.factor(s.lcm([s.denom(s.cancel(x)) for x in w]));v=v.applyfunc(lambda x:s.cancel(x*dv));w=w.applyfunc(lambda x:s.cancel(x*dw));vn=s.zeros(45,1);wn=s.zeros(45,1)
  for j,i in enumerate(co):vn[i]=v[j];wn[i]=w[j]
  cubic=s.Matrix([s.expand(sum(c*(vn[mon[0]]*wn[mon[1]]+vn[mon[1]]*wn[mon[0]]) for mon,c in ff.items() if len(mon)==2)) for ff in F]);z=(cubic-A*U*cubic[pr,:]).applyfunc(s.cancel);den=dv*dw;znum=z;aug=E[:,cols].row_join(znum);nn=aug.subs(sample);rr=list(nn.T.to_DM().rref()[1]);print('which',which,'aug rank',len(rr),'den',den,flush=True)
  if which==0:assert L[:,0]==s.zeros(44,1) and znum==s.zeros(135,1)
  if len(rr)==19:
   dm=aug[rr,:].to_DM();d=s.factor(dm.domain.to_sympy(dm.det()));print('cubicminor',d,flush=True);results.append(dict(which=which,rows=rr,cols=cols,den=str(den),minor=str(d),z_numerator=enc(znum)))
 return dict(beta_kernel=enc(K),full_residual=enc(E),full_kernel=enc(N),full_minor=str(det),full_minor_rows=rows,full_minor_columns=cols,cubic=results,first_pure_cubic_zero_verified=True,seconds=time.time()-st)

def orbit_data(m,Q,K):
 b=s.symbols('b0:10');u,act=beta_actions(m);bv=K*s.Matrix(b);ar=act.subs(dict(zip(u,bv)));two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];q=s.Matrix(coords(m,Q))[two,:];fixed=next(i for i,v in enumerate(q) if v);V=s.Matrix(X[1:]);vectors=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fixed]/q[fixed];vectors.append(dq)
 orbit=s.Matrix.hstack(*vectors);stab=orbit.nullspace();rows=list(K.T.rref()[1]);iv=K[rows,:].inv();actions=[]
 for h in stab:
  v=iv*(ar*h)[rows,:];assert (K*v-ar*h).applyfunc(s.expand)==s.zeros(16,1);actions.append(v)
 D=K.row_join(ar);joint=D.col_join(s.zeros(45,10).row_join(orbit));sample=dict(zip(b,range(1,11)));dm=minor_witness(D.subs(sample));jm=minor_witness(joint.subs(sample));return dict(orbit_dimension=orbit.rank(),stabilizer=enc(s.Matrix.hstack(*stab)),stabilizer_actions=[list(map(str,v)) for v in actions],beta_minor=dm,joint_minor=jm),b,bv,ar,actions

def curve_data(m,a0,K,local=None):
 b=s.symbols('b0:10');d=s.Symbol('d');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,dict(zip(one,K*s.Matrix(b))),True);a=a0+2*b[9]*d/(3*b[7]);q=s.Matrix(coords(m,tensor(a,0,d)))[two,:];assert all(s.cancel(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F);result=dict(parameter='d',shape_a=str(a),shape_b='0',q=list(map(str,q)),all135_highest_zero=True)
 if local is not None:
  qnum=q.subs(dict(zip(b,range(1,11))));q0=qnum.subs(d,0);co=local['chart_columns'];fixed=next(i for i in range(45) if i not in co);dq=qnum.diff(d);dq-=q0*dq[fixed]/q0[fixed];tv=dq[co,:].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.to_DM().rref()[1]);vt=T[rr,:].inv()*tv[rr,:];assert T*vt==tv and vt[0]!=0;result['sample_tangent']=list(map(str,vt))
 return result

def factor_cases(m,K,univ):
 geo,b,bv,ar,actions=orbit_data(m,tensor(1,0,0),K);assert geo['orbit_dimension']==7 and geo['beta_minor']['rank']==16 and geo['joint_minor']['rank']==17;one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];sample=dict(zip(b,range(1,11)));out=[]
 for f,e in s.factor_list(s.sympify(univ['full_minor']))[1]:
  weights=[s.cancel(sum(s.diff(f,bb)*v[j] for j,bb in enumerate(b))/f) for v in actions];assert all(not(z.free_symbols & set(b)) for z in weights);var=next(bb for bb in b if s.degree(f,bb)==1);assert s.gcd(s.diff(f,var),f.subs(var,0))==1;value=s.solve(f,var)[0];bp=bv.subs(var,value);params=[bb for bb in b if bb!=var];D=bp.jacobian(params).row_join(ar.subs(var,value));numeric=D.subs(sample);rank=numeric.to_DM().rank();record=dict(factor=str(f),power=int(e),stabilizer_weights=list(map(str,weights)),solved_variable=str(var),solved_value=str(value),sample_beta_rank=rank)
  if rank<16:
   cert=proper_kernel(K,b,bv,ar,f,var,value);assert cert['image_dimension_upper']==rank;record.update(type='proper_beta_image',certificate=cert)
  elif f==b[0]:record.update(type='reviewed_main_junction',known_k='-1/3',beta_minor=minor_witness(numeric))
  else:
   assert s.expand(f-(28*b[0]*b[9]-45*b[7]**2))==0;point=s.Matrix(coords(m,tensor(1,0,0)))
   for i,v in zip(one,bp.subs(sample)):point[i]=v
   local=calculate(m,True,point);assert local['leading_jacobian_rank']==42 and local['homogenized_jacobian_rank']==43 and local['quadratic_span_rank']==0;initial=initial_data(m,local);record.update(type='cubic_artifact_stratum',beta_minor=minor_witness(numeric),local_data=local,initial=initial,completed_ring='C[[u,e]]/(e^3)',t_zero=True)
  out.append(record)
 return geo,out

def zero_orbit(m):
 Q=tensor(0,0,0);K=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert K.shape==(16,10);geo,b,bv,ar,actions=orbit_data(m,Q,K);assert geo['orbit_dimension']==6 and geo['beta_minor']['rank']==16 and geo['joint_minor']['rank']==16;one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];point=s.Matrix(coords(m,Q))
 for i,v in zip(one,K*s.Matrix(range(1,11))):point[i]=v
 local=calculate(m,True,point);assert local['leading_jacobian_rank']==42 and local['homogenized_jacobian_rank']==43 and local['quadratic_span_rank']==0;curve=curve_data(m,0,K,local);initial=initial_data(m,local,quartic=True)
 a=s.Symbol('a');family_beta=list(bv);family_beta[9]=3*a*b[7];two,F=specialize(m,dict(zip(one,family_beta)),True);q=s.Matrix(coords(m,tensor(a,0,0)))[two,:];assert all(s.expand(sum(c*s.prod(q[i] for i in mon) for mon,c in f.items() if len(mon)==2))==0 for f in F)
 return dict(beta_kernel=enc(K),geometry=geo,local_data=local,same_beta_curve=curve,initial=initial,closure_family_beta=list(map(str,family_beta)),closure_family_q=list(map(str,q)),generic_rank=3,transverse_length=4,t_nilpotent=True,flatness_asserted=False,whole_generic_beta_stratum_removed=True)

def finite_normalization():
 a,b,d=s.symbols('a b d');V=s.Matrix(X[1:]);T1=s.eye(4);T1[3,2]=-2*d;T2=s.eye(4);T2[1,3]=b;T2[2,0]=-b;S=T1*T2;assert S*J[1:,1:]*S.T==J[1:,1:];transformed=S.inv()*tensor(a,b,d).subs(dict(zip(V,S*V)),simultaneous=True)*S.inv().T;assert (transformed-tensor(a-2*b*d,0,0)).applyfunc(s.expand)==s.zeros(4);old=normalization();E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);paired=E.inv()*tensor(a,b,d).subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T
 from s61_minimal_baseline_thickenings_20261101 import tensor as baseline
 target=baseline('M1',b,-a,-1);target[1,2]-=2*d*V[0]*V[3];target[2,1]+=2*d*V[0]*V[3];assert (paired-target).applyfunc(s.expand)==s.zeros(4)
 return dict(finite_shear=enc(S),normalization_invariant='a-2*b*d',nonzero_orbit=old,zero_orbit_included=True,paired_exchange=enc(E),paired_source_c=-1)

def run():
 st=time.time();m=S61();source=ROOT/'reports/S61_MINIMAL_CUBIC_TAIL_20261101.json';univ=uniform(m,source);K=s.Matrix(univ['beta_kernel']).applyfunc(s.sympify);geo,factors=factor_cases(m,K,univ);curve=curve_data(m,1,K);zero=zero_orbit(m);normal=finite_normalization();mainpath=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';main=json.loads(mainpath.read_text())['cases']['M1'];b=s.symbols('b0:10');bv=K*s.Matrix(b);P,R,A,B,U,V,KK,L=s.symbols('P R A B U V K L');sub={P:bv[1],R:bv[3],A:bv[4],B:bv[13],U:bv[6],V:bv[9],KK:bv[14]+5*bv[9],L:5*bv[6]-bv[15]};ma=s.factor(s.sympify(main['solution_a']).subs(sub,simultaneous=True));mb=s.factor(s.sympify(main['solution_b']).subs(sub,simultaneous=True));assert ma.subs(b[0],0)==1 and mb==0;paths=[Path(__file__),source,mainpath,ROOT/'code/s61_minimal_two_variable_initial_20261101.py',ROOT/'code/s61_minimal_cubic_tail_20261101.py',ROOT/'code/s61_minimal_semisimple_tails_removal_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_model.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json']
 return dict(id='S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101',all_passed=True,field='QQ and QQ(b0,...,b9)',universal=univ,geometry=geo,factors=factors,same_beta_curve=curve,zero_orbit=zero,normalization=normal,main_intersection=dict(a=str(ma),b=str(mb),condition='b0=0'),conclusion='The entire finite M1 g2 k=-1/3 and g4 k=-2/3 tail pair is removed over sufficiently general beta.',limitations=['The zero orbit uses rank3/transverse4 and no flatness assertion.','Projective parameter-infinity binary endpoints remain separate.','Other tail families and full affine exhaustion remain open.'],seconds=time.time()-st,proof_note='research/notes/S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
