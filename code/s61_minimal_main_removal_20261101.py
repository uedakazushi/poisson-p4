#!/usr/bin/env python3
"""All finite M0/M1 main-curve points over sufficiently general first jets."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_minimal_baseline_thickenings_20261101 import tensor,coords,enc
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_main_jordan_removal_20261101 import orbit
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_semisimple_tails_removal_20261101 import minor_witness
from s61_minimal_main_uniform_data_20261101 import uniform
from s61_minimal_main_exception_field_20261101 import number_field_case
from s61_minimal_local_initial_bound_20261101 import cubic_bound
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS
ROOT=Path(__file__).resolve().parents[1]

def proper_image(beta,params,act,u):
 base=beta.jacobian(params);ar=act.subs(dict(zip(u,beta)));N=s.Matrix.hstack(*base.T.nullspace());W=ar.T*N;rf,piv=W.to_DM().rref();rf=rf.to_Matrix();free=[i for i in range(N.cols)if i not in piv];C=s.zeros(N.cols,len(free))
 for j,i in enumerate(free):
  C[i,j]=1
  for r,p in enumerate(piv):C[p,j]=-rf[r,i]
 L=(N*C).applyfunc(s.cancel);D=base.row_join(ar);assert(D.T.to_DM()*L.to_DM()).is_zero_matrix;assert L.to_DM().rank()==len(free)
 return dict(image_dimension_upper=16-len(free),beta=list(map(str,beta)),parameters=list(map(str,params)),left_kernel=enc(L),small_matrix=enc(W),small_kernel=enc(C))

def main_branch(m,kind,local,main):
 pars=list(map(s.sympify,main['source_variables']));sb=s.Matrix(main['source_beta']).applyfunc(s.sympify);q=s.Matrix(main['cases'][kind]['q_coordinates']).applyfunc(s.sympify);one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2];p=s.Matrix(local['point61']).applyfunc(s.Rational);beta=p[one,:];q0=p[two,:];sample=s.solve(list(sb-beta),pars,dict=True)[0];assert q.subs(sample)==q0;u,act=beta_actions(m);BD=sb.jacobian(pars).subs(sample).row_join(act.subs(dict(zip(u,beta))));co=local['chart_columns'];fi=next(i for i in range(45)if i not in co);Dq=q.jacobian(pars).subs(sample);Dq-=q0*Dq[fi,:]/q0[fi]
 # Reconstruct the full affine quadratic tensor from its source formula.
 A,B,P,R,U,V,KK,L=s.symbols('A B P R U V K L');x,y,z,w=X[1:];qv=R*z*z+(2*A-B)*z*w-P*w*w;H=2*qv/(3*B)if kind=='M0'else -qv/(6*B)-z*w/2;aa=s.sympify(main['cases'][kind]['solution_a']);bb=s.sympify(main['cases'][kind]['solution_b']);G=aa*x*x*z+bb*x*x*w+x*H;Q=s.zeros(4);Q[0,1]=x*x if kind=='M0'else x*x/2;Q[2,3]=0 if kind=='M0'else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1'else 0);Q[1,3]=s.diff(G,z);Q=(Q-Q.T).subs(sample);O=orbit(m,Q,two,fi);QD=Dq[co,:].row_join(O[co,:]);N=s.Matrix.hstack(*BD.nullspace());C=QD*N;assert C.to_DM().rank()==1;j=next(j for j in range(C.cols)if C[:,j]!=s.zeros(44,1));tv=C[:,j].col_join(s.zeros(1,1));T=s.Matrix(local['tangent_basis']).applyfunc(s.Rational);rr=list(T.T.to_DM().rref()[1]);vt=T[rr,:].inv()*tv[rr,:];assert T*vt==tv;first=next(i for i,v in enumerate(vt)if v);perm=[first]+[i for i in range(vt.rows)if i!=first];local['tangent_basis']=enc(T[:,perm]);vt=vt[perm,:];lower=3 if kind=='M0'else 2;bound=cubic_bound(m,local,[dict(kind=kind,lower=lower,tangent=list(map(str,vt)))],lower);bd=minor_witness(BD);jd=minor_witness(BD.col_join(QD));assert bd['rank']==16 and jd['rank']==17
 return dict(source_point={str(k):str(v)for k,v in sample.items()},source_direction=list(map(str,N[:,j])),beta_minor=bd,joint_minor=jd,tangent_permutation=perm,main_tangent=list(map(str,vt)),initial_bound=bound)

def geometry_and_cases(m,kind,data,main):
 c=s.Symbol('c');b=s.symbols('b0:9');bv=s.Matrix(data['beta']).applyfunc(s.sympify);K=bv.jacobian(b);Q=tensor(kind,0,0,c if kind=='M0'else(c-1)/2);two=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==2];one=[i for i,w in enumerate(m.reduced['parameter_weights'])if w==1];q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q)if v);O=orbit(m,Q,two,fi);S=s.Matrix.hstack(*O.nullspace());assert S.cols==3;u,act=beta_actions(m);ar=act.subs(dict(zip(u,bv)));kr=list(K.T.to_DM().rref()[1]);acts=(K[kr,:].inv()*ar[kr,:]*S).applyfunc(s.cancel);assert(K*acts-ar*S).applyfunc(s.cancel)==s.zeros(16,3);assert acts[6:,:]==s.Matrix([[0,-b[6],0],[0,0,b[7]],[0,0,-b[8]]])
 params=(c,*b);D=bv.jacobian(params).row_join(ar);joint=D.col_join(q.jacobian(params).row_join(O));sample={c:5,**dict(zip(b,range(1,10)))};bd=minor_witness(D.subs(sample));jd=minor_witness(joint.subs(sample));assert bd['rank']==16 and jd['rank']==17
 # Every polynomial point is the established main formula, not a sampled branch.
 pars=list(map(s.sympify,main['source_variables']));sb=s.Matrix(main['source_beta']).applyfunc(s.sympify);mq=s.Matrix(main['cases'][kind]['q_coordinates']).applyfunc(s.sympify);sp=s.solve(list(sb-bv),pars,dict=True)[0];assert(mq.subs(sp)-q).applyfunc(s.cancel)==s.zeros(45,1)
 proper=[];drop=bv.subs(b[6],0);pp=proper_image(drop,[c]+list(b[:6])+list(b[7:]),act,u);assert pp['image_dimension_upper']==13;proper.append(dict(condition='b6=0',certificate=pp))
 for value in ([s.Integer(0),s.Rational(6,5),-s.Rational(6,5)]if kind=='M0'else[s.Integer(0),s.Rational(3,5),-s.Rational(3,5)]):
  pp=proper_image(bv.subs(c,value),b,act,u);assert pp['image_dimension_upper']==(14 if value==0 else 15);proper.append(dict(condition='c='+str(value),certificate=pp))
 cases=[]
 for condition,value in [('b7',s.Integer(5))]+([('c',s.Rational(2,5))]if kind=='M0'else[('c',s.Rational(1,5)),('c',s.Integer(2))]):
  beta=bv.subs(b[7],0)if condition=='b7'else bv.subs(c,value);free=[c]+[x for x in b if x!=b[7]]if condition=='b7'else list(b);Ds=beta.jacobian(free).row_join(act.subs(dict(zip(u,beta))));qs=q if condition=='b7'else q.subs(c,value);Os=O if condition=='b7'else O.subs(c,value);joint=Ds.col_join(qs.jacobian(free).row_join(Os));sample={c:value,**dict(zip(b,range(1,10)))};dom=minor_witness(Ds.subs(sample));jdom=minor_witness(joint.subs(sample));assert dom['rank']==jdom['rank']==16;point=s.Matrix(coords(m,Q.subs(c,value)))
  for i,v in zip(one,beta.subs(sample)):point[i]=v
  local=calculate(m,True,point);branch=main_branch(m,kind,local,main);cases.append(dict(condition=condition+'='+('0'if condition=='b7'else str(value)),sample={str(k):str(v)for k,v in sample.items()},stratum_beta_minor=dom,stratum_joint_minor=jdom,local_data=local,main_branch=branch));print(kind,condition,value,'exception passed',flush=True)
 critical=None
 if kind=='M0':
  numerator=s.fraction(s.sympify(data['normal_determinant']))[0];P=max((f for f,e in s.factor_list(numerator)[1]),key=lambda f:len(s.Poly(f,c,*b).terms()));coef=s.diff(P,b[4]);assert s.degree(P,b[4])==1 and s.gcd(coef,P.subs(b[4],0))==1;der=s.expand(sum(s.diff(P,bb)*acts[i,0]for i,bb in enumerate(b)));assert s.gcd(P,der)==1;critical=dict(polynomial=str(P),primitive_linear_variable='b4',coefficient=str(s.factor(coef)),stabilizer_derivative=str(der),coprime=True,intrinsic_saturation_support_image_dimension_upper=15)
 return dict(stabilizer=enc(S),stabilizer_actions=enc(acts),orbit_dimension=7,beta_minor=bd,joint_minor=jd,universal_main_source_point={str(k):str(v)for k,v in sp.items()},proper_images=proper,exceptional_strata=cases,noninvariant_normal_factor=critical)

def finite_transformations(m,uniforms):
 a,bb,c=s.symbols('a b c');V=s.Matrix(X[1:]);E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);assert E*J[1:,1:]*E.T==J[1:,1:];z,u,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(z[h])for(i,j),h in COORDS]);basis=vec.jacobian(u);rr=list(basis.T.to_DM().rref()[1]);left=basis[rr,:].inv();PB=E.T*Beta.subs(dict(zip(z,E*s.Matrix(z))),simultaneous=True)*E;vv=s.Matrix([PB[i,j].diff(z[h])for(i,j),h in COORDS]);finitebeta=left*vv[rr,:];out={};b=s.symbols('b0:9')
 for kind in ['M0','M1']:
  cc=c if kind=='M0'else(c-1)/2;A=-a/(cc+1)if kind=='M0'else-a/cc;B=-bb/(cc-1)if kind=='M0'else-bb/(cc+1);T=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]]);assert(T*J[1:,1:]*T.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);pull=T.inv()*tensor(kind,a,bb,cc).subs(dict(zip(V,T*V)),simultaneous=True)*T.inv().T;assert(pull-tensor(kind,0,0,cc)).applyfunc(s.cancel)==s.zeros(4);exchange=E.inv()*tensor(kind,0,0,cc).subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T;assert(exchange-tensor(kind,0,0,-c if kind=='M0'else(-c-1)/2)).applyfunc(s.expand)==s.zeros(4)
  bv=s.Matrix(uniforms[kind]['beta']).applyfunc(s.sympify);nb=finitebeta.subs(dict(zip(u,bv)));den=-6 if kind=='M0'else 8;newpars=[nb[i]for i in [5,7,8,10,11,12]]+[nb[13],(5*nb[9]+nb[14])/den,(5*nb[6]-nb[15])/den];assert(nb-bv.subs(c,-c).subs(dict(zip(b,newpars)),simultaneous=True)).applyfunc(s.expand)==s.zeros(16,1);assert s.expand(newpars[8].subs(b[7],0))==0
  out[kind]=dict(normalization=enc(T),denominator='c^2-1',exchange=enc(E),exchanged_beta_parameters=list(map(str,newpars)),b7_zero_maps_to_b8_zero=True)
 return out

def rank_cover(m):
 c=s.Symbol('c');out={}
 for kind in ['M0','M1']:
  Q=tensor(kind,0,0,c if kind=='M0'else(c-1)/2);M=beta_matrix(m,Q);N=M.subs(c,5);cc=list(N.to_DM().rref()[1]);rr=list(N[:,cc].T.to_DM().rref()[1]);dm=M[rr,cc].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));assert len(cc)==7;ex=[]
  for value in s.solve(det,c):
   N=M.subs(c,value);rank=N.to_DM().rank();ex.append(dict(c=str(value),rank=rank,kernel_dimension=16-rank))
  out[kind]=dict(minor=str(det),rows=rr,columns=cc,exceptions=ex)
 return out

def factor_cover(uniforms):
 c=s.Symbol('c');b=s.symbols('b0:9');out={}
 tables={'M0':{s.Integer(0):'proper image14',s.Rational(6,5):'proper image15',s.Integer(2):'reviewed M0 baseline',s.Integer(1):'reviewed M0 Jordan',s.Rational(1,2):'reviewed semisimple tails',s.Rational(4,5):'reviewed semisimple tails',s.Rational(2,3):'reviewed main-tail junction',s.Rational(2,5):'new fixed-c quadratic stratum'},'M1':{s.Integer(0):'proper image14',s.Rational(3,5):'proper image15',s.Integer(3):'reviewed M1 baseline',s.Integer(1):'reviewed triple tails',s.Rational(2,5):'reviewed semisimple tails',s.Rational(1,3):'reviewed main-tail junction',s.Rational(1,5):'new fixed-c quadratic stratum',s.Integer(2):'new fixed-c quadratic stratum'}}
 for kind,r in uniforms.items():
  seen={}
  for key in ['first_minor','jacobian_residual_determinant','all_frame_denominator_factors','normal_determinant']:
   for p in s.fraction(s.sympify(r[key])):
    for f,e in s.factor_list(p)[1]:
     if not f.free_symbols:continue
     mon=s.Poly(f,c,*b).monic().as_expr()
     if mon in seen:continue
     if f.free_symbols<={c}:
      if s.degree(f,c)==1:value=-f.subs(c,0)/s.diff(f,c);label=tables[kind][abs(value)]
      else:assert kind=='M1' and s.cancel(f/s.LC(s.Poly(f,c)))==c*c-s.Rational(1,7);label='new two-embedding quadratic-field stratum'
     elif f in [b[6],b[7],b[8]]:label='proper image13'if f==b[6]else'new main-divisor stratum and Weyl pair'
     else:assert kind=='M0' and s.degree(f,b[4])==1 and s.gcd(s.diff(f,b[4]),f.subs(b[4],0))==1;label='noninvariant primitive-linear factor: intrinsic bad image<=15'
     seen[mon]=dict(factor=str(f),handled_by=label)
  out[kind]=list(seen.values())
 return out

def run():
 st=time.time();m=S61();rpath=ROOT/'reports/S61_MINIMAL_INFINITY_FORMAL_RANKS_20261101.json';mpath=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';ranks=json.loads(rpath.read_text());main=json.loads(mpath.read_text());un={};geo={}
 for kind in ['M0','M1']:
  un[kind]=uniform(m,ranks,kind);geo[kind]=geometry_and_cases(m,kind,un[kind],main)
 field=number_field_case(m,ranks,main);trans=finite_transformations(m,un);cover=rank_cover(m);factors=factor_cover(un)
 paths=[Path(__file__),rpath,mpath,ROOT/'code/s61_minimal_main_uniform_data_20261101.py',ROOT/'code/s61_minimal_main_exception_field_20261101.py',ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_local_initial_bound_20261101.py',ROOT/'code/s61_minimal_main_jordan_removal_20261101.py',ROOT/'reports/S61_MINIMAL_MAIN_JORDAN_REMOVAL_20261101.json',ROOT/'reports/S61_MINIMAL_MAIN_JORDAN_REVIEW_20261101.json',ROOT/'reports/S61_MINIMAL_M0_BASELINE_REMOVAL_20261101.json',ROOT/'reports/S61_MINIMAL_M1_BASELINE_REMOVAL_20261101.json',ROOT/'reports/S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL_20261101.json',ROOT/'reports/S61_MINIMAL_TRIPLE_TAILS_REMOVAL_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md']
 return dict(id='S61_MINIMAL_MAIN_REMOVAL_20261101',all_passed=True,uniform=un,geometry=geo,algebraic_exception=field,finite_transformations=trans,rank_cover=cover,complete_factor_cover=factors,conclusion='All finite semisimple-vertex main-curve points of M0 and M1 are removed over sufficiently general beta, using the reviewed exceptional baseline, tail and Jordan strata.',limitations=['Rank-one vertex main points use the separately reviewed rank-one theorems.','Projective main-pole binary/flag endpoints remain separate.','No full generic affine-fiber exhaustion.'],seconds=time.time()-st,proof_note='research/notes/S61_MINIMAL_MAIN_REMOVAL_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x')as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
