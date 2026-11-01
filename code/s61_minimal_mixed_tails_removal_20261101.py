#!/usr/bin/env python3
"""Remove all four mixed g6/g7 minimal tails, including projective endpoints."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_baseline_thickenings_20261101 import coords,enc
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from s61_minimal_semisimple_tails_removal_20261101 import minor_witness
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]

def tensor(kind,a,b,c,d,h):
 x,y,z,w=X[1:];G=a*x*x*z+b*x*x*w+c*x*z*w+d*x*z*z+h*z*z*w;Q=s.zeros(4);Q[0,1]=x*x if kind=='M0' else x*x/2;Q[2,3]=0 if kind=='M0' else x*x/2;Q[1,2]=-s.diff(G,w)-(x*z if kind=='M1' else 0);Q[1,3]=s.diff(G,z);return Q-Q.T

def polynomial_kernel(kind,B,b):
 n0,n1,n2,n3,n4,n5,lam=b
 if kind=='M0':v=[0,0,0,0,lam,n0,B*lam,n1,n2,lam/4,n3,n4,n5,lam/2,-5*lam/4,(5*B+3)*lam]
 else:v=[0,0,0,0,-8*lam/3,n0,(4*B-1)*lam,n1,n2,lam,n3,n4,n5,-4*lam/3,-5*lam,(20*B-1)*lam]
 return s.Matrix(v).jacobian(b)

def full_kernel_certificate(m,Q,K,B):
 M=beta_matrix(m,Q);assert (M*K).applyfunc(s.expand)==s.zeros(M.rows,7)
 # A constant nonzero seven-minor makes this polynomial frame valid at every B.
 rows=list(K.subs(B,0).T.to_DM().rref()[1]);kd=s.factor(K[rows,:].det());assert kd!=0 and not kd.has(B)
 cert=[];g=s.Integer(0);bez=[];value=s.Integer(0)
 while True:
  N=M.subs(B,value);field=s.QQ if value.is_Rational else s.QQ.algebraic_field(value);cc=list(N.to_DM().convert_to(field).rref()[1]);rr=list(N[:,cc].T.to_DM().convert_to(field).rref()[1]);assert len(cc)==9;dm=M[rr,cc].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));a,d,h=s.gcdex(g,det,B);bez=[s.expand(a*v) for v in bez]+[d];g=h;cert.append(dict(sample_B=str(value),rows=rr,columns=cc,determinant=str(det)))
  if not g.has(B):break
  factor=s.factor_list(g,B)[1][0][0];value=s.solve(factor,B)[0];assert len(cert)<9
 assert g==1 and s.expand(sum(v*s.sympify(c['determinant']) for v,c in zip(bez,cert)))==1
 return dict(matrix=enc(M),kernel=enc(K),kernel_minor_rows=rows,kernel_minor=str(kd),rank=9,rank_minor_cover=cert,bezout_coefficients=list(map(str,bez)))

def universal(m,kind,B,b,K):
 Q=tensor(kind,1,B,0,0,1);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two,F=specialize(m,dict(zip(one,K*s.Matrix(b))),True);q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q) if v);co=[i for i in range(45) if i!=fi];A=s.zeros(135,45);tc=s.zeros(135,1)
 for row,f in enumerate(F):
  val=0
  for mon,c in f.items():
   if len(mon)==2:
    i,j=mon;val+=c*q[i]*q[j];A[row,i]+=c*q[j];A[row,j]+=c*q[i]
   elif len(mon)==1:tc[row]+=c*q[mon[0]]
  assert s.expand(val)==0
 A=A[:,co];A0=A[:70,:];num=A0.subs(B,2);pc=list(num.to_DM().rref()[1]);pr=list(num[:,pc].T.to_DM().rref()[1]);fc=[i for i in range(44) if i not in pc];assert len(pc)==32;iv=A0[pr,pc].inv(method='DM');C=s.zeros(44,len(fc));Cpc=-(iv*A0[pr,fc]).applyfunc(s.cancel)
 for j,i in enumerate(fc):C[i,j]=1
 for j,i in enumerate(pc):C[i,:]=Cpc[j,:]
 assert (A0.to_DM()*C.to_DM()).is_zero_matrix;U=s.zeros(44,len(pc))
 for j,i in enumerate(pc):U[i,:]=iv[j,:]
 dm=A0[pr,pc].to_DM();d0=s.factor(dm.domain.to_sympy(dm.det()));R=(A*C).applyfunc(s.cancel);ct=(tc-A*U*tc[pr,:]).applyfunc(s.cancel);E=R.row_join(ct);sample={B:2,**dict(zip(b,range(1,8)))};nn=E.subs(sample);cp=list(nn.to_DM().rref()[1]);rp=list(nn[:,cp].T.to_DM().rref()[1]);assert len(cp)==12;D=E[rp,cp];rowden=[]
 for i in range(D.rows):
  dd=s.factor(s.lcm([s.denom(v) for v in D[i,:]]));rowden.append(dd)
  for j in range(D.cols):D[i,j]=s.cancel(D[i,j]*dd)
 dm=D.to_DM();det=s.factor(dm.domain.to_sympy(dm.det())/s.prod(rowden));expected=-s.Rational(33480783,4000000)*B**4*b[6]**13 if kind=='M0' else s.Rational(8455716864,15625)*B**2*b[6]**13*(2*B+1)*(4*B-1);assert s.cancel(det-expected)==0
 return dict(q_point=list(map(str,q)),chart_columns=co,first_rank=32,first_minor=str(d0),first_minor_rows=pr,first_minor_columns=pc,first_kernel=enc(C),leading_jacobian=enc(A),t_column=list(map(str,tc)),residual_full_matrix=enc(E),full_minor=str(det),full_minor_rows=rp,full_minor_columns=cp,row_denominators=list(map(str,rowden)),full_rank_lower_bound=44,all135_highest_zero=True,top70_t_compatibility_zero=ct[:70,:]==s.zeros(70,1))

def proper_image(K,act,u):
 b=s.symbols('n0:'+str(K.cols));bv=K*s.Matrix(b);ar=act.subs(dict(zip(u,bv)));N=s.Matrix.hstack(*K.T.nullspace());W=ar.T*N;rr,piv=W.to_DM().rref();rr=rr.to_Matrix();fc=[i for i in range(N.cols) if i not in piv];C=s.zeros(N.cols,len(fc))
 for j,i in enumerate(fc):
  C[i,j]=1
  for r,p in enumerate(piv):C[p,j]=-rr[r,i]
 L=(N*C).applyfunc(s.cancel);assert (K.row_join(ar).T*L).applyfunc(s.cancel)==s.zeros(K.cols+10,len(fc));assert L.subs(dict(zip(b,range(1,K.cols+1)))).to_DM().rank()==len(fc)
 return dict(beta=list(map(str,bv)),image_dimension_upper=16-len(fc),left_kernel=enc(L),small_matrix=enc(W),small_kernel=enc(C))

def geometry(m,kind,B,b,K):
 u,act=beta_actions(m);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];V=s.Matrix(X[1:]);Q=tensor(kind,1,B,0,0,1);bv=K*s.Matrix(b);ar=act.subs(dict(zip(u,bv)));q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q) if v);assert s.diff(q[fi],B)==0;ds=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fi]/q[fi];ds.append(dq)
 O=s.Matrix.hstack(*ds);params=(*b,B);D=bv.jacobian(params).row_join(ar);joint=D.col_join(q.jacobian(params).row_join(O));sample={B:2,**dict(zip(b,range(1,8)))};stab=O.nullspace();assert len(stab)==1;assert (O*stab[0]).applyfunc(s.cancel)==s.zeros(45,1);rows=list(K.T.to_DM().rref()[1]);v=K[rows,:].inv()*(ar*stab[0])[rows,:];assert (K*v-ar*stab[0]).applyfunc(s.cancel)==s.zeros(16,1)
 generic=dict(beta_minor=minor_witness(D.subs(sample)),joint_minor=minor_witness(joint.subs(sample)),stabilizer=enc(stab[0]),stabilizer_action=list(map(str,v)),orbit_minor=minor_witness(O.subs(B,2)));assert generic['beta_minor']['rank']==16 and generic['joint_minor']['rank']==17
 fixed=[]
 for val in ([-1,0] if kind=='M0' else [0,-s.Rational(1,2),s.Rational(1,4)]):
  Kb=K.subs(B,val);Db=Kb.row_join(ar.subs(B,val));num=Db.subs(sample);r=num.to_DM().rank();case=dict(B=str(val),beta_rank=r,orbit_minor=minor_witness(O.subs(B,val)),beta_minor=minor_witness(num));assert case['orbit_minor']['rank']==9
  # The complete fixed-Q kernel is preserved by its entire connected stabilizer.
  ss=O.subs(B,val).nullspace();assert len(ss)==1;db=ar.subs(B,val)*ss[0];rows=list(Kb.T.to_DM().rref()[1]);vv=Kb[rows,:].inv()*db[rows,:];assert (Kb*vv-db).applyfunc(s.expand)==s.zeros(16,1);case.update(stabilizer=enc(ss[0]),stabilizer_action=list(map(str,vv)))
  if r<16:
   case['proper_image']=proper_image(Kb,act,u);assert case['proper_image']['image_dimension_upper']==15
  else:
   point=s.Matrix(coords(m,Q.subs(B,val)))
   for i,x in zip(one,bv.subs(B,val).subs(sample)):point[i]=x
   loc=calculate(m,True,point);assert loc['homogenized_jacobian_rank']==44;case['local_data']=loc;case['full_jacobian_minor']=minor_witness(s.Matrix(loc['jacobian']).applyfunc(s.Rational))
  fixed.append(case)
 lamzero=proper_image(K[:,:6],act,u);assert lamzero['image_dimension_upper']==9
 return dict(generic=generic,fixed_B=fixed,lambda_zero=lamzero)

def normalizations_and_coverage(m):
 V4=s.Matrix(X[1:]);x,y,z,w=V4;a,b,c,d,h,rho=s.symbols('a b c d h rho');K,L,U,V=s.symbols('K L U V');out={};old=json.loads((ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json').read_text());covered=[];k=s.Symbol('k');u,act=beta_actions(m);bb=s.symbols('b0:10');bad={}
 for axis in [6,7]:
  vv=[0,0,0,0,2*bb[7] if axis==6 else -bb[7],bb[0],bb[1],bb[2],bb[3],bb[8],bb[4],bb[5],bb[6],bb[7],bb[8] if axis==6 else bb[9],bb[9] if axis==6 else -bb[1]];mat=s.Matrix(vv).jacobian(bb);bad[str(axis)]=proper_image(mat,act,u);assert bad[str(axis)]['image_dimension_upper']==15
 for kind in ['M0','M1']:
  Bshift=-c/(2*h);Ashift=-d/h;S1=s.eye(4);S1[1,3]=-Bshift;S1[2,0]=Bshift;S2=s.eye(4);S2[1,2]=Ashift;S2[3,0]=Ashift;T=S1*S2;aa=a-d*(c+1 if kind=='M0' else c)/h;beta=b-c*(c-2 if kind=='M0' else c+2)/(4*h);assert T*J[1:,1:]*T.T==J[1:,1:];res=T.inv()*tensor(kind,a,b,c,d,h).subs(dict(zip(V4,T*V4)),simultaneous=True)*T.inv().T;assert (res-tensor(kind,aa,beta,0,0,h)).applyfunc(s.cancel)==s.zeros(4)
  R=s.diag(rho,1/rho,rho/h,h/rho);assert R*J[1:,1:]*R.T==J[1:,1:];res=R.inv()*tensor(kind,aa,beta,0,0,h).subs(dict(zip(V4,R*V4)),simultaneous=True)*R.inv().T/rho**2;assert (res-tensor(kind,aa*rho**2/h,beta*h,0,0,1)).applyfunc(s.cancel)==s.zeros(4)
  E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);paired=-c if kind=='M0' else -c-1;G7=a*x*x*z+b*x*x*w+c*x*z*w+d*x*w*w+h*z*w*w;Q7=s.zeros(4);Q7[0,1]=x*x if kind=='M0' else x*x/2;Q7[2,3]=0 if kind=='M0' else x*x/2;Q7[1,2]=-s.diff(G7,w)-(x*z if kind=='M1' else 0);Q7[1,3]=s.diff(G7,z);Q7-=Q7.T;res=E.inv()*Q7.subs(dict(zip(V4,E*V4)),simultaneous=True)*E.inv().T;assert (res-tensor(kind,b,-a,paired,d,-h)).applyfunc(s.expand)==s.zeros(4)
  out[kind]=dict(shear=enc(T),a_reduced=str(aa),b_reduced=str(beta),torus=enc(R),projective_scale='rho^(-2)',rho_square='h/a_reduced',canonical_B=str(s.factor(h*beta)),exchange=enc(E),paired_c=str(paired))
  fam=next(z for z in old['families'] if z['label']==kind)
  for tail in fam['rational_tails']:
   if tail['axis'] not in [6,7]:continue
   axis=tail['axis'];hh=s.Symbol(tail['parameter']);cf={name:s.sympify(expr) for name,expr in tail['coordinate_formulas'].items()};cf[tail['parameter']]=hh
   for name,expr in cf.items():assert s.cancel(expr.subs(hh,0)-s.sympify(fam['baseline_solution'][name]).subs(k,s.Rational(tail['k'])))==0
   if axis==6:sub={a:cf['g0'],b:cf['g1'],c:cf['g3'],d:cf['g2'],h:hh};badfactor=K-6*V
   else:sub={a:cf['g1'],b:-cf['g0'],c:-cf['g3']-(1 if kind=='M1' else 0),d:cf['g4'],h:-hh};badfactor=L-6*U
   numerator=s.factor((h*aa).subs(sub,simultaneous=True));assert s.cancel(numerator/(hh*badfactor))!=0 and not s.cancel(numerator/(hh*badfactor)).free_symbols;canonB=s.factor((h*beta).subs(sub,simultaneous=True))
   # The top coefficient at parameter infinity is a pure double-root R0 cubic.
   G=sum(cf[f'g{i}']*mon for i,mon in enumerate([x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]));Ginfty=s.expand(s.diff(G,hh));fact=(z-L*x/2)**2*(w+K*x/3) if axis==6 else (w+K*x/2)**2*(z-L*x/3);difference=s.expand(Ginfty-fact);assert s.diff(difference,z)==0 and s.diff(difference,w)==0
   Aend=-K/3 if axis==6 else -K/2;Bend=L/2 if axis==6 else L/3;Tend=s.Matrix([[1,0,0,0],[0,1,Aend,-Bend],[Bend,0,1,0],[Aend,0,0,1]]);assert Tend*J[1:,1:]*Tend.T==J[1:,1:];normalized=s.expand(Ginfty.subs(dict(zip(V4,Tend*V4)),simultaneous=True)-(z*z*w if axis==6 else z*w*w));assert s.diff(normalized,z)==0 and s.diff(normalized,w)==0
   covered.append(dict(kind=kind,axis=axis,k=tail['k'],parameter=str(hh),normalizer_numerator=str(numerator),excluded_factor=str(badfactor),canonical_B=str(canonB),zero_parameter_is_reviewed_main_junction=True,infinity_cubic=str(Ginfty),infinity_factor_mod_x3=str(fact),infinity_symplectic_matrix=enc(Tend),infinity_type='pure R0 z^2*w or z*w^2'))
 assert len(covered)==4
 return dict(normalizations=out,normalizer_zero_exclusions=bad,covered_presentations=covered)

def run():
 st=time.time();m=S61();B=s.Symbol('B');b=s.symbols('b0:7');cases=[]
 for kind in ['M0','M1']:
  K=polynomial_kernel(kind,B,b);Q=tensor(kind,1,B,0,0,1);kernel=full_kernel_certificate(m,Q,K,B);univ=universal(m,kind,B,b,K);geo=geometry(m,kind,B,b,K);cases.append(dict(kind=kind,beta=list(map(str,K*s.Matrix(b))),kernel=kernel,universal=univ,geometry=geo));print(kind,'verified',flush=True)
 coverage=normalizations_and_coverage(m);paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_minimal_semisimple_tails_removal_20261101.py',ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_EXHAUSTION_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',m.source/'symplectic_slice_reduced_equations.json']
 return dict(id='S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101',all_passed=True,field='QQ and QQ(B,b0,...,b6)',cases=cases,coverage=coverage,conclusion='All four mixed g6/g7 tail presentations, including parameter zero and projective parameter infinity, are removed over sufficiently general beta.',limitations=['Other minimal main-curve exceptional loci and other binary/flag endpoints remain separate.','No full generic affine-fiber finiteness or five-point exhaustion.'],proof_note='research/notes/S61_MINIMAL_MIXED_TAILS_REMOVAL_20261101.md',seconds=time.time()-st,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
