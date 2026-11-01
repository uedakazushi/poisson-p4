#!/usr/bin/env python3
"""Remove six finite semisimple g2/g4 tail presentations over general beta."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from s61_nil22_infinity_20261101 import beta_matrix
from s61_flag_transverse_obstruction_20261101 import specialize
from s61_minimal_infinity_local_data_20261101 import calculate
from s61_minimal_baseline_thickenings_20261101 import tensor,coords,enc
from s61_minimal_rankone_A1_removal_20261101 import beta_actions
from weight1_structure import sp_basis
ROOT=Path(__file__).resolve().parents[1]
def minor_witness(matrix):
 columns=list(matrix.to_DM().rref()[1]);rows=list(matrix[:,columns].T.to_DM().rref()[1]);part=matrix[rows,columns];dm=part.to_DM();det=dm.domain.to_sympy(dm.det());assert det
 return dict(rank=len(columns),rows_zero_based=rows,columns_zero_based=columns,determinant=str(det),matrix=enc(part))
def universal(m,kind,c):
    Q=tensor(kind,0,0,c);pt=coords(m,Q);BK=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert BK.shape==(16,10);one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];b=s.symbols('b0:10');two,F=specialize(m,dict(zip(one,BK*s.Matrix(b))),True);q=s.Matrix([pt[i] for i in two]);fixed=next(i for i,v in enumerate(q) if v);cols=[i for i in range(45) if i!=fixed];A=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(F):
        val=0
        for mon,v in f.items():
            if len(mon)==2:
                i,j=mon;val+=v*q[i]*q[j];A[row,i]+=v*q[j];A[row,j]+=v*q[i]
            elif len(mon)==1:tc[row]+=v*q[mon[0]]
        assert s.expand(val)==0
    A=A[:,cols];A0=A[:70,:];cc=list(A0.rref()[1]);rr=list(A0[:,cc].T.rref()[1]);C=s.Matrix.hstack(*A0.nullspace());U=s.zeros(44,len(cc));iv=A0[rr,cc].inv()
    for j,i in enumerate(cc):U[i,:]=iv[j,:]
    R=(A*C).applyfunc(s.cancel);ct=(tc-A*U*tc[rr,:]).applyfunc(s.cancel);full=R.row_join(ct);rf,piv=full.to_DM().rref();rf=rf.to_Matrix();free=[i for i in range(full.cols) if i not in piv];N=s.zeros(full.cols,len(free))
    for j,i in enumerate(free):
        N[i,j]=1
        for r,pc in enumerate(piv):N[pc,j]=-rf[r,i]
    assert (full.to_DM()*N.to_DM()).is_zero_matrix;_,lpiv=R.to_DM().rref();actual=(len(cc)+len(lpiv),len(cc)+len(piv),N[-1,:]==s.zeros(1,N.cols));assert actual==(43,44,True)
    # Retain all top70 compatibility rows of the homogenizer column.
    return dict(beta_kernel=enc(BK),q_point=list(map(str,q)),chart_columns=cols,leading_jacobian=enc(A),t_column=list(map(str,tc)),first_rank=len(cc),first_pivot_rows=rr,first_pivot_columns=cc,first_kernel=enc(C),residual_matrix=enc(R),residual_t_column=list(map(str,ct)),residual_kernel=enc(N),free_columns=free,leading_rank=actual[0],full_rank=actual[1],homogenizer_zero_on_tangent=actual[2],top70_t_compatibility_zero=(ct[:70,:]==s.zeros(70,1)))

def normalizations():
 a,b,c,d=s.symbols('a b c d');V=s.Matrix(X[1:]);out={}
 for kind in ['M0','M1']:
  Q=tensor(kind,a,b,c);Q[1,3]+=2*d*V[0]*V[2];Q[3,1]-=2*d*V[0]*V[2]
  ell=-d/c if kind=='M0' else -2*d/(2*c+1);S=s.eye(4);S[3,2]=ell;newa=a+b*ell;A=-newa/(c+1) if kind=='M0' else -newa/c;B=-b/(c-1) if kind=='M0' else -b/(c+1);T=s.Matrix([[1,0,0,0],[0,1,A,-B],[B,0,1,0],[A,0,0,1]]);C=S*T
  assert (C*J[1:,1:]*C.T-J[1:,1:]).applyfunc(s.cancel)==s.zeros(4);result=C.inv()*Q.subs(dict(zip(V,C*V)),simultaneous=True)*C.inv().T
  assert (result-tensor(kind,0,0,c)).applyfunc(s.cancel)==s.zeros(4)
  E=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]]);paired=-c if kind=='M0' else -c-1;Qc=tensor(kind,0,0,c)
  assert E*J[1:,1:]*E.T==J[1:,1:];assert (E.inv()*Qc.subs(dict(zip(V,E*V)),simultaneous=True)*E.inv().T-tensor(kind,0,0,paired)).applyfunc(s.expand)==s.zeros(4)
  out[kind]=dict(matrix=enc(C),denominator_open='c*(c-1)*(c+1)' if kind=='M0' else 'c*(c+1)*(2*c+1)',paired_c=str(paired),paired_exchange=enc(E))
 return out

def proper_kernel(K,b,bv,ar,f,var,val):
 rows=list(K.T.rref()[1]);iv=K[rows,:].inv();N0=s.Matrix.hstack(*K.T.nullspace());grad=s.Matrix([s.diff(f,bb) for bb in b]);ext=s.zeros(16,1);eg=iv.T*grad
 for i,j in enumerate(rows):ext[j]=eg[i]
 N=N0.row_join(ext).subs(var,val);W=(ar.subs(var,val).T*N).applyfunc(s.cancel);den=s.factor(s.lcm([s.denom(v) for v in W]));WP=W.applyfunc(lambda v:s.cancel(v*den));dm,piv=WP.to_DM().rref(method='FF');rr=dm.to_Matrix();free=[i for i in range(7) if i not in piv];C=s.zeros(7,len(free))
 for j,i in enumerate(free):
  C[i,j]=1
  for h,c in enumerate(piv):C[c,j]=-rr[h,i]
 assert (WP.to_DM()*C.to_DM()).is_zero_matrix;L=(N*C).applyfunc(s.cancel);bp=bv.subs(var,val);params=[bb for bb in b if bb!=var];D=bp.jacobian(params).row_join(ar.subs(var,val));assert (D.T*L).applyfunc(s.cancel)==s.zeros(19,len(free));sample=dict(zip(b,range(1,11)));assert L.subs(sample).to_DM().rank()==len(free)
 return dict(image_dimension_upper=16-len(free),left_kernel=enc(L),small_matrix=enc(WP),small_kernel=enc(C),cleared_denominator=str(den))

def case(m,kind,c,source_main,src):
 b=s.symbols('b0:10');K=s.Matrix(src['beta_kernel']).applyfunc(s.sympify);bv=K*s.Matrix(b);u,act=beta_actions(m);ar=act.subs(dict(zip(u,bv)));one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];sample=dict(zip(b,range(1,11)));E=s.Matrix(src['residual_matrix']).applyfunc(s.sympify).row_join(s.Matrix(src['residual_t_column']).applyfunc(s.sympify));num=E.subs(sample);pc=list(num.to_DM().rref()[1]);pr=list(num[:,pc].T.to_DM().rref()[1]);dm=E[pr,pc].to_DM();det=s.factor(dm.domain.to_sympy(dm.det()));assert len(pc)+src['first_rank']==44;src.update(full_minor_rows=pr,full_minor_columns=pc,full_minor=str(det));Q=tensor(kind,0,0,c);q=s.Matrix(coords(m,Q))[two,:];fi=next(i for i,v in enumerate(q) if v);V=s.Matrix(X[1:]);qd=[]
 for H in sp_basis():
  dQ=sum(((H*V)[i]*Q.diff(V[i]) for i in range(4)),s.zeros(4))-H*Q-Q*H.T;dq=s.Matrix(coords(m,dQ))[two,:];dq-=q*dq[fi]/q[fi];qd.append(dq)
 orbit=s.Matrix.hstack(*qd);stab=orbit.nullspace();assert orbit.to_DM().rank()==7 and len(stab)==3;rows=list(K.T.rref()[1]);iv=K[rows,:].inv();actions=[]
 for st in stab:
  db=ar*st;v=iv*db[rows,:];assert (K*v-db).applyfunc(s.expand)==s.zeros(16,1);actions.append(v)
 D=K.row_join(ar);joint=D.col_join(s.zeros(45,10).row_join(orbit));bd=minor_witness(D.subs(sample));jd=minor_witness(joint.subs(sample));assert bd['rank']==16 and jd['rank']==17;factors=[]
 for f,e in s.factor_list(det)[1]:
  quot=[s.cancel(sum(s.diff(f,bb)*v[i] for i,bb in enumerate(b))/f) for v in actions];assert all(not (z.free_symbols & set(b)) for z in quot);var=min((bb for bb in b if s.degree(f,bb)==1),key=lambda bb:(len(s.Poly(s.diff(f,bb),*b).terms()),len(s.diff(f,bb).free_symbols)));coeff=s.diff(f,var);const=f.subs(var,0);assert s.gcd(coeff,const)==1;value=s.solve(f,var)[0];bp=bv.subs(var,value);params=[bb for bb in b if bb!=var];Ds=bp.jacobian(params).row_join(ar.subs(var,value));numeric=Ds.subs(sample);rank=numeric.to_DM().rank();result=dict(factor=str(f),exponent=int(e),variable=str(var),solved_value=str(value),primitive_linear=True,stabilizer_weights=list(map(str,quot)),sample_beta_rank=rank)
  if rank<16:
   cert=proper_kernel(K,b,bv,ar,f,var,value);assert cert['image_dimension_upper']==rank and rank<=15;result.update(type='proper_beta_image',certificate=cert)
  elif f==b[0]:
   result.update(type='reviewed_main_junction',known_k='1/6' if c==-s.Rational(1,2) else '-11/12',beta_minor=minor_witness(numeric))
  else:
   point=s.Matrix(coords(m,Q));beta=bp.subs(sample)
   for i,v in zip(one,beta):point[i]=v
   loc=calculate(m,True,point);assert loc['leading_jacobian_rank']==43 and loc['homogenized_jacobian_rank']==44;result.update(type='smooth_curve_artifact_divisor',beta_minor=minor_witness(numeric),local_data=loc,full_jacobian_minor=minor_witness(s.Matrix(loc['jacobian']).applyfunc(s.Rational)))
  factors.append(result)
 P,R,A,B,U,V,KK,L=s.symbols('P R A B U V K L');sub={P:bv[1],R:bv[3],A:bv[4],B:bv[13],U:bv[6],V:bv[9],KK:bv[14]+5*bv[9],L:5*bv[6]-bv[15]};ma=s.factor(s.sympify(source_main[kind]['solution_a']).subs(sub,simultaneous=True));mb=s.factor(s.sympify(source_main[kind]['solution_b']).subs(sub,simultaneous=True));assert ma.subs(b[0],0)==0 and mb==0;main_c=2*(2*bv[4]-bv[13])/(3*bv[13]) if kind=='M0' else -(2*bv[4]-bv[13])/(6*bv[13])-s.Rational(1,2);assert s.cancel(main_c-c)==0
 return dict(kind=kind,c=str(c),universal=src,beta=list(map(str,bv)),orbit_dimension=7,stabilizer=enc(s.Matrix.hstack(*stab)),stabilizer_actions=[list(map(str,v)) for v in actions],beta_minor=bd,joint_minor=jd,incidence_dimension=17,factors=factors,main_intersection=dict(a=str(ma),b=str(mb),condition='b0=0',main_c=str(c)))

def run():
 st=time.time();m=S61();main_path=ROOT/'reports/S61_MINIMAL_MAIN_VERTEX_FORMULA_20261101.json';main=json.loads(main_path.read_text())['cases'];cases=[];pre={}
 for kind,c in [('M0',s.Rational(4,5)),('M0',-s.Rational(1,2)),('M1',-s.Rational(7,10))]:
  pre[kind,c]=universal(m,kind,c);print(kind,c,'universal verified',flush=True)
 for kind,c in [('M0',s.Rational(4,5)),('M0',-s.Rational(1,2)),('M1',-s.Rational(7,10))]:
  cases.append(case(m,kind,c,main,pre[kind,c]));print(kind,c,'verified',flush=True)
 normal=normalizations();source=ROOT/'reports/S61_MINIMAL_VERTEX_EXHAUSTION_20261101.json';old=json.loads(source.read_text());covered=[]
 for fam in old['families']:
  kind=fam['label']
  for tail in fam['rational_tails']:
   if tail['axis'] not in [2,4]:continue
   c=s.Rational(tail['coordinate_formulas']['g3']);matched=[z for z in cases if z['kind']==kind and (c==s.Rational(z['c']) or c==(-s.Rational(z['c']) if kind=='M0' else -s.Rational(z['c'])-1))]
   if matched:covered.append(dict(kind=kind,k=tail['k'],axis=tail['axis'],c=str(c)))
 assert len(covered)==6;paths=[Path(__file__),main_path,source,ROOT/'code/s61_model.py',ROOT/'code/s61_minimal_baseline_thickenings_20261101.py',ROOT/'code/s61_minimal_rankone_A1_removal_20261101.py',ROOT/'code/s61_minimal_infinity_local_data_20261101.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'reports/S61_MINIMAL_JUNCTION_GENERIC_STRATA_20261101.json',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
 return dict(id='S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL_20261101',all_passed=True,field='QQ and QQ(b0,...,b9)',cases=cases,normalizations=normal,covered_presentations=covered,conclusion='All finite points of the six specified g2/g4 tail presentations are removed over sufficiently general beta.',limitations=['Their projective parameter-infinity nonsquare binary endpoints remain separate.','The other twelve rational-tail presentations are not covered.','No full generic affine-fiber exhaustion.'],seconds=time.time()-st,proof_note='research/notes/S61_MINIMAL_SEMISIMPLE_TAILS_REMOVAL_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
 with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
 print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
