#!/usr/bin/env python3
"""Exact number-field certificates for three double-root artifact divisors."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from sympy.polys.matrices import DomainMatrix as DM
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def run():
    start=time.time();m=S61();one=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==1];two=[i for i,v in enumerate(m.reduced['parameter_weights']) if v==2];qi={i:j for j,i in enumerate(two)};u=s.symbols('u1:62');a=s.Symbol('a');src=json.loads((ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json').read_text());JM=s.Matrix(src['jacobian']).applyfunc(s.sympify);sub={u[13]:a,u[48]:1,u[25]:2,u[34]:3,u[24]:2,u[26]:3,u[33]:5,u[35]:7,u[36]:11,u[47]:13};JM=JM.xreplace(sub);beta={i:u[i] for i in one};beta.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]});beta={i:s.sympify(v).xreplace(sub) for i,v in beta.items()};F=[]
    for eq in m.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]:
     f={}
     for mon,nu,de in eq:
      mm=tuple(qi[i] for i in mon if i in qi);c=s.Rational(nu,de)*s.prod(beta[i] for i in mon if i in beta);f[mm]=f.get(mm,0)+c
     F.append({mm:s.expand(v) for mm,v in f.items() if v})
    fixed=qi[45];cols=[i for i in range(45) if i!=fixed];x,y,z,w=X[1:];Y=s.Matrix([x,y,z,w]);Q=s.zeros(4);Q[1,2]=z*z;Q[1,3]=-2*z*w;Q-=Q.T
    # Independent orbit and beta action matrices.
    def qcoords(q):
     pp=X[0]**2*J;pp[1:,1:]+=q;return s.Matrix(m.normalize(pp)['parameters'])[two,:]
    q0=qcoords(Q);orbit=[]
    for H in sp_basis():
     dq=qcoords(sum(((H*Y)[j]*Q.diff(Y[j]) for j in range(4)),s.zeros(4))-H*Q-Q*H.T);dq-=q0*dq[fixed]/q0[fixed];orbit.append(dq[cols,:])
    QD=s.Matrix.hstack(*orbit);stab=s.Matrix.hstack(*QD.nullspace());assert QD.rank()==7 and stab.cols==3
    zz,bu,Beta=primitive_beta(m);vec=s.Matrix([Beta[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rr=list(basis.T.rref()[1]);left=basis[rr,:].inv();acts=s.Matrix.hstack(*[left*infinitesimal_action(vec,H)[rr,:] for H in sp_basis()]);Hconds=s.Matrix([bu[one.index(i)] for i in [3,4,11,12]]+[bu[one.index(49)]+5*bu[one.index(34)],bu[one.index(50)]-5*bu[one.index(25)]]).jacobian(bu);HH=s.Matrix.hstack(*Hconds.nullspace());assert HH.cols==10
    result=[]
    for f in [a+1,20*a*a-59*a+65,11*a**3-39*a*a+51*a-7]:
     timer=time.time()
     assert s.Poly(f,a).is_irreducible and s.gcd(f,s.diff(f,a))==1
     if s.degree(f,a)==1:
      root=s.solve(f,a)[0];K=s.QQ;alpha=K.from_sympy(root)
     else:root=s.CRootOf(f,0);K=s.QQ.algebraic_field(root);alpha=K.unit
     zero=K.zero;oneK=K.one
     def cv(expr):
      p=s.Poly(expr,a);return sum((K.from_sympy(c)*alpha**mon[0] for mon,c in p.terms()),zero)
     def dm(M):return DM([[cv(e) for e in row] for row in M.tolist()],M.shape,K)
     def scalar(v):
      if K==s.QQ:return str(K.to_sympy(v))
      coeff=v.to_list();return str(sum(s.QQ.to_sympy(c)*a**(len(coeff)-i-1) for i,c in enumerate(coeff)))
     def enc(D):return [[scalar(v) for v in row] for row in D.to_list()]
     def minor(D):
      cp=list(D.rref()[1]);rp=list(D.extract(range(D.shape[0]),cp).transpose().rref()[1]);value=D.extract(rp,cp).det();assert value
      return dict(rows=rp,columns=cp,determinant=scalar(value))
     def null(D):
      R,pv=D.rref();R=R.to_list();free=[i for i in range(D.shape[1]) if i not in pv];out=[[zero]*len(free) for _ in range(D.shape[1])]
      for j,i in enumerate(free):
       out[i][j]=oneK
       for r,c in enumerate(pv):out[c][j]=-R[r][i]
      return DM(out,(D.shape[1],len(free)),K)
     D=dm(JM);T=null(D);L=null(D.transpose()).transpose();assert T.shape==(45,3);assert (D*T).is_zero_matrix;ET=[[zero]*3 for _ in range(46)];tl=T.to_list()
     for j,i in enumerate(cols):ET[i]=tl[j]
     ET[45]=tl[44];assert ET[45]==[zero]*3
     mon2=list(combinations_with_replacement(range(3),2));raw=[[zero]*6 for _ in range(135)]
     for row,ff in enumerate(F):
      for mm,expr in ff.items():
       if fixed in mm:continue
       i,j=(*mm,*([45]*(2-len(mm))))
       if not any(ET[i]) or not any(ET[j]):continue
       c=cv(expr)
       for k,(r,t) in enumerate(mon2):raw[row][k]+=c*(ET[i][r]*ET[j][t]+(ET[i][t]*ET[j][r] if r!=t else zero))
     obstruction=L*DM(raw,(135,6),K);assert obstruction.rank()==3;hes=[]
     for row in obstruction.to_list():
      if not any(row):continue
      Hess=[[zero]*3 for _ in range(3)]
      for c,(r,t) in zip(row,mon2):
       if r==t:Hess[r][r]+=2*c
       else:Hess[r][t]+=c;Hess[t][r]+=c
      hes+=Hess
     HD=DM(hes,(len(hes),3),K);assert HD.rank()==2
     # H_{a}: A=aB, fixed a. Its stabilizer invariance is polynomial even in a.
     C=s.Matrix([bu[one.index(13)]-a*bu[one.index(48)]]).jacobian(bu);HK=HH*s.Matrix.hstack(*(C*HH).nullspace());bs=s.symbols('b0:9');assert HK.cols==9
     stable=(acts*stab).subs(dict(zip(bu,HK*s.Matrix(bs))),simultaneous=True).applyfunc(s.expand);assert C*stable==s.zeros(1,3) or (C*stable).applyfunc(s.expand)==s.zeros(1,3);assert (Hconds*stable).applyfunc(s.expand)==s.zeros(6,3)
     ac=acts.xreplace(dict(zip(bu,[beta[i] for i in one])));BD=dm(HH.row_join(ac));RB=dm(HK.row_join(ac));joint=RB.vstack(dm(s.zeros(44,9).row_join(QD)));assert BD.rank()==16 and RB.rank()==16 and joint.rank()==16
     fib=dm(s.zeros(44,10).row_join(QD))*null(BD);assert fib.rank()==1
     fc=fib.rref()[1][0];td=fib.extract(range(44),[fc]).vstack(DM([[zero]],(1,1),K));tr=list(T.transpose().rref()[1]);tv=T.extract(tr,range(3)).inv()*td.extract(tr,[0]);assert T*tv==td and (HD*tv).is_zero_matrix and tv.rank()==1
     assert D.extract(range(135),range(44)).rank()==41
     rec={'factor':str(f),'field':'QQ[a]/('+str(f)+')','beta_sample':{str(u[i]):str(v) for i,v in beta.items()},'full_rank':42,'leading_rank':41,'full_minor':minor(D),'tangent_basis':enc(T),'quadratic_rank':3,'quadratic_nonzero_rows':{str(i):list(map(scalar,row)) for i,row in enumerate(obstruction.to_list()) if any(row)},'quadratic_minor':minor(obstruction),'hessian_rank':2,'hessian_minor':minor(HD),'fixed_beta_curve_tangent':enc(tv),'restricted_beta_rank':16,'restricted_beta_minor':minor(RB),'joint_rank':16,'curve_rank':1,'stabilizer_preserves_entire_linear_slice':True,'point_completed_ring':'C[[u,e1,e2]]/(e1,e2)^2, t=0','stratum_dimension':16,'whole_generic_beta_stratum_removed':True,'seconds':time.time()-timer};result.append(rec);print(json.dumps({k:rec[k] for k in ['factor','full_rank','leading_rank','quadratic_rank','restricted_beta_rank','seconds']}),flush=True)
    paths=[Path(__file__),ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.json',ROOT/'reports/S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101.json',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',m.source/'symplectic_slice_reduced_equations.json',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.md',ROOT/'research/notes/S61_NONBINARY_DOUBLE_ROOT_REVIEW_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_ARTIFACTS_20261101',all_passed=True,cases=result,projective_orbit_dimension=7,projective_stabilizer_dimension=3,projective_Q_differential=[[str(v) for v in row] for row in QD.tolist()],stabilizer=[[str(v) for v in row] for row in stab.tolist()],scope='The six complex ratio strata A/B=-1 or roots of P2/P3 only; other double-root divisors and the whole orbit remain separate.',proof_note='research/notes/S61_NONBINARY_DOUBLE_ROOT_ARTIFACTS_20261101.md',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({'all_passed':True,'seconds':r['seconds']}))
