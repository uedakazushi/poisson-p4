#!/usr/bin/env python3
"""Universal square-zero formal ring on the double-root nonbinary open."""
import argparse,json,time,hashlib
from pathlib import Path
from itertools import combinations_with_replacement
import sympy as s
from s61_model import S61,X,J as J0
from s61_flag_transverse_obstruction_20261101 import specialize
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def verify():
    start=time.time();m=S61();u=s.symbols('u1:62');b={i:u[i] for i,w in enumerate(m.reduced['parameter_weights']) if w==1};b.update({3:0,4:0,11:0,12:0,49:-5*u[34],50:5*u[25]})
    two,F=specialize(m,b,True);fixed=two.index(45);cols=[j for j in range(45) if j!=fixed];J=s.zeros(135,45)
    for i,f in enumerate(F):
     for mon,c in f.items():
      if len(mon)==2:
       a,z=mon
       if a!=fixed and z==fixed:J[i,cols.index(a)]-=2*c
       if z!=fixed and a==fixed:J[i,cols.index(z)]-=2*c
      elif len(mon)==1 and mon[0]==fixed:J[i,44]-=2*c

    # Fraction-field RREF directly uses sparse low-degree structure.
    rr,piv=J.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==42
    free=[i for i in range(45) if i not in piv];T=s.zeros(45,len(free))
    for k,j in enumerate(free):
     T[j,k]=1
     for i,c in enumerate(piv):T[c,k]=-rr[i,j]
    assert (J.to_DM()*T.to_DM()).is_zero_matrix and T[-1,:]==s.zeros(1,3)
    # Independent left projection using pivot rows; calculate only three quadratic columns.
    lr,lp=J.T.to_DM().rref();lr=lr.to_Matrix();lf=[j for j in range(135) if j not in lp];L=s.zeros(len(lf),135)
    for k,j in enumerate(lf):
     L[k,j]=1
     for i,c in enumerate(lp):L[k,c]=-lr[i,j]
    assert (L.to_DM()*J.to_DM()).is_zero_matrix
    minor=J[list(lp),list(piv)];dm=minor.to_DM();det=s.factor(dm.domain.to_sympy(dm.det()))
    n=T.cols;mons=list(combinations_with_replacement(range(n),2));Q=s.zeros(135,len(mons))
    for r,f in enumerate(F):
     for mon,c in f.items():
      if fixed in mon:continue
      a,z=[cols.index(v) for v in mon]+[44]*(2-len(mon))
      for k,(i,j) in enumerate(mons):Q[r,k]+=c*(T[a,i]*T[z,j]+(T[a,j]*T[z,i] if i!=j else 0))
    A=(L.to_DM()*Q.to_DM()).to_Matrix().applyfunc(s.factor)



    assert len(J[:,:44].to_DM().rref()[1])==41
    AA,BB,UU,VV=u[13],u[48],u[25],u[34]
    h=(5*AA*u[24]-10*BB*u[24]+12*UU*u[47])/(2*UU*(5*AA-4*BB))
    normal=s.Matrix([[1,0,2*h,0,0,h*h],[0,1,0,0,h,0],[0,0,0,1,0,0]])
    assert (A.to_DM()-A[:,[0,1,3]].to_DM()*normal.to_DM()).is_zero_matrix
    assert [i for i in range(A.rows) if any(A[i,:])]==[41,46,57,73]
    ndm=A[[41,46,57],[0,1,3]].to_DM();normal_det=s.factor(ndm.domain.to_sympy(ndm.det()))
    assert s.cancel(normal_det-UU*VV*(AA-2*BB)**2*(5*AA-4*BB)**2/(4*BB*(2*AA-BB)**3))==0
    allowed={str(s.factor(f)) for f,e in s.factor_list(det)[1]};denominator_factors=set()
    for value in set(T).union(set(L)):
     denominator=s.fraction(s.cancel(value))[1]
     denominator_factors.update(str(s.factor(f)) for f,e in s.factor_list(denominator)[1])
    assert denominator_factors.issubset(allowed)
    assert all(str(s.factor(f)) in allowed for f,e in s.factor_list(s.fraction(normal_det)[0])[1])
    # Finite symplectic symmetry exchanges the two double-root descriptions.
    x,y,z,w=X[1:];vv=s.Matrix([x,y,z,w]);S=s.Matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]])
    def tensor(G):
     Q=s.zeros(4);Q[1,2]=s.diff(G,w);Q[1,3]=-s.diff(G,z);return Q-Q.T
    qz=tensor(z*z*w);qw=tensor(z*w*w)
    assert S.T*J0[1:,1:]*S==J0[1:,1:]
    assert (S.inv()*qz.subs(dict(zip(vv,S*vv)),simultaneous=True)*S.inv().T-qw).applyfunc(s.expand)==s.zeros(4)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_flag_transverse_obstruction_20261101.py',ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json',ROOT/'research/notes/S61_R0_GENERIC_EXHAUSTION_20261101.md',ROOT/'research/notes/S61_R0_UNIFORM_VERTEX_20261101.md',ROOT/'research/notes/NULL_CORRELATION_DUALITY_20261101.md',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101',all_passed=True,field='QQ(beta)',canonical_Q='d_y wedge(z^2 d_z-2zw d_w)',canonical_u46=-2,beta={str(u[i]):str(v) for i,v in b.items()},jacobian=enc(J),tangent=enc(T),left_kernel=enc(L),pivot_rows=list(lp),pivot_columns=list(piv),jacobian_minor=str(det),denominator_factors=sorted(denominator_factors),quadratic=enc(A),normal_form_linear_shift=str(h),normal_quadratic_coefficient_minor=str(normal_det),full_rank=42,leading_rank_upper_bound=41,tangent_t_zero=True,completed_ring_on_minor_open='C[[u,e1,e2]]/(e1,e2)^2',completed_ring_rank=3,homogenizer_zero=True,symmetry_matrix=enc(S),symmetric_zw2_orbit_verified=True,proof_note='research/notes/S61_NONBINARY_DOUBLE_ROOT_OPEN_20261101.md',limitations=['The explicit nonzero-minor open over sufficiently general beta, not all points of this orbit.','The A=2B endpoint meets two minimal g6 tails; other minor factors need separate checks.','No whole infinity saturation or affine exhaustion.'],input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','full_rank','leading_rank_upper_bound','completed_ring_rank','homogenizer_zero','seconds']}))
