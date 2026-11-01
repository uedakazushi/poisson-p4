#!/usr/bin/env python3
"""Three nonbinary resonance divisors have proper first-jet images."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X
from s61_nil22_infinity_20261101 import beta_matrix
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def enc(M):return [[str(v) for v in row] for row in M.tolist()]
def verify():
    start=time.time();m=S61();x,y,z,w=X[1:];Q=s.zeros(4);Q[1,2]=x*z;Q[1,3]=-x*w;Q=Q-Q.T
    one=[i for i,a in enumerate(m.reduced['parameter_weights']) if a==1];H=s.Matrix.hstack(*beta_matrix(m,Q).nullspace());assert H.cols==10
    zz,bu,B=primitive_beta(m);vector=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vector.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv()
    acts=s.Matrix.hstack(*[(left*infinitesimal_action(vector,A)[rows,:]).applyfunc(s.expand) for A in sp_basis()]);out=[]
    for ratio,expected in [(s.Rational(1,5),6),(s.Rational(4,5),6),(s.Rational(1,2),5)]:
        constraint=s.zeros(1,16);constraint[0,one.index(13)]=1;constraint[0,one.index(48)]=-ratio
        K=H*s.Matrix.hstack(*(constraint*H).nullspace());assert K.cols==9
        L=s.Matrix.hstack(*K.T.nullspace()).T;assert L.shape==(7,16) and L*K==s.zeros(7,9)
        b=s.symbols('b0:9');values=K*s.Matrix(b);E=(L*acts).subs(dict(zip(bu,values)),simultaneous=True).applyfunc(s.expand)
        rr,piv=E.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==expected
        free=[j for j in range(10) if j not in piv];N=s.zeros(10,len(free))
        for j,i in enumerate(free):
            N[i,j]=1
            for k,c in enumerate(piv):N[c,j]=-rr[k,i]
        assert (E.to_DM()*N.to_DM()).is_zero_matrix and N[free,:]==s.eye(len(free))
        sample={bb:s.Integer(j*j+j+2) for j,bb in enumerate(b)};Es=E.subs(sample);cp=list(Es.to_DM().rref()[1]);rp=list(Es[:,cp].T.to_DM().rref()[1]);assert len(cp)==expected;minor=Es.extract(rp,cp).det();assert minor!=0
        out.append(dict(ratio=str(ratio),beta_kernel=enc(K),quotient_projection=enc(L),normal_action=enc(E),normal_action_kernel=enc(N),kernel_free=free,normal_rank=expected,first_jet_image_dimension=9+expected,sample={str(k):str(v) for k,v in sample.items()},minor_rows=rp,minor_columns=cp,minor=str(minor)))
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/s61_nil22_infinity_20261101.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101',all_passed=True,canonical_Q='d_y wedge (xz d_z-xw d_w)',cases=out,conclusion='The u14/u49=1/5,4/5,1/2 nonbinary strata have first-jet images of dimensions15,15,14; none occurs over sufficiently general beta.',scope='These three entire Sp4-saturated linear slices, including u49=0 intersections. Not all nonbinary rank-drop loci.',proof_note='research/notes/S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101.md',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,dimensions=[x['first_jet_image_dimension'] for x in r['cases']],seconds=r['seconds'])))
if __name__=='__main__':main()
