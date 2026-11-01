#!/usr/bin/env python3
"""Independent coordinate-quotient audit of three nonbinary resonance images."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61,X,J
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]
def run():
    start=time.time();path=ROOT/'reports/S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101.json';src=json.loads(path.read_text());assert src['all_passed']
    m=S61();one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2]
    x,y,z,w=X[1:];Q=s.zeros(4);Q[1,2]=x*z;Q[1,3]=-x*w;Q-=Q.T;P=X[0]**2*J;P[1:,1:]+=Q;pp=m.normalize(P)['parameters'];assert all(pp[i]==0 for i in one)
    b=s.symbols('b0:16');values=[b[one.index(i)] if i in one else pp[i] for i in range(61)];equations=[]
    for eq in m.reduced['equations']:
        f=s.expand(sum(s.Rational(a,c)*s.prod(values[i] for i in mon) for mon,a,c in eq if sum(i in two for i in mon)==2))
        assert s.expand(f-sum(f.coeff(v)*v for v in b))==0;equations.append(f)
    D=s.Matrix(equations).jacobian(b);assert D.rank()==6
    constraints=s.zeros(6,16)
    for row,i in enumerate((3,4,11,12)):constraints[row,one.index(i)]=1
    constraints[4,one.index(49)]=1;constraints[4,one.index(34)]=5
    constraints[5,one.index(50)]=1;constraints[5,one.index(25)]=-5
    assert constraints.rank()==6 and D.col_join(constraints).rank()==6
    zz,bu,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(zz[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv()
    actions=s.Matrix.hstack(*[left*infinitesimal_action(vec,A)[rows,:] for A in reversed(sp_basis())]);out=[]
    for ratio,expected in [(s.Rational(1,5),6),(s.Rational(4,5),6),(s.Rational(1,2),5)]:
        c=s.zeros(1,16);c[0,one.index(13)]=1;c[0,one.index(48)]=-ratio;C=constraints.col_join(c);assert C.rank()==7
        K=s.Matrix.hstack(*C.nullspace());assert K.shape==(16,9) and C*K==s.zeros(7,9)
        # Explicit coordinate equations C replace the main certificate's arbitrary annihilator basis.
        pars=s.symbols('r0:9');E=(C*actions).subs(dict(zip(bu,K*s.Matrix(pars))),simultaneous=True).applyfunc(s.expand)
        rr,piv=E.to_DM().rref();rr=rr.to_Matrix();assert len(piv)==expected;free=[i for i in range(10) if i not in piv];N=s.zeros(10,len(free))
        for j,i in enumerate(free):
            N[i,j]=1
            for row,col in enumerate(piv):N[col,j]=-rr[row,i]
        assert (E.to_DM()*N.to_DM()).is_zero_matrix and N[free,:]==s.eye(len(free))
        sample=E.subs(dict(zip(pars,[3,5,7,11,13,17,19,23,29])));cp=list(sample.to_DM().rref()[1]);rp=list(sample[:,cp].T.to_DM().rref()[1]);minor=sample.extract(rp,cp).det(method='domain-ge');assert len(cp)==expected and minor!=0
        old=next(r for r in src['cases'] if s.Rational(r['ratio'])==ratio);assert old['first_jet_image_dimension']==9+expected
        oldK=s.Matrix(old['beta_kernel']).applyfunc(s.Rational);assert C*oldK==s.zeros(7,9) and oldK.rank()==9
        out.append(dict(ratio=str(ratio),coordinate_constraints=[[str(v) for v in row] for row in C.tolist()],normal_rank=expected,image_dimension=9+expected,independent_sample_minor=str(minor),minor_rows=rp,minor_columns=cp,reversed_sp_basis=True))
    paths=[Path(__file__),path,ROOT/'code/s61_nonbinary_resonance_exclusions_20261101.py',ROOT/'research/notes/S61_NONBINARY_RESONANCE_EXCLUSIONS_20261101.md',ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',m.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_NONBINARY_RESONANCE_EXCLUSIONS_REVIEW_20261101',all_passed=True,raw_beta_kernel_rank=6,cases=out,accepted_dimensions=[15,15,14],scope='Entire Sp4-saturated linear slices including B=0; no all-resonance exhaustion claim.',seconds=time.time()-start,input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',required=True,type=Path);a=parser.parse_args();r=run()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps({k:r[k] for k in ['all_passed','accepted_dimensions','seconds']}))
