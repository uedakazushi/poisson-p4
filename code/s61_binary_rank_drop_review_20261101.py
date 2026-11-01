#!/usr/bin/env python3
"""Independent rank/finite-pullback audit of the dominant binary rank-drop family."""
import argparse,hashlib,json,time
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis
ROOT=Path(__file__).resolve().parents[1]
def verify():
    start=time.time();source=ROOT/'reports/S61_BINARY_RANK_DROP_FAMILY_20261101.json';d=json.loads(source.read_text());m=S61();u=s.symbols('u1:62');one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1];two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={v:j for j,v in enumerate(two)};fixed=qi[29];sub={u[i-1]:0 for i in [12,13,35,50]}
    extra=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json';F=m.reduced['equations']+[json.loads(extra.read_text())['universal_new_equation']]
    D=s.zeros(135,45);tc=s.zeros(135,1)
    for row,f in enumerate(F):
        value=0
        for mon,n,den in f:
            qmon=[qi[i] for i in mon if i in qi];c=s.Rational(n,den)*s.prod(u[i].subs(sub) for i in mon if i not in qi)
            if len(qmon)==2:
                i,j=qmon
                if i==fixed:D[row,j]+=c
                if j==fixed:D[row,i]+=c
                if i==j==fixed:value+=c
            elif len(qmon)==1 and qmon[0]==fixed:tc[row]+=c
        assert s.expand(value)==0
    K=s.Matrix(d['constant_kernel']).applyfunc(s.sympify);N=s.Matrix(d['residual_kernel']).applyfunc(s.sympify);E=D[70:,:]*K
    assert K.rank()==29 and D[:70,:].rank()==16 and D[:70,:]*K==s.zeros(70,29)
    assert N[d['residual_free_columns'],:]==s.eye(12) and (E.to_DM()*N.to_DM()).is_zero_matrix
    sample={s.Symbol(k):s.Rational(v) for k,v in d['sample'].items()};chart=[i for i in range(45) if i!=fixed]
    for A,w,rank in zip([D,D[:,chart].row_join(tc)],d['sample_minors'],[33,34]):
        A=A.subs(sample);assert A.rank()==rank and A.extract(w['rows'],w['columns']).det(method='domain-ge')==s.Rational(w['value'])
    z,bu,B=primitive_beta(m);vec=s.Matrix([B[i,j].diff(z[k]) for (i,j),k in COORDS]);basis=vec.jacobian(bu);rows=list(basis.T.rref()[1]);left=basis[rows,:].inv();eps=s.Symbol('eps');acts=[];qacts=[]
    Q=s.zeros(4);Q[1,2]=z[0]*z[3];Q[2,1]=-Q[1,2];mons=[z[i]*z[j] for i,j in combinations_with_replacement(range(4),2)]
    def qvector(P):return s.Matrix([s.Poly(s.expand(P[i,j]),*z).coeff_monomial(mon) for i,j in combinations(range(4),2) for mon in mons])
    q0=qvector(Q);at=next(i for i,v in enumerate(q0) if v)
    for A in sp_basis():
        g=s.eye(4)+eps*A;inverse=s.eye(4)-eps*A;change=dict(zip(z,g*s.Matrix(z)))
        pulled=g.T*B.subs(change,simultaneous=True)*g;derivative=pulled.diff(eps).subs(eps,0)
        v=s.Matrix([derivative[i,j].diff(z[k]) for (i,j),k in COORDS]);coeff=left*v[rows,:];assert basis*coeff==v;acts.append(coeff)
        qpull=inverse*Q.subs(change,simultaneous=True)*inverse.T;dq=qvector(qpull.diff(eps).subs(eps,0));qacts.append(dq-q0*dq[at]/q0[at])
    acts=s.Matrix.hstack(*acts);qacts=s.Matrix.hstack(*qacts);assert acts==s.Matrix(d['primitive_action']).applyfunc(s.sympify)
    selected=[bu.index(v) for v in sub];normal=acts[selected,:].subs(sample);assert normal[:,d['beta_dominance_columns']].det()==22680
    stabilizer=s.Matrix.hstack(*qacts.nullspace());assert stabilizer.cols==5
    assert (acts*stabilizer)[selected,:].subs(sub)==s.zeros(4,5)
    free=[i for i in range(16) if i not in selected];joint=s.zeros(76,12)
    for j,i in enumerate(free):joint[i,j]=1
    joint=joint.row_join(acts.subs(sample).col_join(qacts));assert joint.rank()==17
    assert joint.extract(d['joint_minor_rows'],d['joint_minor_columns']).det()==2160
    paths=[Path(__file__),source,ROOT/'code/s61_binary_rank_drop_family_20261101.py',ROOT/'research/notes/S61_BINARY_RANK_DROP_FAMILY_20261101.md',extra]
    return dict(id='S61_BINARY_RANK_DROP_REVIEW_20261101',all_passed=True,independent_full135_extraction=True,independent_first_order_finite_pullback=True,leading_rank=33,full_rank=34,joint_rank=17,beta_normal_minor=22680,stabilizer_dimension=5,scope='The universal rank-drop family and dominance are accepted. Smooth-curve invariants are separately rerun by the parent; no local saturation claim.',input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},seconds=time.time()-start)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=verify()
    with a.output.open('x') as f:json.dump(r,f,indent=2);f.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=r['seconds'])))
if __name__=='__main__':main()
