#!/usr/bin/env python3
"""A dominant 17-dimensional binary rank-drop incidence, over QQ."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from itertools import combinations, combinations_with_replacement
import sympy as s
from s61_model import S61
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS, sp_basis, infinitesimal_action
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta


def verify():
    start=time.time();root=Path(__file__).resolve().parents[1];m=S61();u=s.symbols('u1:62')
    one=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==1]
    two=[i for i,w in enumerate(m.reduced['parameter_weights']) if w==2];qi={v:i for i,v in enumerate(two)}
    constraints={u[i-1]:s.S.Zero for i in [12,13,35,50]}
    beta={i:u[i].subs(constraints) for i in one};fixed=qi[29];q0=s.zeros(45,1);q0[fixed]=1
    extra_path=root/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    equations=m.reduced['equations']+[json.loads(extra_path.read_text())['universal_new_equation']]
    D=s.zeros(135,45);tc=s.zeros(135,1)
    for r,eq in enumerate(equations):
        value=0
        for mon,n,d in eq:
            qmon=[qi[i] for i in mon if i in qi]
            if not qmon:continue
            coefficient=s.Rational(n,d)*s.prod(beta[i] for i in mon if i not in qi)
            if len(qmon)==2:
                i,j=qmon;value+=coefficient*q0[i]*q0[j]
                D[r,i]+=coefficient*q0[j];D[r,j]+=coefficient*q0[i]
            else:tc[r,0]+=coefficient*q0[qmon[0]]
        assert s.expand(value)==0
    K=s.Matrix.hstack(*D[:70,:].nullspace());assert K.shape==(45,29)
    assert D[:70,:]*K==s.zeros(70,29)
    E=(D[70:,:]*K).applyfunc(s.expand);rr,piv=E.to_DM().rref();rr=rr.to_Matrix()
    assert len(piv)==17
    free=[i for i in range(29) if i not in piv];N=s.zeros(29,12)
    for j,col in enumerate(free):
        N[col,j]=1
        for i,p in enumerate(piv):N[p,j]=-rr[i,col]
    assert (E*N).applyfunc(s.cancel)==s.zeros(65,12) and N[free,:]==s.eye(12)
    sample={u[i]:s.Integer(j+2) for j,i in enumerate(one)};sample.update(constraints)
    chart=[i for i in range(45) if i!=fixed];full=D[:,chart].row_join(tc)
    minors=[]
    for A,expected in [(D,33),(full,34)]:
        A=A.subs(sample);cols=list(A.to_DM().rref()[1]);rows=list(A[:,cols].T.to_DM().rref()[1])
        assert len(cols)==expected
        determinant=A.extract(rows,cols).det(method='domain-ge');assert determinant
        minors.append(dict(rows=rows,columns=cols,value=str(determinant)))
    z,uu,form=primitive_beta(m);tensor=s.Matrix([form[i,j].diff(z[a]) for (i,j),a in COORDS])
    B=tensor.jacobian(uu);br=list(B.T.rref()[1]);left=B[br,:].inv();sp=sp_basis()
    acts=s.Matrix.hstack(*[(left*infinitesimal_action(tensor,A)[br,:]).applyfunc(s.expand) for A in sp])
    for j,A in enumerate(sp):
        assert B*acts[:,j]==infinitesimal_action(tensor,A)
    selected=[uu.index(v) for v in constraints];normal_action=acts[selected,:].subs(sample)
    ac=list(normal_action.to_DM().rref()[1]);assert len(ac)==4
    beta_minor=normal_action[:,ac].det();assert beta_minor==22680
    Q=s.zeros(4);Q[1,2]=z[0]*z[3];Q[2,1]=-Q[1,2]
    mons=[z[i]*z[j] for i,j in combinations_with_replacement(range(4),2)]
    def vector(M):return s.Matrix([s.Poly(s.expand(M[i,j]),*z).coeff_monomial(mon) for i,j in combinations(range(4),2) for mon in mons])
    qvec=vector(Q);at=next(i for i,v in enumerate(qvec) if v)
    qacts=[]
    for A in sp:
        # Pull back the bivector, matching the pullback convention for beta.
        v=vector(-A*Q-Q*A.T+sum(((A*s.Matrix(z))[i]*Q.diff(z[i]) for i in range(4)),s.zeros(4)))
        qacts.append(v-qvec*v[at]/qvec[at])
    qacts=s.Matrix.hstack(*qacts);stabilizer=qacts.nullspace();assert len(stabilizer)==5
    for h in stabilizer:
        assert (acts*h)[selected,:].subs(constraints)==s.zeros(4,1)
    # The connected projective stabilizer preserves the twelve-dimensional
    # beta slice, so every joint-image fiber has dimension at least five.
    freebeta=[i for i in range(16) if i not in selected];joint=s.zeros(76,12)
    for j,i in enumerate(freebeta):joint[i,j]=1
    joint=joint.row_join(acts.subs(sample).col_join(qacts))
    jc=list(joint.to_DM().rref()[1]);jr=list(joint[:,jc].T.to_DM().rref()[1]);assert len(jc)==17
    joint_minor=joint.extract(jr,jc).det();assert joint_minor==2160
    ff=form.subs(sample);invariants=curve_invariants_from_beta([ff[i,j].diff(z[k]) for i,j in combinations(range(4),2) for k in range(4)])
    assert invariants[2]!=0
    a,b,d=s.symbols('a b d');v1=s.Matrix([a,1,0,b]);v2=s.Matrix([-b,0,1,d])
    sec=s.Matrix([(v1.T*ff.subs(dict(zip(z,v1)))*v2)[0],(v1.T*ff.subs(dict(zip(z,v2)))*v2)[0]])
    cj=sec.jacobian([a,b,d]).subs({a:0,b:0,d:0});assert cj.rank()==2
    enc=lambda M:[[str(v) for v in row] for row in M.tolist()]
    inputs=[Path(__file__),root/'code/s61_model.py',root/'code/primitive_jet_elliptic_curve_20261101.py',
            root/'code/weight1_structure.py',root/'code/elliptic_curve_invariants_20261101.py',extra_path,
            root/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md']
    return dict(id='S61_BINARY_RANK_DROP_FAMILY_20261101',all_passed=True,
        field='QQ(twelve free first-jet coordinates)',normal_form='Q=xw*d_y wedge d_z',
        first_jet_constraints=list(map(str,constraints)),all135_leading_values_zero=True,
        constant_first70_rank=16,residual_rank=17,constant_kernel=enc(K),residual_kernel=enc(N),residual_free_columns=free,
        generic_leading_jacobian_rank=33,generic_full_homogenized_jacobian_rank=34,
        generic_full_projective_tangent_dimension=11,sample_minors=minors,
        sample={str(k):str(v) for k,v in sample.items()},primitive_action=enc(acts),
        beta_normal_action=enc(normal_action),beta_dominance_columns=ac,beta_dominance_minor=str(beta_minor),
        projective_Q_stabilizer_dimension=5,stabilizer_preserves_beta_slice=True,
        joint_differential_rank=17,joint_minor_rows=jr,joint_minor_columns=jc,joint_minor=str(joint_minor),
        joint_image_dimension=17,first_jet_image_dimension=16,generic_fixed_beta_rank_drop_family_dimension=1,
        canonical_curve_derivative=enc(cj),c4=str(invariants[0]),c6=str(invariants[1]),discriminant=str(invariants[2]),
        proof_note='research/notes/S61_BINARY_RANK_DROP_FAMILY_20261101.md',
        limitations=['One dominant rank-drop stratum, not a full determinantal-locus classification.',
                    'No completed-local-ring or t-saturation conclusion on this stratum.',
                    'This stratum lies within the binary infinity family and is not claimed as a new ambient component.'],
        input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},seconds=time.time()-start)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=verify()
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:result[k] for k in ['all_passed','generic_leading_jacobian_rank','generic_full_homogenized_jacobian_rank','joint_image_dimension','beta_dominance_minor','seconds']}))


if __name__=='__main__':main()
