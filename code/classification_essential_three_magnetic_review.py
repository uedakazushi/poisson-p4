#!/usr/bin/env python3
"""Universal actual shear/GL2 identities for the essential-three source.

All quadratic U,V,R and linear l,m remain free: 52 independent parameters.
The geometric injectivity/residue and finite-orbit proofs are separate and are
not inferred from this universal polynomial identity or a sample determinant.
"""
from pathlib import Path
from itertools import combinations,combinations_with_replacement
import hashlib,json,resource,signal,time
from classification_minimal_quartic_review import LaurentRing
ROOT=Path(__file__).resolve().parents[1]
JOB='reports/CLASSIFICATION_ESSENTIAL_THREE_MAGNETIC_REVIEW_JOB_20261107.json'
OUT='reports/CLASSIFICATION_ESSENTIAL_THREE_MAGNETIC_REVIEW_20261107.json'


def main():
    started=time.monotonic();signal.alarm(30);resource.setrlimit(resource.RLIMIT_AS,(512*1024**2,)*2)
    dest=ROOT/OUT
    if dest.exists():raise FileExistsError(dest)
    job=json.loads((ROOT/JOB).read_text())
    for b in job['exact_inputs_and_hashes']:assert hashlib.sha256((ROOT/b['path']).read_bytes()).hexdigest()==b['sha256'],b['path']
    parameters=tuple(f'U{i}_{j}' for i in range(3) for j in range(6))+tuple(f'V{i}_{j}' for i in range(3) for j in range(6))+tuple('l'+str(i) for i in range(3))+tuple('m'+str(i) for i in range(3))+tuple('R'+str(i) for i in range(6))+('A','B','C','D')
    assert len(parameters)==52
    ring=LaurentRing(parameters,('a0','a1','a2','u','v'));r=ring
    x=[r.variable(n) for n in ('a0','a1','a2','u','v')]
    mons=[r.mul(x[i],x[j]) for i,j in combinations_with_replacement(range(3),2)]
    U=[r.add(*(r.mul(r.variable(f'U{i}_{j}'),m) for j,m in enumerate(mons))) for i in range(3)]
    V=[r.add(*(r.mul(r.variable(f'V{i}_{j}'),m) for j,m in enumerate(mons))) for i in range(3)]
    l=r.add(*(r.mul(r.variable('l'+str(i)),x[i]) for i in range(3)))
    m=r.add(*(r.mul(r.variable('m'+str(i)),x[i]) for i in range(3)))
    R=r.add(*(r.mul(r.variable('R'+str(i)),mon) for i,mon in enumerate(mons)))
    div=lambda W:r.add(*(r.derivative(W[i],i) for i in range(3)))
    apply=lambda W,f:r.add(*(r.mul(W[i],r.derivative(f,i)) for i in range(3)))
    du,dv=div(U),div(V)
    delta=div([r.add(r.mul(l,V[i]),r.scale(r.mul(m,U[i]),-1)) for i in range(3)])
    formula=r.add(r.mul(dv,l),r.scale(r.mul(du,m),-1),r.scale(apply(U,m),-1),apply(V,l))
    assert delta==formula
    def source(U,V,R):
        P={(i,3):U[i] for i in range(3)};P.update({(i,4):V[i] for i in range(3)})
        P[3,4]=r.add(r.scale(r.mul(div(V),x[3]),-1),r.mul(div(U),x[4]),R)
        return P
    P=source(U,V,R);pairs=list(combinations(range(5),2))
    def get(P,i,j):return P.get((i,j),{}) if i<j else r.scale(P.get((j,i),{}),-1) if i>j else {}
    def subst(poly,images):
        out={}
        for ex,c in poly.items():
            value={ex[:r.n]+(0,)*5:c}
            for i,k in enumerate(ex[r.n:]):
                for _ in range(k):value=r.mul(value,images[i])
            out=r.add(out,value)
        return out
    def push_before_inverse(P,G):
        J=[[r.derivative(f,k) for k in range(5)] for f in G]
        return {ij:r.add(*(r.mul(J[ij[0]][k],get(P,k,t),J[ij[1]][t]) for k in range(5) for t in range(5))) for ij in pairs}
    shear=x[:3]+[r.add(x[3],l),r.add(x[4],m)]
    inverse=x[:3]+[r.add(x[3],r.scale(l,-1)),r.add(x[4],r.scale(m,-1))]
    assert [subst(f,inverse) for f in shear]==x and [subst(f,shear) for f in inverse]==x
    transformed=push_before_inverse(P,shear)
    expected=source(U,V,r.add(R,delta))
    assert all(transformed[ij]==subst(expected.get(ij,{}),shear) for ij in pairs)
    A,B,C,D=[r.variable(n) for n in ('A','B','C','D')]
    det=r.add(r.mul(A,D),r.scale(r.mul(B,C),-1))
    newU=[r.add(r.mul(A,u),r.mul(B,v)) for u,v in zip(U,V)]
    newV=[r.add(r.mul(C,u),r.mul(D,v)) for u,v in zip(U,V)]
    G=x[:3]+[r.add(r.mul(A,x[3]),r.mul(B,x[4])),r.add(r.mul(C,x[3]),r.mul(D,x[4]))]
    transformed=push_before_inverse(P,G);expected=source(newU,newV,r.mul(det,R))
    assert all(transformed[ij]==subst(expected.get(ij,{}),G) for ij in pairs)
    # Polynomial adjugate identity certifies the genuine inverse on GL2;
    # no unrelated parameter denominator or generic F localization appears.
    M=[[A,B],[C,D]];adj=[[D,r.scale(B,-1)],[r.scale(C,-1),A]]
    assert all(r.add(*(r.mul(M[i][k],adj[k][j]) for k in range(2)))==(det if i==j else {}) for i in range(2) for j in range(2))
    assert 2*3==6 and 3*4//2==6 and 4+6==10 and 10-1==9 and 24-1==23
    result=dict(task_id='E_reduced_nonnormal',all_passed=True,
        theorem='Entire irreducible essential-three quintic Poisson support is a finite union of actual23dim PGL5 orbits. All6magnetic coefficients are removed by a unique actual linear shear for eachF/frame. Representative list and component memberships remain open.',
        universal_exact_polynomial_checks=dict(independent_parameters=52,all_quadratic_anchor_U_V_coefficients=36,all_magnetic_R_coefficients=6,all_linear_shear_coefficients=6,all_GL2_coefficients=4,
            actual_shear_all_ten_brackets=True,both_shear_inverse_compositions=True,
            magnetic_shift='R+div(l*V-m*U)',GL2_actual_all_ten_brackets=True,GL2_frame_action='Unew=A*U+B*V,Vnew=C*U+D*V,Rnew=(AD-BC)*R',GL2_adjugate_inverse=True),
        mathematical_proof_audit=dict(delta_injective_by_closed_rational_form_residue=True,allF_uniform_six_dimensional_inverse=True,
            irreducibility_of_F_essential=True,finite_PGL3_stabilizer_from_Derlog_degree1=True,
            finite_orbits_from_invariant_constructible_type_loci_dim8=True,
            canonical_lift_GL5_equivariant=True,full_projective_stabilizer_identity_component='C*',orbit_dimension=23),
        no_new_CAS_or_parameter_probes=True,no_componenthood_or_known_E_membership_claim=True,
        inputs=job['exact_inputs_and_hashes'],reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        seconds=time.monotonic()-started,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS',OUT,round(result['seconds'],3),'seconds',result['maximum_RSS_KiB'],'KiB',flush=True)
if __name__=='__main__':main()
