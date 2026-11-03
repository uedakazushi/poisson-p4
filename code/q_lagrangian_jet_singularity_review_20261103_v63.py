#!/usr/bin/env python3
"""Independent full12-jet classification, Pfaffian singularity and six sources."""
import argparse
import itertools
import json
from pathlib import Path
import sympy as s
from additional_rank4_families_20261103_v62 import families
from elliptic_curve_invariants_20261101 import curve_invariants_from_beta
from s61_model import S61, X, J
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def skew(entries):
    M=s.zeros(4)
    for (i,j),v in entries.items():M[i,j]=v;M[j,i]=-v
    return M


def pf(M):
    return s.expand(M[0,1]*M[2,3]-M[0,2]*M[1,3]+M[0,3]*M[1,2])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'reports/Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63'
    d=read(directory/'report.json');fresh=read(args.replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    assert (directory/'certificate.json.gz').read_bytes()==(args.replay/'certificate.json.gz').read_bytes()
    for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    c=read(directory/'certificate.json.gz');assert sha(directory/'certificate.json.gz')==d['certificate_sha256']
    y=s.Matrix(s.symbols('z0:4'));x,yy,z,w=y;j=J[1:5,1:5]
    pairs=list(itertools.combinations(range(4),2));v=s.symbols('v0:12')
    # Start with ALL12 coefficients satisfying the geometric shape, before
    # imposing compatibility; this checks completeness of the ten-jet ansatz.
    raw=skew({(0,1):v[0]*x+v[1]*z,(0,3):v[2]*x+v[3]*z,
              (1,2):v[4]*x+v[5]*z,(1,3):v[6]*x+v[7]*yy+v[8]*z+v[9]*w,
              (2,3):v[10]*x+v[11]*z})
    conditions=[s.expand(sum(j[i,h]*s.diff(raw[k,l],y[h])+j[k,h]*s.diff(raw[l,i],y[h])+j[l,h]*s.diff(raw[i,k],y[h])for h in range(4)))for i,k,l in itertools.combinations(range(4),3)]
    compatibility=s.Matrix(conditions).jacobian(v)
    assert compatibility.rank()==2
    solutions=s.solve(conditions,[v[7],v[9]])
    assert solutions=={v[7]:v[1]-v[2],v[9]:-v[5]-v[10]}
    compatible=raw.subs(solutions,simultaneous=True)
    free=[a for a in v if a not in solutions]
    alpha=s.Matrix([s.diff(compatible[0,1]+compatible[2,3],a)/5 for a in y])
    tangent=j*alpha
    primitive=(compatible-2*alpha.dot(y)*j+y*tangent.T-tangent*y.T).applyfunc(s.expand)
    a=s.symbols('a b c d e f g h');aa,b,cc,dd,e,f,g,h=a
    universal=skew({(0,1):aa*x+b*z,(0,3):cc*x+dd*z,(1,2):e*x+f*z,
                   (1,3):g*x+h*z+(b-cc)*yy+(aa-f)*w,(2,3):-aa*x-b*z})
    def encode(M):return s.Matrix([s.diff(M[i,k],t)for i,k in pairs for t in y])
    image=encode(primitive).jacobian(free);basis=encode(universal).jacobian(a)
    assert image.rank()==basis.rank()==image.row_join(basis).rank()==8
    assert universal==s.Matrix(c['universal_primitive_linear_jet'])
    beta=-j*universal*j
    pp=s.symbols('p0:5');p0,p1,p2,p3,p4=pp
    Pl=skew({(0,1):p0,(0,2):p1,(0,3):p2,(1,2):p3,(1,3):p4,(2,3):-p0})
    Star=skew({(0,1):-p0,(0,2):-p4,(0,3):p3,(1,2):p2,(1,3):-p1,(2,3):p0})
    cov=s.Matrix([sum(s.diff(beta[i,k],t)*Pl[i,k]for i,k in pairs)for t in y])
    Phi=Star.row_join(cov).col_join((-cov.T).row_join(s.zeros(1)))
    F=s.Matrix([(-1)**i*pf(Phi.minor_submatrix(i,i))for i in range(5)])
    assert (F-(Pl*cov).col_join(s.Matrix([pf(Pl)]))).applyfunc(s.expand)==s.zeros(5,1)
    point=dict(zip(pp,[0,0,0,0,1]))
    assert F.subs(point)==s.zeros(5,1)
    jac=F.jacobian(pp).subs(point)
    assert jac==s.Matrix(c['entire_Jacobian_at_singular_point'])
    assert jac[:,1][4]==-1 and all(jac[:,i]==s.zeros(5,1)for i in [0,2,3,4])
    assert s.expand(pf(Pl).diff(p1).subs(point))==-1
    # General value normalization, preserving the vertical plane, with ALL5
    # value parameters rather than just one rational value witness.
    q=s.symbols('t0:5');A=s.Matrix([[q[0],q[1]],[q[2],q[3]]]);B=s.Matrix([[0,q[4]],[-q[4],0]])
    C=s.zeros(2).row_join(A).col_join((-A.T).row_join(B))
    M=s.eye(2).row_join(s.zeros(2)).col_join((-A.inv().T*B*A.inv()/2).row_join(A.inv().T))
    standard=s.zeros(2).row_join(s.eye(2)).col_join((-s.eye(2)).row_join(s.zeros(2)))
    assert (M*C*M.T-standard).applyfunc(s.cancel)==s.zeros(4)
    assert s.cancel(M.det()-1/A.det())==0
    names=d['six_actual_component_firstjet_image_closures_satisfy_Delta_identically_zero']
    for record in families():
        name=record['component_id']
        if name not in names:continue
        P=record['P']
        affine=s.Matrix(4,4,lambda i,k:s.expand((P[i+1,k+1]-X[i+1]*P[0,k+1]+X[k+1]*P[0,i+1]).subs({X[0]:1}|dict(zip(X[1:],y)))))
        symbols=set().union(*(entry.free_symbols for entry in affine))
        local={str(symbol):symbol for symbol in symbols}
        saved=s.Matrix([[s.sympify(entry,locals=local)for entry in row]for row in c['families'][name]['full_symbolic_affine_P']])
        assert affine==saved
        assert affine[0,2]==0
        assert all(s.diff(affine[i,k],t)==0 for i,k in [(0,1),(0,3),(2,1),(2,3)]for t in [yy,w])
    # Actual positive beta from physical normalized reconstruction; do not
    # rely on a coordinate-name substitution between z1..4 and z0..3.
    ray=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz')
    eta=s.Symbol('eta');u=[s.sympify(t).subs(eta,1)for t in ray['S61']]
    model=S61();P=model.reconstruct(u)
    affine=s.Matrix(4,4,lambda i,k:s.expand((P[i+1,k+1]-X[i+1]*P[0,k+1]+X[k+1]*P[0,i+1]).subs({X[0]:1}|dict(zip(X[1:],y)))))
    L=s.Matrix(4,4,lambda i,k:sum(s.Poly(affine[i,k],*y).coeff_monomial(t)*t for t in y))
    assert L[0,1]+L[2,3]==0
    positive=(-j*L*j).applyfunc(s.expand)
    assert positive==skew({(0,2):s.Rational(27,8)*z,(0,3):s.Rational(3,2)*x,(1,3):-s.Rational(2048,27)*yy})
    c4,c6,delta,invariant=curve_invariants_from_beta(encode(positive))
    assert (c4,c6,delta,invariant)==(0,-949978046398464,-522255954073519803473068032,0)
    assert delta==(c4**3-c6**2)/1728!=0
    catalog=read(ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json')
    for name in names:assert next(r for r in catalog['records']if r['component_id']==name)['full_s61_component_dominance']
    orbit=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json')
    assert orbit['nonzero_eta_one_complex_projective_orbit_and_closure_Qa8over9']
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json.gz',
           ROOT/'code/additional_rank4_families_20261103_v62.py',ROOT/'code/elliptic_curve_invariants_20261101.py',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.json',
           ROOT/'reports/KNOWN51_COORDINATE_MAPS_20261103_v62_REVIEW.json',
           ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/certificate.json.gz',
           ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/s61_model.py',
           ROOT/'research/notes/Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63_REVIEW.md']
    result=dict(id='Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_and_complete_symbolic_certificate_match=True,
        independent_FULL12_shape_jet_compatibility_rank2_and_kernel10=True,
        primitive_projection_image_exactly_full8_parameter_family=True,
        independent_full5_Pfaffians_and_universal_singular_point_Jac_rank1=True,
        general5parameter_value_frame_preserves_vertical_plane_and_normalizes=True,
        all_six_FULL_source_affine_projections_and_eight_partial_identities_checked=True,
        universal_singularity_plus_Pfaffian_discriminant_criterion_and_dense_closure_accepted=True,
        independent_physical_positive_beta_and_exact_nonzero_discriminant=str(delta),
        all_nonzero_ray_outside_six_actual_components=names,
        universal_MPS_vanishing_or_unknown_component_count_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as output:
        json.dump(result,output,indent=2);output.write('\n')
    print('PASS universal complete Lagrangian jet singularity and six actual component exclusions')


if __name__=='__main__':main()
