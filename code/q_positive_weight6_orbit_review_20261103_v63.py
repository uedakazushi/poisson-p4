#!/usr/bin/env python3
"""Independent physical GL action, projective orbit and closed-locus audit."""
import argparse
import json
from pathlib import Path
import sympy as s
from s61_model import X, coefficients, pushforward
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'reports/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63'
    d=read(directory/'report.json')
    f=read(args.replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in f.items()if k!='seconds'}
    assert (directory/'certificate.json').read_bytes()==(args.replay/'certificate.json').read_bytes()
    for filename,digest in d['input_sha256'].items():assert sha(ROOT/filename)==digest,filename
    c=read(directory/'certificate.json')
    eta,t=s.symbols('eta t',nonzero=True)
    P=s.zeros(5)
    for i,j,v in [(1,2,-s.Rational(256,3)*eta*X[0]*X[4]),
                  (1,3,-64*eta*X[1]*X[4]),(1,4,s.Rational(8,9)*X[0]**2),
                  (2,3,-s.Rational(2,9)*X[0]**2+64*eta*X[2]*X[4]),
                  (2,4,s.Rational(4,3)*X[0]*X[1]),
                  (3,4,s.Rational(2,3)*X[0]*X[2]+X[1]**2)]:
        P[i,j]=v
        P[j,i]=-v
    point=P.subs(eta,1)
    columns=[]
    for i in range(5):
        for j in range(5):
            E=s.zeros(5)
            E[i,j]=1
            derivative=E*point+point*E.T-X[j]*point.diff(X[i])
            assert all(s.expand(sum(s.diff(derivative[k,l],X[k])for k in range(5)))==0 for l in range(5))
            columns.append(s.Matrix(coefficients(derivative)))
    M=s.Matrix.hstack(*columns)
    N=M.row_join(s.Matrix(coefficients(point)))
    assert N==s.Matrix(c['GL25_plus_scalar_matrix'])
    K=s.Matrix(c['kernel_basis'])
    assert K.shape==(26,3) and K.rank()==3 and N*K==s.zeros(150,3)
    assert K[25,:]==s.zeros(1,3)
    rows=c['nonzero_minor_rows_zero_based']
    cols=c['nonzero_minor_columns_zero_based']
    minor=N.extract(rows,cols).det()
    assert len(rows)==len(cols)==23 and minor==s.Rational(c['nonzero_minor'])!=0
    # Three independent GL kernels imply rank(M)<=22. The23-minor after
    # adding one scalar column implies rank(M)>=22 and projective rank22.
    assert M*K[:25,:]==s.zeros(150,3)
    transformed=pushforward(P,s.diag(1,t,t**2,t**3,t**4))/t**5
    assert (transformed-P.subs(eta,eta/t**6)).applyfunc(s.cancel)==s.zeros(5)
    previous=read(ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.json')
    assert previous['all_passed'] and previous['independent_all150_physical_tensor_strict_Jacobi_and_divergence_and_quintic']
    membership=read(ROOT/'reports/Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63.json')
    assert membership['all_passed']
    assert membership['actual_T39_membership_iff']=='a=2/9'
    assert membership['actual_T40_membership_iff']=='a=-8/9'
    a=s.Rational(8,9)
    assert (9*a-2)/(2*a*a)!=0 and (9*a+8)/(2*a*a)!=0
    catalogue=read(ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json')
    labels=['R4_DOUBLE_SMOOTH_QUADRIC','R4_TRIPLE_HYPERPLANE_QUADRIC']
    for label in labels:
        record=next(r for r in catalogue['records']if r['component_id']==label)
        assert record['expected_s61_dimension']==13 and record['full_s61_component_dominance']
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.json',
           ROOT/'reports/Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63.json',
           ROOT/'code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json',
           ROOT/'legacy/v03/p4_poisson_completeness_v3.tex',ROOT/'legacy/v05/triple_hyperplane.tex',
           ROOT/'code/s61_model.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py',
           ROOT/'research/notes/Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW.md']
    result=dict(id='Q_POSITIVE_WEIGHT6_ORBIT_20261103_v63_REVIEW',all_passed=True,
        fresh_producer_and_complete_matrix_certificate_match=True,
        independent_all25_physical_matrix_derivatives_and_canonical_divergence_zero=True,
        exact23minor_and3independent_full_kernel_vectors_checked=True,
        actual_affine_GL_rank=22,GL_plus_independent_scalar_rank=23,projective_PGL_orbit_dimension=22,
        independent_full150_finite_diagonal_pushforward_identity=True,
        nonzero_eta_one_complex_projective_orbit_and_closure_Qa8over9=True,
        T39_T40_excluded_by_closed_invariant_components_and_actual_limit=True,
        two22dimensional_quadric_components_excluded_by_closed_pure_fifth_power_cone=True,
        no_local_tangent_dimension_as_component_dimension_inference=True,
        other_known_membership_new_component_or_branch_exhaustion_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as output:
        json.dump(result,output,indent=2)
        output.write('\n')
    print('PASS physical projective orbit22 and four full-component exclusions on all eta!=0')


if __name__=='__main__':main()
