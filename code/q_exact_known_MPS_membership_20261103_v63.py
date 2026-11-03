#!/usr/bin/env python3
"""Exact known-MPS40 membership of physical Qa on the entire a!=0 locus."""
import argparse
import json
from pathlib import Path
import sympy as s
from s61_model import S61, X, J, coefficients, pushforward
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    a=s.Symbol('a')
    P=s.zeros(5)
    for i,j,v in [(1,4,a*X[0]**2),(2,3,-s.Rational(2,9)*X[0]**2),
                  (2,4,s.Rational(4,3)*X[0]*X[1]),
                  (3,4,s.Rational(2,3)*X[0]*X[2]+X[1]**2)]:
        P[i,j]=v
        P[j,i]=-v
    M=s.zeros(5)
    M[0,0]=M[1,1]=M[3,2]=1
    M[2,4]=1/a
    M[4,3]=-s.Rational(9,2)
    assert s.factor(M.det())==-s.Rational(9,2)/a
    C=P.subs(dict(zip(X,[1,0,0,0,0])))[1:5,1:5]
    assert (M[1:5,1:5]*C*M[1:5,1:5].T-J[1:5,1:5]).applyfunc(s.cancel)==s.zeros(4)
    model=S61()
    normal=model.normalize(pushforward(P,M))
    assert normal['perspective']==[0]*4
    values=list(map(s.cancel,normal['parameters']))
    expected=[s.S.Zero]*61
    expected[24]=-4/(3*a)
    expected[35]=3/a
    expected[37]=9/(2*a)
    assert values==expected
    checks=model.verify_point(normal)
    for v in [*M,*values]:
        den=s.Poly(s.denom(s.cancel(v)),a)
        assert len(den.terms())==1
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/q_smith_global_model_review_20261103_v63.py',
           ROOT/'reports/Q_ALL_TORIC_MPS_COMPONENTS_20261103_v63.json']
    all38=read(paths[-1])
    assert all38['all_passed'] and all38['actual_known_MPS_components_through_ALL_Qa']==[f'T{i:02d}'for i in range(1,39)]
    residuals={}
    u=s.symbols('u1:62')
    for label,root in [('Q51',s.Rational(2,9)),('Q52',-s.Rational(8,9))]:
        p=ROOT/f'reports/{label}_UNIVERSAL_SECTION_20261101.json'
        review=ROOT/f'reports/{label}_UNIVERSAL_SECTION_REVIEW_20261101.json'
        d=read(p)
        r=read(review)
        assert d['all_passed'] and r['all_passed']
        paths.extend([p,review,ROOT/f'research/notes/{label}_UNIVERSAL_SECTION_20261101.md',
                      ROOT/f'research/notes/{label}_UNIVERSAL_SECTION_REVIEW_20261101.md'])
        for dd in [d,r]:
            for filename,digest in dd.get('input_sha256',dd.get('source_sha256',{})).items():
                assert sha(ROOT/filename)==digest,filename
        beta=d['weight1_u_indices_zero_based']
        q=d['weight2_u_indices_zero_based']
        assert len(beta)==16 and len(q)==45 and sorted(beta+q)==list(range(61))
        sub={u[i]:values[i] for i in beta}
        relation=[s.factor(values[j]-s.sympify(f).subs(sub))for j,f in zip(q,d['quadratic_section_q'])]
        target=(9*a-2)/(2*a*a) if label=='Q51' else (9*a+8)/(2*a*a)
        assert all(s.cancel(v-(target if j==25 else 0))==0 for j,v in enumerate(relation))
        assert s.cancel(target.subs(a,root))==0
        assert all(s.cancel(v.subs(a,root))==0 for v in relation)
        assert s.Poly(s.numer(target),a).degree()==1
        residuals[label]={'only_nonzero_ordinal_zero_based':25,'relation':str(target),
                          'unique_nonzero_parameter_root':str(root),'all45_residuals':list(map(str,relation))}
    paths.append(ROOT/'research/notes/Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63.md')
    report=dict(id='Q_EXACT_KNOWN_MPS_MEMBERSHIP_20261103_v63',all_passed=True,
        scope='Physical Qa for every nonzero finite complex a; no negative conclusion at a=0',
        direct_physical_Qa_normalization_without_finite35_or_chosen_group_slice=True,
        matrix_determinant=str(s.factor(M.det())),only_variable_denominator='a',
        all61_normalized_parameters=list(map(str,values)),physical_model_checks=checks,
        global_reduced_T39_T40_S61_graph_identification_from_existing_universal_proofs=True,
        all45_graph_residuals_independently_evaluated=residuals,
        actual_T39_membership_iff='a=2/9',actual_T40_membership_iff='a=-8/9',
        exact_known_MPS40_membership='T01--T38, plus T39 only at2/9 or T40 only at-8/9',
        unknown_components_other_catalogue_families_or_local_branch_exhaustion_claimed=False,
        finite35_group_minor_open_or_discriminant_open_used=False,
        global_known_component_lower_bound_remains=57,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:
        json.dump(report,f,indent=2)
        f.write('\n')
    print('PASS ALL nonzero physical Qa: exact MPS40 membership, T39 iff2/9 and T40 iff-8/9')


if __name__=='__main__':
    main()
