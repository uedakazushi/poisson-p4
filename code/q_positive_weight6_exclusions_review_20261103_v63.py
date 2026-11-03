#!/usr/bin/env python3
"""Independent physical tensor, closed-jet polynomial and whole-flag exclusions."""
import argparse
from fractions import Fraction
import gzip
import itertools
import json
from pathlib import Path
import sympy as s
from s61_model import S61, X, J, coefficients, pushforward
from toric_closed_jet_coordinates_20261102 import coords
from q_smith_global_model_review_20261103_v63 import ROOT, read, sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--replay',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    directory=ROOT/'reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63'
    d=read(directory/'report.json')
    fresh=read(args.replay/'report.json')
    assert {k:v for k,v in d.items()if k!='seconds'}=={k:v for k,v in fresh.items()if k!='seconds'}
    assert (directory/'certificate.json.gz').read_bytes()==(args.replay/'certificate.json.gz').read_bytes()
    assert sha(directory/'certificate.json.gz')==d['certificate_sha256']
    for filename,digest in d['input_sha256'].items():
        assert sha(ROOT/filename)==digest,filename
    c=read(directory/'certificate.json.gz')
    eta=s.Symbol('eta')
    P=s.zeros(5)
    for i,j,v in [(1,2,-s.Rational(256,3)*eta*X[0]*X[4]),
                  (1,3,-64*eta*X[1]*X[4]),(1,4,s.Rational(8,9)*X[0]**2),
                  (2,3,-s.Rational(2,9)*X[0]**2+64*eta*X[2]*X[4]),
                  (2,4,s.Rational(4,3)*X[0]*X[1]),
                  (3,4,s.Rational(2,3)*X[0]*X[2]+X[1]**2)]:
        P[i,j]=v
        P[j,i]=-v
    assert list(map(str,coefficients(P)))==c['raw150']
    assert all(s.expand(sum(s.diff(P[i,j],X[i])for i in range(5)))==0 for j in range(5))
    for i,j,k in itertools.combinations(range(5),3):
        assert s.expand(sum(P[i,h]*s.diff(P[j,k],X[h])+P[j,h]*s.diff(P[k,i],X[h])+P[k,h]*s.diff(P[i,j],X[h])for h in range(5)))==0
    F=0
    for e in range(5):
        i,j,k,l=[v for v in range(5)if v!=e]
        F+=(-1)**e*X[e]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
    assert s.expand(F)==-s.Rational(16,81)*X[0]**5
    frame=s.zeros(5)
    frame[0,0]=frame[1,1]=frame[3,2]=1
    frame[2,4]=s.Rational(9,8)
    frame[4,3]=-s.Rational(9,2)
    model=S61()
    normal=model.normalize(pushforward(P,frame))
    assert normal['perspective']==[0]*4
    u=s.Matrix(normal['parameters'])
    assert list(map(str,u))==c['S61']
    assert model.verify_point(normal)['reduced_equations']==134
    y=s.Matrix(s.symbols('y0:4'))
    raw=model.reconstruct(u)
    affine=s.Matrix(4,4,lambda i,j:s.expand((raw[i+1,j+1]-X[i+1]*raw[0,j+1]+X[j+1]*raw[0,i+1]).subs({X[0]:1}|dict(zip(X[1:],y)))))
    L=s.Matrix(4,4,lambda i,j:sum(s.Poly(affine[i,j],*y).coeff_monomial(v)*v for v in y))
    alpha=[s.diff(L[0,1]+L[2,3],v)/5 for v in y]
    assert alpha==[0]*4
    pairs=list(itertools.combinations(range(4),2))
    m2=list(itertools.combinations_with_replacement(range(4),2))
    D=coords()
    linear=[s.diff(L[i,j],v)for i,j in pairs for v in y]
    quadratic=s.Matrix([s.Poly(affine[i,j],*y).coeff_monomial(y[k]*y[l])for i,j in pairs for k,l in m2])
    # Independently check that Pi really implements -J d(i_E(-JQJ))/4 J
    # on this physical quadratic jet, rather than the scaled producer's4Pi.
    Q=s.zeros(4)
    for h,(i,j)in enumerate(pairs):
        Q[i,j]=sum(quadratic[10*h+k]*y[a]*y[b]for k,(a,b)in enumerate(m2))
        Q[j,i]=-Q[i,j]
    j=J[1:5,1:5]
    form=-j*Q*j
    contraction=-form*y
    closed=s.Matrix(4,4,lambda i,k:s.expand((s.diff(contraction[k],y[i])-s.diff(contraction[i],y[k]))/4))
    projected=-j*closed*j
    encoded=s.Matrix([s.Poly(projected[i,k],*y).coeff_monomial(y[a]*y[b])for i,k in pairs for a,b in m2])
    assert encoded==D['E']*D['Pi']*quadratic
    values=[linear[i]for i in D['li']]+list(D['Pi']*quadratic)
    assert list(map(str,values))==c['closedjet61']
    stored=read(ROOT/'reports/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.json')
    selected={tuple(mon):Fraction(coef)for mon,coef in stored['highest_lifts_in_projected_coordinates']['(1, 1)'][0]}
    assert selected=={tuple(mon):Fraction(coef)for mon,coef in c['selected_global_T01_polynomial']}
    assert len(selected)==403
    # Evaluate coefficient dictionaries without the producer's evaluation/orbit code.
    result=0
    surviving=[]
    for mon,coefficient in selected.items():
        term=s.Rational(coefficient.numerator,coefficient.denominator)
        for index in mon:
            term*=values[index]
        if term:
            surviving.append(mon)
            result+=term
    assert len(surviving)==4 and s.expand(result)==-s.Rational(1163264,135)*eta**2
    beta=-j*L*j
    K=s.Matrix([[s.diff(beta[i,k],v)for v in y]for i,k in pairs if (i,k)!=(2,3)])
    assert K==s.Matrix(c['K_beta'])
    assert K[:,3]==s.zeros(5,1)
    assert s.factor(K[[1,2,4],[0,1,2]].det())==-384*eta
    transport=s.eye(5)
    transport[1:5,1:5]=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    assert transport*J*transport.T==J
    assert transport[1:5,1:5]*s.Matrix([0,0,0,1])==s.Matrix([0,1,0,0])
    off=s.Matrix(model.normalize(pushforward(raw,transport))['parameters'])
    bundle=read(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json')
    review=read(ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json')
    assert review['proper_closed_incidence_and_rank3_kernel_inverse_accepted']
    for name,family in bundle['families'].items():
        A=s.Matrix(family['beta_embedding'])
        C=s.Matrix(family['beta_closed_constraints'])
        assert A.rank()==6 and C*A==s.zeros(C.rows,6)
        residual=C*off[bundle['beta_indices_zero_based'],0]
        assert residual==s.Matrix(c['all_CA_CB_beta_constraints_at_sole_flag'][name])
        assert s.Rational(3,2)in residual
    paths=[Path(__file__),directory/'report.json',directory/'certificate.json.gz',
           ROOT/'code/q_smith_global_model_review_20261103_v63.py',ROOT/'code/s61_model.py',
           ROOT/'code/toric_closed_jet_coordinates_20261102.py',
           ROOT/'reports/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.json',
           ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json',
           ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json',
           ROOT/'research/notes/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW.md']
    report=dict(id='Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63_REVIEW',all_passed=True,
        fresh_QQ_274_root_stability_and_dense_coframe_producer_matches=True,
        independent_all150_physical_tensor_strict_Jacobi_and_divergence_and_quintic=True,
        independent_actual_S61_normalization_and134_and150=True,
        independent_physical_closed_projection_and_Pi_not_fourPi=True,
        global_T01_polynomial_extension_by_root_stability_and_dense_coframe_reviewed=True,
        selected403_term_global_relation_independently_evaluated='-1163264*eta**2/135',
        independent_intrinsic_Kbeta_rank3_minor='-384*eta',
        entire_projective_flag_unique_for_eta_nonzero='e_w',
        finite_Sp4_transport_and_both_complete_beta_constraints_checked=True,
        all_nonzero_ray_outside_actual_T01_CA_CB_N=True,
        other_known_component_membership_new_component_branch_count_or_nonflatness_claimed=False,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in dict.fromkeys(paths)})
    with args.output.open('x')as f:
        json.dump(report,f,indent=2)
        f.write('\n')
    print('PASS independent actual +8/9 ray exclusions T01/CA/CB/N; all eta!=0')


if __name__=='__main__':
    main()
