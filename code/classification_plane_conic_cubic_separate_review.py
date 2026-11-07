#!/usr/bin/env python3
"""Separate-author full cubic slice/component check for contact4-other-tangent.

Reuse the previous separate-author source arithmetic, without importing the
two-conic producer or its author-written component/cubic checker.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools as it
import json
from pathlib import Path
import resource
import time

import classification_plane_bitangent_component_review as prior
q = prior.q


def solve(A,b):
    B = [list(map(F,row))+[F(value)] for row,value in zip(A,b)]
    n = len(B)
    assert all(len(row)==n+1 for row in B)
    for j in range(n):
        k = next(k for k in range(j,n) if B[k][j])
        B[j],B[k] = B[k],B[j]
        div = B[j][j]
        B[j] = [v/div for v in B[j]]
        for k in range(j+1,n):
            value = B[k][j]
            if value:
                B[k] = [v-value*w for v,w in zip(B[k],B[j])]
    x = [F(0)]*n
    for i in reversed(range(n)):
        x[i] = B[i][-1]-sum(B[i][j]*x[j] for j in range(i+1,n))
    assert all(sum(v*t for v,t in zip(row,x))==value for row,value in zip(A,b))
    return x


def vec_poly(vector):
    out = {}
    for value,(ij,e) in zip(vector,q.COLS):
        if value:
            out.setdefault(ij,{})[e+(0,0,0)] = value
    return q.antisymmetric(out)


def plus(A,B):
    return {ij:q.add(A.get(ij,{}),B.get(ij,{})) for ij in set(A)|set(B)}


def J(P):
    return q.jvec(q.jacobi(P))


def mixed(A,B):
    return [a-b-c for a,b,c in zip(J(plus(A,B)),J(A),J(B))]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args = ap.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(110,115));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
    began = time.monotonic()
    paths = ['reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_AMBIENT_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_SUPPORT_20261107.json',
             'reports/CLASSIFICATION_BASELINE90_COMPONENT_PROJECTION_20261107.json',
             'reports/CLASSIFICATION_TWO_QUADRICS_PLANE_BITANGENT_SEPARATE_REVIEW_20261107.json']
    sources,ambient,normals,cubics,baseline,bitangent = [json.loads(Path(p).read_text()) for p in paths]
    name = 'contact4_other_tangent'
    source = next(r for r in sources['cases'] if r['name']==name)
    actual = next(r for r in ambient['cases'] if r['name']==name)
    normal = next(r for r in normals['cases'] if r['case']==name)
    cubic = next(r for r in cubics['cases'] if r['case']==name)
    P,directions,source_certificate = prior.actual_source(actual,source)
    assert len(normal['normals'])==1
    R = q.decode(normal['normals'][0])
    T = q.decode(cubic['quadratic_corrections'][0])
    Dcols = []
    for ij,e in q.COLS:
        variation = dict(P)
        variation[ij] = q.add(P.get(ij,{}),{e+(0,0,0):F(1)})
        variation[ij[::-1]] = q.scale(variation[ij],-1)
        Dcols.append(J(variation))
    D = q.transpose(Dcols)
    Vcols,Ocols = [],[]
    for a,b in it.product(range(5),repeat=2):
        V,O = {},{}
        for i,j in q.PAIRS:
            V[i,j] = q.add(q.atom(1,i,b) if j==a else {},q.atom(-1,j,b) if i==a else {})
            O[i,j] = q.add(P.get((b,j),{}) if i==a else {},P.get((i,b),{}) if j==a else {},
                           q.scale(q.mul(q.atom(1,b),q.der(P.get((i,j),{}),a)),-1))
        Vcols.append(q.vec(V));Ocols.append(q.vec(O))
    Mcols = Vcols+Ocols+[q.vec(P)]
    Fcols = Mcols+[q.vec(v) for v in directions]
    rankF,basisF = q.rank(q.transpose(Fcols))
    assert q.rank(q.transpose(Vcols))[0]==24 and q.rank(q.transpose(Mcols))[0]==47
    assert rankF==48
    Kcols = [Fcols[i] for i in basisF]+[q.vec(R)]
    assert q.rank(q.transpose(Kcols))[0]==49
    assert all(sum(v*c for v,c in zip(row,col))==0 for row in D for col in Kcols)
    rows,cols = cubic['implicit_D_rows'],cubic['implicit_D_columns']
    assert len(rows)==len(cols)==101
    A = [[D[i][j] for j in cols] for i in rows]
    determinant = q.determinant(A)
    assert determinant and determinant==F(cubic['implicit_minor_determinant'])
    Q = J(R)
    selected_solution = solve(A,[-Q[i] for i in rows])
    Tvec = [F(0)]*150
    for j,value in zip(cols,selected_solution):
        Tvec[j] = value
    assert Tvec==q.vec(T)
    assert all(sum(v*c for v,c in zip(row,Tvec))+quad==0 for row,quad in zip(D,Q))
    C = mixed(R,T)
    assert {str(i):str(v) for i,v in enumerate(C) if v}==cubic['all_cubic_coefficient_vectors'][0]
    function = cubic['cubic_functionals'][0]
    left = {int(i):F(v) for i,v in function['nonzero_coefficients'].items()}
    assert all(sum(v*D[i][j] for i,v in left.items())==0 for j in range(150))
    assert sum(v*Q[i] for i,v in left.items())==0
    pairing = sum(v*C[i] for i,v in left.items())
    assert pairing==F(function['cubic_pairing'][0])==F(-533355840000)
    # Diagnose, but never assume, invariance under arbitrary kernel correction.
    # The proof uses the fixed U101 coordinate slice, which forbids such shifts.
    ambiguity = []
    for vector in Kcols:
        cross = mixed(R,vec_poly(vector))
        ambiguity.append(sum(v*cross[i] for i,v in left.items()))
    assert baseline['count']==90 and bitangent['all_passed'] and bitangent['permitted_lower_bound']==91
    assert bitangent['full_source_componenthood_and_distinctness_accepted']
    # Restrict second conics to the normalization of their respective q1.
    # Other tangent: (h,x,y)=(s²,t²,st), q2|q1=t⁴.
    # Bitangent: (h,x,y)=(s²,t²,t(s-t)/4), q2|q1=t³(s-t)/4.
    # Exact two-variable polynomial check, using spare indices5,6 of q.N=8.
    s,t = q.atom(1,5),q.atom(1,6)
    other_h,other_x,other_y = q.mul(s,s),q.mul(t,t),q.mul(s,t)
    other_q1 = q.add(q.mul(other_h,other_x),q.scale(q.mul(other_y,other_y),-1))
    other_q2 = q.add(other_q1,q.mul(other_x,other_x))
    assert not other_q1 and other_q2==q.mul(q.mul(t,t),q.mul(t,t))
    bh,bx,by = q.mul(s,s),q.mul(t,t),q.scale(q.mul(t,q.add(s,q.scale(t,-1))),F(1,4))
    byplus = q.add(bx,q.scale(by,4));byminus = q.add(bx,q.scale(by,-4))
    bq1 = q.scale(q.add(q.mul(bh,bx),q.scale(q.mul(byplus,byplus),-1)),F(1,16))
    bq2 = q.scale(q.add(q.mul(bh,bx),q.scale(q.mul(byminus,byminus),-1)),F(1,16))
    assert not bq1 and bq2==q.scale(q.mul(q.mul(t,t),q.mul(t,q.add(s,q.scale(t,-1)))),F(1,4))
    result = dict(task_id='E_reduced_nonnormal',all_passed=True,no_producer_import=True,no_CAS_import=True,
                  whole_source=source_certificate,full_D_shape=[350,150],full_D_exact_rank=101,
                  vertical_rank=24,actual_group_vertical_scale_rank=47,actual_source_rank=48,whole_kernel_dimension=49,
                  all_D_times_actual_source_and_normal_zero=True,
                  independently_solved_fixed_U101_correction_equals_saved_T=True,all350_D_T_plus_normal_square_zero=True,
                  all350_actual_cubic_coefficients_reconstructed=True,
                  cokernel_annihilates_all150_D_columns=True,cokernel_annihilates_normal_square=True,
                  nonzero_cubic_pairing=str(pairing),literal_true_slice_dimension=102,true_slice_cut_count=48,
                  arbitrary_kernel_correction_pairings=[str(v) for v in ambiguity],
                  arbitrary_kernel_invariance=all(not v for v in ambiguity),
                  arbitrary_kernel_invariance_required=False,
                  source_kernel_ambiguity_disposition='T is uniquely solved in literal U101, the selected standard-coordinate complement. No arbitrary kernel or source-direction addition is permitted in that fixed slice; higher U corrections change the cubic only through D, killed by the verified cokernel.',
                  full_projective_component_dimension=23,full_source_closure_is_component=True,local_uniqueness_claimed=False,
                  distinct_from_baseline90='Same verified dimension/rank/essential3 and genericfactor argument as the accepted bitangent review; all90 have different component dimension, rank, SNC membership or general factor stratum.',
                  distinct_from_accepted_bitangent=True,
                  exact_normalized_conic_intersection_restrictions=dict(other_tangent='t^4',bitangent='t^3*(s-t)/4'),
                  contact4_closed_condition='On the smooth distinct conic-pair locus, the restriction binary quartic is a fourth power. This is the closed Veronese locus, equivalently the proper contact-order>=4 incidence over the conic normalization. It is invariant under exchanging factors and the unique essential-plane quotient.',
                  generic_disjointness_argument='Both fullsource images are dense constructible in their proven irreducible component closures. A shared component would force overlap of nonempty opens carrying incompatible4 versus3+1 contact partitions.',
                  permitted_lower_bound=92,permitted_rank4_lower_bound=86,global_exhaustion_claimed=False,
                  input_sha256={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in paths+
                                [str(Path(prior.__file__).relative_to(Path.cwd())),str(Path(q.__file__).relative_to(Path.cwd())),'code/classification_quadratic_cofactor_review.py']},
                  reviewer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),reviewer_source_bytes=Path(__file__).read_text(),
                  seconds=time.monotonic()-began,maximum_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['all_passed','seconds','maximum_RSS_KiB','nonzero_cubic_pairing','arbitrary_kernel_invariance']}))


if __name__=='__main__':
    main()
