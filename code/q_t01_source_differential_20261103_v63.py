#!/usr/bin/env python3
"""Exact differential and actual source-fiber audit for Qa in actual T01.

This is separate from the frozen universal-membership certificate.  All
matrices use the full canonical homogeneous 150 coefficients, with no
projective quotient and no jet truncation of the source map.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as sp
from sympy.polys.matrices import DomainMatrix

from q_t01_membership_20261103_v63 import (
    PAIRS, MONS, matrix_zero, identity, multiply, plus, scale,
    wedge_coefficients, combine,
)

ROOT=Path(__file__).resolve().parents[1]
POWER_PAIRS=list(itertools.combinations(range(1,5),2))
BASIS=list(itertools.product(PAIRS,MONS))
INDEX={value:i for i,value in enumerate(BASIS)}


def as_matrix(columns,K):
    return DomainMatrix([[col[i]for col in columns]for i in range(len(columns[0]))],
                        (len(columns[0]),len(columns)),K)


def dm_rows(rows,K):
    return DomainMatrix(rows,(len(rows),len(rows[0])if rows else 0),K)


def a_from_s(s,K):
    return (4*s+6)/(3*s*(s+3))


def build(K,s):
    """Return raw150, the 150x31 Frechet differential, N and powers."""
    Z,O=K.zero,K.one;I=identity(Z,O);N=matrix_zero(Z)
    for i,c in enumerate([O,s,(s+3)/2,-9*O/(s+3)]):N[i+1][i]=c
    powers=[I]
    for _ in range(7):powers.append(multiply(powers[-1],N,Z))
    a=a_from_s(s,K);lam=-2*a/(9*s);mu=-4*O/(9*s*s*(s+3))
    raw=combine(wedge_coefficients(powers[1],powers[4]),
                wedge_coefficients(powers[2],powers[3]),lam,mu)
    target=[Z]*150
    for key,value in {((1,4),(0,0)):a,((2,3),(0,0)):-2*O/9,
                      ((2,4),(0,1)):4*O/3,((3,4),(0,2)):2*O/3,
                      ((3,4),(1,1)):O}.items():target[INDEX[key]]=value
    assert raw==target
    cols=[]
    for h in range(5):
        for k in range(5):
            E=matrix_zero(Z);E[h][k]=O;dW=[None]
            for j in range(1,5):
                D=matrix_zero(Z)
                for r in range(j):D=plus(D,multiply(multiply(powers[r],E,Z),powers[j-1-r],Z))
                tr=sum((D[i][i]for i in range(5)),Z)/5
                dW.append(plus(D,scale(I,-tr)))
            col=[Z]*150
            for A,B,c in [(dW[1],powers[4],lam),(powers[1],dW[4],lam),
                          (dW[2],powers[3],mu),(powers[2],dW[3],mu)]:
                col=combine(col,wedge_coefficients(A,B),O,c)
            cols.append(col)
    for i,j in POWER_PAIRS:cols.append(wedge_coefficients(powers[i],powers[j]))
    return dict(a=a,lam=lam,mu=mu,N=N,powers=powers,raw=raw,columns=cols,J=as_matrix(cols,K))


def physical_lie_columns(raw,K):
    """Direct differentiation of g P(g^-1 x) g^T in all150 coefficients."""
    Z,O=K.zero,K.one;columns=[]
    for h in range(5):
        for k in range(5):
            col=[Z]*150
            def put(i,j,p,q,c):
                if i==j:return
                if i>j:i,j,c=j,i,-c
                col[INDEX[((i,j),tuple(sorted((p,q))))]]+=c
            for ((i,j),(p,q)),c in zip(BASIS,raw):
                if not c:continue
                if i==k:put(h,j,p,q,c)
                if j==k:put(i,h,p,q,c)
                if p==h:put(i,j,k,q,-c)
                if q==h:put(i,j,p,k,-c)
            columns.append(col)
    return columns


def polynomial_reparameterizations(model,K):
    """Five integrable generator changes M -> p0 I+...+p4 M^4."""
    Z=K.zero;columns=[]
    for r in range(5):
        col=[c for row in model['powers'][r]for c in row]+[Z]*6
        dc=[Z]*6
        def add(i,j,c):
            if i==j or not(1<=i<=4 and 1<=j<=4):return
            if i>j:i,j,c=j,i,-c
            dc[POWER_PAIRS.index((i,j))]-=c
        for i,j,c in [(1,4,model['lam']),(2,3,model['mu'])]:
            add(i-1+r,j,i*c);add(i,j-1+r,j*c)
        col[25:]=dc;columns.append(col)
    return columns


def commutator_columns(model,K):
    Z,O=K.zero,K.one;N=model['N'];cols=[]
    for h in range(5):
        for k in range(5):
            E=matrix_zero(Z);E[h][k]=O
            C=plus(multiply(E,N,Z),scale(multiply(N,E,Z),-O))
            cols.append([c for row in C for c in row]+[Z]*6)
    return as_matrix(cols,K)


def rank_receipt(M,K,include_minor=True):
    _,cp=M.rref();rank=len(cp)
    result={'shape':list(M.shape),'rank':rank,'pivot_columns':list(cp)}
    if rank and include_minor:
        _,rp=M.extract(list(range(M.shape[0])),list(cp)).transpose().rref()
        rp=list(rp);minor=M.extract(rp,list(cp));det=minor.det()
        assert det
        result.update(pivot_rows=rp,nonzero_minor=str(sp.factor(K.to_sympy(det))))
    return result


def sparse(M,K):
    return {'shape':list(M.shape),'nonzero':[[i,j,str(K.to_sympy(c))]
            for i,row in enumerate(M.to_list())for j,c in enumerate(row)if c]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    symbol=sp.Symbol('s');K=sp.QQ.frac_field(symbol);s=K.gens[0]
    model=build(K,s);J=model['J'];jrec=rank_receipt(J,K)
    assert jrec['rank']==24
    lie=as_matrix(physical_lie_columns(model['raw'],K),K)
    comm=commutator_columns(model,K)
    assert J*comm==lie
    lrec=rank_receipt(lie,K);assert lrec['rank']==17
    stabilizer=lie.nullspace().transpose()
    assert stabilizer.shape==(25,8) and (lie*stabilizer).is_zero_matrix
    stabsource=comm*stabilizer
    srec=rank_receipt(stabsource,K);assert srec['rank']==3
    poly=as_matrix(polynomial_reparameterizations(model,K),K)
    assert (J*poly).is_zero_matrix
    prec=rank_receipt(poly,K);assert prec['rank']==5
    jointfiber=poly.hstack(stabsource)
    frec=rank_receipt(jointfiber,K);assert frec['rank']==7
    assert (J*jointfiber).is_zero_matrix and 31-jrec['rank']==frec['rank']
    # The entire 5-dimensional regular-nilpotent centralizer lies in strict
    # stabilizer.  This identifies H orbit dimension 8-5=3, not a conjecture
    # inferred from the source differential alone.
    centralizer=as_matrix([[c for row in power for c in row]for power in model['powers'][:5]],K)
    assert (lie*centralizer).is_zero_matrix
    assert rank_receipt(centralizer,K,False)['rank']==5
    assert (comm*centralizer).is_zero_matrix
    sigma=-3*(s+3)/(2*s+3)
    other=build(K,sigma)
    assert other['a']==model['a'] and other['raw']==model['raw']
    joint=J.hstack(other['J']);jointrec=rank_receipt(joint,K)
    assert other['J']*commutator_columns(other,K)==lie
    qdot=[K.zero]*150;qdot[INDEX[((1,4),(0,0))]]=K.one
    common=lie.hstack(as_matrix([qdot],K));commonrec=rank_receipt(common,K)
    assert commonrec['rank']==18
    curve_derivative=[]
    for value in [c for row in model['N']for c in row]+[
            model['lam']if pair==(1,4)else model['mu']if pair==(2,3)else K.zero
            for pair in POWER_PAIRS]:
        curve_derivative.append(K.from_sympy(sp.diff(K.to_sympy(value),symbol)))
    da=K.from_sympy(sp.diff(K.to_sympy(model['a']),symbol));assert da
    assert J*as_matrix([curve_derivative],K)==as_matrix([[da*c for c in qdot]],K)
    assert rank_receipt(other['J'].hstack(common),K,False)['rank']==24
    assert 24+24-jointrec['rank']==commonrec['rank']
    # Exact small closed-point checks are supplementary; the generic joint
    # rank and closedness of Sing(T01) cover ALL finite a if joint rank>26.
    closed={}
    for label,d,ss,a0 in [('a1',97,(-5+sp.sqrt(97))/6,sp.S.One),
                         ('a_two_ninths',5,(3+3*sp.sqrt(5))/2,sp.Rational(2,9))]:
        L=sp.QQ.algebraic_field(sp.sqrt(d));u=L.from_sympy(ss)
        left=build(L,u);right=build(L,-3*(u+3)/(2*u+3))
        assert left['a']==L.from_sympy(a0)and left['raw']==right['raw']
        r1=rank_receipt(left['J'],L,False);r2=rank_receipt(right['J'],L,False)
        rj=rank_receipt(left['J'].hstack(right['J']),L,False)
        closed[label]={'field':str(L),'s':str(ss),'left':r1,'right':r2,'joint':rj,
                       'intersection_dimension':r1['rank']+r2['rank']-rj['rank']}
    # Persist full generic matrices, not merely sampled rank claims.
    payload={'coefficient_field':'QQ(s)','a_of_s':str(K.to_sympy(model['a'])),
             'sigma_s':str(K.to_sympy(sigma)),'source_columns':'M00,...,M44,c12,c13,c14,c23,c24,c34',
             'raw_basis':'lex i<j (0..4), then lex p<=q (0..4)',
             'jacobian':sparse(J,K),'other_jacobian':sparse(other['J'],K),
             'strict_GL5_action':sparse(lie,K),'strict_stabilizer_basis':sparse(stabilizer,K),
             'polynomial_generator_kernel':sparse(poly,K),
             'stabilizer_conjugation_kernel':sparse(stabsource,K),
             'source_fiber_spanning_directions':sparse(jointfiber,K),
             'source_curve_s_derivative':sparse(as_matrix([curve_derivative],K),K),
             'common_GL_orbit_and_Qa_curve':sparse(common,K)}
    packed=args.output/'matrices.json.gz'
    with gzip.open(packed,'wt')as handle:json.dump(payload,handle,separators=(',',':'))
    inputs=[Path(__file__),ROOT/'code/q_t01_membership_20261103_v63.py',
            ROOT/'reports/Q_T01_MEMBERSHIP_20261103_v63/report.json',
            ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json']
    report={'id':'Q_T01_SOURCE_DIFFERENTIAL_20261103_V63','all_passed':True,
            'coefficient_field':'QQ(s)','full_canonical_raw150_used':True,
            'source_jacobian':jrec,'strict_GL5_action':lrec,
            'strict_GL5_stabilizer_dimension':8,'regular_N_centralizer_dimension':5,
            'stabilizer_source_conjugation':srec,'generator_reparameterizations':prec,
            'combined_actual_source_fiber_directions':frec,
            'actual_source_fiber_regular_dimension_generic_Qa':7,
            'source_fiber_regular_proof':'Integrable H-conjugations plus invertible polynomial generator changes give dimension>=7; full150 differential gives embedding dimension7.',
            'conjugate_source_joint':jointrec,
            'generic_two_image_intersection_dimension':48-jointrec['rank'],
            'generic_two_image_intersection_equals_GL_orbit_plus_Qa_curve':True,
            'common_orbit_and_Qa_curve':commonrec,
            'closed_point_checks':closed,
            'generic_Qa_T01_tangent_lower_bound':jointrec['rank'],
            'entire_finite_Qa_affine_T01_tangent_lower_bound_by_closedness':jointrec['rank'],
            'entire_finite_Qa_projective_T01_tangent_lower_bound_by_closedness':jointrec['rank']-1,
            'entire_finite_Qa_T01_singular_by_closedness':jointrec['rank']>26,
            'local_branch_count_or_normalization_fiber_claimed':False,
            'target_smoothness_inferred_from_source_fiber':False,
            'seconds':time.monotonic()-started,
            'matrices_sha256':hashlib.sha256(packed.read_bytes()).hexdigest(),
            'input_sha256':{str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in inputs}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:report[key]for key in ['all_passed','source_jacobian','strict_GL5_action',
          'combined_actual_source_fiber_directions','conjugate_source_joint','closed_point_checks','seconds']},indent=2))


if __name__=='__main__':main()
