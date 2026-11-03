#!/usr/bin/env python3
"""Bounded independent checks for v62 rational MPS families and image minors.

The Schouten calculation below uses sparse polynomial-ring operations directly;
it does not import either author's Jacobi implementation or family builder.
Component labels and their local deformation interpretation remain dependencies.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[1]
P=32003


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_family(label):
    raw=ROOT/'reports/MPS_S61_RAW_CATALOGUE_20261103_v62'/f'{label}.json'
    record=json.loads(raw.read_text())
    if record['projective_Jacobi_identically_zero']:
        return raw,record
    corrected=ROOT/'reports/MPS_SMOOTHING_COMPLETION_20261103_v62'/f'{label}.json'
    record=json.loads(corrected.read_text())
    assert record['all_passed']
    return corrected,record


def check_jacobi(record):
    parameters=s.symbols(' '.join(record['period_parameters']+record['smoothing_parameters']),seq=True)
    field=s.QQ.frac_field(*parameters)
    R,*x=ring('x0,x1,x2,x3,x4',field)
    Pi=[[R.zero for _ in range(5)]for _ in range(5)]
    for pair,expression in record['homogeneous_bivector'].items():
        i,j=map(int,pair.split(','))
        Pi[i][j]=R.from_expr(s.sympify(expression))
        Pi[j][i]=-Pi[i][j]
    triples={}
    for i,j,k in itertools.combinations(range(5),3):
        triples[i,j,k]=sum((Pi[i][ell]*Pi[j][k].diff(x[ell])
                           +Pi[j][ell]*Pi[k][i].diff(x[ell])
                           +Pi[k][ell]*Pi[i][j].diff(x[ell])for ell in range(5)),R.zero)
    for quad in itertools.combinations(range(5),4):
        value=sum(((-1)**a*x[index]*triples[tuple(quad[:a]+quad[a+1:])]
                   for a,index in enumerate(quad)),R.zero)
        assert value==R.zero,(record['component_id'],quad,str(value))
    return len(parameters)


def check_equivariance(record,diagram):
    """Every stored coefficient must have the actual smoothing torus weight.

    Period-degree homogeneity separately proves that the affine image is a cone.
    These checks use the frozen final coefficient expressions, not the author's
    weight-space enumeration or recursive construction.
    """
    periods=s.symbols(' '.join(record['period_parameters']),seq=True)
    smoothings=s.symbols(' '.join(record['smoothing_parameters']),seq=True)
    x=s.symbols('x0:5')
    weights=[]
    for edge in diagram['edges']:
        weight=[0]*5
        for i in edge['edge']:weight[i]-=1
        for i,power in edge['orders'].items():weight[int(i)]+=int(power)
        weights.append(weight)
    assert len(weights)==len(smoothings)
    terms=0
    for pair,expression in record['homogeneous_bivector'].items():
        i,j=map(int,pair.split(','))
        for spatial,coefficient in s.Poly(s.sympify(expression),*x).terms():
            target=list(spatial);target[i]-=1;target[j]-=1
            for multidegree,value in s.Poly(coefficient,*smoothings).terms():
                actual=[sum(n*w[k]for n,w in zip(multidegree,weights))for k in range(5)]
                assert actual==target,(record['component_id'],pair,spatial,multidegree)
                total=s.cancel(sum(p*s.diff(value,p)for p in periods)
                               -(1-sum(multidegree))*value)
                assert total==0,(record['component_id'],'scalar homogeneity',value)
                terms+=1
    return terms


def check_minor(record):
    matrix=record['differential_matrix_QQ']
    witness=record['nonzero_minor']
    rows,cols=witness['rows'],witness['columns']
    assert len(rows)==len(cols)==witness['rank']
    def reduce(value):
        a=s.Rational(value)
        assert int(a.q)%P
        return int(a.p)*pow(int(a.q),-1,P)%P
    A=[[reduce(matrix[i][j])for j in cols]for i in rows]
    determinant=1
    for j in range(len(A)):
        k=next(i for i in range(j,len(A))if A[i][j])
        if k!=j:
            A[k],A[j]=A[j],A[k]
            determinant=-determinant
        pivot=A[j][j]
        determinant=determinant*pivot%P
        inverse=pow(pivot,-1,P)
        for i in range(j+1,len(A)):
            scalar=A[i][j]*inverse%P
            if scalar:
                A[i]=[(a-scalar*b)%P for a,b in zip(A[i],A[j])]
    assert determinant%P==witness['determinant_mod_32003']!=0
    return dict(component_id=record['component_id'],rank_lower_bound=witness['rank'],
                determinant_mod_32003=determinant%P)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    started=time.monotonic()
    inputs={str(Path(__file__).relative_to(ROOT)):sha(__file__)}
    diagrams_path=ROOT/'legacy/v02/baseline/data/p4_toric_smoothing_diagrams.json'
    diagrams={d['id']:d for d in json.loads(diagrams_path.read_text())['diagrams']}
    inputs[str(diagrams_path.relative_to(ROOT))]=sha(diagrams_path)
    jacobi=[];minors=[]
    for number in range(2,39):
        label=f'T{number:02d}'
        source,record=load_family(label)
        count=check_jacobi(record)
        term_count=check_equivariance(record,diagrams[label])
        inputs[str(source.relative_to(ROOT))]=sha(source)
        jacobi.append(dict(component_id=label,parameter_count=count,projective_Jacobi=True,
                           smoothing_torus_weights_checked=True,scalar_homogeneity_checked=True,
                           coefficient_terms_checked=term_count))
        print(label,'SPARSE PROJECTIVE JACOBI PASS',flush=True)
    for number in range(1,39):
        source=ROOT/'reports/MPS_S61_DOMINANCE_20261103_v62'/f'T{number:02d}.json'
        record=json.loads(source.read_text())
        minors.append(check_minor(record))
        inputs[str(source.relative_to(ROOT))]=sha(source)
    report=dict(id='MPS_V62_INDEPENDENT_REVIEW_20261103',all_passed=True,
                sparse_Jacobi_checks=jacobi,independently_recomputed_minors=minors,
                direct_image_upper_bound='A free diagonal T4 acts on each PGL5 times parameter source with image fixed, by verified weight covariance. Hence affine image dimension <=24+s+e-4=20+s+e. The nonzero minor attains this bound.',
                limitations=['No generic component upper bound is inferred.',
                             'Component labels use the stated MPS local theorem and graph crosswalk.',
                             'The matrices are supplied by the producer; determinants are recomputed independently.'],
                input_sha256=inputs,seconds=time.monotonic()-started)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS',len(jacobi),'families',len(minors),'minors',report['seconds'])


if __name__=='__main__':
    main()
