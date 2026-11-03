#!/usr/bin/env python3
"""Independent exact QQ review of the six-relation global cubic row basis."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import tarfile
import time

import sympy as sp

from s61_nonsquare_sparse_kernel_review_20261103_v62 import (
    R, aa, ee, dd, parse_matrix, parse_poly, multiply, known_unit, to_ring,
)

ROOT = Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def integer_polynomials(matrix):
    denominator = math.lcm(*(q.denominator for row in matrix for p in row for q in p.values()))
    return [[{u:q.numerator*(denominator//q.denominator) for u,q in p.items()} for p in row]
            for row in matrix], denominator


def exact_relation(coefficients, integer_matrix):
    denominator = math.lcm(*(q.denominator for p in coefficients.values() for q in p.values()))
    integral = {i:{u:q.numerator*(denominator//q.denominator) for u,q in p.items()}
                for i,p in coefficients.items()}
    for column in range(91):
        total = {}
        for i,p in integral.items():
            for u,c in p.items():
                for v,d in integer_matrix[i][column].items():
                    monomial = (u[0]+v[0],u[1]+v[1],u[2]+v[2])
                    total[monomial] = total.get(monomial,0)+c*d
        assert not any(total.values())
    return denominator


def determinant_mod(matrix, prime):
    matrix = [list(row) for row in matrix]
    value = 1
    for j in range(len(matrix)):
        pivot = next(i for i in range(j,len(matrix)) if matrix[i][j]%prime)
        if pivot != j:
            matrix[j],matrix[pivot] = matrix[pivot],matrix[j]
            value = -value
        p = matrix[j][j]%prime
        value = value*p%prime
        inv = pow(p,-1,prime)
        for i in range(j+1,len(matrix)):
            c = matrix[i][j]*inv%prime
            for k in range(j+1,len(matrix)):
                matrix[i][k] = (matrix[i][k]-c*matrix[j][k])%prime
            matrix[i][j] = 0
    return value%prime


def evaluate(polynomial,point,prime):
    answer = 0
    for u,q in polynomial.items():
        assert q.denominator%prime
        term = q.numerator*pow(q.denominator,-1,prime)%prime
        for x,power in zip(point,u):
            term = term*pow(x,power,prime)%prime
        answer = (answer+term)%prime
    return answer


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    begin=time.monotonic()
    directory=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62'
    report=json.loads((directory/'report.json').read_text())
    assert report['all_passed'] and not report['all_point_rank_lower57_claimed'] and not report['free_cokernel_claimed']
    for name,value in report['input_sha256'].items():assert sha((ROOT/name).read_bytes())==value
    for name,value in report['generated_sha256'].items():assert sha((directory/name).read_bytes())==value
    data_path=ROOT/'code/data/s61_binary_nonsquare_cubic_row_basis_20261103_v62.json'
    data=json.loads(data_path.read_text())
    assert data['known_units']==['a','500*a+90*e^2+3*e-24','(50*a+18*e^2-9*e-2)*d+60*e-40','100*a+18*e^2-9*e-2']
    units=[aa,500*aa+90*ee**2+3*ee-24,(50*aa+18*ee**2-9*ee-2)*dd+60*ee-40,100*aa+18*ee**2-9*ee-2]
    previous=ROOT/'reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62'
    with tarfile.open(previous/'matrices.tar.gz')as archive:full_bytes=archive.extractfile('full').read()
    assert sha(full_bytes)==data['original_full_matrix_sha256']
    full,target=parse_matrix(full_bytes);assert target==57
    v3path=ROOT/'research/handoff/20261103_v62/nonsquare_core_wip_v3_manifest.json'
    v3=json.loads(v3path.read_text())
    name='work/nonsquare-v62/fresh-cubic/syzygies-normalized-degree1.normalized-input'
    path=ROOT/name
    if path.exists()and sha(path.read_bytes())==v3['files'][name]['sha256']:normalized=path.read_bytes()
    else:
        with tarfile.open(v3path.parent/v3['archive'])as archive:normalized=archive.extractfile(name).read()
    assert sha(normalized)==v3['files'][name]['sha256']
    A,_=parse_matrix(normalized);assert len(A)==63 and len(A[0])==91
    multipliers=[parse_poly(line)for line in(directory/'normalization.txt').read_text().splitlines()[1:]]
    assert len(multipliers)==63
    for i,multiplier in enumerate(multipliers):
        known_unit(multiplier,units)
        for column in range(91):assert multiply(A[i][column],multiplier)==full[i][column]
    integer_A,denominator=integer_polynomials(A)
    deleted=[];relations=[];records=[]
    expected=[{(0,0,0):Fraction(1)},{(0,0,0):Fraction(1)},{(0,0,0):Fraction(1)},
              {(1,0,0):Fraction(1)},{(1,0,0):Fraction(-1)},{(1,0,0):Fraction(12,5)}]
    for r,item in enumerate(data['relations']):
        coefficients={entry['original_row']:parse_poly(str(len(entry['terms']))+''.join(' '+' '.join(map(str,term))for term in entry['terms']))for entry in item['coefficients']}
        drop=item['deleted_original_row'];bound=item['coefficient_degree_bound']
        assert drop==[48,55,59,57,61,14][r]and coefficients[drop]==expected[r]
        assert all(not coefficients.get(i)for i in deleted)
        assert all(sum(u)<=bound for p in coefficients.values()for u in p)
        known_unit(coefficients[drop],units)
        relation_denominator=exact_relation(coefficients,integer_A)
        records.append(dict(deleted_original_row=drop,coefficient=str(to_ring(coefficients[drop])),
                            degree_bound=bound,exact_QQ_columns_verified=91,
                            cleared_relation_denominator=str(relation_denominator)))
        relations.append(coefficients);deleted.append(drop)
    primary_matrix=sp.Matrix([[to_ring(row.get(column,{})).as_expr()for column in deleted]for row in relations])
    assert all(primary_matrix[i,j]==0 for i in range(6)for j in range(i))
    primary=sp.prod(primary_matrix[i,i]for i in range(6))
    assert sp.expand(primary)==-sp.Rational(12,5)*sp.Symbol('a')**3
    alternative_order=deleted[:-1]+[22]
    alternative_matrix=sp.Matrix([[to_ring(row.get(column,{})).as_expr()for column in alternative_order]for row in relations])
    assert all(alternative_matrix[i,j]==0 for i in range(6)for j in range(i))
    assert sp.expand(sp.prod(alternative_matrix[i,i]for i in range(6)))==sp.Symbol('a')**3*sp.Symbol('d')/10
    assert relations[-1][22]=={(1,0,1):Fraction(-1,10)}
    assert all(i<42 for i in relations[-1])
    for filename,order,key in [('remaining.txt.gz',deleted,'remaining_matrix_sha256'),('alternative-ad.txt.gz',alternative_order,'alternative_ad_matrix_sha256')]:
        contents=gzip.decompress((directory/filename).read_bytes())
        matrix,_=parse_matrix(contents)
        assert sha(contents)==report[key]and matrix==[row for i,row in enumerate(A)if i not in order]
    witness=report['fresh_generic_point_minor'];prime=witness['prime'];point=witness['point'];columns=witness['columns']
    assert prime==32003 and point==[1,2,3]and len(set(columns))==57
    basis=[row for i,row in enumerate(A)if i not in deleted]
    evaluated=[[evaluate(row[j],point,prime)for j in columns]for row in basis]
    det=determinant_mod(evaluated,prime)
    assert det==witness['determinant']==27632
    assert all(evaluate({u:Fraction(str(c))for u,c in unit.items()},point,prime)for unit in units)
    labels=json.loads((previous/'source-report.json').read_text())['labels']
    assert len(labels)==63 and all(label[0]=='product'for label in labels[:42])and all(label[0]=='cubic'for label in labels[42:])
    assert len([i for i in range(42)if i not in deleted])==41
    assert len([i for i in range(42,63)if i not in deleted])==16
    (args.output/'independent-checks.json').write_text(json.dumps(dict(relations=records,
        global_matrix_denominator=str(denominator),primary_triangular_minor=str(primary),
        alternative_minor='a^3*d/10',independently_recomputed_point_determinant=det,
        primary_deleted_order=deleted,alternative_deleted_order=alternative_order),indent=2)+'\n')
    inputs=[Path(__file__),ROOT/'code/s61_nonsquare_sparse_kernel_review_20261103_v62.py',
            directory/'report.json',data_path,previous/'source-report.json',previous/'matrices.tar.gz',v3path,
            ROOT/'research/notes/S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW.md']
    result=dict(id='S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62_REVIEW',all_passed=True,
        exact_all63_known_unit_reconstructions=True,normalization_uses_only_a_F_B_C=True,
        exact_all6_QQ_relations_all91_columns=True,independent_integer_arithmetic_after_exact_fraction_denominator_clearing=True,
        primary_unit_triangular_minor='-12*a^3/5',alternative_ad_unit_triangular_minor='a^3*d/10',
        independent_point_minor_mod32003=det,global_free_row_module_rank=57,
        full_localized_row_kernel_generated_by_six_relations=True,
        quadric_product_row_submodule_free_rank=41,abstract_row_module_quotient_free_rank=16,
        all_point_image_rank_lower57_claimed=False,all_point_product_image_rank41_claimed=False,
        free_cokernel_in_rank91_ambient_claimed=False,quartic_or_ambient_removal_claimed=False,
        scope='Global abstract row-module and kernel bases on U=aFBC; specialization image ranks can drop.',
        heavy_search_repeated=False,input_sha256={str(p.relative_to(ROOT)):sha(p.read_bytes())for p in inputs},
        evidence_sha256={p.name:sha(p.read_bytes())for p in args.output.iterdir()if p.is_file()},seconds=time.monotonic()-begin)
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS FULL6 global row basis; exact6x91 QQ identities; independent57minor27632;',result['seconds'])


if __name__=='__main__':main()
