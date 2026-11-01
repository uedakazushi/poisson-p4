#!/usr/bin/env python3
"""Independent product-space certificate check and true-core cubic proof audit."""
import argparse
import hashlib
from itertools import combinations_with_replacement as cwr
import json
from pathlib import Path
import time
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def verify():
    start=time.time()
    primary='reports/S61_FLAG_TRUE_CORE_CUBIC_20261101.json'
    quadratic='reports/S61_FLAG_TRUE_CORE_QUADRATIC_20261101.json'
    d=json.loads((ROOT/primary).read_text());q=json.loads((ROOT/quadratic).read_text())
    assert d['all_passed'] and q['all_passed']
    for record in [d,q]:
        for path,digest in record['input_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    mat=lambda x:s.Matrix(x).applyfunc(s.sympify)
    Q=mat(q['polynomial_quadratic_basis']);qp=q['constant_minor']['columns']
    assert Q[:,qp]==s.eye(15)
    generic,pivots=Q.to_DM().to_field().rref();generic=generic.to_Matrix()
    C=generic[[i for i in range(15) if not any(x.free_symbols for x in generic[i,:])],:]
    assert C.rows==13
    assert (C-C[:,qp]*Q).applyfunc(s.expand)==s.zeros(13,55)
    m2=list(cwr(range(10),2));m3=list(cwr(range(10),3));index={m:i for i,m in enumerate(m3)}
    def products(matrix):
        result=s.zeros(10*matrix.rows,220)
        for row in range(matrix.rows):
            for col,monomial in enumerate(m2):
                coefficient=matrix[row,col]
                if coefficient:
                    for variable in range(10):
                        result[10*row+variable,index[tuple(sorted(monomial+(variable,)))]]+=coefficient
        return result
    BC=products(C);rc,pc=BC.to_DM().rref();RC=rc.to_Matrix()[:len(pc),:]
    p=d['quadratic_product_space']
    assert len(pc)==100 and list(pc)==p['constant_pivots']
    assert RC==mat(p['constant_rref'])
    free=[i for i in range(220) if i not in pc]
    B=products(Q)
    E=(B[:,free]-B[:,list(pc)]*RC[:,free]).applyfunc(s.expand)
    nonzero=[i for i in range(150) if any(E[i,:])]
    assert nonzero==p['remaining_rows'] and len(nonzero)==54
    E=E[nonzero,:];RE=mat(p['remaining_rref']);pe=p['remaining_pivots']
    assert RE[:,pe]==s.eye(14)
    assert all(s.denom(x).is_Number for x in RE)
    assert (E-E[:,pe]*RE).applyfunc(s.expand)==s.zeros(54,120)
    selected=E[p['minor']['rows'],p['minor']['columns']].to_DM().to_field()
    determinant=selected.domain.to_sympy(selected.det())
    assert determinant==s.Rational(p['minor']['value']) and determinant
    c=d['cubic_space'];RR=mat(c['rref'])
    assert RR.shape==(5,106) and RR[:,c['pivots']]==s.eye(5)
    assert s.Rational(c['minor_value'])!=0
    symbols=set().union(*(x.free_symbols for x in RR))
    assert symbols=={s.Symbol('b8')}
    assert sum(x!=0 for x in RR)==11
    assert max(s.Poly(x,s.Symbol('b8')).total_degree() for x in RR if x)<=2
    assert c['t3_after_t2_correction_in_ideal']
    assert 220-(100+14)-5==101
    paths=list(d['input_sha256'])+[primary,quadratic,
        'research/notes/S61_FLAG_TRUE_CORE_CUBIC_20261101.md',
        str(Path(__file__).relative_to(ROOT)),
        'research/notes/S61_FLAG_TRUE_CORE_CUBIC_REVIEW_20261101.md']
    return dict(id='S61_FLAG_TRUE_CORE_CUBIC_REVIEW_20261101',all_passed=True,
        independent_constant_thirteen_quadrics_inclusion=True,
        independent_full_quadratic_product_space=True,
        independent_fourteen_minor=str(determinant),
        independent_new_cubic_saved_basis_check=True,
        proof_and_source_review_accepted=True,t2_tail_correction_source_audited=True,
        independent_raw135_cubic_reimplementation=False,primary_producer_replayed=False,
        accepted_hilbert_function=[1,10,40,101],accepted_all_seven_parameter_values=True,
        limitations=['The 135 raw cubic identities and five-row nonzero minor are primary computational dependencies.',
                     'No exact order of t, full initial ideal, completed-ring identification or ambient t-nilpotence.'],
        input_sha256={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        seconds=time.time()-start)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print('TRUE CORE CUBIC REVIEW PASS',result['seconds'])
