#!/usr/bin/env python3
"""Independently evaluate the actual mixed-slice quadratic certificate."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import sympy as s
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM

from s61_model import S61

ROOT = Path(__file__).resolve().parents[1]


def dm(rows, K):
    return DomainMatrix.from_rep(DDM(rows, (len(rows), len(rows[0])), K))


def verify(output):
    folder = ROOT/'reports/S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62_FINAL'
    report = json.loads((folder/'report.json').read_text())
    certificate = json.loads((folder/'certificate.json').read_text())
    case = report['cases'][0]
    K = s.QQ.algebraic_field(s.sqrt(-2), alias='h')
    a = s.Symbol('a')
    parse = lambda value: K.from_sympy(s.sympify(value.replace('^','**')).subs(a,s.sqrt(-2)))
    point = list(map(parse, certificate['center']))
    directions = [list(map(parse, row)) for row in certificate['directions']]
    assert all(row[3] == -row[19] for row in directions)
    model = S61()
    two = {i for i,w in enumerate(model.reduced['parameter_weights']) if w == 2}
    equations = model.reduced['equations']+[json.loads((ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json').read_text())['universal_new_equation']]
    zero = (0,0,0)
    mons = [(2,0,0),(1,1,0),(1,0,1),(0,2,0),(0,1,1),(0,0,2)]
    variables = []
    for i in range(62):
        value = {zero:point[i]} if point[i] else {}
        for j in range(3):
            if directions[j][i]:
                value[tuple(int(t == j) for t in range(3))] = directions[j][i]
        variables.append(value)
    def multiply(left, right):
        result = {}
        for m,c in left.items():
            for n,d in right.items():
                exponent = tuple(x+y for x,y in zip(m,n))
                if sum(exponent) <= 2:
                    result[exponent] = result.get(exponent,K.zero)+c*d
        return {m:c for m,c in result.items() if c}
    full = [[K.zero]*62 for _ in range(135)]
    raw = []
    for row, equation in enumerate(equations):
        value = {}
        for mon,n,d in equation:
            ids = list(mon)+[61]*(2-sum(i in two for i in mon))
            coefficient = K.convert(n)/K.convert(d)
            term = {zero:coefficient}
            for i in ids:
                term = multiply(term,variables[i])
            for m,c in term.items():
                value[m] = value.get(m,K.zero)+c
            for i,count in Counter(ids).items():
                rest = list(ids);rest.remove(i)
                linear = coefficient*count
                for j in rest:
                    linear *= point[j]
                full[row][i] += linear
        assert not value.get(zero,K.zero)
        assert all(not value.get(tuple(int(t==j) for t in range(3)),K.zero) for j in range(3))
        raw.append([value.get(m,K.zero) for m in mons])
    assert raw == [[parse(v) for v in row] for row in certificate['raw_quadratic_matrix']]
    indices = certificate['global_variables']
    J = [[row[i]-(row[3] if i == 19 else K.zero) for i in indices] for row in full]
    rows, columns = certificate['linear_pivot_rows'], certificate['linear_pivot_columns']
    JD = dm(J,K)
    block = JD.extract(rows,columns)
    assert block.det() == parse(case['jacobian_minor'])
    rawD = dm(raw,K)
    corrected = rawD-JD.extract(range(135),columns)*block.inv()*rawD.extract(rows,range(6))
    assert corrected.to_list() == [[parse(v) for v in row] for row in certificate['corrected_quadratic_matrix']]
    basis_rows = certificate['quadratic_basis_rows']
    basis = corrected.extract(basis_rows,range(6))
    minor = basis.extract(range(5),certificate['quadratic_minor_columns']).det()
    assert minor == parse(case['quadratic_minor'])
    for matrix in [full,J,raw,corrected.to_list(),directions]:
        assert all(int(c.denominator)%32003 for row in matrix for v in row for c in v.to_list())
    residues = []
    for root in [3083,28920]:
        def residue(v):
            result = 0
            for c in v.to_list():
                result = (result*root+int(c.numerator)*pow(int(c.denominator),-1,32003))%32003
            return result
        values = [residue(block.det()),residue(minor)]
        assert all(values)
        residues.append(dict(h_root=root,jacobian_and_quadratic_minor=values))
        assert 'HS=(1+3*T+T^2)/(1)' in case['outputs'][f'root{root}']
    assert 'HS=(1+3*T+T^2)/(1)' in case['outputs']['QQ']
    # Independently check the compact ideal normal form's five standard terms.
    x,y,z = s.symbols('x y z')
    for root, domain in [(s.sqrt(-2),K),(3083,s.GF(32003)),(28920,s.GF(32003))]:
        polys = []
        for row in basis.to_list():
            if domain is K:
                coefficients = list(map(K.to_sympy,row))
            else:
                def reduce_coefficient(v):
                    total = 0
                    for c in v.to_list():
                        total=(total*root+int(c.numerator)*pow(int(c.denominator),-1,32003))%32003
                    return total
                coefficients=list(map(reduce_coefficient,row))
            polys.append(sum(c*x**m[0]*y**m[1]*z**m[2] for c,m in zip(coefficients,mons)))
        gb=s.groebner(polys,x,y,z,domain=domain)
        leading=[p.monoms()[0] for p in gb.polys]
        standard=[]
        for degree in range(5):
            for i in range(degree+1):
                for j in range(degree-i+1):
                    m=(i,j,degree-i-j)
                    if not any(all(a>=b for a,b in zip(m,n)) for n in leading):standard.append(m)
        assert len(standard)==5 and all(sum(m)<=2 for m in standard)
    for path,digest in report['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    paths=[Path(__file__),folder/'report.json',folder/'certificate.json',
           ROOT/'code/s61_marked_sigma_old_double_mixed_certificate_20261103_v62.py']
    result=dict(id='S61_MARKED_SIGMA_OLD_DOUBLE_MIXED_CERTIFICATE_20261103_v62_REVIEW',all_passed=True,
                exact_mixed_affine_slice_in_all62_directions=True,
                independent_original135_multivariate_Taylor_coefficients_match=True,
                independent_Jacobian_and_implicit_quadratic_correction_match=True,
                all_coefficients_p_integral=True,unit_witness_residues=residues,
                independent_QQ_and_two_Fp_Groebner_length=5,hilbert_function=[1,3,1],
                input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()for path in paths})
    with output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print('PASS FULL135 MIXED SLICE; QQ + BOTH GF LENGTH5',residues)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();verify(args.output)
