#!/usr/bin/env python3
"""First-jet generation and a square-free degree-24 discriminant restriction.

The geometric incidence proof is in the accompanying note. Finite-field gcd
is used only to certify square-freeness of an explicit rational polynomial.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s
from weight1_structure import constraints, form_matrices, minor_witness

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    basis=s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    b=s.symbols('b0:16');matrices=form_matrices(basis*s.Matrix(b))
    a,c,d=s.symbols('a c d');u=s.Matrix([1,a,0,c]);v=s.Matrix([0,c,1,d])
    section=[s.expand((u.T*sum((w[k]*matrices[k] for k in range(4)),s.zeros(4))*v)[0]) for w in (u,v)]
    first=s.Matrix([*section,*[s.diff(f,t) for f in section for t in (a,c,d)]]).subs({a:0,c:0,d:0}).jacobian(b)
    witness=minor_witness(first);assert witness['rank']==8
    source=ROOT/'reports/PRIMITIVE_CURVE_DOMINANT_SLICE_20261101.json'
    data=json.loads(source.read_text());t=s.Symbol('t')
    base=[1,2,3,4,5,6];direction=[2,-1,3,5,7,11];invariants=[]
    for key in ('c4_sparse','c6_sparse'):
        poly=s.Poly(0,t,domain=s.QQ)
        for mon,coefficient in data[key]:
            term=s.Poly(s.Rational(coefficient),t,domain=s.QQ)
            for i,e in enumerate(mon):
                if e:term*=s.Poly(base[i]+direction[i]*t,t,domain=s.QQ)**e
            poly+=term
        invariants.append(poly)
    delta=(invariants[0]**3-invariants[1]**2).quo_ground(1728)
    assert [p.degree() for p in invariants]==[8,12] and delta.degree()==24
    denominator,integral=delta.clear_denoms(convert=True)
    prime=32003;assert int(denominator)%prime and int(integral.LC())%prime
    reduced=s.Poly(integral.as_expr(),t,modulus=prime)
    left,right,gcd=s.gcdex(reduced,reduced.diff());assert gcd==s.Poly(1,t,modulus=prime)
    assert left*reduced+right*reduced.diff()==gcd
    assert s.gcd(delta,delta.diff()).degree()==0
    # The invariant line really comes from the prescribed beta-space slice.
    assert s.Matrix(data['weight_basis_24_by16'])==basis
    report=dict(id='PRIMITIVE_CURVE_DISCRIMINANT_20261101',all_passed=True,
        field='QQ with an explicit finite-field Bezout square-freeness certificate',
        lagrangian_chart_columns=[list(map(str,u)),list(map(str,v))],
        section_components=list(map(str,section)),
        first_jet_matrix=[[str(x) for x in row] for row in first.tolist()],
        first_jet_rank_witness=witness,
        slice_line_base=base,slice_line_direction=direction,
        c4=str(invariants[0].as_expr()),c6=str(invariants[1].as_expr()),delta=str(delta.as_expr()),
        discriminant_degree=24,denominator=str(denominator),integral_discriminant=str(integral.as_expr()),
        square_free_witness=dict(prime=prime,leading_coefficient_mod_prime=int(integral.LC())%prime,
            bezout_for_delta=str(left.as_expr()),bezout_for_derivative=str(right.as_expr()),gcd='1'),
        rational_gcd_degree=0,
        geometric_conclusion='The canonical curve discriminant in W16 is a geometrically irreducible reduced homogeneous hypersurface of degree24, by the critical-zero incidence proof.',
        limitations=['The geometric proof uses the previously established section/Pfaffian identification and Fisher-Sadek Theorem2.1(ii).',
            'No claim about all full-S61 fibers or regular extension of C9.',
            'Irreducibility over C is proved by incidence, not by factoring a rational or finite-field line polynomial.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'code/weight1_structure.py',source)},
        proof_note='research/notes/PRIMITIVE_CURVE_DISCRIMINANT_20261101.md')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,first_jet_rank=8,discriminant_degree=24,line_square_free=True)))

if __name__=='__main__':main()
