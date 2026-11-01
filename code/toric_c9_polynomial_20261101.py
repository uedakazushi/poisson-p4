#!/usr/bin/env python3
"""A pencil certificate that the rational C9 covariant is polynomial globally.

The key algebraic fact is that the rank-drop locus of B_beta has codimension
at least two. The accompanying note supplies the normality extension argument.
"""
import argparse,hashlib,json,random
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from weight1_structure import constraints,minor_witness
from toric_algebra_fiber import TRIPLES,multiplication_from_jet,associativity_quadrics
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    basis=s.Matrix.hstack(*s.Matrix.vstack(*constraints()).nullspace())
    cubic_zero=dict.fromkeys(TRIPLES,0)
    points=[list(range(1,17)),[(-1)**i*(i*i+3*i+5) for i in range(16)]]
    cubic_matrices=[];pure=[]
    for triple in TRIPLES:
        coeff=dict(cubic_zero);coeff[triple]=1
        mc=multiplication_from_jet(s.zeros(24,1),coeff);cubic_matrices.append(mc)
        pure.append(s.Matrix(associativity_quadrics(mc)[1]))
    matrices=[]
    for point in points:
        m0=multiplication_from_jet(basis*s.Matrix(point),cubic_zero)
        constant=s.Matrix(associativity_quadrics(m0)[1]);columns=[]
        for mc,h in zip(cubic_matrices,pure,strict=True):
            mixed=[a+b for a,b in zip(m0,mc,strict=True)]
            columns.append(s.Matrix(associativity_quadrics(mixed)[1])-constant-h)
        matrices.append(s.Matrix.hstack(*columns))
    B0,B1=matrices;t=s.Symbol('t');domain=s.QQ.poly_ring(t);minors=[]
    for seed in (0,1):
        order=list(range(96));random.Random(seed).shuffle(order)
        pivots=DomainMatrix.from_Matrix(B0[order,:].T).convert_to(s.QQ).rref()[1]
        assert len(pivots)==20
        rows=[order[i] for i in pivots]
        minor=DomainMatrix.from_Matrix((B0+t*B1)[rows,:]).convert_to(domain).det()
        poly=s.Poly(minor.as_expr(),t,domain=s.QQ);assert poly.degree()==20
        minors.append((rows,poly))
    left,right,gcd=s.gcdex(minors[0][1],minors[1][1]);assert gcd==s.Poly(1,t,domain=s.QQ)
    assert left*minors[0][1]+right*minors[1][1]==gcd
    # A degree20 minor has nonzero leading coefficient, which proves rank20
    # at the point at infinity of this SAME projective line.
    determinant_infinity=B1[minors[0][0],:].det(method='domain-ge')
    assert determinant_infinity==minors[0][1].LC()!=0
    report=dict(id='TORIC_C9_POLYNOMIAL_20261101',all_passed=True,field='QQ[t]',
        base_weight_coordinates=points[0],direction_weight_coordinates=points[1],
        construction='Exact polarization E(m0+mc)-E(m0)-E(mc) of the homogeneous quadratic96 associativity equations',
        B_base=[[str(v) for v in row] for row in B0.tolist()],B_direction=[[str(v) for v in row] for row in B1.tolist()],
        maximal_minors=[dict(rows=rows,polynomial=str(p.as_expr()),degree=p.degree()) for rows,p in minors],
        bezout_coefficients=[str(left.as_expr()),str(right.as_expr())],gcd='1',
        rank_at_infinity=20,minor_determinant_at_infinity=str(determinant_infinity),
        rank_drop_codimension_at_least_two=True,
        global_C9_is_homogeneous_polynomial_degree9=True,
        proof_note='research/notes/TORIC_C9_POLYNOMIAL_20261101.md',
        limitations=['The full16-variable C9 expansion is not generated; global regularity follows from uniqueness and extension across codimension two.',
            'No fullS61 generic-fiber exhaustion.'],
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'code/toric_algebra_fiber.py',ROOT/'code/weight1_structure.py',ROOT/'reports/TORIC_JACOBIAN_TWO_TORSION_20261101.json')})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as out:json.dump(report,out,ensure_ascii=False,indent=2);out.write('\n')
    print(json.dumps(dict(all_passed=True,pencil_has_no_rank_drop=True,global_C9_polynomial=True)))

if __name__=='__main__':main()
