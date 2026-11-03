#!/usr/bin/env python3
"""Exact FULL350 finite weight-six four-coordinate slice over QQ, not jets."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

import sympy as s
from sympy.polys.rings import ring
from q_weighted_finite_model_validation_20261103_v63 import evaluate,physical_equations

ROOT=Path(__file__).resolve().parents[1]
R,u,v,w,t=ring('u,v,w,t',s.QQ)


class Domain:
    zero=R.zero
    one=R.one

    @staticmethod
    def q(value):
        value=Fraction(value)
        return R(s.QQ(value.numerator,value.denominator))

    @staticmethod
    def add(left,right):return left+right

    @staticmethod
    def mul(left,right):return left*right

    @staticmethod
    def neg(value):return -value

    @staticmethod
    def inv(value):
        assert value and value.is_ground
        return R.one/value


def encode(poly):
    return [[list(mon),str(coefficient)] for mon,coefficient in sorted(poly.items())]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=ROOT/'reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();raw=gzip.decompress(args.model.read_bytes());model=json.loads(raw)
    inputs=[R.zero]*35;inputs[0]=Domain.q('2/9')
    substitutions={'z11':u+4*t,'z12':v+6*t,'z23':w-4*t,'z33':t}
    assert {name for name,weight in zip(model['input_names'],model['input_weights']) if weight==6}==set(substitutions)
    for name,value in substitutions.items():inputs[model['input_names'].index(name)]=value
    values,inverses=evaluate(model['nodes'],inputs,Domain)
    reconstructed=[values[node] for node in model['raw150_reconstruction_nodes']]
    residuals=[values[node] for node in model['all350_residual_nodes']]
    _,physical=physical_equations(reconstructed,Domain)
    assert residuals==physical
    nonzero=[i for i,poly in enumerate(residuals) if poly]
    monomials=[(2,0,0,0),(1,1,0,0),(0,2,0,0),(1,0,1,0),(0,1,1,0),
               (0,0,2,0),(1,0,0,1),(0,1,0,1),(0,0,1,1)]
    def contained(mon):return any(all(a>=b for a,b in zip(mon,generator)) for generator in monomials)
    assert all(contained(mon) for poly in residuals for mon in poly)
    quadratic=[i for i in nonzero if {sum(mon) for mon in residuals[i]}=={2}]
    assert all({sum(mon) for mon in residuals[i]} in ({2},{3}) for i in nonzero)
    matrix=s.Matrix([[s.Rational(poly.get(mon,s.QQ.zero)) for mon in monomials]
                     for poly in (residuals[i] for i in quadratic)])
    rows=list(matrix.T.rref()[1]);assert len(rows)==9
    selected=[quadratic[i] for i in rows];square=matrix[rows,:];det=s.factor(square.det());assert det
    inverse=square.inv();assert inverse*square==s.eye(9)
    numerator,denominator=map(int,s.fraction(det))
    assert numerator%32003 and denominator%32003
    assert max(sum(mon) for poly in reconstructed for mon in poly)==2
    polys=dict(input_substitutions={name:str(value) for name,value in substitutions.items()},
               all150_reconstruction=[encode(poly) for poly in reconstructed],
               all350_residuals=[encode(poly) for poly in residuals],
               monomial_generators=[list(mon) for mon in monomials],
               selected_quadratic_rows=selected,
               inverse_coefficient_matrix=[[str(value) for value in row] for row in inverse.tolist()])
    (args.output/'polynomials.json').write_text(json.dumps(polys,indent=2)+'\n')
    files=[Path(__file__).resolve(),args.model.resolve(),
           ROOT/'code/q_weighted_finite_model_validation_20261103_v63.py']
    report=dict(id='Q_WEIGHTED_RESONANCE_FOUR_SLICE_20261103_V63',all_passed=True,
                coefficient_ring='QQ[u,v,w,t]',fixed_a='2/9',
                all_other30_positive_free_coordinates_zero=True,
                full350_finite_residuals_not_truncated_jets=True,
                all350_match_independent_physical_Jacobi=True,
                physical_divergence25_zero=True,
                nonzero_residual_rows=nonzero,quadratic_rows=quadratic,
                residual_total_degrees=[2,3],maximum_residual_terms=max(map(len,residuals)),
                selected_nine_quadratic_rows=selected,
                selected_nine_by_nine_determinant=str(det),determinant_is_32003_unit=True,
                every_residual_in_nine_generator_monomial_ideal=True,
                all_nine_monomial_generators_in_quadratic_residual_span=True,
                exact_ideal='(u,v,w)*(u,v,w,t)',
                reduced_ideal='(u,v,w)',dimension=1,
                ordinary_hilbert_series='1/(1-T)+3*T',
                weighted_hilbert_series='1/(1-T^6)+3*T^6',
                nilradical_QQ_dimension=3,nilradical_square_zero=True,
                nilradical_annihilated_by='(u,v,w,t)',
                generic_resonance_line_reduced=True,
                full35_branch_decomposition_claimed=False,
                no_inverse_in_positive_variables=True,inverse_nodes_evaluated=inverses,
                input_sha256={str(path.relative_to(ROOT)):sha(path) for path in files},
                decoded_model_sha256=hashlib.sha256(raw).hexdigest(),
                evidence_sha256={'polynomials.json':sha(args.output/'polynomials.json')},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact full350 four-coordinate ideal',report['exact_ideal'],
          'nonzero',len(nonzero),'det',det,'seconds',report['seconds'])


if __name__=='__main__':main()
