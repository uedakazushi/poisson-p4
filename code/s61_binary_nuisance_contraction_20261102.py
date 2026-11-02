#!/usr/bin/env python3
"""Exact contracting nuisance action and its finite-type support lemma."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61, BASIS, X

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'research/notes/S61_BINARY_NUISANCE_CONTRACTION_20261102.md'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify():
    start = time.time()
    model = S61()
    base_weights = model.reduced['parameter_weights']
    torus = [0, 1, -1, -1, 1]
    raw = []
    for index in model.certificate['raw_free_indices']:
        pair, monomial = BASIS[index]
        exponents = s.Poly(monomial, *X).monoms()[0]
        raw.append(sum(torus[i] for i in pair)-sum(a*b for a, b in zip(torus, exponents)))
    weights = [a-b for a, b in zip(base_weights, raw)]
    extra_path = ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    extra = json.loads(extra_path.read_text())
    assert extra['all_passed']
    equations = model.reduced['equations']+[extra['universal_new_equation']]
    equation_weights = []
    qdegrees = []
    for polynomial in equations:
        values = {sum(weights[i] for i in monomial) for monomial, _, _ in polynomial}
        assert len(values) == 1
        equation_weights.append(next(iter(values)))
        degree = max(sum(base_weights[i] == 2 for i in monomial) for monomial, _, _ in polynomial)
        assert degree <= 2
        qdegrees.append(degree)
    assert len(equations) == 135
    assert weights[29] == 6 and base_weights[29] == 2
    a,e,d,p,q,r,h,k = s.symbols('a e d p q r h k')
    lam = s.Symbol('lam', nonzero=True)
    source = {i:s.S.Zero for i,w in enumerate(base_weights) if w == 1}
    source.update({3:a,13:(e+1)/5,34:s.S.One,35:d,48:e,49:-5,
                   24:p,25:q,26:r,36:h,47:k})
    source = {i:s.sympify(value) for i,value in source.items()}
    nuisance = [p,q,r,h,k]
    nuisance_weights = [weights[i] for i in [24,25,26,36,47]]
    assert nuisance_weights == [4,2,4,2,2]
    substitution = {v:lam**w*v for v,w in zip(nuisance,nuisance_weights)}
    for i,value in source.items():
        assert s.expand(lam**weights[i]*value-value.subs(substitution, simultaneous=True)) == 0
    assert all(weights[i] == 0 for i in [3,13,34,35,48,49])
    assert all(w > 0 for w in nuisance_weights)
    dependencies = [Path(__file__), NOTE, ROOT/'code/s61_model.py', extra_path,
                    model.source/'symplectic_slice_certificate.json',
                    model.source/'symplectic_slice_reduced_equations.json']
    return dict(
        id='S61_BINARY_NUISANCE_CONTRACTION_20261102', all_passed=True,
        parameter_weights=weights, raw_torus_weights=raw,
        equations_checked=135, equation_weights=equation_weights,
        equation_weight_counts=dict(Counter(equation_weights)),
        all_equations_individually_weight_homogeneous=True,
        homogenization_q_degree=2, homogenizing_coordinate_global_weight=0,
        projective_section_parameter_index_zero_based=29,
        section_coordinate_weight=6, normalized_chart_t_weight=-6,
        nuisance_variables=list(map(str,nuisance)), nuisance_weights=nuisance_weights,
        core_variables=['a','e','d'], core_weights=[0,0,0],
        source_beta={f'u{i+1}':str(value) for i,value in source.items()},
        source_action_identity_checked=True,
        contraction_support_lemma='If the full relative finite-type scheme has t-nilpotent local rings along the zero-nuisance section over a core locus, then its t-saturation has no section points over any nuisance values on that same core locus.',
        proof_note=str(NOTE.relative_to(ROOT)),
        limitations=['Hypothesis is t-nilpotence in the full relative scheme, not just a transverse or fixed-beta quotient.',
                     'No claim that the action extends to lambda=0 on the whole scheme.',
                     'No new core rank coverage or ambient S61 exhaustion by this lemma alone.'],
        input_sha256={str(path.relative_to(ROOT)):sha(path) for path in dependencies},
        seconds=time.time()-start)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();assert not args.output.exists()
    result=verify()
    with args.output.open('x') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({key:result[key] for key in ['all_passed','equations_checked','nuisance_weights','seconds']}))
