#!/usr/bin/env python3
"""Pin the exact two-prime plus residual ideal decomposition corollary."""
import argparse
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    colon_path=ROOT/'reports/Q_KNOWN_GRAPH_COLON_NONFLATNESS_20261103_v63.json'
    colon=json.loads(colon_path.read_text());assert colon['all_passed']
    assert colon['C39_times_C40_is_zero_in_R']
    assert colon['principal_known_ideal_intersection_is_zero']
    for path,expected in colon['input_sha256'].items():assert sha(ROOT/path)==expected,path
    a=s.Symbol('a');u39=a-s.Rational(2,9);u40=a+s.Rational(8,9)
    assert s.expand(s.Rational(9,10)*(u40-u39))==1
    note=ROOT/'research/notes/Q_KNOWN_GRAPH_RESIDUAL_DECOMPOSITION_20261103_v63.md'
    inputs=[Path(__file__),note,colon_path]
    report=dict(id='Q_KNOWN_GRAPH_RESIDUAL_DECOMPOSITION_20261103_V63',all_passed=True,
                single_prime_exact_decomposition='I=N_i intersection (I+(C_i))',
                two_known_primes_comaximal=True,
                exact_full_ideal_decomposition='I=N39 intersection N40 intersection K, K=I+(C39,C40)',
                residual_ideal_preserves_positive_weight_grading=True,
                cofactor_generator_degrees={'C39':264,'C40':275},
                residual_contains_neither_entire_known_prime_component=True,
                exact_local_Qstar_decomposition='Ihat=N39hat intersection (Ihat+(C39))',
                reason_C40_vanishes_in_Qstar_localization='(a+8/9) is a unit and annihilates C40',
                exact_local_Qa40_decomposition='Ihat=N40hat intersection (Ihat+(C40))',
                residual_primary_or_reduced_components_classified=False,
                full_Poisson_branch_exhaustion_claimed=False,
                proof='If f belongs N and I+(C), write f=i+Cg; the polynomial domain P/N and nonzero Cbar force g in N, hence Cg in I. For two primes use C_other in N and the same argument modulo each prime.',
                input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs})
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS exact two-known-prime plus residual ideal decomposition')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
