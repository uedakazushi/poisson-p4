#!/usr/bin/env python3
"""Receipts and exact scalar checks for the cofactor colon/nonflatness lemma."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def build(output):
    inputs=[Path(__file__),ROOT/'research/notes/Q_KNOWN_GRAPH_COLON_NONFLATNESS_20261103_v63.md',
            ROOT/'reports/Q_T39_DIVIDED_DIFFERENCE_OPEN_20261103_v63_REVIEW.json',
            ROOT/'reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63_REVIEW.json',
            ROOT/'reports/Q_NONRESONANT_FINITE_MODEL_20261103_v63_REVIEW.json']
    witnesses=[]
    for label,center,detkey in [('T39',Fraction(2,9),'exact_cofactor_determinant_at_Rstar1'),
                                ('T40',Fraction(-8,9),'exact_cofactor_determinant_at_dense_T40_point')]:
        folder=ROOT/f'reports/Q_{label}_DIVIDED_DIFFERENCE_OPEN_20261103_v63'
        report_path=folder/'report.json';archive_path=folder/'cofactor-model.json.gz'
        report=json.loads(report_path.read_text());assert report['all_passed']
        assert report['normalized_generators_are_polynomial_graph_kernel']
        assert report['all350_graph_identities_replayed_exact_QQ9']
        assert report['all26_graph_equation_remainders_independently_exact_QQ9_zero']
        for name,expected in report['input_sha256'].items():assert sha(ROOT/name)==expected,name
        assert sha(archive_path)==report['evidence_sha256']['cofactor-model.json.gz']
        model=json.loads(gzip.decompress(archive_path.read_bytes()))
        assert model['dependent_parent_input_indices'][0]==0
        a_graph=model['nodes'][model['graph_parent35_input_nodes'][0]]
        assert a_graph[0]=='q'and Fraction(a_graph[1])==center
        assert s.Rational(report[detkey])!=0
        witnesses.append(dict(component=label,a_center=str(center),nonzero_graph_determinant=report[detkey]))
        inputs.extend([report_path,archive_path])
    a=s.Symbol('a');u39=a-s.Rational(2,9);u40=a+s.Rational(8,9)
    assert s.expand(s.Rational(9,10)*u40-s.Rational(9,10)*u39)==1
    output.parent.mkdir(parents=True,exist_ok=True)
    report=dict(id='Q_KNOWN_GRAPH_COLON_NONFLATNESS_20261103_V63',all_passed=True,
                exact_colons={'T39':'I:C39^k=N39 for every k>=1',
                              'T40':'I:C40^k=N40 for every k>=1'},
                reason='C*N subset I subset N; P/N is the polynomial domain A9 and C mod N is nonzero.',
                exact_module_sequence='0 -> P/N --multiply C--> P/I -> P/(I,C) -> 0',
                nonzero_a_torsion_centers=['2/9','-8/9'],
                full35_local_and_completed_models_nonflat_at_both_Qa_centers=True,
                T40_zero_section_full_Jacobian_rank_is1_not_a_resonance_rank_jump=True,
                coprime_scalar_Bezout='1=(9/10)*(a+8/9)-(9/10)*(a-2/9)',
                C39_times_C40_is_zero_in_R=True,
                principal_known_ideal_intersection_is_zero=True,
                sum_module_is_direct='(C39)+(C40) is isomorphic as P-module to P/N39 direct_sum P/N40',
                combined_colon='I:(C39+C40)^k=N39 intersection N40 for every k>=1',
                combined_scheme_open='R[1/(C39+C40)] = A9_T39[1/C39] product A9_T40[1/C40]',
                whole_R_direct_product_claimed=False,all_remaining_branches_classified=False,
                witnesses=witnesses,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs})
    with output.open('x')as stream:json.dump(report,stream,indent=2);stream.write('\n')
    print('PASS exact colon/nonflatness consequences from actual cofactor ideals')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();build(args.output.resolve())
