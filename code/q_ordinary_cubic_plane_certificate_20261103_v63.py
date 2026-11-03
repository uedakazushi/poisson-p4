#!/usr/bin/env python3
"""Replay the constant syzygies and certify the explicit false-plane divisor."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def subtract(row,other,scale):
    for mon,c in other.items():
        value=row.get(mon,Fraction(0))-scale*c
        if value:row[mon]=value
        else:row.pop(mon,None)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cubic',type=Path,default=ROOT/'reports/Q_ORDINARY_CUBIC_INITIAL_20261103_v63')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    archive=args.cubic/'cubic-initials.json.gz';data=json.load(gzip.open(archive,'rt'))
    oldpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz'
    old=json.load(gzip.open(oldpath,'rt'));names=data['variables']
    rawQ=[{(i,j):Fraction(c) for i,j,c in row} for row in old['all350_quadratic_rows']]
    rawC=[{tuple(mon):Fraction(c) for mon,c in row} for row in data['raw350_cubic_terms']]
    expected={row['source']:{tuple(mon):Fraction(c) for mon,c in row['terms']}
              for row in data['cubic_initial_relations']}
    bases={};carried={};seen=set()
    for operation in data['exact_constant_row_elimination']:
        source=operation['source'];assert source not in seen;seen.add(source)
        q=rawQ[source].copy();c=rawC[source].copy()
        for pivottext,coefficient in operation['reductions']:
            pivot=tuple(pivottext);coefficient=Fraction(coefficient)
            assert min(q)==pivot and q[pivot]==coefficient
            subtract(q,bases[pivot],coefficient);subtract(c,carried[pivot],coefficient)
        if 'new_pivot' in operation:
            pivot=tuple(operation['new_pivot']);normalizer=Fraction(operation['normalizer'])
            assert min(q)==pivot and q[pivot]==normalizer
            assert normalizer.numerator%32003 and normalizer.denominator%32003
            bases[pivot]={mon:value/normalizer for mon,value in q.items()}
            carried[pivot]={mon:value/normalizer for mon,value in c.items()}
        else:
            assert not q and c==expected[source]
    assert seen==set(range(350)) and len(bases)==57 and len(expected)==293
    active={names.index(name) for name in data['plane_active_variables']}
    assert all(not all(i in active for i in mon) for row in rawQ for mon in row)
    restricted={source:{mon:c for mon,c in poly.items() if all(i in active for i in mon)}
                for source,poly in expected.items()}
    f={tuple(sorted([names.index('z8'),names.index('z8'),names.index(name)])):Fraction(coefficient)
       for name,coefficient in [('z14',12),('z17',1),('z3',6)]}
    for source,poly in restricted.items():
        scalar=Fraction(3,16) if source==234 else Fraction(-1,16) if source==349 else Fraction(0)
        assert poly==({mon:scalar*c for mon,c in f.items()} if scalar else {})
    oldplane=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/coordinate-plane.json'
    oldproof=json.loads(oldplane.read_text());assert oldproof['all_passed']
    assert oldproof['maximum_coordinate_subspace_dimension']==13
    assert len(oldproof['all_maximal_dimension_coordinate_subspaces'])==1
    assert set(oldproof['all_maximal_dimension_coordinate_subspaces'][0])=={names[i] for i in active}
    report=dict(id='Q_ORDINARY_CUBIC_PLANE_CERTIFICATE_20261103_V63',all_passed=True,coefficient_field='QQ',
                all350_recorded_constant_row_operations_independently_replayed=True,
                quadratic_basis_rank=57,true_cubic_constant_relations=293,
                all_quadratic_restrictions_to_plane_zero=True,
                exactly_two_nonzero_cubic_restrictions=True,
                source234='(3/16)*z8^2*(12*z14+z17+6*z3)',
                source349='(-1/16)*z8^2*(12*z14+z17+6*z3)',
                restricted_ideal_generated_by_degree2_and3_initial_parts='(z8^2*(12*z14+z17+6*z3))',
                restricted_scheme_dimension=12,restricted_scheme_degree=3,
                restricted_scheme_components_with_multiplicity=[dict(equation='z8=0',multiplicity=2),
                    dict(equation='12*z14+z17+6*z3=0',multiplicity=1)],
                maximum_coordinate_linear_subspace_dimension_in_through3_initial_scheme=12,
                global_initial_cone_or_actual_local_dimension_not_inferred=True,
                input_sha256={str(path.resolve().relative_to(ROOT)):sha(path) for path in [Path(__file__),archive,oldpath,oldplane]},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact constant row replay and plane divisor; seconds',report['seconds'])


if __name__=='__main__':main()
