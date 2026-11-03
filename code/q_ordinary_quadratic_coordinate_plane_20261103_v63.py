#!/usr/bin/env python3
"""Exact QQ 13-plane contained in the complete ordinary quadratic obstruction cone."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cone',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);start=time.monotonic()
    archive=args.cone/'quadrics.json.gz';data=json.load(gzip.open(archive,'rt'))
    metadata=json.loads((args.cone/'extraction.json').read_text())
    assert sha(archive)==metadata['quadratic_archive_sha256']
    names=data['variables'];assert len(names)==35
    active_names=['u','z3','z5','z8','z13','z14','z15','z16','z17','z24','z25','z26','z28']
    active={names.index(name) for name in active_names};assert len(active)==13
    rows=data['all350_quadratic_rows'];assert len(rows)==350
    supports=[];loops=set()
    for row in rows:
        seen=set()
        for i,j,coefficient in row:
            assert Fraction(coefficient) and i<=j
            assert (i,j) not in seen;seen.add((i,j))
            assert not (i in active and j in active), ('nonzero restricted quadratic',i,j,coefficient)
            supports.append((i,j))
            if i==j:loops.add(i)
    # Optional stronger combinatorial statement: no larger coordinate linear
    # space can be contained.  There are only 15 eligible axes after excluding
    # squared monomials; all 2^15 subsets are explicitly checked.
    eligible=sorted(set(range(35))-loops);assert len(eligible)==15
    edges={(i,j) for i,j in supports if i!=j and i in eligible and j in eligible}
    position={axis:i for i,axis in enumerate(eligible)}
    masks=[(1<<position[i])|(1<<position[j]) for i,j in sorted(edges)]
    maximum=0;maximizers=[]
    for subset in range(1<<len(eligible)):
        size=subset.bit_count()
        if size<maximum:continue
        if any(subset&edge==edge for edge in masks):continue
        if size>maximum:maximum=size;maximizers=[]
        maximizers.append([eligible[i] for i in range(len(eligible)) if subset>>i&1])
    assert maximum==13
    report=dict(id='Q_ORDINARY_QUADRATIC_COORDINATE_PLANE_20261103_V63',all_passed=True,
                coefficient_field='QQ',all350_quadrics_restricted_identically_zero=True,
                active_variables=active_names,active_variable_indices=sorted(active),
                zero_variables=[name for i,name in enumerate(names) if i not in active],
                closed_coordinate_linear_subspace_dimension=13,
                ideal_inclusion='I_quadratic subset (all22 inactive variables)',
                quadratic_cone_dimension_lower_bound_over_QQ_and_C=13,
                squared_axes_excluded=len(loops),eligible_axes=len(eligible),
                subsets_exhaustively_checked=1<<len(eligible),
                maximum_coordinate_subspace_dimension=maximum,
                all_maximal_dimension_coordinate_subspaces=[[names[i] for i in subset] for subset in maximizers],
                ordinary_quadratic_cone_cannot_give_dimension_upper12=True,
                actual_formal_model_dimension13_or_component_claimed=False,
                higher_initial_equations_still_needed=True,
                input_sha256={str(p.resolve().relative_to(ROOT)):sha(p) for p in
                              [Path(__file__),archive,args.cone/'extraction.json']},
                seconds=time.monotonic()-start)
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS exact QQ coordinate13-plane; maximum coordinate plane',maximum,
          'number of maximizers',len(maximizers),'seconds',report['seconds'])


if __name__=='__main__':main()
