#!/usr/bin/env python3
"""Restrict only certified true quartic initials to the two frozen15-planes."""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path
import time

from q_true_initial_eta_open_20261103_v63 import compose,pack,unpack
from q_ordinary_cubic_initial_20261103_v63 import rank_QQ

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quartic',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    ppath=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz'
    planes=json.loads(gzip.decompress(ppath.read_bytes()))
    tailpath=args.quartic/'quartic-tails.json.gz';tails=json.loads(gzip.decompress(tailpath.read_bytes()))
    names=tails['variables'];assert names==planes['variables'];h=names.index('z22')
    true=[row for row in tails['all293_constant_quadratic_cancellation_rows']if not row['cubic_terms']]
    assert len(true)==158
    branches=[]
    for branch in planes['branches']:
        graph={i:{(h,)if not m else tuple(m):F(c)for m,c in data}for i,data in enumerate(branch['all35_graph'])}
        assert all(len(m)==1 for p in graph.values()for m in p)
        restricted=[{'source':row['source'],'terms':pack(compose(unpack(row['quartic_terms']),graph))}
                    for row in true if row['quartic_terms']]
        polys=[unpack(row['terms'])for row in restricted];basis,piv=rank_QQ(polys)
        branches.append({'branch':branch['branch'],'homogeneous_graph':{names[i]:pack(p)for i,p in graph.items()},
                         'all_true_quartic_restrictions':restricted,
                         'nonzero_restricted_quartics':sum(bool(p)for p in polys),
                         'constant_QQ_restricted_span_rank':len(basis),
                         'restricted_QQ_basis':[pack(p)for _,p in sorted(basis.items())],
                         'selected_true_source_rows':[restricted[i]['source']for i in piv]})
    cp=args.output/'restrictions.json.gz'
    with gzip.open(cp,'wt')as handle:json.dump({'variables':names,'branches':branches},handle,separators=(',',':'))
    report={'id':'Q_ORDINARY_QUARTIC_PLANE_RESTRICTIONS_20261103_V63','all_passed':True,'coefficient_field':'QQ',
            'source_rows_have_exact_quadratic_and_cubic_cancellations':True,
            'true_quartic_source_zeroC_count':158,
            'branch_results':[{key:branch[key]for key in ['branch','nonzero_restricted_quartics','constant_QQ_restricted_span_rank','selected_true_source_rows']}for branch in branches],
            'full_quartic_initial_space_or_local_dimension_claimed':False,
            'seconds':time.monotonic()-started,'evidence_sha256':hashlib.sha256(cp.read_bytes()).hexdigest(),
            'input_sha256':{str(p.resolve().relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),ppath,tailpath,
                ROOT/'code/q_ordinary_cubic_initial_20261103_v63.py',ROOT/'code/q_true_initial_eta_open_20261103_v63.py']}}
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
