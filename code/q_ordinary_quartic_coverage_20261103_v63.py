#!/usr/bin/env python3
"""Assemble the complete degree4 cancellation coverage, no further jet/GB run."""
import argparse,gzip,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import time
from q_true_initial_eta_open_20261103_v63 import compose,unpack
ROOT=Path(__file__).resolve().parents[1]

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--packet',type=Path,required=True);args=ap.parse_args();p=args.packet.resolve();t=time.monotonic()
    files=[p/'report.json',p/'plane-restrictions/report.json',p/'constant-kernel/report.json',p/'linear-syzygies/report.json',p/'mixed-syzygies/report.json']
    reports=[json.loads(f.read_text())for f in files]
    for r in reports:
        assert r['all_passed']
        for name,h in r['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    assert reports[0]['identically_zero_cubic_rows']==158 and reports[4]['exact_full_nonzero_cubic_space_dependencies']==219
    assert reports[4]['all219_actual350_lower2_and_lower3_polynomial_cancellations_replayed']
    assert all(r['constant_QQ_restricted_span_rank']==0 for r in reports[1]['branch_results'])
    assert all(r['QQ_rank']==0 for r in reports[4]['branches'])
    tailpath=p/'quartic-tails.json.gz';tails=json.loads(gzip.decompress(tailpath.read_bytes()))
    qpath=ROOT/'reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63/quadrics.json.gz';q=json.loads(gzip.decompress(qpath.read_bytes()))
    plane_path=ROOT/'reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz';planes=json.loads(gzip.decompress(plane_path.read_bytes()))
    h=q['variables'].index('z22')
    for b in planes['branches']:
        graph={i:{(h,)if not m else tuple(m):F(c)for m,c in data}for i,data in enumerate(b['all35_graph'])}
        assert all(not compose({(i,j):F(c)for i,j,c in r},graph)for r in q['all350_quadratic_rows'])
        assert all(not compose(unpack(r['cubic_terms']),graph)for r in tails['all293_constant_quadratic_cancellation_rows'])
    files +=[Path(__file__),tailpath,qpath,plane_path]
    report={'id':'Q_ORDINARY_QUARTIC_COVERAGE_20261103_V63','all_passed':True,'coefficient_field':'QQ',
            'degree3_source_rows_including_zero_cubics':1995+293,'degree3_image_rank_QQ':1911,
            'degree3_kernel_dimension':219+158,'nonzero_cubic_row_kernel':219,'zero_cubic_row_kernel':158,
            'all377_quartic_carries_vanish_identically_on_both15planes':True,
            'fresh_all350_Q_and_all293_cubic_tails_vanish_on_both15planes':True,
            'entire_true_ordinary_initial_components_degrees_le4_vanish_on_both15planes':True,
            'completeness_formula':'I*_4 = span{S_2 Q57, S_1 C293, quartic carry of ker(S_1 Q57 + C293 -> S_3)}',
            'through_degree4_homogeneous_ideal_dimension_lower_bound':15,
            'any12_linear_section_through_degree4_dimension_lower_bound':3,
            'actual_full_initial_or_completed_dimension_lower15_claimed':False,
            'actual15dim_branch_claimed':False,'no_degree5_jet_or_large_GB_launched':True,
            'input_sha256':{str(f.resolve().relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()for f in files},'seconds':time.monotonic()-t}
    out=p/'coverage.json';assert not out.exists();out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
