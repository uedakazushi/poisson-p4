#!/usr/bin/env python3
"""Freeze the bounded fifth computation and exactQQ completion receipts."""
import argparse
import json
from pathlib import Path
import shutil
from q_positive_weight6_quadratic_model_20261103_v63 import ROOT,read,sha


def main(source,output):
    construction=read(source/'report.json');assert construction['all_passed']
    for name,digest in construction['input_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in construction['evidence_sha256'].items():assert sha(source/name)==digest,name
    log=(source/'complete-fifth-cone.log').read_text()
    assert all(marker in log for marker in ['ALL377_EXACT_QQ_REDUCTIONS_ZERO','HF5=171','A_DIM=2','B_DIM=3','COMPLETE_FIFTH_AND_RADICAL_SPLIT_DONE'])
    assert 'A_MINIMAL_PRIMES={ideal(u6,u4,u3), ideal(u6,u4,u2)}'in log
    assert 'B_MINIMAL_PRIMES={ideal(v6,v2), ideal(v6,v5,v3)}'in log
    selected=(source/'selected-fifth.log').read_text()
    assert 'DIM5=3'in selected and 'DEG5=5'in selected
    output.mkdir(parents=True,exist_ok=False)
    for path in source.iterdir():
        if path.is_file():shutil.copyfile(path,output/('construction.json'if path.name=='report.json'else path.name))
    paths=[Path(__file__),ROOT/'code/q_positive_weight6_initial5_20261103_v63.py',
           ROOT/'code/q_positive_weight6_initial5_rank_20261103_v63.cpp',ROOT/'code/q_positive_weight6_initial5_cone_20261103_v63.py',
           ROOT/'research/notes/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63.md',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_INITIAL4_20261103_v63/report.json',
           ROOT/'reports/Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL_20261103_v63/report.json']
    report=dict(id='Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63',all_passed=True,
        bounded_construction_completed_seconds=construction['seconds'],memory_cap_GiB=2.2,
        all350_fifth_QQ_formal_residuals_and95_implicit_coefficients_saved=True,
        complete_initial_degree5_spanning_rows=1329,lower_multiple_rows=952,
        all329_order5_generators_and48_fourth_syzygy_carries_accounted_for=True,
        all377_new_QQ_rows_reduce_to_selected15_plus_prior21=True,
        exact_QQ_degree5_initial_rank=621,actual_transverse_Hilbert_function_through5=[1,8,31,72,126,171],
        through5_ideal_dimension=3,through5_ideal_degree=5,
        radical_minimal_primes_in_normal_z_coordinates=[[1,2,4,6,7],[1,3,4,5,6,7]],
        radical_unique3plane_free_coordinates=['z0','z3','z5'],radical2plane_free_coordinates=['z0','z2'],
        actual_transverse_local_dimension_upper_bound=3,actual_projective_Poisson_local_dimension_upper_bound=25,
        actual3dimensional_formal_family_or_new_component_proved=False,
        projective_dimension_at_most24_or_full_initial_generation_claimed=False,
        no_sixth_jet_or_additional_large_GB_launched=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
        evidence_sha256={p.name:sha(p)for p in output.iterdir()if p.is_file()})
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print('FROZEN complete exact initial5; HF5=171; remaining3plane+2plane, no actuallift claim')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--input',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();main(args.input,args.output)
