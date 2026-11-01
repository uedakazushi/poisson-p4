#!/usr/bin/env python3
"""Universal T5 normal length lower bound23 and the specified free rank23 local open."""
import argparse,hashlib,json,time
from pathlib import Path
from s61_model import S61
from s61_r0_t5_jet_engine_20261101 import compute
from s61_r0_t5_initial_engine_20261101 import verify as verify_initial
from s61_r0_t5_slice_20261101 import verify as verify_slice

ROOT=Path(__file__).resolve().parents[1]

def run():
    start=time.time();model=S61()
    frame_path=ROOT/'reports/S61_R0_T5_NORMAL_FRAME_INPUT_20261101.json'
    orbit_path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json'
    old_path=ROOT/'reports/S61_R0_MIXED_SURFACE_T5_OPEN_20261101.json'
    g6_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    frame=json.loads(frame_path.read_text());orbit=json.loads(orbit_path.read_text())['orbits']['T5']
    old=json.loads(old_path.read_text());g6=json.loads(g6_path.read_text())['universal_new_equation']
    for source in [frame,old]:
        for name,digest in source['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    assert old['all_passed'] and old['point']['actual_transverse_length']==23
    assert old['point']['stats'][-1]['transverse_counts']==[1,3,5,6,5,3,0]
    assert old['point']['point_data']['parameter_rows']==[0,1]
    geometry=verify_slice(model,orbit,frame);print('finite slice and actual surface PASS',time.time()-start,flush=True)
    R,F,statistics=compute(model,orbit,frame,g6)
    initial=verify_initial(R,F);assert initial['truncated_length']==23
    paths=[Path(__file__),frame_path,orbit_path,old_path,g6_path,
           ROOT/'code/s61_r0_t5_jet_engine_20261101.py',ROOT/'code/s61_r0_t5_initial_engine_20261101.py',ROOT/'code/s61_r0_t5_slice_20261101.py',
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_R0_MIXED_SURFACE_T5_OPEN_20261101.md',ROOT/'research/notes/S61_R0_MIXED_SURFACE_T5_OPEN_REVIEW_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',model.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_R0_T5_FORMAL_RANK23_20261101',all_passed=True,field='QQ and the rational function field of an exact eight-variable stabilizer slice',
                geometry=geometry,jet_statistics=statistics,normal_initial_certificate=initial,
                specified_point=dict(finite_free_rank=23,base_dimension=2,reduced_support='the same actual smooth T5 surface',homogenizer_power_zero=23),
                relative=dict(base_dimension=14,finite_free_rank=23,homogenizer_power_zero=23,smooth_equivariant_descent_to_ambient=True),
                conclusion='The specified T5 completion is finite free rank23 over C[[u1,u2]] with t^23=0; a containing ambient open is removed by t-saturation.',
                limitations=['No presentation of the rank23 Artin multiplication algebra.','No assertion t=0 in the full local ring.','The entire sufficiently general-beta T5 surface is NOT removed.','No full generic affine-fiber exhaustion.'],
                prior_transverse_length23_reused=True,prior_sixth_jet_upper_bound_not_rerun=True,
                proof_note='research/notes/S61_R0_T5_FORMAL_RANK23_20261101.md',seconds=time.time()-start,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
