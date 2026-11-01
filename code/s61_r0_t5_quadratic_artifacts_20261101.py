#!/usr/bin/env python3
"""Remove both T5 quadratic pivot artifacts by a number-field normal change."""
import argparse,hashlib,json,time
from pathlib import Path
from s61_model import S61
from s61_r0_t5_jet_engine_20261101 import compute
from s61_r0_t5_quadratic_artifact_engine_20261101 import verify
from s61_r0_t5_slice_20261101 import verify as verify_slice
ROOT=Path(__file__).resolve().parents[1]


def run():
    start=time.time();model=S61()
    frame_path=ROOT/'reports/S61_R0_T5_NORMAL_FRAME_INPUT_20261101.json'
    orbit_path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json'
    lower_path=ROOT/'reports/S61_R0_T5_FORMAL_RANK23_20261101.json'
    open_path=ROOT/'reports/S61_R0_T5_OPEN_RANK23_20261101.json'
    g6_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    frame=json.loads(frame_path.read_text());orbit=json.loads(orbit_path.read_text())['orbits']['T5']
    g6=json.loads(g6_path.read_text())['universal_new_equation']
    for path in [lower_path,open_path]:
        report=json.loads(path.read_text());assert report['all_passed']
        for name,digest in report['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    weights=model.reduced['parameter_weights'];equations=model.reduced['equations']+[g6]
    for i,eq in enumerate(equations):
        expected=4 if i<70 else 5 if i<134 else 6
        assert {sum(weights[j] for j in mon) for mon,n,d in eq}=={expected}
        # With qhat fixed, beta -> c beta, t -> c^-2 t scales this homogenized equation by c^(expected-4).
        assert {sum(weights[j]==1 for j in mon)-2*(2-sum(weights[j]==2 for j in mon)) for mon,n,d in eq}=={expected-4}
    geometry=verify_slice(model,orbit,frame)
    print('weighted scaling and finite slice PASS',time.time()-start,flush=True)
    R,F,statistics=compute(model,orbit,frame,g6)
    certificate=verify(R,F)
    print('both quadratic artifacts upper23 PASS',time.time()-start,flush=True)
    paths=[Path(__file__),frame_path,orbit_path,lower_path,open_path,g6_path,
           ROOT/'code/s61_r0_t5_quadratic_artifact_engine_20261101.py',ROOT/'code/s61_r0_t5_jet_engine_20261101.py',ROOT/'code/s61_r0_t5_slice_20261101.py',
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_R0_T5_FORMAL_RANK23_20261101.md',ROOT/'research/notes/S61_R0_T5_OPEN_RANK23_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',model.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_R0_T5_QUADRATIC_ARTIFACTS_20261101',all_passed=True,weighted_scaling_all135_verified=True,
                geometry=geometry,jet_statistics=statistics,quadratic_artifact_certificate=certificate,
                conclusion='Both Aq-root strata with b9*L*P != 0 have transverse length23, free rank23 completion over the same smooth surface and t^23=0. Combined with OPEN_RANK23 this covers EVERY b9*L*P != 0 point.',
                local=dict(base_dimension=2,finite_free_rank=23,homogenizer_power_zero=23),relative_ambient_removal=True,
                limitations=['The L=0 chart remains separate.','No full t=0 assertion.','No entire T5 removal or affine exhaustion.'],
                proof_note='research/notes/S61_R0_T5_QUADRATIC_ARTIFACTS_20261101.md',seconds=time.time()-start,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
