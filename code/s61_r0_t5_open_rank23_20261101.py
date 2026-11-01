#!/usr/bin/env python3
"""Every T5 point in an explicit invariantly normalized open is free rank23."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_r0_t5_jet_engine_20261101 import compute
from s61_r0_t5_open_sextic_engine_20261101 import verify as verify_upper
from s61_r0_t5_slice_20261101 import verify as verify_slice
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action

ROOT=Path(__file__).resolve().parents[1]


def proper_images(model,orbit,geometry):
    b=s.symbols('b0:12');H=s.Matrix(orbit['beta_kernel']).applyfunc(s.Rational)
    N=s.Matrix(geometry['beta_normal_action']).applyfunc(s.sympify)
    assert N[:,[0,1,2,3,8,9]]==s.zeros(4,6)
    assert s.factor(N[:,[4,5,6,7]].det())==s.sympify(geometry['beta_normal_minor'])
    z,u,beta=primitive_beta(model)
    v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis=v.jacobian(u);rows=list(basis.T.rref()[1]);inverse=basis[rows,:].inv()
    act=s.Matrix.hstack(*[inverse*infinitesimal_action(v,A)[rows,:] for A in sp_basis()])
    ar=act.subs(dict(zip(u,H*s.Matrix(b))));hr=list(H.T.rref()[1]);free=[0,1,2,3,8,9]
    DA=(H[hr,:].inv()*ar[hr,free]).applyfunc(s.expand)
    assert H*DA==ar[:,free]
    char=25*b[5]+7*b[10]
    dchar=(s.Matrix([char]).jacobian(b)*DA).applyfunc(s.expand)
    restricted=dchar.subs(b[10],-25*b[5]/7).applyfunc(s.expand)
    assert restricted==s.zeros(1,6)
    return dict(beta_normal_rank_at_most3_when_b9P_zero=True,
                beta_image_dimension_b9P_zero_at_most=15,
                additional_proper_hypersurface=str(char),
                normal_free_columns=free,character_derivatives=list(map(str,dchar)),
                character_restricted_derivatives_zero=True,
                beta_image_dimension_character_zero_at_most=15,
                note='The character exclusion is preparation for a future alternative chart; it does not by itself handle L=0 or Aq=0.')


def run():
    start=time.time();model=S61()
    frame_path=ROOT/'reports/S61_R0_T5_NORMAL_FRAME_INPUT_20261101.json'
    orbit_path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json'
    lower_path=ROOT/'reports/S61_R0_T5_FORMAL_RANK23_20261101.json'
    g6_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    frame=json.loads(frame_path.read_text());orbit=json.loads(orbit_path.read_text())['orbits']['T5']
    lower=json.loads(lower_path.read_text());g6=json.loads(g6_path.read_text())['universal_new_equation']
    assert lower['all_passed']
    for source in [frame,lower]:
        for name,digest in source['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    geometry=verify_slice(model,orbit,frame)
    proper=proper_images(model,orbit,geometry)
    print('finite slice, same surface and proper images PASS',time.time()-start,flush=True)
    R,F,statistics=compute(model,orbit,frame,g6)
    initial=verify_upper(R,F)
    print('universal sextic upper23 PASS',time.time()-start,flush=True)
    paths=[Path(__file__),frame_path,orbit_path,lower_path,g6_path,
           ROOT/'code/s61_r0_t5_jet_engine_20261101.py',ROOT/'code/s61_r0_t5_open_sextic_engine_20261101.py',ROOT/'code/s61_r0_t5_slice_20261101.py',
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_R0_T5_FORMAL_RANK23_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',model.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_R0_T5_OPEN_RANK23_20261101',all_passed=True,
                field='QQ and the exact eight-variable stabilizer slice',geometry=geometry,proper_images=proper,
                jet_statistics=statistics,normal_initial_certificate=initial,
                conclusion='At EVERY point with b9*L*P != 0 and Aq != 0 after the stated finite normalization, the fixed-beta completion is free rank23 over its same smooth surface, and t^23=0.',
                local=dict(base_dimension=2,finite_free_rank=23,homogenizer_power_zero=23),
                relative=dict(base_dimension=14,finite_free_rank=23,homogenizer_power_zero=23,smooth_equivariant_ambient_descent=True),
                limitations=['The frame artifact L=0 and two normalized Aq roots remain separate.','No assertion t=0.','No entire T5 orbit removal or affine exhaustion.'],
                proof_note='research/notes/S61_R0_T5_OPEN_RANK23_20261101.md',seconds=time.time()-start,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
