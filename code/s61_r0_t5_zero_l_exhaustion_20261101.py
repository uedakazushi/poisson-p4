#!/usr/bin/env python3
"""Finish the L=0 chart and remove the entire general-beta T5 infinity orbit."""
import argparse,hashlib,json,time
from pathlib import Path
import sympy as s
from s61_model import S61
from s61_r0_t5_zero_l_jet_engine_20261101 import compute
from s61_r0_t5_zero_l_sextic_engine_20261101 import verify as verify_upper
from s61_r0_t5_slice_20261101 import verify as verify_slice
from primitive_jet_elliptic_curve_20261101 import primitive_beta
from weight1_structure import COORDS,sp_basis,infinitesimal_action
ROOT=Path(__file__).resolve().parents[1]


def geometry_audit(model,orbit,old_geometry,frame):
    b=s.symbols('b0:12');L=68*b[0]-23*b[10]-121*b[5];N=25*b[5]+7*b[10]
    P=1425*b[0]*b[10]*b[5]+1375*b[0]*b[5]**2-1050*b[0]*b[6]*b[9]+105*b[10]*b[2]*b[9]-180*b[11]*b[5]*b[9]+275*b[2]*b[5]*b[9]+210*b[7]*b[9]**2+84*b[8]*b[9]**2
    substitutions={s.Symbol(k):s.sympify(v) for k,v in frame['slice'].items()}
    assert s.cancel(L.subs(substitutions,simultaneous=True))==0
    assert s.cancel(N.subs(substitutions,simultaneous=True))==1
    assert s.cancel(P.subs(substitutions,simultaneous=True))==s.Symbol('p')
    assert substitutions[b[9]]==substitutions[b[5]]==1 and substitutions[b[2]]==substitutions[b[7]]==0
    for label,record in old_geometry['finite_actions'].items():
        bn=list(map(s.sympify,record['beta_action']));new=dict(zip(b,bn));mult=s.Symbol('lam') if label=='torus' else 1
        assert s.expand(L.subs(new,simultaneous=True)-mult*L)==0
        assert s.expand(N.subs(new,simultaneous=True)-mult*N)==0
        assert s.expand(bn[9]-(mult**3)*b[9])==0
    H=s.Matrix(orbit['beta_kernel']).applyfunc(s.Rational)
    z,u,beta=primitive_beta(model);v=s.Matrix([beta[i,j].diff(z[k]) for (i,j),k in COORDS])
    basis=v.jacobian(u);rows=list(basis.T.rref()[1]);inv=basis[rows,:].inv()
    action=s.Matrix.hstack(*[inv*infinitesimal_action(v,A)[rows,:] for A in sp_basis()])
    ar=action.subs(dict(zip(u,H*s.Matrix(b))));hr=list(H.T.rref()[1]);free=[0,1,2,3,8,9]
    DA=(H[hr,:].inv()*ar[hr,free]).applyfunc(s.expand)
    assert H*DA==ar[:,free]
    derivative=-(s.Matrix([L]).jacobian(b)*DA).applyfunc(s.expand)
    assert derivative[2]==136*b[9]/5
    restricted=derivative.subs(b[10],(68*b[0]-121*b[5])/23).applyfunc(s.expand)
    surface=s.Matrix(old_geometry['actual_surface_tangent']).applyfunc(s.sympify)
    assert surface[:2,:].rank()==2 and surface[2:,:]==s.zeros(3,6)
    assert frame['free_columns'][:2]==[19,24]
    return dict(slice_conditions_verified=True,N_and_L_relative_invariance_verified=True,
                weighted_and_torus_normalization='lambda^2=N/b9, c=1/(lambda*N), followed by the three verified unipotents',
                fixed_beta_compensating_b_derivative=[list(map(str,row)) for row in (-DA).tolist()],
                fixed_beta_L_derivative=list(map(str,derivative)),L_zero_restriction=list(map(str,restricted)),
                nonzero_crossing_direction='136*b9/5',same_actual_surface_crosses_L_zero_at_every_b9P_nonzero_point=True,
                actual_surface_parameter_rows=[0,1],source_smooth_on='b9*P != 0')


def run():
    start=time.time();model=S61()
    frame_path=ROOT/'reports/S61_R0_T5_ZERO_L_FRAME_INPUT_20261101.json'
    old_frame_path=ROOT/'reports/S61_R0_T5_NORMAL_FRAME_INPUT_20261101.json'
    orbit_path=ROOT/'reports/S61_R0_MIXED_ORBIT_REDUCTION_20261101.json'
    open_path=ROOT/'reports/S61_R0_T5_OPEN_RANK23_20261101.json'
    aq_path=ROOT/'reports/S61_R0_T5_QUADRATIC_ARTIFACTS_20261101.json'
    g6_path=ROOT/'reports/S61_UNIVERSAL_WEIGHT6_20261031.json'
    frame=json.loads(frame_path.read_text());old_frame=json.loads(old_frame_path.read_text())
    orbit=json.loads(orbit_path.read_text())['orbits']['T5'];g6=json.loads(g6_path.read_text())['universal_new_equation']
    for path in [frame_path,open_path,aq_path]:
        report=json.loads(path.read_text())
        if path!=frame_path:assert report['all_passed']
        for name,digest in report['input_sha256'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    old_geometry=verify_slice(model,orbit,old_frame)
    geometry=geometry_audit(model,orbit,old_geometry,frame)
    print('same smooth surface crosses L=0: PASS',time.time()-start,flush=True)
    R,F,statistics=compute(model,orbit,frame,g6)
    # compute() verifies a polynomial inverse over the whole slice. Thus the
    # block determinant is a nonzero constant; one exact evaluation fixes it.
    at={s.Symbol(v):s.S.Zero for v in frame['variables']};at[s.Symbol('p')]=1
    inverse=s.Matrix(frame['inverse']).applyfunc(lambda v:s.sympify(v).subs(at))
    assert inverse.det()*s.Rational(frame['minor'])==1
    upper=verify_upper(R,F)
    print('L=0 all-point upper23: PASS',time.time()-start,flush=True)
    paths=[Path(__file__),frame_path,old_frame_path,orbit_path,open_path,aq_path,g6_path,
           ROOT/'code/s61_r0_t5_zero_l_jet_engine_20261101.py',ROOT/'code/s61_r0_t5_zero_l_sextic_engine_20261101.py',ROOT/'code/s61_r0_t5_slice_20261101.py',
           ROOT/'code/s61_model.py',ROOT/'code/primitive_jet_elliptic_curve_20261101.py',ROOT/'code/weight1_structure.py',
           ROOT/'research/notes/S61_R0_T5_FORMAL_RANK23_20261101.md',ROOT/'research/notes/S61_R0_T5_OPEN_RANK23_20261101.md',ROOT/'research/notes/S61_R0_T5_QUADRATIC_ARTIFACTS_20261101.md',
           ROOT/'research/notes/S61_HOMOGENIZED_EQUIVARIANCE_20261101.md',model.source/'symplectic_slice_reduced_equations.json']
    return dict(id='S61_R0_T5_ZERO_L_EXHAUSTION_20261101',all_passed=True,
                constant_linear_minor=frame['minor'],linear_rank=40,geometry=geometry,jet_statistics=statistics,initial_certificate=upper,
                local=dict(base_dimension=2,finite_free_rank=23,homogenizer_power_zero=23,same_actual_surface=True),
                relative_ambient_removal=True,entire_sufficiently_general_beta_T5_orbit_removed=True,
                coverage=['b9*P=0 has beta image dimension at most15.','N=25b5+7b10=0 has beta image dimension at most15.','L!=0 is covered by OPEN_RANK23 and both QUADRATIC_ARTIFACTS.','L=0,N*b9*P!=0 is covered at every point by the present constant-frame calculation.'],
                conclusion='The entire T5 projective infinity orbit is removed by t-saturation over sufficiently general beta. Every remaining orbit point has free rank23 over its same smooth surface and t^23=0.',
                limitations=['No assertion t=0 in the full ring.','No conclusion about remaining flag/binary infinity strata or full affine five-point exhaustion.'],
                proof_note='research/notes/S61_R0_T5_ZERO_L_EXHAUSTION_20261101.md',seconds=time.time()-start,
                input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();result=run()
    with args.output.open('x') as stream:json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps(dict(all_passed=True,seconds=result['seconds'])))
