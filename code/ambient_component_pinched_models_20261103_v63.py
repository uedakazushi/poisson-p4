#!/usr/bin/env python3
"""Exact completed local models of reduced CA/CB/N along rank2 critical beta.

The base slice supplies p²,p³,sp as actual target coordinates. Therefore after
choosing target q fiber coordinates, all higher p terms already belong to the
image; the first odd coefficient computes the WHOLE missing module, not only
a truncated obstruction.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61,X,pushforward
from additional_rank4_families_20261103_v62 import four
from ambient_base_pinched_surface_20261103_v63 import finite_beta
from ambient_weight1_image import common_plane

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    model=S61();h,x,y,z,w=X;p,sv,rv=s.symbols('p s r');k=s.symbols('k0:9')
    mon=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3]
    G=sum(a*b for a,b in zip(k,mon));h2=x*w+z*z-p*z*w/2+sv*w*w;potential=h*h2+G
    group=s.eye(5);group[1:5,1:5]=s.Matrix([[1,0,-p,0],[0,1,0,0],[0,0,1,0],[0,p,0,1]])
    weights=[a['weight'] for a in model.certificate['free_variable_basis']]
    bi=[i for i,a in enumerate(weights) if a==1];qi=[i for i,a in enumerate(weights) if a==2]
    adapter=s.Matrix(json.loads((ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json').read_text())['beta16_to_intrinsic24_adapter'])
    bp,b=common_plane();t,A,B,C,D,F=bp
    expected_beta=finite_beta(b.subs({t:p,A:0,B:1,C:1,D:-p/2,F:sv}),0,0,p)
    records={}
    for name in ['N','CA','CB']:
        source_r=s.S.Zero if name=='N' else rv;pars=list(k) if name=='N' else [rv,*k]
        f=h*h+(s.Rational(5,2) if name!='CB' else -5)*p*h*x+source_r*x*x
        radial=s.S.Zero if name!='CB' else -s.diff(f,x)/2
        raw=four([f,0,0,radial*z-s.diff(potential,w),radial*w+s.diff(potential,z),h*h if name!='CB' else f])
        normalized=model.normalize(pushforward(raw,group));u=s.Matrix(normalized['parameters'])
        assert adapter*u[bi,0]==expected_beta
        verification=model.verify_point(normalized,check_equations=True)
        q=u[qi,0];M=q.jacobian(pars)
        assert not any(entry.free_symbols & set(pars) for entry in M)
        assert q==M*s.Matrix(pars)+q.subs(dict.fromkeys(pars,0))
        q0=q.subs({p:0,sv:0});q1=q.diff(p).subs({p:0,sv:0});M0=q0.jacobian(pars)
        selected=list(M0.T.rref()[1]);L=M0[selected,:]
        assert len(selected)==len(pars) and L.det()!=0
        # q_selected is the actual target K. At p=s=0 it is L*source+b.
        # Differentiating at fixed K corrects the source fiber by -L^-1*q1sel.
        first_odd=(q1-M0*L.inv()*q1[selected,0]).applyfunc(s.expand)
        odd_matrix=first_odd.jacobian(pars);constant=first_odd.subs(dict.fromkeys(pars,0))
        assert constant==s.zeros(45,1)
        desired=[k[i] for i in [1,3,4,6,7,8]]
        if name!='N':desired=[rv,*desired]
        expected_matrix,_=s.linear_eq_to_matrix(desired,pars)
        assert odd_matrix.rank()==len(desired)
        assert odd_matrix.col_join(expected_matrix).rank()==len(desired)
        assert s.linsolve(list(first_odd),pars)==s.linsolve(desired,pars)
        codim=len(desired)+1;dimension=18 if name=='N' else 19
        free_dimension=dimension-1-codim
        assert free_dimension==10
        records[name]=dict(source_fiber_parameters=list(map(str,pars)),full_61_coordinates=list(map(str,u)),
            exact_verification=verification,actual_q_selected_rows=selected,
            selected_q_coordinate_functions=list(map(str,q[selected,0])),fiber_linear_block_at_origin=L.tolist(),
            fiber_linear_block_determinant=str(L.det()),target_K_constant_shift=list(map(str,q0[selected,0].subs(dict.fromkeys(pars,0)))),
            first_odd_coefficient_at_fixed_target_K=list(map(str,first_odd)),
            first_odd_coefficient_matrix=odd_matrix.tolist(),first_odd_matrix_rank=odd_matrix.rank(),
            missing_module_ideal_in_target_residue_coordinates=['s',*map(str,desired)],
            pinching_codimension=codim,local_dimension=dimension,regular_factor_dimension=free_dimension,
            embedding_dimension=dimension+1+codim,depth=dimension-codim,multiplicity=2,
            normalization_quotient_dimension=free_dimension,
            normalization_ring=f'C[[u1,...,u10,s1,...,s{codim},p]]',
            image_ring=f'C[[u1,...,u10,s1,...,s{codim},p^2,p^3,s1*p,...,s{codim}*p]]',
            conductor_in_normalization=f'(s1,...,s{codim},p^2)')
        print(name,'134+150+PROJECTIVE JACOBI PASS','ODD RANK',odd_matrix.rank(),'c',codim,'seconds',time.monotonic()-started,flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/additional_rank4_families_20261103_v62.py',
        ROOT/'code/ambient_base_pinched_surface_20261103_v63.py',ROOT/'code/ambient_weight1_image.py',
        ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_BASE_PINCHED_SURFACE_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_RANKDROP_FLAG_FIBERS_V2_20261103_v63/report.json']
    result=dict(id='AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63',all_passed=True,field='QQ',
        base_slice='t=p,A=0,B=1,C=1,D=-p/2,F=s; flag=(0:1:0:p)',
        cubic_monomials=list(map(str,mon)),beta24_coordinates=list(map(str,expected_beta)),families=records,
        geometric_scope='Full completed local rings of reduced CA/CB/N at their NONNORMAL points above the7-dimensional critical rank2 beta orbit, including distinguished v24 Q(a=1) by equivariance. Not full S61 or whole intersection ideal.',
        exactness_reason='Actual beta outputs contain s,p²,p³,sp. Chosen actual q outputs give all remaining regular source coordinates. Thus p² times the whole normalization is already in the target image; the computed first odd coefficient ideal determines the entire missing normalization module.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    def serial(v):
        if isinstance(v,s.Rational):return str(v)
        raise TypeError(type(v))
    (output/'report.json').write_text(json.dumps(result,indent=2,default=serial)+'\n');print('COMPONENT PINCHED MODELS SAVED',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
