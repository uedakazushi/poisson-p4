#!/usr/bin/env python3
"""Actual complete rank1-beta fibers used to close reduced CA/CB intersection.

This proves no scheme-theoretic intersection ideal or full S61 exhaustion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time
import sympy as s
from s61_model import S61,X,pushforward
from additional_rank4_families_20261103_v62 import four
from ambient_weight1_image import common_plane
from ambient_zero_firstjet_fibers_20261103_v63 import divergence,wedge

ROOT=Path(__file__).resolve().parents[1]

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    model=S61();h,x,y,z,w=X;l,r=s.symbols('lambda r');k=s.symbols('k0:9')
    mons=[x*x*z,x*x*w,x*z*z,x*z*w,x*w*w,z**3,z*z*w,z*w*w,w**3];G=sum(a*b for a,b in zip(k,mons))
    report=json.loads((ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json').read_text())
    bi=report['beta_indices_zero_based'];qi=report['q_indices_zero_based'];adapter=s.Matrix(report['beta16_to_intrinsic24_adapter'])
    params,b=common_plane();expected=b.subs(dict(zip(params,[0,1,0,0,0,0])))
    chart=s.Matrix([[1,0,-l,0],[0,1,0,0],[0,0,1,0],[0,l,0,1]])
    swap=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    records={}
    for chart_name,spatial,h2,line,linear in [('y_chart',chart,x*z-l*z*z/2,s.Matrix([0,1,0,l]),x+l*z),
            ('off_point',swap,-z*z/2,s.Matrix([0,0,0,1]),z)]:
        full=s.eye(5);full[1:5,1:5]=spatial
        for family in ['CA','CB']:
            potential=h*h2+G;f=h*h+r*x*x;radial=0 if family=='CA' else -r*x
            raw=four([f,0,0,radial*z-s.diff(potential,w),radial*w+s.diff(potential,z),h*h if family=='CA' else f])
            transformed=pushforward(raw,full);normalized=model.normalize(transformed);u=s.Matrix(normalized['parameters'])
            assert normalized['perspective']==[0]*4
            assert adapter*u[bi,0]==expected
            verification=model.verify_point(normalized,check_equations=True)
            Q=transformed[1:5,1:5].subs(h,0)
            Y=s.Matrix(s.symbols('x y z w'))
            Q4=Q.subs(dict(zip(X[1:],Y)),simultaneous=True);vv=line.subs(dict(zip(X[1:],Y)),simultaneous=True)
            linear4=linear.subs(dict(zip(X[1:],Y)),simultaneous=True)
            assert (divergence(Q4)-(2 if family=='CA' else 4)*r*linear4*vv).applyfunc(s.expand)==s.zeros(4,1)
            assert wedge(vv,Q4.subs(r,0))==[0]*4
            q=u[qi,0];M=q.jacobian([r,*k]);assert (q-M*s.Matrix([r,*k])).applyfunc(s.expand)==s.zeros(45,1)
            assert not any(v.free_symbols & set([r,*k]) for v in M)
            records[chart_name+'_'+family]=dict(normalized_61_coordinates=list(map(str,u)),zero_perspective=True,
                exact_verification=verification,actual_quadratic_matrix=[[str(v) for v in row] for row in Q.tolist()],
                divergence=list(map(str,divergence(Q4))),flag=list(map(str,line)),
                quadratic_q_direction_matrix=[[str(v) for v in row] for row in M.tolist()])
            print(chart_name,family,'134+150+JACOBI+DIV+ZERO SHIFT PASS',time.monotonic()-started,flush=True)
    paths=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/ambient_zero_firstjet_fibers_20261103_v63.py',
        ROOT/'code/ambient_weight1_image.py',ROOT/'code/additional_rank4_families_20261103_v62.py',
        ROOT/'reports/AMBIENT_ZERO_FIRSTJET_FIBERS_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_RANK2_INTERSECTION_UNIQUE_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json']
    body=dict(id='AMBIENT_WHOLE_REDUCED_INTERSECTION_20261103_v63',all_passed=True,field='QQ',
        canonical_rank1_beta=list(map(str,expected)),families=records,
        conclusion='GLOBAL reduced intersection (CA intersection CB)_red=N in S61; hence N is its sole irreducible component. The scheme-theoretic intersection ideal/thickness remains open.',
        proof_dependency='All nonzero rank1 closed primitive beta are Sp4-equivalent to x dx wedge dz. Full proper resolution fibers above it have flags in Pspan(e_y,e_w). Actual q has no quadratic beta shift in this frame. Nonzero divergence recovers the flag and contradicts fixedflag intersection; zero divergence gives a common rank9 family. Rank>=2 and beta0 cases were separately certified.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(body,indent=2)+'\n');print('WHOLE REDUCED INTERSECTION SAVED',body['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
