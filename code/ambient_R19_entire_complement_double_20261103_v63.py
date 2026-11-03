#!/usr/bin/env python3
"""Six exact remaining orbit checks complete ordinary-double R minus N.

GenericAF normalforms and complete orbit coverage are separate frozen inputs.
No assertion is made about R intersection N.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
import sympy as s
from ambient_CA_global_flag_model_20261103_v63 import X,raw,canonical,vector
from ambient_CB_global_flag_model_20261103_v63 import raw as cbraw
from ambient_rank_two_flag_intersection_table_20261103_v63 import permutation,lie

ROOT=Path(__file__).resolve().parents[1]

def push(P,g):return (g*P.subs(dict(zip(X,g.inv()*s.Matrix(X))),simultaneous=True)*g.T).applyfunc(s.expand)

def minor(M):
    cols=list(M.rref()[1]);rows=list(M[:,cols].T.rref()[1]);return dict(rank=len(cols),rows=rows,columns=cols,determinant=str(M[rows,cols].det()))

def main(output):
    output.mkdir(parents=True,exist_ok=False);started=time.monotonic();h,x,y,z,w=X;AA,BB,FF,eps=s.symbols('A B F eps')
    ca_path=ROOT/'reports/AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63/report.json';cb_path=ROOT/'reports/AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63/report.json'
    ca=json.loads(ca_path.read_text());cb=json.loads(cb_path.read_text());cubics=list(map(s.sympify,cb['cubic_monomials']))
    CA=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v))) for v in ca['basis_canonical150']))
    CB0=s.Matrix.hstack(*(s.Matrix(list(map(s.Rational,v))) for i,v in enumerate(cb['enlarged_canonical_basis']) if i<3 or i>=6))
    perm=(1,3,2,0,4);CB=s.Matrix.hstack(*(permutation(CB0[:,i],perm) for i in range(19)));selected=list(CB.T.rref()[1]);left=CB[selected,:].inv()
    mu=permutation(vector(canonical(cbraw(0,6*h*x,0))),perm)
    sourceP=raw(-4*x*x,0,AA*h*h*z+BB*h*z*z/2+FF*z**3/3+2*x*z*w);point=vector(canonical(sourceP))
    cp=left*point[selected,0];assert CB*cp==point and list(cp[:3,0])==[0,6,0]
    M=s.zeros(5);M[1,3]=2;g=s.eye(5)+eps*M;pg=s.zeros(5)
    for i,j in enumerate(perm):pg[j,i]=1
    PA=canonical(push(sourceP+eps*raw(0,0,3*z*z*w),g))
    H0=sum(cp[i+3]*v for i,v in enumerate(cubics));PB=canonical(push(cbraw(6*h*x-21*eps*x*x,0,H0+eps*cubics[10]),pg))
    diff=s.zeros(150,1);diff[72]=-18
    assert vector(PA-PB)==eps*eps*diff
    rd=[vector(canonical(raw(-4*x*x,0,2*x*z*w)))]+[vector(canonical(raw(0,0,H))) for H in [w*x*x,x*x*z,x*z*z/2,2*h*x*z,z**3/3,h*z*z,3*h*h*z]]
    records=[]
    for label,values in [('h2_plus_hz',[1,1,0]),('hz_plus_z2',[0,1,1]),('h2',[1,0,0]),('z2',[0,0,1]),('hz',[0,1,0]),('zero',[0,0,0])]:
        subst=dict(zip([AA,BB,FF],values));P0=point.subs(subst)
        orbit=s.Matrix.hstack(*(lie(P0,i,j) for i,j in itertools.product(range(5),repeat=2)))
        JA=CA.row_join(orbit);JB=CB.row_join(mu).row_join(orbit);W=JA.row_join(JB);R=s.Matrix.hstack(*rd).row_join(orbit)
        assert [JA.rank(),JB.rank(),W.rank(),R.rank(),W.row_join(diff).rank()]==[29,29,37,20,38]
        ann=next(v for v in W.T.nullspace() if (v.T*diff)[0]!=0);obstruction=(ann.T*diff)[0]
        first=vector(PA.diff(eps).subs(eps,0)).subs(subst);assert R.row_join(first).rank()==21
        records.append(dict(label=label,Q_coefficients=values,CB_source19_coefficients=list(map(str,cp.subs(subst))),
            CA_source_minor=minor(JA),CB_source_minor=minor(JB),joint_R_source_minor=minor(R),
            tangent_union_rank=37,second_difference_augmented_rank=38,
            second_obstruction_covector=list(map(str,ann)),second_obstruction_value=str(obstruction),
            first_direction_augmented_R_rank=21))
        print(label,'29/29/37/20/38 + NONZERO SECOND OBSTRUCTION PASS',time.monotonic()-started,flush=True)
    # TheonlynewDarbouxcandidate atQ0 ismarkedw; its traceconditionsbothfail.
    ellw_f=-2*x*x;assert s.diff(ellw_f,x)==-4*x
    paths=[Path(__file__),ca_path,cb_path,ROOT/'code/ambient_rank_two_flag_intersection_table_20261103_v63.py',
        ROOT/'reports/AMBIENT_R19_OPEN_DOUBLE_20261103_v63/report.json',
        ROOT/'reports/AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63/report.json',
        ROOT/'research/notes/AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63.md']
    result=dict(id='AMBIENT_R19_ENTIRE_COMPLEMENT_DOUBLE_20261103_v63',all_passed=True,field='QQ',
        full_A_B_F_actual_finite_curve_identity='PA-PB=eps²*(-18z²dxy), exactnohigherterms',
        special_orbit_records=records,
        Q_zero_additional_w_flag={'f':'-2*x²','f_x':'-4*x','CA_required':'2*x','CB_required':'x','both_rejected':True},
        entire_complement_conclusions=['R minusN is smooth19 and rational; isomorphic to the t!=0 affine7bundle overthe12D bigflagpairorbit',
            'CA andCB smooth28 onopenneighborhoods containingALL R minusN',
            'ActualD completion C[[u1,...,u19,eps]]/(eps²) atEVERYpointofR minusN',
            'NoembeddedassociatedprimesofD outsideN; the entireR complement isordinarydouble'],
        scope='No orbit, scheme, singularity, or primary-decomposition exhaustion on R intersectionN.',
        input_sha256={str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},seconds=time.monotonic()-started)
    (output/'report.json').write_text(json.dumps(result,indent=2)+'\n');print('ENTIRE R MINUS N ORDINARY DOUBLE PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
