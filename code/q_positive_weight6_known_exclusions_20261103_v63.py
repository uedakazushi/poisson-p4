#!/usr/bin/env python3
"""Exact +8/9 weight-six ray exclusions from actual T01, CA, CB (and N).

No new full component or exhaustive branch claim is made. The T01 relation
is globally valid on its reduced component, including the discriminant-zero
boundary. CA/CB exclusion uses their whole proper flag-source models.
"""
import argparse
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import time
import numpy as np
import sympy as s
from s61_model import S61,X,pushforward,bivector_from_coefficients
from q51_universal_section_20261101 import affine
from toric_closed_jet_coordinates_20261102 import coords,pairs,m2
import toric_zero_fiber_coframe_quadrics_20261102 as lifts

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
def same(A,B):return (A-B).applyfunc(s.expand)==s.zeros(A.rows,A.cols)

def build(out):
    start=time.monotonic();eta=s.Symbol('eta');axis_path=ROOT/'reports/Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63/certificate.json.gz'
    axis=json.loads(gzip.decompress(axis_path.read_bytes()));raw=list(map(s.sympify,axis['plus_raw150']))
    expected={64:-s.Rational(256,3)*eta,83:-64*eta,90:s.Rational(8,9),105:-s.Rational(2,9),116:64*eta,121:s.Rational(4,3),137:s.Rational(2,3),140:s.Integer(1)}
    assert {i:v for i,v in enumerate(raw)if v}==expected
    P=bivector_from_coefficients(raw);div=[s.expand(sum(s.diff(P[i,j],X[i])for i in range(5)))for j in range(5)];assert not any(div)
    all350=[]
    for quad in itertools.combinations(range(5),4):
        value=0
        for pos,e in enumerate(quad):
            i,j,k=quad[:pos]+quad[pos+1:]
            jac=sum(P[i,b]*s.diff(P[j,k],X[b])+P[j,b]*s.diff(P[k,i],X[b])+P[k,b]*s.diff(P[i,j],X[b])for b in range(5))
            value+=(-1)**pos*X[e]*jac
        poly=s.Poly(s.expand(value),*X)
        all350 += [poly.coeff_monomial(s.prod(X[i]for i in mon))for mon in itertools.combinations_with_replacement(range(5),4)]
    assert len(all350)==350 and not any(all350)
    quintic=0
    for e in range(5):
        i,j,k,l=[i for i in range(5)if i!=e]
        quintic+=(-1)**e*X[e]*(P[i,j]*P[k,l]-P[i,k]*P[j,l]+P[i,l]*P[j,k])
    assert s.expand(quintic)==-s.Rational(16,81)*X[0]**5
    # This finite physical frame is nonsingular since a=8/9 is fixed nonzero.
    aa=s.Rational(8,9);frame=s.eye(5);frame[1:,1:]=s.Matrix([[1,0,0,0],[0,0,0,1/aa],[0,1,0,0],[0,0,-s.Rational(9,2),0]])
    model=S61();normalized=model.normalize(pushforward(P,frame));u=s.Matrix(normalized['parameters'])
    model.verify_point(normalized)
    assert {i:v for i,v in enumerate(u)if v}=={3:-s.Rational(2048,27)*eta,15:256*eta,24:-s.Rational(3,2),35:s.Rational(27,8),37:s.Rational(81,16),56:-256*eta}

    # Correct primitive closed-jet coordinates; Pi is NOT the producer's 4Pi.
    z=s.symbols('y0:4');D=coords();Paff=affine(model.reconstruct(u),s.Matrix(z))
    linear=s.Matrix(4,4,lambda i,j:sum(s.Poly(Paff[i,j],*z).coeff_monomial(v)*v for v in z))
    alpha=[s.diff(linear[0,1]+linear[2,3],v)/5 for v in z];assert alpha==[0]*4
    def closedjet(P):
        ll=s.Matrix([s.Poly(P[i,j],*z).coeff_monomial(v)for i,j in pairs for v in z])
        qq=s.Matrix([s.Poly(P[i,j],*z).coeff_monomial(z[k]*z[l])for i,j in pairs for k,l in m2])
        return list(ll[D['li'],:])+list(D['Pi']*qq)
    values=closedjet(Paff)
    assert {i:s.factor(v)for i,v in enumerate(values)if v}=={5:-s.Rational(2048,27)*eta,11:-s.Rational(3,2),15:s.Rational(27,8),26:192*eta,33:-64*eta,41:s.Rational(81,16),45:-64*eta,56:-192*eta}
    lifts.W=D['weights'];lifts.ACT=[lifts.action_lists(A)for A in D['acts']];lifts.sel=D['li'];lifts.P4=np.array((4*D['Pi']).tolist(),dtype=object)
    from toric_coframe_fiber import affine_bivector,Y
    eye=np.eye(4,dtype=int).astype(object)
    calibrated=closedjet(affine_bivector(s.eye(4)).subs(dict(zip(Y,z)),simultaneous=True))
    assert list(lifts.jets(eye,eye,1))==[10*v for v in calibrated[:16]]+[100*v for v in calibrated[16:]]
    relation_path=ROOT/'reports/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.json';stored=load(relation_path);assert stored['all_passed']
    for name,digest in stored['input_sha256'].items():assert sha(ROOT/name)==digest,name
    highest={weight:[{tuple(mon):s.Rational(c)for mon,c in poly}for poly in polys]for weight,polys in stored['highest_lifts_in_projected_coordinates'].items()}
    relations=lifts.orbit(highest);assert len(relations)==274
    # Replay the full QQ root-stability and all six-variable dense coframe
    # coefficient identities, so the selected polynomial is GLOBAL on T01.
    stable=lifts.stability(relations);assert stable==814
    cached=lifts.chart_check(relations);assert cached==stored['cached_monomial_images']
    selected=highest['(1, 1)'][0];assert len(selected)==403
    evaluations=[s.factor(sum(c*s.prod(values[i]for i in mon)for mon,c in f.items()if all(values[i]for i in mon)))for f in relations]
    selected_value=s.factor(sum(c*s.prod(values[i]for i in mon)for mon,c in selected.items()if all(values[i]for i in mon)))
    assert selected_value==-s.Rational(1163264,135)*eta**2
    assert sum(v!=0 for v in evaluations)==39
    surviving=[[list(mon),str(c)]for mon,c in selected.items()if all(values[i]for i in mon)]
    assert len(surviving)==4
    print('Actual T01 GLOBAL witness and fresh274 coframe/stability proof PASS',flush=True)

    # All CA/CB source flags must lie in P(ker K_beta); its sole point
    # when eta!=0 is e_w, which fails both closed beta-image constraints.
    bundle_path=ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63/report.json';bundle=load(bundle_path);bi=bundle['beta_indices_zero_based']
    beta=s.Matrix(bundle['beta16_to_intrinsic24_adapter'])*u[bi,0]
    K=s.Matrix([[beta[4*i+j]for j in range(4)]for i,ij in enumerate(pairs)if ij!=(2,3)])
    assert K==s.Matrix([[0,0,0,0],[0,0,s.Rational(27,8),0],[s.Rational(3,2),0,0,0],[0,0,0,0],[0,-s.Rational(2048,27)*eta,0,0]])
    assert K.nullspace()==[s.Matrix([0,0,0,1])]
    minor=s.factor(K[[1,2,4],[0,1,2]].det());assert minor!=0 and s.cancel(minor/eta).is_Rational
    off=s.eye(5);off[1:,1:]=s.Matrix([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])
    uoff=s.Matrix(model.normalize(pushforward(model.reconstruct(u),off))['parameters']);constraints={}
    for name,fa in bundle['families'].items():
        witness=s.Matrix(fa['beta_closed_constraints'])*uoff[bi,0]
        assert [(i,v)for i,v in enumerate(witness)if v]==[(1,s.Rational(2048,27)*eta),(3,s.Rational(3,2))]
        constraints[name]=list(map(str,witness))
    print('ENTIRE CA/CB proper-source unique flag excluded exactly PASS',flush=True)

    out.mkdir(parents=True,exist_ok=False)
    certificate=dict(raw150=list(map(str,raw)),all350=list(map(str,all350)),S61=list(map(str,u)),primitive_alpha=list(map(str,alpha)),closedjet61=list(map(str,values)),
       selected_global_T01_polynomial=[ [list(mon),str(c)]for mon,c in sorted(selected.items())],selected_value=str(selected_value),four_surviving_terms=surviving,
       all274_evaluations=list(map(str,evaluations)),K_beta=[[str(v)for v in row]for row in K.tolist()],K_beta_unit_minor_on_eta_open=str(minor),
       flag_transport=[[str(v)for v in row]for row in off.tolist()],all_CA_CB_beta_constraints_at_sole_flag=constraints)
    archive=out/'certificate.json.gz'
    with archive.open('wb')as stream:
        with gzip.GzipFile(fileobj=stream,mode='wb',filename='',mtime=0)as zipped:zipped.write((json.dumps(certificate,separators=(',',':'))+'\n').encode())
    inputs=[Path(__file__),ROOT/'research/notes/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63.md',axis_path,
       ROOT/'reports/Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63/report.json',relation_path,bundle_path,
       ROOT/'reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json',ROOT/'reports/Q_CA_CB_N_ALL_A_20261103_v63_REVIEW.json',
       ROOT/'code/toric_zero_fiber_coframe_quadrics_20261102.py',ROOT/'code/toric_closed_jet_coordinates_20261102.py',
       ROOT/'code/toric_coframe_fiber.py',ROOT/'code/s61_model.py',ROOT/'code/q51_universal_section_20261101.py',
       ROOT/'research/notes/TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102.md']
    result=dict(id='Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63',all_passed=True,field='QQ[eta], nonzero ray over C',
       physical_raw150_divergence_and_all350_zero=True,actual_quintic='-16*x0^5/81',full_S61_normalization_134_and150_verified=True,
       primitive_alpha_zero=True,Pi_versus_four_Pi_scaling_independently_calibrated=True,
       all274_global_T01_lifts_QQ_root_stability_and_six_variable_chart_replayed=True,nonzero_root_images=stable,
       selected_global_T01_witness='highest (1,1), index0,403terms',selected_value=str(selected_value),all274_nonzero_evaluation_count=39,
       all_nonzero_ray_points_outside_actual_T01=True,
       K_beta_rank3_exact_eta_open_minor=str(minor),unique_entire_projective_flag='e_w',
       sole_flag_fails_both_actual_CA_CB_base_constraints=True,all_nonzero_ray_points_outside_entire_CA_CB_and_N=True,
       eta_zero_is_known_Qa_eight_over_nine=True,
       positive_ray_membership_in_other_known_components='unresolved',new_component_or_exhaustion_claimed=False,
       actual_full_local_branch_or_nonflatness_consequence_claimed=False,
       input_sha256={str(path.relative_to(ROOT)):sha(path)for path in inputs},certificate_sha256=sha(archive),seconds=time.monotonic()-start)
    (out/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print('POSITIVE WEIGHT6 RAY ACTUAL T01/CA/CB EXCLUSIONS PASS',result['seconds'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    build(parser.parse_args().output)
