#!/usr/bin/env python3
"""ONE actual relative tangent/source gate at the existing accepted V1.

No new beta/source point, source jet, implicit solve or Groebner basis.
"""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
from ambient_CA_global_flag_model_20261103_v63 import X,vector,canonical,action
from ambient_CB_global_flag_model_20261103_v63 import raw as rawB
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63'
PAIRS=[(i,j) for i in range(5) for j in range(i+1,5)]
MONS=[(i,j) for i in range(5) for j in range(i,5)]
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def matrix(x):return s.Matrix([[s.Rational(c) for c in row] for row in x])
def em(M):return [[str(c) for c in row] for row in M.tolist()]
def fromvector(v):
    P=s.zeros(5)
    for k,(i,j) in enumerate(PAIRS):
        P[i,j]=sum(v[15*k+l]*X[a]*X[b] for l,(a,b) in enumerate(MONS));P[j,i]=-P[i,j]
    return P
def push(P,G):
    Gi=G.inv();subs=dict(zip(X,Gi*s.Matrix(X)));return (G*P.subs(subs,simultaneous=True)*G.T).applyfunc(s.expand)
def dm(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ)
def rank_certificate(M):
    independent=list(dm(M).rref()[1]);B=M[:,independent]
    rows=list(dm(B.T).rref()[1]);minor=B[rows,:];det=dm(minor).det();assert det
    factor=dm(minor).inv().to_Matrix()*M[rows,:]
    assert B*factor==M
    return {'rank':len(independent),'independent_columns':independent,'minor_rows':rows,
            'minor_determinant':str(s.Rational(det)),'all_column_factorization':em(factor)}
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    files={
      'V1checkpoint':'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz',
      'V1setup':'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/setup.json.gz',
      'V1odd':'reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/certificate.json.gz',
      'proper19':'reports/AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63/literal-proper-source14-preflight.json',
      'CBsourceflags':'reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/certificate.json.gz'}
    C,S,B,P,CB=map(read,[ROOT/files[q] for q in files])
    point=s.Matrix(list(map(s.Rational,C['point150'])));Pold=fromvector(point)
    assert vector(Pold)==point and canonical(Pold)==Pold
    assert C['exact_GL_orbit_rank']==20 and C['target_source_differential_rank']==28
    # V1 is the already accepted critical-V source point, not a new point.
    beta=s.Rational(1,10);F=s.eye(5);F[1,0]=-beta;assert F.det()==1
    Pprime=push(Pold,F);vprime=vector(Pprime)
    B19=matrix(P['all150_CB_fixed_flag19']);rows36=P['zero_coordinate_rows36']
    br=P['CB_constant_left_inverse_rows19'];bi=matrix(P['CB_left_inverse19'])
    coefficients=bi*vprime[[rows36[i] for i in br],0]
    assert B19*coefficients==vprime
    assert push(Pprime,F.inv())==Pold
    # Full nonlinear source-ratio derivative, not a copied CA column.
    h,x,y,z,w=X
    f=coefficients[0]*h*h+coefficients[1]*h*x+coefficients[2]*x*x
    swap=s.eye(5);swap[0,0]=swap[1,1]=0;swap[0,1]=swap[1,0]=1
    qcolumn=vector(push(canonical(rawB(0,-f,0)),swap))
    fiber=B19[:,:3].row_join(qcolumn).row_join(B19[:,3:])
    assert fiber.shape==(150,20)
    flags=CB['flag_positions'];flagcols=[]
    for i,j in flags:
        M=s.zeros(5);M[i,j]=1;flagcols.append(vector(action(M,Pprime)))
    Jprime=fiber.row_join(s.Matrix.hstack(*flagcols));assert Jprime.shape==(150,29)
    JCB=s.Matrix.hstack(*[vector(push(fromvector(Jprime[:,i]),F.inv())) for i in range(29)])
    # Complete CA tangent:28source image + actual branch-even direction +9odd
    # target directions. The old finite-image theorem proves completeness.
    JCA=matrix(C['full150_by29_target_Jacobian']);cols=C['complement28_columns'];assert len(cols)==28
    wc=s.Matrix(list(map(s.Rational,S['corrected_quadratic150'])))
    odd=matrix(B['all150_actual_odd_linear']);oc=B['nine_columns'];assert len(oc)==9
    RCA=S['actual_target_Y28']+[S['actual_branch29th_row']]
    assert RCA==[1,5,15,16,18,19,20,23,31,46,61,105,106,108,110,112,113,117,120,125,109,118,119,121,123,127,132,140,22]
    assert odd[RCA,:]==s.zeros(29,29)
    tangentCA=JCA[:,cols].row_join(wc).row_join(odd[:,oc]);assert tangentCA.shape==(150,38)
    certA=rank_certificate(tangentCA);assert certA['rank']==38
    certB=rank_certificate(JCB)
    union=tangentCA.row_join(JCB);certU=rank_certificate(union)
    projection=JCB[RCA,:];certP=rank_certificate(projection)
    exact_source_identities={'existing_V1_point_reconstructed_all150':True,
      'proper_Npair_frame5':em(F),'frame_determinant':'1','inverse_frame':em(F.inv()),
      'CB_fixed19_coefficients':list(map(str,coefficients)),
      'CB_fixed19_reconstruction_all150_exact':True,'frame_inverse_returns_exact_old_V1':True,
      'CB_binary_f_coefficients':list(map(str,coefficients[:3,0])),
      'ratio_q_derivative_rule':'-sum f_i*G_i; rawB(0,-f,0), at q=0',
      'CB_full20fiber_columns_prime_frame':em(fiber),'CB_flag9_positions':flags,
      'CB_normalizer_J150_by29_original_frame':em(JCB),
      'CA_complete_target_tangent150_by38_original_frame':em(tangentCA),
      'CA_tangent_columns_rule':{'source28_indices':cols,'even_branch_raw_row':22,'actual_odd9_base_columns':oc},
      'CA_even29_actual_raw_rows':RCA,'CB_to_CA_even_projection29':em(projection),
      'CA38_certificate':certA,'CB29_certificate':certB,'actual_union_certificate':certU,
      'CB_to_CA_even_projection_certificate':certP,
      'actual_tangent_intersection_dimension':38+certB['rank']-certU['rank'],
      'affine_cone_source_counts':[29,29],'no_extra_projective_or_scalar_dimension':True}
    terms=sum(1 for M in [tangentCA,JCB,fiber,projection] for c in M if c)
    terms+=sum(c!='0' for rec in [certA,certB,certU,certP] for row in rec['all_column_factorization'] for c in row)
    assert terms<30000 and terms<150000
    ep=out/'actual-V1-source-relative-tangent.json';ep.write_text(json.dumps(exact_source_identities,separators=(',',':'))+'\n');assert ep.stat().st_size<32*1024*1024
    valid=certB['rank']==29 and certU['rank'] in [38,39]
    clean=valid and certU['rank']==39
    gates=[ROOT/'reports'/x for x in [
      'AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW/report.json',
      'AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63/report.json',
      'AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_20261104_v63_REVIEW/report.json',
      'AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63_REVIEW.json',
      'AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63_REVIEW.json',
      'AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW/report.json',
      'AMBIENT_B18_T1_T2_QUADRATIC_PHI_20261104_v63_REVIEW/report.json',
      'AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63_REVIEW/report.json']]
    assert all(read(p)['all_passed'] for p in gates)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in files.values()]+gates+[
      ROOT/'code/ambient_CA_global_flag_model_20261103_v63.py',ROOT/'code/ambient_CB_global_flag_model_20261103_v63.py']
    report={'id':STEM,'all_passed':valid,'status':('PROVED_clean_actual_D_and_kappa_zero_at_existing_dense_V1' if clean else 'actual_tangent_D29_only_kappa_unresolved' if valid else 'STOP_unexpected_source_or_tangent_rank'),
      'CA_complete_tangent_rank':certA['rank'],'CB_actual_source_rank':certB['rank'],
      'actual_union_rank':certU['rank'],'actual_D_tangent_dimension':38+certB['rank']-certU['rank'],
      'CB_to_CA_even29_projection_rank':certP['rank'],
      'CB_to_CA_even_projection_unit':certP['rank']==29,
      'actual_D_and_N_smooth28_proved':clean,'generic_CAcritical_kappa_zero_proved':clean,
      'if_rank38_no_kappa_or_higher_order_conclusion':True,
      'all_source_frames_and_counts_same_existing_V1_affine_cone':True,
      'no_new_source_point_beta_jet_IFT_or_GB':True,'D0_origin_socle_or_depth_claimed':False,
      'stored_nonzero_QQ_coefficient_count':terms,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ep.name:sha(ep)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
      'independent_review_pending':True}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','status','actual_union_rank','actual_D_tangent_dimension','CB_to_CA_even29_projection_rank','seconds','max_rss_kib']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
