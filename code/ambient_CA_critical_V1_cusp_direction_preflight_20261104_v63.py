#!/usr/bin/env python3
"""Exact cusp-tangent direction from existing V1 matrices only.

No source action, source solve, new jet or field/point query.
"""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_CUSP_DIRECTION_PREFLIGHT_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(a):return s.Matrix([[s.Rational(c) for c in row] for row in a])
def em(M):return [[str(c) for c in r] for r in M.tolist()]
def inverse(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ).inv().to_Matrix()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    files={
      'relative':'reports/AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63/actual-V1-source-relative-tangent.json',
      'checkpoint':'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz',
      'setup':'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/setup.json.gz',
      'odd':'reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/certificate.json.gz',
      'oldjets':'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/jets.json.gz',
      'weights':'reports/AMBIENT_B18_POSITIVE_GRADING_OBSTRUCTION_20261104_v63/certificate.json'}
    R,C,S,O,J,W=map(read,[ROOT/p for p in files.values()])
    A=mat(O['all150_actual_odd_linear']);oldJ=mat(C['full150_by29_target_Jacobian'])
    jr=oldJ[:,C['complement28_columns']];Y=S['actual_target_Y28'];branch=S['actual_branch29th_row'];assert branch==22
    invY=inverse(jr[Y,:])
    L=jr[branch,:]*invY;wc=s.Matrix(list(map(s.Rational,S['corrected_quadratic150'])));weight=wc[branch];assert weight
    nu_linear=(-L/weight).row_join(s.Matrix([[1/weight]]))
    # Full even differential da/dR: source tangent28 plus the actual fold
    # square supplies all29 columns. Subtract da_i when using raw odd rows.
    even=(jr*invY-wc/weight*L).row_join(wc/weight)
    assert even[Y+[branch],:]==s.eye(29)
    assert even*oldJ[Y+[branch],:]==oldJ
    odd_target_rows=s.zeros(8,150)
    for j,row in enumerate(S['actual_odd_selected8_rows']):
        odd_target_rows[j,row]=1
        for k,base_row in enumerate(Y+[branch]):odd_target_rows[j,base_row]-=even[row,k]
    eight=S['actual_odd_selected8_rows'];assert eight==[27,28,29,35,60,65,67,68]
    I9=A[eight,:].col_join(nu_linear);selected=O['nine_columns'];minor=I9[:,selected];assert minor.det()!=0
    embedding=A[:,selected]*inverse(minor);assert embedding*I9==A
    assert embedding[eight,:]==s.eye(9)[:8,:]
    # Calibrate Z=nu*tau=p^3 in the old exact paired source, not raw row63.
    cubic=s.Matrix([s.Rational(row[3]) for row in J['stages'][0]['all150_output_series']])
    assert embedding[:,8]==cubic and embedding[63,8]==s.Rational(800,27)
    cb=mat(R['CB_normalizer_J150_by29_original_frame']);projection=mat(R['CB_to_CA_even_projection29'])
    cert=R['CB_to_CA_even_projection_certificate'];assert cert['rank']==28
    piv=cert['independent_columns'];coeff=mat(cert['all_column_factorization']);free=next(i for i in range(29) if i not in piv)
    kernel=s.zeros(29,1);kernel[free]=1
    for k,i in enumerate(piv):kernel[i]-=coeff[k,free]
    assert projection*kernel==s.zeros(29,1)
    normal=cb*kernel;assert any(normal) and normal[Y+[branch],:]==s.zeros(29,1)
    y=s.Matrix([normal[i] for i in eight])
    z=(normal[63]-sum(embedding[63,i]*y[i] for i in range(8)))/embedding[63,8]
    normalcoords=y.col_join(s.Matrix([z]));assert embedding*normalcoords==normal
    first=next(v for v in normalcoords if v);kernel/=first;normal/=first;normalcoords/=first
    frame=mat(C['normalizer_source29_frame']);fold=s.Matrix(list(map(s.Rational,C['kernel29'])))
    Nnormal=(frame[2,:]*fold)[0];assert Nnormal
    # True N adaptation. Unknown second derivatives of nu(R) disappear in
    # the mixed direction because J*fold=0; tau(fold)=1 is calibrated above.
    T=fold.row_join(s.eye(29)[:,C['complement28_columns']])
    mixed=mat(S['full150_by29_mixed_Hessian'])
    tau_source=(nu_linear*mixed[Y+[branch],:]/2)*inverse(T)
    assert (tau_source*fold)[0]==1
    h_source=frame[2,:]/Nnormal
    u_source=A[eight,:]*oldJ[Y+[branch],:]
    source_parameters=oldJ[S['orbit20_raw_rows'],:].col_join(u_source).col_join(tau_source)
    source_parameters_det=DomainMatrix.from_Matrix(source_parameters).convert_to(s.QQ).det()
    assert source_parameters_det
    columns=list(DomainMatrix.from_Matrix(u_source).convert_to(s.QQ).rref()[1]);assert len(columns)==8
    lam=(tau_source-h_source)[:,columns]*inverse(u_source[:,columns])
    assert lam*u_source==tau_source-h_source
    gamma=(lam*normalcoords[:8,0])[0]
    nx=-2*gamma
    # The actual derivative of x_N=(tau-beta_N)^2 is dnu-2sum lambda_i dYi.
    adapted_x_row=s.zeros(1,150)
    for j,row in enumerate(Y+[branch]):adapted_x_row[0,row]+=nu_linear[0,j]
    adapted_x_row-=2*lam*odd_target_rows
    assert (adapted_x_row*normal)[0]==nx
    adapted_projection=cb[S['orbit20_raw_rows'],:].col_join(A[eight,:]*projection).col_join(adapted_x_row*cb)
    adapted_projection_det=DomainMatrix.from_Matrix(adapted_projection).convert_to(s.QQ).det()
    assert adapted_projection_det
    weights=W['physical150_torus_weights'];characters=sorted({weights[i] for i,v in enumerate(normal) if v})
    assert len(characters)==1
    result={'actual_CB_source_kernel_of_CA_even_projection29':list(map(str,kernel)),
      'actual_CB_normal_direction_physical150':list(map(str,normal)),
      'cusp_even_nu_linear_in_old_R29':list(map(str,nu_linear)),
      'actual_even_tangent_embedding150_by29':em(even),
      'actual_odd_Y8_linear_physical150_after_even_subtraction':em(odd_target_rows),
      'actual_regular_I9_rows8_plus_nu':em(I9),'actual_cusp_odd_embedding150_by9':em(embedding),
      'cusp_odd_coordinate_order':[f'Y_row{i}=b_row{i}*tau' for i in eight]+['Z=nu*tau=p^3'],
      'cusp_direction_even_nu':'0','cusp_direction_even_u8':['0']*8,
      'cusp_direction_Y8':list(map(str,normalcoords[:8,0])),'cusp_direction_Z':str(normalcoords[8]),
      'calibration_Z_all150_equals_saved_exact_paired_cubic':True,'calibration_raw63_Z':'800/27',
      'source_N_binary_third_coefficient_on_CA_kernel':str(Nnormal),
      'N_source_divisor_transverse_to_fold_and_differential_immersed':True,
      'actual_trace_tau_linear_on_source29':list(map(str,tau_source)),
      'actual_N_divisor_linear_on_source29_normalized':list(map(str,h_source)),
      'actual_even_u8_linear_on_source29':em(u_source),
      'actual_source_parameter_order':['actual_orbit20_outputs','actual_even_b_selected8','trace_tau'],
      'actual_source_parameter_Jacobian29':em(source_parameters),
      'actual_source_parameter_Jacobian_determinant':str(s.Rational(source_parameters_det)),
      'actual_N_trace_beta_linear_coefficients_lambda8':list(map(str,lam)),
      'all29_source_beta_factorization_verified':True,
      'N_adapted_source_p':'p_N=tau-s_N(g20,u8); solve the actual Ndivisor uniquely for tau in regular source parameters(g20,u8,tau), then factor the defining function by a unit',
      'N_root_independent_of_old_nu':'nu=tau^2 is a source function, not a separate source variable; solving Ndivisor for tau yields s_N(g,u) with no nu argument. It vanishes when u8=0.',
      'N_adapted_target_x_linear_physical150':list(map(str,adapted_x_row)),
      'N_adapted_target_coordinate_conventions':{'x':'p_N^2','Yi':'u_i*p_N','Z':'p_N^3','N':'p_N=0'},
      'N_adapted_direction_nx':str(nx),'N_adapted_direction_nY8':list(map(str,normalcoords[:8,0])),
      'N_adapted_direction_nZ':str(normalcoords[8]),
      'N_adapted_x_projection_nonzero':bool(nx),
      'N_adapted_CB_even29_projection_order':['old_actual_orbit20_outputs','actual_b_selected8','x_N=p_N^2'],
      'N_adapted_CB_even29_projection_matrix':em(adapted_projection),
      'N_adapted_CB_even29_projection_determinant':str(s.Rational(adapted_projection_det)),
      'normal_direction_torus_character':characters[0],
      'no_new_source_or_jet_or_fieldpoint':True,
      'kappa_or_actual_D_type_not_computed':True,
      'next_method_gate':'Use exact normal Y/Z/nu orientation to derive a bounded formal cusp-core test; all additional ambient graph relations must still be included. No automatic source2 or higher jet.'}
    ev=out/'existing-V1-exact-cusp-direction.json';ev.write_text(json.dumps(result,indent=2)+'\n')
    gates=[ROOT/'reports'/p for p in [
      'AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63_REVIEW/report.json',
      'AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW/report.json',
      'AMBIENT_B18_ODD_FIRSTJETS_N_PROJECTION_20261103_v63_REVIEW.json']]
    assert all(read(p)['all_passed'] for p in gates)
    paths=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in files.values()]+gates+[
      ROOT/'code/ambient_B18_CA_critical_scheme_20261103_v63.py',
      ROOT/'code/ambient_B18_odd_firstjets_N_projection_20261103_v63.py',
      ROOT/'research/notes/AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63.md']
    report={'id':STEM,'all_passed':True,'actual_CB_normal_cusp_direction_identified_from_old_matrices':True,
      'n_nu':'0','n_Y8':result['cusp_direction_Y8'],'n_Z':str(normalcoords[8]),
      'old_n_nu_is_trace_frame_not_N_adapted':True,
      'N_beta_lambda8':list(map(str,lam)),'N_adapted_nx':str(nx),
      'N_adapted_nY8':list(map(str,normalcoords[:8,0])),'N_adapted_nZ':str(normalcoords[8]),
      'N_adapted_CB_even29_unit_minor':str(s.Rational(adapted_projection_det)),
      'normal_character':characters[0],'N_source_fold_transversality':str(Nnormal),
      'kappa_actual_D_type_or_new_jet_claimed':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','N_adapted_nx','N_adapted_nY8','N_adapted_nZ','N_beta_lambda8','normal_character','N_source_fold_transversality','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
