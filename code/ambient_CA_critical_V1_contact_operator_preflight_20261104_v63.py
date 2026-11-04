#!/usr/bin/env python3
"""Fixed-matrix inputs for ONE relative contact vector. No new source jet."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
from sympy.polys.matrices import DomainMatrix
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_CA_CRITICAL_V1_CONTACT_OPERATOR_PREFLIGHT_20261104_v63'
def read(p):return json.loads(gzip.decompress(p.read_bytes()) if p.suffix=='.gz' else p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mat(a):return s.Matrix([[s.Rational(c) for c in row] for row in a])
def em(M):return [[str(c) for c in row] for row in M.tolist()]
def inv(M):return DomainMatrix.from_Matrix(M).convert_to(s.QQ).inv().to_Matrix()
def main(out):
    start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out.mkdir(parents=True,exist_ok=False)
    names={
      'direction':'reports/AMBIENT_CA_CRITICAL_V1_CUSP_DIRECTION_PREFLIGHT_20261104_v63/existing-V1-exact-cusp-direction.json',
      'relative':'reports/AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63/actual-V1-source-relative-tangent.json',
      'checkpoint':'reports/AMBIENT_CA_V1_ODD_ORDER_CHECKPOINT_20261103_v63/checkpoint.json.gz',
      'setup':'reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/setup.json.gz'}
    D,R,C,S=map(read,[ROOT/p for p in names.values()]);J=mat(C['full150_by29_target_Jacobian']);B=mat(R['CB_normalizer_J150_by29_original_frame'])
    k=s.Matrix(list(map(s.Rational,C['kernel29'])));cols=C['complement28_columns'];Jr=J[:,cols];basis=s.eye(29)[:,cols]
    acoords=mat(D['actual_source_parameter_Jacobian29']);h=s.Matrix([list(map(s.Rational,D['actual_N_divisor_linear_on_source29_normalized']))])
    frame=mat(C['normalizer_source29_frame']);assert h==frame[2,:] and list(h)==[s.Integer(i==22) for i in range(29)]
    # Replace the last trace-tau row by the literal N equation h_N=k_2.
    params=acoords.copy();params[28,:]=h;assert params.det()!=0
    nx=s.Rational(D['N_adapted_direction_nx']);by=s.Matrix(list(map(s.Rational,D['N_adapted_direction_nY8'])))/nx
    rhs=s.zeros(29,1);rhs[20:28,0]=by
    a=inv(params)*rhs;assert (h*a)[0]==0 and (h*k)[0]==1
    cbproj=mat(D['N_adapted_CB_even29_projection_matrix']);b=inv(cbproj)*rhs
    xr=s.zeros(29,1);xr[28]=1;kb=inv(cbproj)*xr
    V=B*kb;assert B*b==J*a
    assert V==s.Matrix(list(map(s.Rational,D['actual_CB_normal_direction_physical150'])))/nx
    T=k.row_join(basis);H=mat(S['full150_by29_mixed_Hessian']);W=H[:,0]/2;M=H*inv(T)*a
    residual=V-M-W;Y=S['actual_target_Y28'];r=inv(Jr[Y,:])*residual[Y,:]
    assert Jr*r==residual
    rr=basis*r
    assert V==M+W+J*rr
    odd=mat(D['actual_odd_Y8_linear_physical150_after_even_subtraction']);even=mat(D['actual_even_tangent_embedding150_by29'])
    embed=mat(D['actual_cusp_odd_embedding150_by9']);Rrows=Y+[S['actual_branch29th_row']]
    z=s.zeros(1,150);z[63]=1
    for j,row in enumerate(Rrows):z[0,row]-=even[63,j]
    z=(z-embed[63,:8]*odd)/embed[63,8]
    assert z*J==s.zeros(1,29) and (z*W)[0]==0 and (z*M)[0]==0
    assert z*embed==s.Matrix([[0]*8+[1]])
    assert z*B==s.zeros(1,29)
    CA=mat(R['CA_complete_target_tangent150_by38_original_frame']);selected=R['CA38_certificate']['minor_rows'];left=inv(CA[selected,:])
    # Saved constant quotient map Pi(x)=x-CA*left*x[selected].
    assert CA*left*J[selected,:]==J and CA*left*V[selected,:]==V
    inputs={
      'CA_source_order':'ordered exp(group20 GL generators), followed by the fixed9 linear fiber complement',
      'CA_group_GL_indices':C['group20_GL_indices'],'CA_extra_fiber_indices':C['extra9_fiber_indices'],
      'CA_actual_N_equation':'h_N=source_coordinate22; fixed binary third coefficient; p=h_N',
      'CA_source_kappa29':list(map(str,k)),'CA_source_v29':list(map(str,a)),'CA_source_r29':list(map(str,rr)),
      'CB_source_order':'fixed20=(f0,f1,f2,q,H16), then ordered9 flag exponentials in saved prime frame; pull back by saved fixed frame inverse',
      'CB_fixed19_coefficients':R['CB_fixed19_coefficients'],'CB_fixed_frame':R['proper_Npair_frame5'],'CB_flag9_positions':R['CB_flag9_positions'],
      'CB_source_v29':list(map(str,b)),'CB_source_x29':list(map(str,kb)),
      'CB_mixed_fixed_fiber_required':'q*(-f), including all three q*delta_f_i columns, not just the center ratio derivative',
      'actual_v_even_u8_direction':list(map(str,by)),
      'actual_CAnormalizer_v150':list(map(str,J*a)),'actual_CB_x150':list(map(str,V)),
      'cached_CA_fold_square150':list(map(str,W)),'cached_CA_fold_v_mixed150':list(map(str,M)),
      'exact_decomposition':'CB_x=CA_H(kappa,v)+CA_coeff(p^2)+CA_J*r',
      'CA_target_tangent150_by38':em(CA),'CA_tangent_left_rows38':selected,'CA_tangent_left_inverse38':em(left),
      'canonical_Z_linear150':list(map(str,z)),
      'future_source_jobs':[{'family':'CA','directions':['kappa','v'],'precision':3,'monomials':10,'outputs':150},
                            {'family':'CA','directions':['v','r'],'precision':2,'monomials':6,'outputs':150},
                            {'family':'CB','directions':['v','x'],'precision':2,'monomials':6,'outputs':150}],
      'future_vector_formula':'Delta=CB_coeff(v*x)-CA_coeff(p*v^2)-CA_coeff(p^2*v)-CA_coeff(v*r)',
      'future_c_v_formula':'canonical_Z_linear150*Delta',
      'future_all_ambient_graph_bits':'Pi_CA*Delta, Pi_CA=Id-CA_tangent*left_inverse*select38',
      'future_mononomial_slots_output_bound':3300,'future_coefficient_hard_budget_seconds':60,'future_coefficient_hard_budget_bytes':2000000000,
      'source_coefficient_evaluation_executed':False,'no_new_source_jet_or_point_or_field':True,
      'conditional_completion':'If c_v!=1, PiDelta!=0 gives kappa=0; PiDelta=0 gives kappa=1. If c_v=1 this method stops inconclusive; no automatic higher order.'}
    ev=out/'literal-V1-contact-vector-input.json';ev.write_text(json.dumps(inputs,indent=2)+'\n')
    gates=[ROOT/'reports'/n/'report.json' for n in ['AMBIENT_CA_CRITICAL_V1_CUSP_DIRECTION_PREFLIGHT_20261104_v63','AMBIENT_CA_CRITICAL_V1_RELATIVE_TANGENT_20261104_v63_REVIEW','AMBIENT_CA_V1_ODD_ORDER_20261103_v63_REVIEW','AMBIENT_B18_COMPARISON_C_GENERIC_DVR_RANK_20261104_v63_REVIEW','AMBIENT_CA_GENERIC_V_FULL_COMPLETION_20261103_v63']]
    assert all(read(p)['all_passed'] for p in gates)
    files=[Path(__file__),ROOT/'research/notes'/f'{STEM}.md']+[ROOT/p for p in names.values()]+gates
    report={'id':STEM,'all_passed':True,'fixed_existing_source_directions_and_all150_decomposition_verified':True,
      'N_equation_exact_source_coordinate22':True,'CA_actual_u8_and_g20_functions_used':True,
      'CA_source_parameter_minor':str(params.det()),'CB_unit_even29_minor':str(cbproj.det()),
      'source_coefficient_job_executed':False,'kappa_decided':False,
      'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'evidence_sha256':{ev.name:sha(ev)},
      'seconds':time.monotonic()-start,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['id','all_passed','source_coefficient_job_executed','seconds']}))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
