#!/usr/bin/env python3
"""Small exact checks and finite full-tensor operator recipe; no expansion."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 paths={
  'six':'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json',
  'kernel':'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/kernel-certificate.json.gz',
  'constants':'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json',
  'family':'reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz',
  'geometry':'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json',
  'dense':'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63/report.json',
  'source_review':'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json',
  'centralizer_review':'reports/Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW/report.json',
  'rational15':'reports/Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63/report.json',
  'orbit_review':'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW/report.json',
  'group_convention':'reports/Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63/report.json',
  'target_review':'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW/report.json',
 }
 d={k:read(ROOT/v)for k,v in paths.items()}
 for k in ['geometry','dense','source_review','centralizer_review','orbit_review','group_convention','target_review']:assert d[k]['all_passed'],k
 assert d['rational15']['birational_Sp4_times_A5_to_actual_B_NEW15_over_QQ_and_C']
 assert d['kernel']['compensating_perspective4']==['0']*4
 u,v,c,X,W,B,t,z=s.symbols('u v c X W B t z');tau={u:-u,v:-v,c:-c,X:-X,W:W,B:B}
 R=s.Matrix(d['six']['rational_Darboux_frame4']);alpha=s.Matrix(list(map(s.sympify,d['six']['primitive_perspective_alpha'])))
 D=s.diag(-1,1,-1,1);H=s.diag(-1,-1,1,1);gu=s.diag(u,1/u,1,1)
 assert(R.subs(tau,simultaneous=True)*D-H*R).applyfunc(s.cancel)==s.zeros(4)
 assert(gu.subs(u,-u)-H*gu)==s.zeros(4)
 assert(alpha.subs(tau,simultaneous=True)-H*alpha).applyfunc(s.cancel)==s.zeros(4,1)
 assert alpha[3]==-4*c*u/15
 xi,r,ss,q=s.symbols('xi r s q')
 abar=(gu.T*alpha).subs({X:ss/(48*u),c:xi*B*u/40},simultaneous=True)
 abar=abar.applyfunc(lambda x:s.cancel(x).subs(u*u,1600*r))
 assert abar==s.Matrix([0,0,(ss+1)/(240*B),-s.Rational(32,3)*xi*B*r])
 zz=u*(t-2)/t;tt=2*u/(u-z);tau_t=t/(t-1)
 assert s.cancel(zz.subs({u:-u,t:tau_t},simultaneous=True)-zz)==0
 assert s.cancel(zz.subs(t,tt)-z)==0
 assert s.cancel(tt.subs(z,zz)-t)==0
 assert s.cancel(tt.subs(u,-u)-tt/(tt-1))==0
 # Exact finite parameter identity was proved for FULL P in the pinned geometry
 # packet, not inferred from a jet calculation. Specialize its literal recipe.
 eq=next(e for e in d['geometry']['parameter_equivalences']if e['sign']==-1)
 a=s.Symbol('a')
 assert s.Matrix(eq['coordinate_linear']).subs(a,-1)==D
 assert s.sympify(eq['z_translation']).subs(a,-1)==0
 assert s.sympify(eq['tensor_scalar']).subs(a,-1)==-1
 J=s.zeros(4);J[0,1]=J[2,3]=1;J[1,0]=J[3,2]=-1
 assert gu.T*J*gu==J and H.T*J*H==J
 # Verify rawpush/pullback conversion symbolically for a general symplectic
 # diagonal; the global non-diagonal convention is a pinned independent input.
 k,l=s.symbols('k l',nonzero=True);A=s.diag(k,1/k,l,1/l)
 assert J*A==A.T.inv()*J and A.T*J==J*A.inv()
 operators=[
  dict(name='constant_source',field='F5=QQ(U,K,Iv,Ic,Ib)=E5(u)',E5='QQ(xi,r,s,q,B)',extension='u^2=1600*r',
       exact_coordinate_recipe='use embedded first_integral and t2 records, not the flow r,s as shape r,s'),
  dict(name='flow_source',field='F6=F5(t)',inverse='embedded constants inverse_r,inverse_s,inverse_v,inverse_c,inverse_B and inverse_X,inverse_W',
       WARNING='flow r,s in constants record are NOT birational-compression shape r,s'),
  dict(name='full_tensor',input='embedded actual full_affine_P with n=(u,1,v)',point='(X,0,Z,W)',
       Z='-3/(4*u)*(B+W/1536-W^2+X^2/64+2/9)',
       homogenization='homogenize each degree<=2 affine coefficient to degree2 in h, keep row/col h zero; scalar c',
       translation='physical (h,x) -> (h,x-point*h)',Darboux='physical block diag(1,R_t)',
       pushforward='M*P(M^-1*(h,x))*M^T, acting on ALL150 homogeneous coefficients'),
  dict(name='reference_rawjet',parameter='t=2',raw_normalized_linear='from embedded six record before perspective',
       alpha='(0,0,(48*X*u+1)/(240*B),-4*c*u/15)',
       generic_free_reason='primitive beta reference is dense in actual B15 quotient and has identity stabilizer'),
  dict(name='finite_matching_frame',inputs='Ct,Dt,C2,D2 and alpha_t,alpha2',
       linear_system='C2*P=P*Ct, D2*P=P*Dt',
       solve='lex-first nonzero15x15 minor with one remaining entry normalized to1; Cramer cofactors',
       kernel_dimension='1, by scalar common GL4 commutant and same geometric orbit',
       scalar='lambda=(P*alpha_t)[3]/alpha2[3]',Q='P/lambda',
       PHYSICAL_rawpush_frame='S=Q^(-T)',
       convention='rawpush S => beta_new=(S^-1)^*beta_old, Cnew=S^-T*Cold*S^T, alpha_new=S^-T*alpha_old',
       proof_of_orbit_existence='for normalized RAW pair(beta,alpha), compensator deltaR*R^-1-R*A*R^-1 is sp4 since valueJ fixed; pair quotient field is E5 because beta field E5 is contained, all pair invariants lie in kerdelta F5, and FULLtau makes them fixed E5; rational t2 pair section is free and represents generic pair orbit'),
  dict(name='full_match',physical_pipeline='block(g_u^-1) * perspective(alpha2) * block(S_t) * block(R_t) * translation(-point_t)',
       perspective='5x5 identity with row0 columns1..4=alpha2',
       resulting_primitive_beta='root beta_bar16(xi,r,s,q,B), independent of t',
       after_pipeline='Euler quotient and existing S61.normalize over the coefficient field'),
  dict(name='full_tau_descent',source='tau(u,v,c,X,W,B)=(-u,-v,-c,-X,W,B),tau(t)=t/(t-1)',
       FULL_source_identity='c_tau*P_tau=D_*(c*P), by exact eps=a=-1 family identity',
       identities=['R_tau*D=H*R','alpha_tau=H*alpha','S_tau=H*S*H','g_-u=H*g_u'],
       S_identity_reason='unique rawjet matching frame, not just infinitesimal integration',
       conclusion='ALL full projective tensor coefficients after pipeline are tau-fixed'),
  dict(name='rational_six_field',z='u*(t-2)/t',t='2*u/(u-z)',
       fixed_field='F6^tau=E5(z)',
       implementation='perform all finite operations in E5(z)[u]/(u^2-1600*r); odd coefficients vanish by proven full equivariance; take even coefficient'),
 ]
 recipe=dict(operators=operators,full_family=d['family'],six_source_record=d['six'],
             first_integrals_and_inverse=d['constants'],raw26_kernel=d['kernel'],
             full_parameter_tau_record=eq,
             no_expanded_intertwiner_full150_or_ODE_integration_performed=True)
 with gzip.open(out/'operator-recipe.json.gz','wt')as f:json.dump(recipe,f,separators=(',',':'));f.write('\n')
 cert=dict(compensating_perspective_identically_zero=True,
  raw_trace_covector_fourth_entry='-4*c*u/15',R_tau_D_equals_H_R=True,
  rational_reference_pair_covector=['0','0','(s+1)/(240*B)','-32*xi*B*r/3'],
  physical_group_convention='matching P proportional to S^-T; scalar fixed by P*alpha_t=lambda*alpha2',
  beta_positive_weight_g_u_uses_PHYSICAL_g_u_inverse=True,
  full_tau_equivalence_from_parameter_record=True,full_fixed_field='E5(z), z=u*(t-2)/t',
  all_small_rational_matrix_and_coordinate_identities_checked=True,
  proposed_conceptual_conclusion='generic S61 fibre is K-rational by embedding in K(z) and Luroth; independent proof audit pending',
  no_full_component25_or_global_fibre_or_t_degree_one_claim=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 inputs=[Path(__file__),ROOT/'code/s61_model.py',ROOT/'code/r4_positive_weight6_geometry_20261103_v63.py',*[ROOT/v for v in paths.values()]]
 report=dict(id=ID,all_passed=True,finite_full_tensor_operator_recipe=True,
  saved_raw_kernel_has_exactly_zero_perspective=True,correct_dual_physical_frame_and_full_tau_descent_proof_recorded=True,
  no_new_point_or_prime_full_tensor_coefficient_expansion_or_source26=True,
  generic_curve_Luroth_and_S61_rational16_proof_pending_independent_conceptual_review=True,
  full_projective25_rationality_claimed=False,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in inputs},
  evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','finite_full_tensor_operator_recipe','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
