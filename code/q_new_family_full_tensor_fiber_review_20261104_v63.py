#!/usr/bin/env python3
"""Independent conceptual S61 rational16 review plus tiny universal identities."""
import argparse,gzip,hashlib,json,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main(out):
 start=time.monotonic();out.mkdir(parents=True,exist_ok=False)
 stem='Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63';base=ROOT/'reports'/stem
 primary=read(base/'report.json');assert primary['all_passed']
 paths=[base/'report.json',base/'certificate.json',base/'operator-recipe.json.gz',
  ROOT/f'research/notes/{stem}.md',ROOT/'code/q_new_family_full_tensor_fiber_preflight_20261104_v63.py',
  ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/report.json',
  ROOT/'reports/Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW/report.json',
  ROOT/'reports/Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW/report.json',
  ROOT/'reports/Q_NEW_FAMILY_DENSE_S61_MAP_20261103_v63_REVIEW.json',
  ROOT/'reports/Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW/report.json',
  ROOT/'reports/Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW/report.json',
  ROOT/'reports/Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW/report.json',
  ROOT/'reports/Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63_REVIEW/report.json',
  ROOT/'code/s61_model.py',ROOT/'code/r4_positive_weight6_geometry_20261103_v63.py']
 for rel,h in primary['input_sha256'].items():assert sha(ROOT/rel)==h,rel
 for n,h in primary['evidence_sha256'].items():assert sha(base/n)==h,n
 oldrank=read(paths[5]);assert oldrank['full_exact_function_field_rank']==25 and oldrank['source_part_nonzero']
 assert all(read(p)['all_passed']for p in paths[6:13])
 recipe=read(base/'operator-recipe.json.gz');assert recipe['raw26_kernel']['compensating_perspective4']==['0']*4
 assert recipe['raw26_kernel']['all30_original_polynomial_dot_products_exactly_zero']
 assert len(recipe['full_family']['full_affine_P'])==4
 # Universal perspective naturality with arbitrary16-entry linear frame.
 S=s.Matrix(4,4,s.symbols('m0:16'));a2=s.Matrix(s.symbols('a0:4'));at=S.T*a2
 block=s.diag(1,1,1,1,1);block[1:5,1:5]=S
 def persp(v):
  p=s.eye(5)
  for j in range(4):p[0,j+1]=v[j]
  return p
 assert(persp(a2)*block-block*persp(at)).applyfunc(s.expand)==s.zeros(5)
 # Rational source-field tau and actual reference-pair alpha, no point.
 u,t,z,rf,sf,A=s.symbols('u t z rf sf A');tt=(rf/sf-1)/(rf/A-1)
 assert s.cancel(tt.subs({rf:sf,sf:rf},simultaneous=True)-tt/(tt-1))==0
 zz=u*(t-2)/t;inv=2*u/(u-z)
 assert s.cancel(zz.subs({u:-u,t:t/(t-1)},simultaneous=True)-zz)==0
 assert s.cancel(zz.subs(t,inv)-z)==0 and s.cancel(inv.subs(z,zz)-t)==0
 xi,r,ss,q,B=s.symbols('xi r ss q B');h=ss**2-25*r*q**2
 X2=ss/(48*u);c2=xi*B*u/40;v2=-32*u*ss/h;W2=(q+1)/3072
 alpha3=s.cancel((48*X2*u+1)/(240*B));alpha4=s.cancel(-4*c2*u/15).subs(u*u,1600*r)
 assert alpha3==(ss+1)/(240*B)and alpha4==-s.Rational(32,3)*xi*B*r
 sectioneq=s.together(v2*((384*X2)**2-(3072*W2-1)**2)+256*384*X2)
 num=s.fraction(sectioneq)[0];assert s.rem(s.Poly(num,u),s.Poly(u*u-1600*r,u)).is_zero
 proofs=dict(independent_normalized_pair_argument=[
  'Unique raw26 kernel has perspective4 zero, so source6+GL4 22columns rank21.',
  'Value orbit has dimension6 and rational Darboux frame: normalized RAWpair imageY has dimension15.',
  'Actual t2 pairsection overE5 has beta_bar and displayed rationalalpha_bar; itsSp4saturation projects densely toB15, so isdenseY.',
  'GxS5 ->Y ->B15 has birational composite from acceptedB15 theorem; functionfield sandwich makesY->B15 birational.',
  'Therefore genericprimitivebeta matching also matchesalpha; no primitive-only assumption orODEintegration.'],
  finite_frame='P line from C2P=PCt,D2P=PDt is S^-T up toscalar. Nonzeroalpha2[3] fixesQ=P/lambda; physicalS=Q^-T isrational.',
  dominance='Persp(alpha2)*blockS=blockS*Persp(alphat); newpipeline isphysicalSp4 transform ofoldnormalization. Historicaln isperspectivecross-section, not assumedLee0.',
  same_base='Use physicalh^-1 push for group variableh; primitivebeta isrho(h)beta_bar, exactly the acceptedB15 birationalproduct.',
  full_tau='All-degreefamilytau and p_tau=D*p; RtauD=HR, Stau=HSH, alphatau=Halpha,gu_tau=Hgu cancel onFULLprojectivetensor beforeEuler/S61 normalization.',
  fixed_field='F6=F5(t)=F5(z), tau(z)=z, F5=E5(u),u2=1600r ->F6tau=E5(z).',
  source_seed=dict(u='u,u²=1600r',v='-32us/(s²-25rq²)',c='xi*B*u/40',X='s/(48u)',W='(q+1)/3072',B='B',t='2'),
  general_source='Evaluate oldfirstintegrals atseed, retainfiveconstants; set t=2u/(u-z), useoldrationalinverse. No unspecified recovery.',
  Luroth='DominantGxSpecE5(z)->C16 overbirationalGxSpecE5=B15 gives SAMEK inclusionK<QQ(C16)<K(z). trdeg1; Luroth QQ(C16)=K(w).',
  source_z_degree_or_fullscheme_reducedness_not_claimed=True)
 (out/'independent-proof.json').write_text(json.dumps(proofs,indent=2)+'\n')
 note=ROOT/f'research/notes/{ID}.md'
 result=dict(id=ID,all_passed=True,method='Independent conceptual proof; gate hashes and tiny universal symbolic identities only',
  perspective0_rank21_and_normalized_RAWpair_dimension15_reviewed=True,
  independent_pair_to_beta_birational_field_sandwich=True,
  rawalpha_matching_not_inferred_from_primitive_beta_alone=True,
  rational_finite_matching_frame_not_ODE_integration=True,
  correct_dual_physical_S_and_gu_inverse_reviewed=True,
  independent_universal_perspective_matrix_identity=True,
  full_tensor_tau_descent_from_all_degree_family_identity_reviewed=True,
  independent_flow_tau_z_inverse_and_reference_alpha_identities=True,
  actual_full_S61_Sp4_dominance_and_historical_perspective_quotient_reviewed=True,
  same_K_B15_field_inclusion_into_K_z_reviewed=True,
  generic_reduced_fiber_K_rational_by_Luroth=True,
  actual_reduced_S61_C16_QQ_and_C_rational_dimension16=True,
  source_z_map_degree1_or_entire_component25_rationality_claimed=False,
  no_new_source26_matrix_full150_expansion_point_prime_GB_or_author_rerun=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in [Path(__file__),note,*paths]},
  independent_proof_sha256=sha(out/'independent-proof.json'),seconds=time.monotonic()-start,frozen=True)
 (out/'report.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items()if k!='input_sha256'},indent=2),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
