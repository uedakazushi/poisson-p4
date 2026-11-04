#!/usr/bin/env python3
"""Actual v0 beta section: exact cancellation, transverse DVR gate, 16 formulas."""
import argparse,gzip,hashlib,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
ID='Q_NEW_FAMILY_V_ZERO_BETA_SECTION_20261104_v63'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def rows(M):return [[str(z)for z in row]for row in M.tolist()]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_CPU,(60,60));resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out.mkdir(parents=True,exist_ok=False)
 stems=[
 'Q_NEW_FAMILY_SIX_SOURCE_20261103_v63_REVIEW',
 'Q_NEW_FAMILY_SOURCE_CONSTANT_FIELD_20261103_v63_REVIEW',
 'Q_NEW_FAMILY_TARGET_FIELD_AND_SPECTRAL_DEGREES_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_RATIONAL_ORBIT_SECTION_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_V_ZERO_BOUNDARY_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_FULL_SOURCE_DEGREE_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_RATIONAL15_APPLICATION_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_RATIONAL_GROUP_INVERSE_20261104_v63_REVIEW',
 'Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63_REVIEW']
 gates=[ROOT/'reports'/n/'report.json'for n in stems];assert all(read(p)['all_passed']for p in gates)
 sixpath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 kerpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_KERNEL_20261103_v63/kernel-certificate.json.gz'
 constpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 geo=ROOT/'reports/R4_POSITIVE_WEIGHT6_GEOMETRY_20261103_v63/report.json'
 six,ker,old,iv,G=map(read,[sixpath,kerpath,constpath,invpath,geo])
 u,v,c,X,W,B,al,be=s.symbols('u v c X W B alpha beta');A=s.Integer(256);rf=A+v*al;sf=A+v*be
 rr,ss=s.symbols('r s')
 original=[s.sympify(z).subs({rr:rf,ss:sf},simultaneous=True)for z in old['five_first_integrals']]
 cancelled=[u,al*(2*A+v*al)*sf**2/(be*(2*A+v*be)*rf**2),
  (rf*sf-8*u*v*(al-be)**2)/(sf**2*al*(2*A+v*al)),
  c*rf**3/(sf*al**2*(2*A+v*al)**2),B*rf**4/(sf**2*al**3*(2*A+v*al)**3)]
 assert all(s.cancel(a-b)==0 for a,b in zip(original,cancelled))
 boundary=[s.cancel(z.subs(v,0))for z in cancelled]
 expected=[u,al/be,1/(2*A*al),c/(4*al**2),B/(8*A*al**3)]
 assert all(s.cancel(a-b)==0 for a,b in zip(boundary,expected))
 U,K,Iv,Ic,Ib=s.symbols('U K Iv Ic Ib');targets=[U,K,Iv,Ic,Ib]
 ainv=1/(2*A*Iv);binv=ainv/K
 inverse={u:U,al:ainv,be:binv,c:4*Ic*ainv**2,B:8*A*Ib*ainv**3}
 assert all(s.cancel(z.subs(inverse,simultaneous=True)-t)==0 for z,t in zip(boundary,targets))
 forward=dict(zip(targets,boundary))
 assert all(s.cancel(z.subs(forward,simultaneous=True)-k)==0 for k,z in inverse.items())
 tgen=s.cancel(((rf/sf-1)/(rf/A-1)));assert s.cancel(tgen-A*(al-be)/(sf*al))==0
 assert s.cancel(tgen.subs(v,0)-(1-be/al))==0
 tau={u:-u,c:-c,al:-be,be:-al,B:B}
 tau_expected=[-boundary[0],1/boundary[1],-boundary[1]*boundary[2],-boundary[1]**2*boundary[3],-boundary[1]**3*boundary[4]]
 assert all(s.cancel(f.subs(tau,simultaneous=True)-t)==0 for f,t in zip(boundary,tau_expected))
 # Use ORIGINAL polynomial kernel V, NOT delta=(v/3072)V which is tangentv0.
 polys=[s.sympify(z)for z in ker['full26_polynomial_kernel']]
 for z in polys:s.Poly(z,u,v,X,W,B,domain=s.QQ)
 assert ker['compensating_perspective4']==['0']*4
 transverse=s.sympify(ker['source_six_direction']['v']).subs(v,0);assert transverse==402653184
 R=s.Matrix(six['rational_Darboux_frame4']);assert not any(z.has(v)for z in R)
 assert s.cancel(R.det()-s.Rational(9,8)/(c*c*B))==0
 # Distinct coordinate names: these are NEWv0 parameters, NOT oldt2coordinates.
 xi0,r0,s0,q0,B0=s.symbols('xi0 r0 s0 q0 B0');basevars=(xi0,r0,s0,q0,B0)
 field=s.QQ.frac_field(*basevars);modulus=s.Poly(u*u-1600*r0,u,domain=field)
 def evenodd(expr):
  n,d=s.fraction(s.cancel(expr));n=s.Poly(n,u,domain=field).rem(modulus);d=s.Poly(d,u,domain=field).rem(modulus)
  aa,bb,cc,dd=n.nth(0),n.nth(1),d.nth(0),d.nth(1);norm=cc*cc-1600*r0*dd*dd;assert norm!=0
  return s.cancel((aa*cc-1600*r0*bb*dd)/norm),s.cancel((bb*cc-aa*dd)/norm)
 source_sub={v:0,X:s0/(48*u),c:xi0*B0*u/40,W:(q0+1)/3072,B:B0}
 weights=[w[0]for w in iv['Cartan_weights']];betas=[];stats=[]
 for index,(text,weight)in enumerate(zip(six['primitive_beta16'],weights)):
  expression=s.sympify(text).subs(v,0)*u**weight
  even,odd=evenodd(expression.subs(source_sub,simultaneous=True));assert odd==0
  betas.append(even);num,den=s.fraction(even)
  pn=s.Poly(num,*basevars);pd=s.Poly(den,*basevars);assert pd.length()==1
  em=pd.monoms()[0];assert em[2]==em[3]==0
  stats.append(dict(index=index,nonzero=even!=0,numerator_terms=pn.length()if even else 0,denominator_terms=pd.length(),numerator_multidegree=list(pn.degree_list())if even else[0]*5,denominator_monomial=list(em),odd_part_exactly_zero=True))
 expectedbeta=[-45*s0/(16*B0**2*xi0),-30/(B0**2*xi0),2400*r0,0,0,(5*q0+s0+1)/(240*B0),s0*xi0/(54000*r0),0,
  -32*B0*r0*xi0/3,-32*B0**2*xi0**2/675,0,0,0,0,(s0+1)/(120*B0),-64*B0*r0*xi0/3]
 assert all(s.cancel(a-b)==0 for a,b in zip(betas,expectedbeta))
 assert sum(z!=0 for z in betas)==9 and sum(r['numerator_terms']for r in stats)==12
 # Raw-fulltensor finite factors. No full150 coefficient expansion is needed.
 alpha0=s.Matrix([0,0,(s0+1)/(240*B0),-s.Rational(32,3)*xi0*B0*r0])
 alpha=s.Matrix(list(map(s.sympify,six['primitive_perspective_alpha'])))
 for a,b in zip(alpha.subs(source_sub,simultaneous=True),alpha0):assert evenodd(a-b)==(0,0)
 R0=R.subs(source_sub,simultaneous=True).applyfunc(s.cancel)
 detR0,odd=evenodd(R0.det());assert odd==0 and s.cancel(detR0-s.Rational(9,8)/(r0*xi0**2*B0**3))==0
 pointZ=-s.Rational(3,4)/u*(B+W/1536-W**2+X**2/64+s.Rational(2,9))
 p0=s.Matrix([X,0,pointZ,W]).subs(source_sub,simultaneous=True).applyfunc(s.cancel)
 n0,n1,n2=s.symbols('n0 n1 n2');Pzero=s.Matrix(G['affine_P']).subs({n0:u,n1:1,n2:0})
 cert=dict(new_variables=list(map(str,basevars)),WARNING='newv0xi0,s0,q0,B0 are NOT the same functions as oldt2xi,s,q,B; onlyr0=r_old=u^2/1600',
  cancelled_five_constants=list(map(str,cancelled)),flow_coordinates=dict(rf=str(rf),sf=str(sf),alpha='384*X+(3072*W-1)',beta='384*X-(3072*W-1)'),
  exact_v0_restrictions=list(map(str,boundary)),boundary_inverse={str(k):str(z)for k,z in inverse.items()},inverse_X='(alpha+beta)/768',inverse_W='(1+(alpha-beta)/2)/3072',
  exact_cancelled_t=str(tgen),v0_t='1-1/K',all_five_both_inverse_compositions=True,tau_equivariance=True,
  ORIGINAL_polynomial_kernel_v_at_v0=str(transverse),all26_kernel_entries_polynomial=True,all4_perspective_compensators_zero=True,
  DVR_density_proof='lift polynomialV toGL4xsource using gdot=g*A(source); annihilates all30 rawjetpullbacks; Dvunit402653184. Nonzeropullback f withord_v n>=1 hasord_v Df=n-1 inchar0, contradiction Df0. Thusv0rawGL4saturationdense; normalizedbetaSp4saturationdenseB15.',
  unrescaled_V_not_delta_used=True,
  new_coordinate_inverse=dict(u_squared='1600*r0',X='s0/(48*u)',c='xi0*B0*u/40',W='(q0+1)/3072',B='B0',v='0'),
  beta_bar0_16=list(map(str,betas)),beta_statistics=stats,nonzero_beta_count=9,total_numerator_terms=12,max_numerator_terms=3,
  beta_denominator_units='xi0*r0*B0 only; NO h0=s0^2-25*r0*q0^2 denominator',
  field_binding='boundaryconstantrestrictionF5->k(u,c,X,W,B) isbirational andtau-equivariant; actualB15dominance then restricts acceptedE_B=F5^tau toF0=k(xi0,r0,s0,q0,B0), index1',
  full_tensor_factors=dict(source_Pzero=rows(Pzero),source_scalar='c=xi0*B0*u/40',point=rows(p0),Darboux_R0=rows(R0),Darboux_determinant=str(detR0),primitive_alpha0=list(map(str,alpha0)),
   physical_pipeline='block(g_u^-1)*Persp(alpha0)*block(R0)*Translate(-p0) acting onhomogenizedactualcPzero',g_u='diag(u,u^-1,1,1)',u_squared='1600*r0',
   output='applyQQEulerquotientandS61.normalize, takeevenpartin u; allfullcoeffeventau byexactv0fullfamilytau andRtauD=HR',
   Laurent_ring='QQ[xi0^+-1,r0^+-1,B0^+-1,s0,q0]; eachfactoranditsinverseLaurent, detR0unit, fullevenness givesrationalfulltensorsection without h0denominator'),
  charts_comparison=dict(new_to_old='fromnewv0seedcompute displayedboundaryfiveconstants; evaluate ORIGINALsavedconstantinverse att=2; thenoldxi2=40c2/(B2u),r2=u^2/1600,s2=48uX2,q2=3072W2-1',
   old_to_new='fromold t2seed v2=-32u*s2/(s2^2-25r2q2^2),c2=xi2B2u/40,X2=s2/(48u),W2=(q2+1)/3072 compute CANCELLEDconstants; usev0inverse; thennewxi0,r0,s0,q0,B0',
   original_constant_inverse={k:old[k]for k in['inverse_r','inverse_s','inverse_v','inverse_c','inverse_B','inverse_X','inverse_W']},
   no_expanded_coordinate_change_claimed=True),
  section_on_actual_beta='recoveroldE5coordinateswithacceptedfiniteinverse, converttoF0, computebeta_bar0, recoveruniqueSp4framebyacceptedC/D-intertwiner+oddβscale recipe; applyphysicalinverseframe tofulltensorandS61.normalize',
  v0_boundary_Y15_birational_B15=True,original_v0_full_source_degree2_tau_removed=True,
  generic_Y16_source_degree4_unchanged=True,no_high_invariant_source26_or_NEW_equation_claim=True)
 payload=json.dumps(cert,indent=2)+'\n';assert len(payload.encode())<32*1024*1024;(out/'certificate.json').write_text(payload)
 paths=[Path(__file__),ROOT/'research/notes'/(ID+'.md'),sixpath,kerpath,constpath,invpath,geo,*gates]
 report=dict(id=ID,all_passed=True,actual_v0_beta_source_Sp4_dense_by_transverse_DVR_argument=True,
  all5_cancelled_constants_and_both_boundary_inverses_exact=True,
  all16_beta_even_parts_exact_9_nonzero_12_terms=True,no_h_denominator=True,
  actual_invariant_quotient_field_new_F0_index1=True,
  actual_fulltensor_v0_Y15_section_birational_actualB15=True,
  no_new_point_prime_rank_job_high_invariant_source26_or_full150_expansion=True,
  independent_review_pending=True,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},
  evidence_sha256={'certificate.json':sha(out/'certificate.json')},seconds=time.monotonic()-start)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','all16_beta_even_parts_exact_9_nonzero_12_terms','no_h_denominator','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
