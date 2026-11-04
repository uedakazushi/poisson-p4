#!/usr/bin/env python3
"""Counts/support bounds and exact input binding only; no new rank/expansion job."""
import argparse,gzip,hashlib,itertools,json,resource,signal,time
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1]
STEM='Q_NEW_FAMILY_FULL_DEGREE24_PREFLIGHT_20261104_v63'
def read(p):
 with(gzip.open(p,'rt')if p.suffix=='.gz'else p.open())as f:return json.load(f)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def necklaces(n):return [w for w in itertools.product(range(2),repeat=n)if w==min(w[i:]+w[:i]for i in range(n))]
def main(out):
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False)
 sixpath=ROOT/'reports/Q_NEW_FAMILY_SIX_PARAMETER_COVARIANTS_20261103_v63/certificate.json'
 sectionpath=ROOT/'reports/Q_NEW_FAMILY_SOURCE_RATIONAL_QUOTIENT_20261103_v63/certificate.json'
 spinpath=ROOT/'reports/Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63/certificate.json.gz'
 invpath=ROOT/'reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/invariant-certificates.json'
 ambientpath=ROOT/'reports/Q_NEW_FAMILY_DEGREE16_INVARIANT_PROBE_20261103_v63/report.json'
 six,section,spin,info,ambient=map(read,[sixpath,sectionpath,spinpath,invpath,ambientpath])
 oldvectors=ambient['ambient_beta16'][:3];assert len(oldvectors)==3 and all(len(v)==16 for v in oldvectors)
 assert all(isinstance(x,int)and 0<=x<32003 for v in oldvectors for x in v)
 Bsym,c,u,v,X,W=s.symbols('B c u v X W');weights=info['Cartan_weights'];split=[]
 for k,expr in enumerate(six['primitive_beta16']):
  z=s.cancel(s.sympify(expr))
  if z==0:continue
  power=s.Poly(s.denom(z),Bsym).degree();assert power in[0,1]
  assert weights[k][1]>=1 if power else weights[k][1]==-1
  split.append(dict(coordinate=k,denominator_B_power=power,second_Cartan_weight=weights[k][1]))
 D=147456*X**2-(3072*W-1)**2;vs=s.sympify(section['section_v'])
 pp=[s.Poly(s.cancel(Bsym*D*s.sympify(z).subs({c:1,v:vs},simultaneous=True)),u,X,W,Bsym)for z in six['primitive_beta16']]
 denominator=s.ilcm(1,1,*(s.denom(z)for p in pp for z in p.coeffs()));assert denominator==1080
 norms=[sum(abs(int(z*denominator))for z in p.coeffs())for p in pp]
 mon=spin['quadratic_monomials']
 def norm_matrix(key):
  vals=[int(sum(abs(4*s.Rational(z))*norms[mon[k][0]]*norms[mon[k][1]]for k,z in enumerate(row)))for row in spin[key]]
  return [vals[4*i:4*i+4]for i in range(4)]
 C,Dn=map(norm_matrix,['spin_adjoint_C_coefficients16x136','vector_adjoint_D_coefficients16x136'])
 def mm(A,B):return [[sum(A[i][k]*B[k][j]for k in range(4))for j in range(4)]for i in range(4)]
 def count(deg):
  cancel=deg//2;ans=0
  for aa in range(deg+1):
   for YY in range(2*deg+1):
    if not 0<=2*YY-aa<=3*deg:continue
    for ww in range(3*deg+1):
     for bb in range(deg-cancel+1):
      if 2*YY+ww+bb<=4*deg-cancel and 2*(aa+YY+ww+bb)<=9*deg-2*cancel:ans+=1
  return ans
 bounds=[]
 for n in[6,8]:
  numbers=[]
  for word in necklaces(n):
   M=[[int(i==j)for j in range(4)]for i in range(4)]
   for letter in word:M=mm(M,[C,Dn][letter])
   numbers.append(sum(M[i][i]for i in range(4)))
  d=2*n;maxnorm=max(numbers);terms=count(d);bits=maxnorm.bit_length();digits=len(str(maxnorm))
  bounds.append(dict(trace_word_length=n,necklace_count=len(numbers),beta_degree=d,
   guaranteed_B_factors_cancelled_at_least=d//2,remaining_B_denominator_power_at_most=d//2,
   remaining_H_denominator_power=d,source_polynomial_total_degree_after_parity_at_most=4*d,
   rigorous_monomial_bound_after_guaranteed_B_cancel=terms,integer_coefficient_bit_bound=bits,
   compact_JSON_output_byte_bound=terms*(digits+40)))
 size3=sum([3*bounds[0]['compact_JSON_output_byte_bound'],bounds[1]['compact_JSON_output_byte_bound']])
 size4=sum([4*bounds[0]['compact_JSON_output_byte_bound'],bounds[1]['compact_JSON_output_byte_bound']])
 report=dict(id=STEM,all_passed=True,preflight_only_no_rank_or_invariant_expansion_job=True,
  fixed_ambient_vectors_from_saved_probe=oldvectors,formal_plane='beta0+s*beta1+t*beta2 over F32003, coefficient rows only',
  formal_coefficient_row_counts={'degree12':91,'degree16':153,'degree24':325},
  intended_rank_gates=[{'degree12':9},{'product_degree16':19},{'degree16':20},{'degree24':67}],
  trace6_necklaces=necklaces(6),trace8_necklaces=necklaces(8),
  force_existing_N12_first_primitive12=True,new_low_source_output_count_if_success=4,
  source_B_weight_split=split,all_invariant_B_denominator_bound='d/2 by zero second Cartan weight',
  full70_weight24_candidate_counts=dict(spectral=35,primitive12_times_weight12_spectral=20,primitive12_symmetric_square=10,primitive16_times_Inv8=5),
  universal_rank67_required_before_completeness=True,
  source_low_word_output_bounds=bounds,total_new_low_output_upper_bytes_3plus1=size3,
  total_new_low_output_upper_bytes_if_N12_not_reused_4plus1=size4,
  future_degree24_source_column_rigorous_support_bound=count(24),
  common_future_source_denominator='B^12*H^24 after guaranteed B cancellations',
  fixed_hard_caps=dict(cpu_seconds=60,address_space_bytes=2_000_000_000,output_bytes=32*1024*1024),
  no_alternative_plane_or_growth_on_rank_failure=True,no_source_expansion_authorized_by_this_preflight=True,
  input_sha256={str(p.relative_to(ROOT)):sha(p)for p in[Path(__file__),sixpath,sectionpath,spinpath,invpath,ambientpath]},seconds=time.monotonic()-started)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0)
 print(json.dumps({k:report[k]for k in['all_passed','source_low_word_output_bounds','total_new_low_output_upper_bytes_3plus1','future_degree24_source_column_rigorous_support_bound','seconds']},indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output.resolve())
