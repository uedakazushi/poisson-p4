#!/usr/bin/env python3
"""Mixed Mi/G source and multiplication INPUT only; no norm/solve/GCD job."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_MIXED_NORM16_INPUT_20261104_v63'
PEAK=0
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def poly(p):
 z={}
 for m,c in p:
  assert len(m)==3 and m[2]==0 and min(m)>=0 and tuple(m[:2])not in z;q=Q(c);assert q;z[tuple(m[:2])]=q
 return z
def encode(p):return [[[a,e,0],str(c)]for(a,e),c in sorted(p.items(),reverse=True)]
def add(a,b,c=Q(1)):
 global PEAK
 z=dict(a)
 for m,v in b.items():
  z[m]=z.get(m,Q(0))+c*v
  if not z[m]:del z[m]
 PEAK=max(PEAK,len(z));assert len(z)<150000;return z
def shift(p,a=0,e=0):return {(i+a,j+e):v for(i,j),v in p.items()}
def scalar(p,c):return {m:v*c for m,v in p.items()if v*c}
def mul(a,b):
 global PEAK
 z={}
 for(i,j),c in a.items():
  for(k,l),d in b.items():z[i+k,j+l]=z.get((i+k,j+l),Q(0))+c*d
 z={m:v for m,v in z.items()if v};PEAK=max(PEAK,len(z));assert len(z)<150000;return z
def check_product(a,b,c):
 """Independent numeric-denominator-cleared integer product check."""
 def clear(p):
  d=math.lcm(*(v.denominator for v in p.values()));return {m:v.numerator*(d//v.denominator)for m,v in p.items()},d
 aa,da=clear(a);bb,db=clear(b);cc,dc=clear(c);D=math.lcm(da*db,dc);out={m:-v*(D//dc)for m,v in cc.items()};fac=D//(da*db)
 for(i,j),x in aa.items():
  for(k,l),y in bb.items():out[i+k,j+l]=out.get((i+k,j+l),0)+fac*x*y
 assert not any(out.values())
def divide_e_monic(p,H):
 r=dict(p);q={};steps=0
 while r and max(e for a,e in r)>=16:
  e=max(e for a,e in r);lead={(a,e-16):v for(a,j),v in r.items()if j==e};q=add(q,lead);r=add(r,mul(lead,H),Q(-1));steps+=1;assert steps<64
 return q,r,steps
def profile(p):return {'terms':len(p),'a_degree':max((a for a,e in p),default=-1),'e_degree':max((e for a,e in p),default=-1),'total_degree':max((a+e for a,e in p),default=-1)}
def assignment_degree(degrees):
 n=len(degrees);dp=[-10**6]*(1<<n);dp[0]=0
 for mask in range(1<<n):
  i=mask.bit_count()
  if i==n or dp[mask]<0:continue
  for j in range(n):
   if not(mask>>j&1)and degrees[i][j]>=0:dp[mask|1<<j]=max(dp[mask|1<<j],dp[mask]+degrees[i][j])
 return dp[-1]
def main(out):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.alarm(60);out=out.resolve();out.mkdir(parents=True,exist_ok=False)
 pb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_NORM_PENCIL_20261104_v63';rb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_PROJECTION_CURVE_RESIDUAL_INPUT_20261104_v63';ib=ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_NORM16_INPUT_20261104_v63';review=ROOT/'reports/S61_NONSQUARE_RAW_H_B_ADDITIONAL_NORM_PENCIL_20261104_v63_REVIEW'
 for b in [pb,rb,ib]:
  r=read(b/'report.json');filename='native-certificate.json'if b==pb else'certificate.json';assert sha(b/filename)==r['evidence_sha256'][filename]
 p=read(pb/'native-certificate.json');r=read(rb/'certificate.json');old=read(ib/'certificate.json');G=poly(p['final_common_curve_gcd']);assert G==poly(r['common_G']);assert poly(p['original_G_quotient_by_final'])=={(0,0):Q(1)}
 mixed=[poly(x['remaining_h_coprime_U'])for x in p['mixed3_records']];mi=[poly(x)for x in p['mixed3_quotients_by_final']];assert len(mi)==3
 for M,m in zip(mixed,mi):check_product(G,m,M)
 H=poly(old['h1_e16_monic']);assert {(a,e):v for(a,e),v in H.items()if e==16}=={(0,16):Q(1)}
 profiles=[]
 for m in mi:
  degree=max(e for a,e in m);leading={a:v for(a,e),v in m.items()if e==degree};assert set(leading)=={0};profiles.append(dict(profile(m),e_leading_QQunit=str(leading[0])))
 chosen=min(range(3),key=lambda j:(len(mi[j]),profiles[j]['e_degree'],profiles[j]['a_degree']));assert chosen==0
 leading=Q(profiles[chosen]['e_leading_QQunit']);K=scalar(mi[chosen],1/leading);assert profile(K)['e_degree']==19
 qs=[];rs=[];steps=[];M=[[{}for _ in range(16)]for _ in range(16)]
 for j in range(16):
  target=shift(K,e=j);q,rem,nsteps=divide_e_monic(target,H);assert add(mul(H,q),rem)==target;qs.append(q);rs.append(rem);steps.append(nsteps)
  for(a,e),v in rem.items():assert e<16;M[e][j][a,0]=v
  rebuilt={}
  for i in range(16):rebuilt=add(rebuilt,shift(M[i][j],e=i))
  assert rebuilt==rem
 terms=sum(len(z)for row in M for z in row);assert terms<150000;allp=[z for row in M for z in row];D=math.lcm(*(v.denominator for z in allp for v in z.values()));mx=max(abs(int(v*D))for z in allp for v in z.values());s=max(map(len,allp));height=math.factorial(16)*(s*mx)**16;degrees=[[profile(z)['a_degree']for z in row]for row in M];detdeg=assignment_degree(degrees);assert detdeg>=0
 future_det_bytes=(detdeg+1)*(100+(height.bit_length()*30103+99999)//100000)
 cofactor_deg=15*max(max(row)for row in degrees);cofactor_bits=(math.factorial(15)*(s*mx)**15).bit_length();future16cofactor_bytes=16*(cofactor_deg+1)*(100+(cofactor_bits*30103+99999)//100000)
 cert=dict(source_variable_order=['a','e','d_unused'],G=encode(G),mixed3_stripped=[encode(m)for m in mixed],all3_quotients_by_G=[encode(m)for m in mi],all3_exact_source_products_checked=True,quotient_profiles=profiles,chosen_mixed_lambda_coefficient_index=chosen+1,chosen_mi_QQ_leading_unit=str(leading),chosen_mi_monic=encode(K),H_h1_monic16=encode(H),matrix_rows_source_order=[[encode(z)for z in row]for row in M],all16_column_quotients_source_order=[encode(q)for q in qs],all16_column_remainders_source_order=[encode(q)for q in rs],all16_source_divisions_checked=True,necessity_scope='ForactualBbad onD(G) only: h=G*h1,Kold=G*k1,Mi=G*mi andstrippingunitsareFULLUunits implyh1=k1=all3mi=0. G=0 isretained inexistingQ/Lcasehandling; G isNOTinvertedglobally.',algebra_model='QQ[a][e]/(h1monic16), multiplicationbym1monic19. Allentriespolynomialina; onlyQQleadingunit normalized.',norm_nonzero_not_yet_proved=True,no_norm_determinant_solve_univariate_gcd_or_point_job=True)
 (out/'certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
 costs=dict(matrix_shape=[16,16],matrix_terms=terms,max_entry_terms=s,matrix_a_degrees=degrees,matrix_uniform_QQ_denominator=str(D),cleared_integer_max_coefficient_bits=mx.bit_length(),determinant_a_degree_assignment_upper=detdeg,determinant_support_upper=detdeg+1,determinant_integer_coefficient_bits_upper=height.bit_length(),determinant_serialized_bytes_upper=future_det_bytes,future16cofactor_degree_uniform_upper=cofactor_deg,future16cofactor_bits_upper=cofactor_bits,future16cofactor_bytes_upper=future16cofactor_bytes,future_Bezout_total_artifact_size_NOT_bounded_by_det_only=True,previous_34MB_guard_failure_requires_complete_future_storage_preflight=True,no_new_determinant_or_solve_allocation=True,division_steps=steps,maximum_actual_input_arithmetic_support=PEAK)
 (out/'cost-preflight.json').write_text(json.dumps(costs,indent=2)+'\n')
 paths=[Path(__file__),ROOT/f'research/notes/{ID}.md',pb/'report.json',pb/'native-certificate.json',rb/'report.json',rb/'certificate.json',ib/'report.json',ib/'certificate.json',review/'report.json',ROOT/'reports/S61_NONSQUARE_RAW_H_B_OLD_RESIDUAL_LENGTH_20261104_v63_REVIEW/report.json']
 report=dict(id=ID,all_passed=True,all3_exact_Mi_Gquotient_source_products=True,cheapest_m1_e19_input=True,all16_monic_divisions_input_only=True,all_a_coefficients_polynomial_no_unknown_LC_inversion=True,necessary_only_outsideG0=True,G0_cases_retained=True,norm_nonzero_and_annihilator_gcd_not_computed=True,determinant_solve_GB_point_prime_jobs=0,matrix_terms=terms,max_entry_terms=s,input_certificate_bytes=(out/'certificate.json').stat().st_size,input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(z.stat().st_size for z in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(costs,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
