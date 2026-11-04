#!/usr/bin/env python3
"""ONE U coefficients only; source plan for H4 saturation and Bminor alternative."""
import argparse,hashlib,json,math,resource,signal,subprocess,sys,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_RAW_H_B_K13_U_INPUT_20261104_v63';sys.set_int_max_str_digits(0)
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ip(p):return {int(a):int(v)for a,v in p}
def add(a,b,c=1):
 z=dict(a)
 for a,v in b.items():
  z[a]=z.get(a,0)+c*v
  if not z[a]:del z[a]
 return z
def mul(a,b):
 z={}
 for a,c in a.items():
  for b,d in b.items():z[a+b]=z.get(a+b,0)+c*d
 return{k:v for k,v in z.items()if v}
def powers(p,n):
 out=[{0:1}]
 for j in range(n):out.append(mul(out[-1],p))
 return out
def main(out,work):
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda *_:(_ for _ in()).throw(TimeoutError('60s total cap')));signal.alarm(60);out=out.resolve();work=work.resolve();out.mkdir(parents=True,exist_ok=False);work.mkdir(parents=True,exist_ok=False)
 hb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_H4_INPUT_20261104_v63';bb=ROOT/'reports/S61_NONSQUARE_RAW_H_B_35_MINORS_20261104_v63';bp=ROOT/'reports/S61_NONSQUARE_RAW_H_B_K13_MIXED23_SCALARS_20261104_v63';f=read(hb/'H4-monic-field-coefficients.json');g={a:int(v)for a,v in f['g13']};tn=ip(f['e_numerator']);td=int(f['e_denominator']);U=read(bb/'certificate.json')['FULLU_product'];DU=math.lcm(*(Q(v).denominator for m,v in U));S=sum(abs(int(Q(v)*DU))for m,v in U);TM=max(td,sum(abs(v)for v in tn.values()));growth=abs(g[13])+sum(abs(v)for v in g.values());N=S*TM**7*growth**77;D=DU*td**7*abs(g[13])**77;bound=52*(150+((N.bit_length()+D.bit_length())*30103+99999)//100000);assert bound<1024*1024
 (out/'preflight.json').write_text(json.dumps(dict(source49terms=True,raw_a_degree89=89,pseudo_reduction_steps77=77,numerator_bits_upper=N.bit_length(),denominator_bits_upper=D.bit_length(),output52QQ_RAW_bytes_upper=bound,no_H4gcd_or_fullmatrix_authorized=True),indent=2)+'\n')
 inp=out/'native-input.txt'
 with inp.open('w')as s:
  for p in[g,tn]:
   s.write(str(len(p))+'\n')
   for a,v in sorted(p.items(),reverse=True):s.write(str(a)+' '+str(v)+'\n')
  s.write(str(td)+'\n'+str(DU)+'\n49\n')
  for m,v in U:s.write(' '.join(map(str,m))+' '+str(int(Q(v)*DU))+'\n')
 cpp=ROOT/'code/s61_nonsquare_raw_H_B_K13_U_input_20261104_v63.cpp';exe=work/'U-input'
 with(out/'compile.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nc++ -O2 -std=c++17 "$1" -o "$2" $(pkg-config --cflags --libs flint)','compile',str(cpp),str(exe)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 with(out/'native.log').open('w')as log:subprocess.run(['bash','-lc','source /workspace/cloud-setup/activate.sh\nexec "$@"','U-input',str(exe),str(inp),str(out/'native-certificate.json')],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=max(1,58-(time.monotonic()-start)))
 n=read(out/'native-certificate.json');tp=powers(tn,7);raw=[{}for _ in range(4)];common=DU*td**7;assert common==int(n['common_source_denominator'])
 for m,v in U:raw[m[2]]=add(raw[m[2]],{a+m[0]:c for a,c in tp[m[1]].items()},int(Q(v)*DU)*td**(7-m[1]))
 fields=[];num_bits=den_bits=terms=0
 for i,r in enumerate(n['U_coefficients']):
  assert raw[i]==ip(r['raw_common_numerator']);R=ip(r['remainder']);lp=r['L_power'];assert {a:v*g[13]**lp for a,v in raw[i].items()}==add(mul(g,ip(r['g13_quotient'])),R);num=ip(r['numerator']);den=int(r['denominator']);assert {a:v*common*g[13]**lp for a,v in num.items()}=={a:v*den for a,v in R.items()};fields.append([[a,str(Q(v,den))]for a,v in sorted(num.items(),reverse=True)]);terms+=len(num);num_bits=max(num_bits,max((abs(v).bit_length()for v in num.values()),default=0));den_bits=max(den_bits,den.bit_length())
 (out/'U-field-coefficients.json').write_text(json.dumps(dict(g13=f['g13'],U_coefficients_d0_to3=fields,all_U_zero=n['all_U_coefficients_zero'],actual_QQ_terms=terms,max_numerator_bits=num_bits,max_denominator_bits=den_bits),indent=2)+'\n')
 # Plan only; preserve multiplicities and explicit coefficient-unit gates.
 plan=dict(H4_source=str((hb/'H4-monic-field-coefficients.json').relative_to(ROOT)),U_source='U-field-coefficients.json',next_target='hU=H4/gcd_d(H4,U^4) onfiniteetaleKa13branches; no gcdcomputedhere',Ka13_etale_gate='Bindgcd(g13,g13prime)=1 byexplicitQQBezoutoracceptedfactorreceipt. Noirreducibilityguess.',coefficient_inverse_gate='For anyneeded p(a), literalxgcd(p,g13)=1 beforeinverse; propergcdretain/splitcomponents, notdrop.',multiplicities='EachfieldbranchH4length4, U^4 killsallU-supportedprimaryfactors; quotientkeepsoriginalmultiplicitiesofU-nonzero factors.',if_etale_gate_unavailable='QQfinitecoverdim52 givesU^52 instead ofU^4; no suchpowercomputedhere.',max_d_degree_of_U4=12,Euclid_polynomial_degrees_only=[4,12],live_field_coefficients_upper=3*13*13,no_d_Euclid_no12x28_solve_no_newnorm=True)
 (out/'H4-saturation-source-plan.json').write_text(json.dumps(plan,indent=2)+'\n')
 # Read-only one-Bminor cost inventory, no determinant or reduction.
 actualB=read(bp/'native-certificate.json.B21.json')['entries'];rows=[]
 for triple in __import__('itertools').combinations(range(7),3):
  score=sum(sum(len(z['numerator'])for z in actualB[3*r+c])for r in triple for c in range(3));rows.append((score,triple))
 score,triple=min(rows);alt=dict(cheapest_row_triple_by_sum_input_support=list(triple),selected9entry_QQ_terms=score,raw_minor_d_degree_upper=9,raw_minor_a_degree_upper=36,no_minor_or_norm4_computed=True,source='acceptedB21canonicalcoeffs, not22expandedminors',next_required='FirsthU/coprimalityscope; thenactualminorcoefficient/powerheightpreflight beforeoneNorm4test.')
 (out/'one-Bminor-alternative-preflight.json').write_text(json.dumps(alt,indent=2)+'\n')
 paths=[Path(__file__),cpp,ROOT/f'research/notes/{ID}.md',hb/'report.json',hb/'H4-monic-field-coefficients.json',bb/'report.json',bb/'certificate.json',bp/'report.json',bp/'native-certificate.json.B21.json']
 report=dict(id=ID,all_passed=True,all49source_terms_and4_reductions_verified=True,U_all_zero=n['all_U_coefficients_zero'],U_terms=terms,U_numerator_bits=num_bits,U_denominator_bits=den_bits,no_H4gcd_minor_norm4_or_Gaussian=True,no_eEuclid_newnorm_point_prime=True,native_seconds=n['native_seconds'],input_sha256={str(p.relative_to(ROOT)):sha(p)for p in paths},evidence_sha256={p.name:sha(p)for p in out.iterdir()if p.is_file()},seconds=time.monotonic()-start)
 assert sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));print(json.dumps(alt,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
 try:main(a.output,a.work)
 except Exception as e:
  if a.output.exists():(a.output/'failure.json').write_text(json.dumps({'all_passed':False,'error':str(e),'no_auto_retry':True},indent=2)+'\n')
  raise
