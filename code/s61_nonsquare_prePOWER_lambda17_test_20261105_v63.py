#!/usr/bin/env python3
"""Prepare, then ONE fixed17 scalar polynomial-ratio test; all-coefficient check."""
import argparse,gc,hashlib,json,math,resource,signal,struct,time
from fractions import Fraction as Q
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];ID='S61_NONSQUARE_PREPOWER_TARGET17_LAMBDA17_TEST_20261105_v63'
CAP=36714;BITS=2048;RAWCAP=6820000
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def pack(m):
 assert len(m)==3 and all(0<=a<1024 for a in m);return (m[0]<<20)|(m[1]<<10)|m[2]
def unpack(k):return [k>>20,(k>>10)&1023,k&1023]
def prepare(out):
 out.mkdir(parents=True,exist_ok=False);pins={}
 def load(st,name='report.json'):
  rel=f'reports/{st}/{name}';p=ROOT/rel;pins[rel]=sha(p);return read(p)
 old='S61_NONSQUARE_PREPOWER_TARGET17_ANCHOR_QUOTIENT_SOURCE_20261105_v63';s=load(old,'finite-lambda-ratio-gate.json');assert s['QQscalar']=='-1/3'
 hg=load('S61_NONSQUARE_PREPOWER_TARGET17_LAMBDA17_COMMON_DEN_SOURCE_20261105_v63');assert hg['all_passed']and hg['all_explicit_integer_temporaries_proved_under2048']
 peer='S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63_REVIEW';pr=load(peer);assert pr['all_passed'];assert pins[f'reports/{peer}/report.json']=='20f835365d200da04b01e14265937a248a2fedde7b728e26e614e7b7496b3b1b'
 kp=load('S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63','root237.json')['polynomial'];yp=load('S61_NONSQUARE_PREPOWER_TARGET17_THREE_CRAMER_ENTRY_20261105_v63','Cramer2.json')
 DK=int(s['K_common_D']);GK=int(s['K_integer_content']);DY=int(s['Y_common_D']);GY=int(s['Y_integer_content']);K=[];H=[]
 for m,c in kp:
  v=Q(c)*DK;assert v.denominator==1 and int(v)%GK==0;z=int(v)//GK;assert Q(GK*z,DK)==Q(c);K.append([m,str(z)])
 for m,c in yp['integer_coefficients']:
  v=int(c);assert v%GY==0;z=v//GY;assert GY*z==v;H.append([m,str(z)])
 assert len(K)==1174 and len(H)==1768 and -Q(GY*DK,DY*GK)==Q(-1,3)
 assert max(K,key=lambda v:(sum(v[0]),v[0]))==[s['leading_K_monomial'],s['leading_K_integer']]
 data=dict(K_primitive_integer=K,H_primitive_integer=H,candidate_monomials=s['candidate_monomials'],K_LM=s['leading_K_monomial'],K_LC=s['leading_K_integer'],
  common_D=s['common_denominator_upper'],all_R_bound=s['numerator_upper'],all_temporary_bound=s['integer_verification_L1_upper'],actual_lambda_scalar='-1/3',
  source_primitive_units=dict(DK=DK,GK=GK,DY=DY,GY=GY),ratio_unknowns_NOT_computed=True)
 dump(out/'literal-input.json',data)
 code=Path(__file__);pins[str(code.relative_to(ROOT))]=sha(code)
 source=dict(height_completed_bound_bits=1746,implementation_STOP_bits=2048,no_R_times_R=True,
  numeric_ownership=dict(input_sparse_keys_coefficients=5884,verification_keys_coefficients=15540,R_and_quotient_buffer=34,ALL_static_numeric_metadata=15000,extra_scalar_cursor_key_buffers=256,TOTAL=36714),
  preparation_ownership_upper=32908,readback_ownership_upper=16000,
  prefix=dict(version='LAM17P01',all_active_state='inputpath/hash,phase/j/k/current_numerator,commonD,LC,Rprefix,residualmap; K/H regenerated only from immutable literal input',binary_RAW_bound=2097152,inspect_command='python THIS_CODE inspect --state PREFIX'),
  no_candidate_enlargement=True,no_polyGCD_factor_point_prime_or_other_target=True,
  root_authorization='conditional ONE fixed17 common-leading-denominator test, preserve prefix on2048/36714/60s/2GB/6.82MB stop; no rerun',
  source_normalization_all2942_coefficients_exact=True)
 dump(out/'source-and-custody-proof.json',source)
 gate=dict(id=ID+'_WORKER_GATE',all_passed=True,frozen=True,mathematical_test_NOT_run=True,input_sha256=pins,literal_input_sha256=sha(out/'literal-input.json'),
  required_three_Cramer_peer_hash=pins[f'reports/{peer}/report.json'],candidate_count=17,verification_positions=7770,owned_cap=CAP,integer_guard_bits=BITS,proved_max_integer_bits=1746,TOTAL_RAW_cap=RAWCAP,
  preparation_owned_upper=32908,ONE_test_only=True,evidence_sha256={p.name:sha(p)for p in out.iterdir()})
 dump(out/'report.json',gate);print(json.dumps({k:v for k,v in gate.items()if k not in['input_sha256','evidence_sha256']},indent=2))
class StopGuard(Exception):pass
def execute(gatepath,out):
 start=time.monotonic();signal.alarm(59);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));g=read(gatepath);assert g['all_passed']and g['mathematical_test_NOT_run'];pins={str(gatepath.relative_to(ROOT)):sha(gatepath)}
 for rel,h in g['input_sha256'].items():assert sha(ROOT/rel)==h;pins[rel]=h
 inp=gatepath.parent/'literal-input.json';assert sha(inp)==g['literal_input_sha256'];pins[str(inp.relative_to(ROOT))]=sha(inp)
 out.mkdir(parents=True,exist_ok=False);x=read(inp);M=x['candidate_monomials'];lm=x['K_LM'];lc=int(x['K_LC']);D=int(x['common_D']);RB=int(x['all_R_bound']);VB=int(x['all_temporary_bound'])
 K={pack(m):int(c)for m,c in x.pop('K_primitive_integer')};H={pack(m):int(c)for m,c in x.pop('H_primitive_integer')};gc.collect();assert len(K)==1174 and len(H)==1768 and len(M)==17
 R=[];res={};phase='triangular';j=k=0;cur=0;products=divisions=0;maxbits=peak=0;inputref=str(inp.relative_to(ROOT))
 def guard():
  nonlocal peak
  owned=2*(len(K)+len(H)+len(res))+2*len(R)+15000+256;peak=max(peak,owned)
  if owned>CAP:raise StopGuard('owned_before_commit')
  if time.monotonic()-start>48:raise StopGuard('time_before_commit')
 def check(z):
  nonlocal maxbits
  maxbits=max(maxbits,abs(z).bit_length())
  if abs(z).bit_length()>BITS:raise StopGuard('integer_bits_before_commit')
  return z
 def mul(a,b):
  nonlocal products
  if abs(a).bit_length()+abs(b).bit_length()>BITS:raise StopGuard('multiply_upper_before_operation')
  products+=1;return check(a*b)
 def add(a,b):
  if max(abs(a).bit_length(),abs(b).bit_length())+1>BITS:raise StopGuard('addition_upper_before_operation')
  return check(a+b)
 def prefix(reason):
  p=out/'prefix-LAM17P01.bin';meta=dict(version='LAM17P01',reason=reason,input_reference=inputref,input_sha256=sha(inp),phase=phase,j=j,k=k,next_operation_NOT_committed=True,products=products,divisions=divisions)
  b=json.dumps(meta).encode();assert len(b)<8192
  with p.open('wb')as f:
   f.write(b'LAM17P01');f.write(struct.pack('<I',len(b)));f.write(b)
   def z(v):
    b=abs(v).to_bytes((abs(v).bit_length()+7)//8,'big');assert len(b)<=256;f.write(struct.pack('<BH',v<0,len(b)));f.write(b)
   z(D);z(lc);z(cur);f.write(struct.pack('<I',len(R)))
   for v in R:z(v)
   f.write(struct.pack('<I',len(res)))
   for key,v in res.items():f.write(struct.pack('<I',key));z(v)
  assert p.stat().st_size<=2097152
  return p
 trace=(out/'triangular-exact-divisions.jsonl').open('w')
 try:
  for j,m in enumerate(M):
   guard();target=[a+b for a,b in zip(m,lm)];cur=mul(D,H.get(pack(target),0))
   for k in range(j):
    delta=[target[t]-M[k][t]for t in range(3)];kv=K.get(pack(delta),0)if all(v>=0for v in delta)else 0
    cur=add(cur,-mul(R[k],kv))
   assert cur%lc==0;q=check(cur//lc);assert abs(q)<=RB and mul(q,lc)==cur;divisions+=1;R.append(q);guard()
   trace.write(json.dumps(dict(step=j,monomial=m,predivision_numerator=str(cur),divisor=str(lc),integer_R=str(q),remainder_ZERO=True))+'\n');trace.flush()
  phase='verify_initialize'
  for k,(key,z)in enumerate(H.items()):guard();v=mul(D,z);res[key]=v;guard()
  phase='verify_ALL_coefficients'
  for j,q in enumerate(R):
   if not q:continue
   for k,(key,z)in enumerate(K.items()):
    guard();a=unpack(key);m=[a[t]+M[j][t]for t in range(3)];assert sum(m)<=34;newkey=pack(m);v=add(res.get(newkey,0),-mul(q,z));assert abs(v)<=VB
    if v:res[newkey]=v
    else:res.pop(newkey,None)
    guard()
  phase='canonical_readback';poly=[]
  for j,q in enumerate(R):
   if not q:continue
   den=check(3*D);gg=math.gcd(abs(q),den);num=check(-q//gg);dd=check(den//gg);assert mul(num,gg)==-q and mul(dd,gg)==den and math.gcd(abs(num),dd)==1
   poly.append([M[j],str(num)if dd==1 else f'{num}/{dd}'])
  success=not res;candidate=dict(actual_lambda_polynomial=success,candidate_only_if_not_success=not success,QQ_polynomial=poly,common_R_denominator=str(D),R=[str(v)for v in R],exact_identity='D*Hprimitive-Kprimitive*Rpoly=0'if success else 'nonzero verification remainder; original rational lambda retained')
  if res:
   key,z=next(iter(res.items()));candidate['countermonomial']=unpack(key);candidate['integer_countercoefficient']=str(z)
  dump(out/'polynomial-ratio-result.json',candidate);back=read(out/'polynomial-ratio-result.json');assert back==candidate
  trace.close();new=sum(p.stat().st_size for p in out.iterdir());external=sum((ROOT/p).stat().st_size for p in pins);assert new+external<=RAWCAP
  rep=dict(id=ID,all_passed=True,frozen=True,ONE_fixed17_scalar_ratio_test=True,actual_lambda_is_polynomial=success,nonzero_polynomial_terms=len(poly),
   scalar_exact_divisions=divisions,all7770_coefficient_positions_verified=True,full_identity_ZERO=success,
   no_candidate_enlargement_polyGCD_factor_point_prime_or_lastrow_value=True,original_rational_lambda_retained_if_failure=True,
   peak_owned_numeric_metadata=peak,observed_integer_bits=maxbits,integer_products=products,new_RAW=new,external_read_RAW=external,TOTAL_RAW=new+external,
   independent_ratio_peer_pending=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in out.iterdir()},seconds=time.monotonic()-start)
  dump(out/'report.json',rep);print(json.dumps({k:v for k,v in rep.items()if k not in['input_sha256','evidence_sha256']},indent=2))
 except StopGuard as e:
  trace.close();p=prefix(str(e));dump(out/'report.json',dict(id=ID,all_passed=False,frozen=True,guarded_prefix_only=True,no_ratio_or_lastrow_result_claim=True,reason=str(e),prefix_sha256=sha(p),restore_parser='THIS_CODE inspect --state prefix-LAM17P01.bin; no auto resume',input_sha256=pins))
 finally:signal.alarm(0)
def inspect(p):
 with p.open('rb')as f:
  assert f.read(8)==b'LAM17P01';n=struct.unpack('<I',f.read(4))[0];assert n<8192;meta=json.loads(f.read(n))
  def z():
   s,n=struct.unpack('<BH',f.read(3));assert s<2 and n<=256;b=f.read(n);assert len(b)==n;return (-1 if s else 1)*int.from_bytes(b,'big')
  for _ in range(3):z()
  nr=struct.unpack('<I',f.read(4))[0];assert nr<=17
  for _ in range(nr):z()
  nc=struct.unpack('<I',f.read(4))[0];assert nc<=7770
  for _ in range(nc):f.read(4);z()
  assert not f.read(1)
 print(json.dumps(dict(metadata=meta,Rprefix=nr,residual_coefficients=nc,parse_complete=True,no_resume_executed=True),indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='mode',required=True);a=s.add_parser('prepare');a.add_argument('--output',type=Path,required=True);a=s.add_parser('execute');a.add_argument('--gate',type=Path,required=True);a.add_argument('--output',type=Path,required=True);a=s.add_parser('inspect');a.add_argument('--state',type=Path,required=True);a=p.parse_args()
 if a.mode=='prepare':prepare(a.output.resolve())
 elif a.mode=='execute':execute(a.gate.resolve(),a.output.resolve())
 else:inspect(a.state.resolve())
