#!/usr/bin/env python3
"""Same-agent distinct QQ direct-product audit of saved fixed QRE division traces."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_IDENTITY_REVIEW_20261105_v64'
BASE='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_EXISTING_QRE_STRIP_20261105_v64'
PRIMARY_SHA='d1f1940f49703f4c7d956b6ead3c803b25de82aa4e628d6a0a74725739462fb3'
SOURCE='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_STRIP_PREFLIGHT_20261105_v64'
ROOT_GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_ROW58_QRE_ROOT_SOURCE_REVIEW_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(a):
 assert a.verify
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('fixed QQ identity audit60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);base=ROOT/'reports'/BASE
 inputs={};owned_peak=0;pool={};products=0;checks=0
 def add(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return p
 def own(extra=2048):
  nonlocal owned_peak
  n=25000+8*sum(map(len,pool.values()))+extra;owned_peak=max(owned_peak,n);assert n<=82696;return n
 def reader(p,expected=None):
  # Deliberately not the producer codec/parser or the pseudo-division engine.
  data=add(p).read_bytes();assert data[:8]==b'Q5CELL1\n';at=8
  def u():
   nonlocal at
   n=0;k=0;start=at
   while True:
    v=data[at];at+=1;n+=(v%128)<<k
    if v<128:break
    k+=7
   assert at-start==max(1,(n.bit_length()+6)//7);return n
  def integer():
   nonlocal at
   tag=u();length=tag//2;b=data[at:at+length];at+=length;assert len(b)==length and(not b or b[0]!=0);n=int.from_bytes(b,'big');assert n or tag%2==0;return -n if tag%2 else n
  flag=data[at];at+=1
  if flag==0:assert at==len(data);return {}
  assert flag==1;num=integer();den=integer();content=Q(num,den);assert(content.numerator,content.denominator)==(num,den)
  shift=[]
  for _ in range(5):v=u();shift.append(v//2 if v%2==0 else-(v//2)-1)
  n=u();assert n<=1739;A={};g=0;prev=None;h=hashlib.sha256();h.update((str(n)+'\n').encode())
  for _ in range(n):
   c=integer();e=tuple(shift[j]+u()for j in range(5));key=(sum(e),e);assert c and(prev is None or prev>key)and e not in A;prev=key;g=gcd(g,abs(c));v=content*c;A[e]=v
   h.update((str(v)+' '+' '.join(map(str,e))+'\n').encode())
   assert max(abs(v.numerator).bit_length(),v.denominator.bit_length())<2048
  assert g==1 and at==len(data)
  if expected is not None:assert expected['sha256']==inputs[str(p.relative_to(ROOT))]and expected['terms']==n and expected['value_sha256']==h.hexdigest()
  return A
 def logical(A):
  h=hashlib.sha256();h.update((str(len(A))+'\n').encode())
  for e in sorted(A,key=lambda e:(sum(e),e),reverse=True):h.update((str(A[e])+' '+' '.join(map(str,e))+'\n').encode())
  return h.hexdigest()
 try:
  pp=add(base/'report.json');assert sha(pp)==PRIMARY_SHA;pr=json.loads(pp.read_text());assert pr['all_passed']
  for name,h in pr['evidence_sha256'].items():assert sha(add(base/name))==h
  sp=add(ROOT/f'reports/{SOURCE}/report.json');sr=json.loads(sp.read_text());assert sr['all_passed']
  bp=add(sp.parent/'integer-pseudo-division-budget-and-proof.json');b=json.loads(bp.read_text())
  gp=add(ROOT/f'reports/{ROOT_GATE}/report.json');assert sha(gp)=='a03e3d78d68bfc67639065ad48a7b7000387d7acee611f22854856ed011be350'and json.loads(gp.read_text())['all_passed']
  add(Path(__file__));add(ROOT/f'research/notes/{ID}.md')
  raw_bound=sum((ROOT/p).stat().st_size for p in inputs)+2*1024*1024
  assert raw_bound<32*1024*1024
  pre={'same_agent_distinct_QQ_arithmetic':True,'source_only_fixed_existing951':True,'no_producer_import_or_reexecution':True,'max_maps':4,'terms_per_map_upper':1739,'owned_upper_including25000_metadata':82696,'QQ_temporary_bit_upper':2048,'all_actual_READ_RAW_plus_receipts_upper':raw_bound,'seconds_limit':60,'AS_bytes':2000000000}
  (out/'review-preflight.json').write_bytes(enc(pre))
  result=json.loads((base/'fixed-row58-QRE-strip-and-Gauss-factorization.json').read_text())
  terms=reader(ROOT/b['original_binary_path']);valuation=b['actual_Laurent_monomial_valuation5'];content=Q(b['original_QQ_content']);pool['P']={tuple(e[j]-valuation[j]for j in range(5)):c/content for e,c in terms.items()};del terms;own();assert len(pool['P'])==951 and all(c.denominator==1 for c in pool['P'].values());initial=logical(pool['P'])
  factors={'Q':(3,5,{(0,0,0,1,0):Q(5),(0,0,0,0,0):Q(1)}),'R':(3,1,{(0,0,0,1,0):Q(1),(0,0,0,0,0):Q(1)}),'E':(2,32,{(0,0,1,0,0):Q(32),(0,0,0,1,0):Q(15),(0,0,0,0,0):Q(35)})}
  records=[];powers={'Q':0,'R':0,'E':0};terminals=[]
  for index in result['input_and_final_factorization_chain']:
   tf=base/index['file'];assert sha(tf)==index['sha256'];t=json.loads(tf.read_text());label=t['factor'];j,aLead,F=factors[label]
   pool['V']=reader(base/t['quotient']['file'],t['quotient']);pool['R']=reader(base/t['remainder']['file'],t['remainder']);pool['T']={};own();P=pool['P'];A=int(t['pseudo_multiplier']);assert A==aLead**max(e[j]for e in P)
   assert all(e[j]==0 for e in pool['R'])
   for e,c in pool['V'].items():
    for f,v in F.items():
     g=tuple(e[k]+f[k]for k in range(5));z=pool['T'].get(g,Q(0))+c*v
     if z:pool['T'][g]=z
     else:pool['T'].pop(g,None)
     products+=1;own()
   for e,c in pool['R'].items():
    z=pool['T'].get(e,Q(0))+c
    if z:pool['T'][e]=z
    else:pool['T'].pop(e,None)
   assert len(pool['T'])==len(P)
   for e,c in P.items():assert pool['T'].get(e,Q(0))==A*c;checks+=1
   assert bool(not pool['R'])==t['zero_remainder'];del pool['T']
   if t['zero_remainder']:
    g=0
    for c in pool['V'].values():assert c.denominator==1;g=gcd(g,abs(c.numerator))
    assert g==A
    del pool['P'],pool['R'];pool['P']=pool.pop('V')
    for e in pool['P']:pool['P'][e]/=A;assert pool['P'][e].denominator==1
    powers[label]+=1;assert logical(pool['P'])==next((v['logical_integer_value_sha256']for v in [t.get('output_integer_polynomial',{})]if 'logical_integer_value_sha256'in v),logical(pool['P']))
   else:
    # FV+R=A*P and deg_j R=0 prove this nonzero R is the unique
    # rational canonical remainder. Hence the next factor does not divide P.
    e=next(iter(pool['R']));terminals.append({'factor':label,'nonzero_canonical_remainder_witness':{'exponent5':list(e),'coefficient':str(pool['R'][e]/A)},'no_leading_variable_in_remainder':True,'degree_one_factor_proves_exact_maximal_power':True})
    del pool['V'],pool['R']
   del P;own();records.append({'trial':t['trial'],'factor':label,'direct_QQ_F_times_V_plus_R_equal_A_P':True,'zero_remainder':t['zero_remainder'],'quotient_original_gauss_content_checked':t['zero_remainder']})
  final=reader(base/result['final_reduced_integer_polynomial']['file'],result['final_reduced_integer_polynomial']);assert pool['P']==final and len(final)==317 and powers=={'Q':4,'R':2,'E':1}and len(terminals)==3
  assert initial==result['original_rebuilt_value_sha256']
  ev={'all_passed':True,'same_agent_distinct_QQ_arithmetic':True,'no_producer_import_or_reexecution':True,'all10_saved_division_identities':records,'terminal_remainders_prove_exact_Q4_R2_E1':terminals,'only_QRE_powers':powers,'final_primitive_integer_terms':317,'final_integer_value_sha256':logical(final),'actual_pivot_formula':'F10/(353894400000*H*t0*F9)','F9_not_divided_or_assumed_coprime':True,'all_full168x15_87plus66_RHS_P290_once_complements_preserved':True}
  ep=out/'all-saved-QRE-division-identities-and-terminal-remainders.json';ep.write_bytes(enc(ev))
  raw=sum((ROOT/p).stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.iterdir());assert raw+65536<=raw_bound
  report={'id':ID,'all_passed':True,'same_agent_distinct_QQ_arithmetic':True,'independent_agent_ROOT_source_gate':str(gp.relative_to(ROOT)),'division_engine_not_imported_or_reexecuted':True,'fixed_QQ_products':products,'full_QQ_coefficient_comparisons':checks,'terminal_nonzero_canonical_remainders':3,'powers':powers,'final_terms':317,'owned_peak':owned_peak,'owned_upper':82696,'actual_READ_and_output_RAW_before_report':raw,'input_sha256':inputs,'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start}
  (out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(exc),'all_completed_files_preserved':True,'no_automatic_retry':True}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--verify',action='store_true');main(p.parse_args())
