#!/usr/bin/env python3
"""Independent QQ[a,a^-1] axis peer of ONE T1 and its four columns."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text())
def numeric(x):
 if type(x)is int:return 1
 if isinstance(x,str):
  try:Q(x);return 1
  except ValueError:return 0
 if isinstance(x,list):return sum(map(numeric,x))
 if isinstance(x,dict):return sum(map(numeric,x.values()))
 return 0
def check(p):
 assert all(q and abs(q.numerator).bit_length()<=359 and q.denominator.bit_length()<=395 for q in p.values())
 return p
def add(p,q):
 r=dict(p)
 for k,v in q.items():
  r[k]=r.get(k,Q(0))+v
  if not r[k]:del r[k]
 return check(r)
def mul(p,q):
 r={}
 for k,v in p.items():
  for l,w in q.items():r[k+l]=r.get(k+l,Q(0))+v*w
 return check({k:v for k,v in r.items()if v})
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_GUARD_CONTROL_V2_20261105_v63';p=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_CACHED_KAPPA_ENTRY_20261105_v63';rp=p/'report.json'
 assert sha(rp)=='8f43ae03981eed91b74bbdb151be62010f3f69ee976b56ce5b368687c0b02fc8';r=read(rp);assert r['all_passed'] and r['only_one_tail1_value']
 for n,h in r['input_sha256'].items():assert sha(ROOT/n)==h
 for n,h in r['evidence_sha256'].items():assert sha(p/n)==h
 plan=read(b/'selected-plan.json');N=plan['nodes'];V={};peak=0;axis_atom_terms=0
 # These bounds hold BEFORE every QQ operation. A primitive numerator is
 # bounded by the fixed integer l1 majorant; unreduced Fraction add/mul
 # intermediates fit in 2*359+2*395+2=1510 bits, below2048.
 assert 2*359+2*395+2<2048
 with (b/'literal-atoms58.jsonl').open()as f:
  for line in f:
   a=json.loads(line);i=a['node'];v={};q=Q(a['QQ_divisor'])
   for ex,c in a['original_num_den'][0]:
    if ex[1]==ex[2]==0:
     k=ex[0]+a['a_shift'];assert k not in v;v[k]=Q(c)/q
   V[i]=check(v);axis_atom_terms+=len(v)
   peak=max(peak,20000+numeric(plan)+numeric(a)+4*sum(len(x)for x in V.values()));assert peak<80000
 needed=set();todo=plan['column_roots'][:]
 while todo:
  i=todo.pop()
  if i not in needed:needed.add(i);todo+=N[i]['args']
 for i in sorted(needed):
  n=N[i]
  if not n['args']:continue
  x,y=n['args']
  if n['op']=='mul_bound_only':V[i]=mul(V[x],V[y])
  else:V[i]=add(V[x],V[y] if n['op']=='add_bound_only' else{k:-q for k,q in V[y].items()})
  D=int(n['QQ_D']);L=int(n['integer_l1_upper']);assert all((q*D).denominator==1 and abs((q*D).numerator)<=L for q in V[i].values())
  peak=max(peak,24000+numeric(plan)+4*sum(len(x)for x in V.values()));assert peak<80000
 result={}
 for k,c in zip(plan['cached_kappa_nodes'],plan['column_roots']):result=add(result,mul(V[k],V[c]))
 assert result and result.get(3,Q(0))!=0
 # Compare only after independent substitution/evaluation. Read the producer
 # all-state output once, in a stream; no producer engine is imported.
 matched=[]
 with (p/'integer-column4-and-result.txt').open()as f:
  assert f.readline().strip()=='J1KV001'
  for expected in plan['column_roots']+[155]:
   node,D=map(int,f.readline().split());assert node==expected;m=int(f.readline());axis={}
   for _ in range(m):
    key,v=map(int,f.readline().split());aa=key//1048576;rest=key-aa*1048576;ee,dd=divmod(rest,1024);assert abs(v).bit_length()<=359
    if ee==dd==0:axis[aa]=Q(v,D)
   want=result if node==155 else V[node];assert axis==want
   matched.append(dict(node=node,axis_terms=len(axis)))
  assert not f.read().strip()
 inputs=[rp,b/'selected-plan.json',b/'literal-atoms58.jsonl',p/'integer-column4-and-result.txt',Path(__file__)];raw=sum(x.stat().st_size for x in inputs)+65536;assert raw<9517864
 out=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_AXIS_REVIEW_20261105_v63';out.mkdir(exist_ok=False)
 evidence=dict(independent_field='QQ[a,a^-1] after e=d=0',all149_DAG_nodes_from_original_literals=True,four_cached_kappa_axes=True,all_four_column_and_T1_axes_match=matched,axis_T1_coefficients=[[k,str(q)]for k,q in sorted(result.items())],nonzero_a3_e0_d0=str(result[3]),no_producer_engine_or_full3D_convolution_used=True,no_new_source_point_prime_or_full_rank_claim=True)
 ep=out/'independent-axis-witness.json';ep.write_text(json.dumps(evidence,indent=2)+'\n')
 report=dict(id=out.name,all_passed=True,T1_polynomial_nonzero_independently_proved=True,only_axis_subset_of_all7595_coefficients_peered=True,four_column_axes_and_complete_T1_axis_match=True,coefficient_witness=str(result[3]),preoperation_QQ_temporary_bound=1510,source_owned_peak=peak,hard_source_owned_limit=80000,total_RAW_upper=raw,axis_atom_terms=axis_atom_terms,input_sha256={str(x.relative_to(ROOT)):sha(x)for x in inputs},evidence_sha256={ep.name:sha(ep)},seconds=time.monotonic()-t)
 (out/'report.json').write_text(json.dumps(report,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if 'sha256'not in k}));print(sha(out/'report.json'))
if __name__=='__main__':main()
