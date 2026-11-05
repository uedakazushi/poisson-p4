#!/usr/bin/env python3
"""Root independent pre-emission review of fixed j1 source, heights and custody."""
import hashlib,json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
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
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 b=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL1_WORKER_CUSTODY_V2_20261105_v63';rp=b/'report.json';r=read(rp);assert sha(rp)=='7b836b3d4a1a88ce1665d8377ee0138eb48c6ec6496c16de0c593bf3e20376c4'
 for name,h in r['evidence_sha256'].items():assert sha(b/name)==h
 for name,h in r['input_sha256'].items():assert sha(ROOT/name)==h
 plan=read(b/'selected-plan.json');sp=ROOT/'reports/S61_NONSQUARE_PREPOWER_TAIL8_COMPLEMENT_SOURCE_20261105_v63/typed-family-and-source-plans.json';family=read(sp);assert plan==family['source_plans'][0] and family['other8'][0]==1
 assert family['rows']==[2,4,5,8] and family['anchor']==22 and family['pair_H2']==[27,31]
 del family
 N=plan['nodes'];cols=plan['column_roots'];kap=plan['cached_kappa_nodes'];assert cols==[66,94,122,150] and kap==[151,152,153,154]
 need=set();todo=cols[:]
 while todo:
  i=todo.pop()
  if i not in need:need.add(i);todo.extend(N[i]['args'])
 assert len(need)==149 and all(x<i for i in need for x in N[i]['args'])
 with (b/'native-input.txt').open()as f,(b/'literal-atoms58.jsonl').open()as atoms:
  assert f.readline().strip()=='J1K001 155 149'
  assert list(map(int,f.readline().split()))==cols and list(map(int,f.readline().split()))==kap
  Dfinal=int(f.readline());assert Dfinal==int(N[161]['QQ_D']);assert list(map(int,f.readline().split()))==sorted(need)
  for i,n in enumerate(N[:155]):
   op={'external_atom':0,'add_bound_only':1,'sub_bound_only':2,'mul_bound_only':3}[n['op']];args=n['args']or[-1,-1];sup=n['bound']['support']if n['bound'] else 0
   assert list(map(int,f.readline().split()))==[i,op,*args,sup,int(n['integer_l1_upper']).bit_length(),int(n['QQ_D'])]
  assert int(f.readline())==57;total=0;atompeak=0;kappa_checked=0
  for line in atoms:
   a=json.loads(line);i=a['node'];n=N[i];u=Q(a['QQ_divisor']);D=int(n['QQ_D']);p,den=a['original_num_den'];assert a['reference']==n['reference'] and a['a_shift']==n['a_shift'] and a['QQ_divisor']==n['QQ_divisor']
   if 'cached_kappa_root'in n['reference']:
    cached=read(ROOT/n['reference']['path']);assert cached['polynomial']==p;kappa_checked+=1;del cached
   assert list(map(int,f.readline().split()))==[i,len(p)];L=0;seen=set()
   for m,c in p:
    key,v=map(int,f.readline().split());expected_key=(m[0]+n['a_shift'])*1048576+m[1]*1024+m[2];assert key==expected_key and key not in seen;seen.add(key)
    qq=Q(c)/u;assert qq*D==v;L+=abs(v)
   assert L==int(n['integer_l1_upper']);total+=len(p);atompeak=max(atompeak,numeric(a))
  assert not f.read().strip()
 assert total==14188 and kappa_checked==4
 pairs=0
 for i in sorted(need):
  n=N[i]
  if not n['args']:continue
  x,y=[N[j]for j in n['args']];D,dx,dy=map(int,[n['QQ_D'],x['QQ_D'],y['QQ_D']]);L,lx,ly=map(int,[n['integer_l1_upper'],x['integer_l1_upper'],y['integer_l1_upper']])
  if n['op']=='mul_bound_only':assert D==dx*dy and L==lx*ly;pairs+=x['bound']['support']*y['bound']['support']
  else:assert D%dx==D%dy==0 and L==lx*(D//dx)+ly*(D//dy)
  assert L.bit_length()<=359 and max(D.bit_length(),dx.bit_length()+dy.bit_length())<=395
  assert n['bound']['box'][1][0]>=0 and n['bound']['box'][1][1]<1024 and n['bound']['box'][2][0]>=0 and n['bound']['box'][2][1]<1024
 # Fixed fused dot: one common denominator, no product polynomial maps.
 M=0;dotpairs=0
 for k,c in zip(kap,cols):
  dk,dc=int(N[k]['QQ_D']),int(N[c]['QQ_D']);assert Dfinal%(dk*dc)==0
  M+=int(N[k]['integer_l1_upper'])*int(N[c]['integer_l1_upper'])*(Dfinal//(dk*dc));dotpairs+=N[k]['bound']['support']*N[c]['bound']['support']
 assert M==int(N[161]['integer_l1_upper']) and M.bit_length()==359 and dotpairs==14283624 and pairs==1313080
 dot_terms=sum(N[i]['bound']['support']for i in kap+cols)+N[161]['bound']['support'];assert 2*dot_terms+15100==91676
 # All source operands coexist until the current result is finished, then
 # last-use release. Check the forecast at the exact C++ allocation points.
 uses={i:0 for i in range(155)}
 for i in sorted(need):
  for j in N[i]['args']:uses[j]+=1
 for i in kap+cols:uses[i]+=1
 live={i:N[i]['bound']['support']for i in need|set(kap) if N[i]['op']=='external_atom'};peak_terms=sum(live.values())
 for i in sorted(need):
  if not N[i]['args']:continue
  live[i]=N[i]['bound']['support'];peak_terms=max(peak_terms,sum(live.values()))
  for j in N[i]['args']:
   uses[j]-=1
   if uses[j]==0:del live[j]
 source_bound=2*peak_terms+15100
 # A guard miss is an explicit stop, never permission to raise91676.
 result=dict(id='S61_NONSQUARE_PREPOWER_TAIL1_ROOT_SOURCE_REVIEW_20261105_v63',all_passed=True,source_only=True,literal14188_normalizations_and_native155_headers_exact=True,source149_DAG_and_cached4_kappa_bound=True,source_mul_pairs= pairs,dot_pairs_each_pass=dotpairs,three_pass_dot_pairs=3*dotpairs,coefficient_l1_bit_bound=359,denominator_temporary_bound=395,dot_owned_bound=91676,source_owned_forecast=source_bound,source_guard_stop_required_if_forecast_exceeds91676=source_bound>91676,source_review_owned_bound=24000+numeric(plan)+2*atompeak,no_coefficient_expression_evaluated=True,control_companion_review_still_required=True,input_sha256={str(rp.relative_to(ROOT)):sha(rp),str(sp.relative_to(ROOT)):sha(sp),str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))},seconds=time.monotonic()-start)
 assert result['source_review_owned_bound']<91676
 out=ROOT/'reports'/result['id'];out.mkdir(exist_ok=False);(out/'report.json').write_text(json.dumps(result,indent=2)+'\n');signal.alarm(0);print(json.dumps({k:v for k,v in result.items()if k not in ['input_sha256']}));print(sha(out/'report.json'))
if __name__=='__main__':main()
