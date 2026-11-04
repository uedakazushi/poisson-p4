#!/usr/bin/env python3
"""Independent exact QQ source62/cross3 check, no kappa18 node or other column."""
import hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE='S61_NONSQUARE_PREPOWER_RAWCROSS18_ENTRY_20261104_v63';GATE='S61_NONSQUARE_PREPOWER_K92_T18_SOURCE_PREFLIGHT_20261104_v63';ID=BASE+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dec(p):return {tuple(m):F(c)for m,c in p if F(c)}
def write(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(base,name='report.json'):
  p=ROOT/'reports'/base/name;rel=str(p.relative_to(ROOT));pins[rel]=sha(p)
  if rel not in cache:cache[rel]=json.loads(p.read_text())
  return cache[rel]
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='dc6c3ee7d3fea2c8d750ccea47ef653f2dac95853845fccd8dc19a73eca69e13'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/BASE/p)==h
 op=load(GATE,'operator-and-cost.json');N=op['nodes'];profile=load(GATE,'profiles.json')[1];binding=load(BASE,'source-binding.json');need=set(binding['source_ops62']);assert len(N)==76 and len(need)==62 and not need.intersection([52,53,54,55,60,61,62,63,64,65,66])
 assert sorted(need)==profile['needed_nodes']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];polys={}
 for i in sorted(need):
  n=N[i];assert n['id']==i and all(j<i and j in need for j in n['args'])
  if n['op']!='atom':continue
  x=n['reference']
  if 'constant'in x:p=[[[0,0,0],x['constant']]]
  elif x['path']=='accepted_Q2L_source_atom':
   a=x['definition'];assert a==q['scalar_DAG']['nodes'][x['old_DAG_root']]
   p,d=cp['matrix25x32'][a[1]][a[2]]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   assert len(d)==1 and d[0][0][1:]==[0,0]and n['a_shift']==-d[0][0][0]and F(n['QQ_unit'])==F(d[0][1])
  else:
   assert x['key']=='bottom4x6';p=bot
   for j in x['index']:p=p[j]
  polys[i]=p
 ids={}
 def key(x):
  if x not in ids:ids[x]=len(ids)
  return ids[x]
 def atom(p):return key(('p',tuple(sorted((m,str(c))for m,c in p.items()if c))))
 z=atom({});one=atom({(0,0,0):F(1)});minus=atom({(0,0,0):F(-1)})
 def mul(a,b):return z if z in[a,b]else b if a==one else a if b==one else key(('m',*sorted((a,b))))
 def add(a,b,s=1):return a if b==z else b if a==z and s==1 else mul(minus,b)if a==z else z if s==-1 and a==b else key(('s'if s<0 else'a',a,b))
 E={}
 for i in sorted(need):
  n=N[i]
  if n['op']=='atom':v=atom({(m[0]+n['a_shift'],m[1],m[2]):c/F(n['QQ_unit'])for m,c in dec(polys[i]).items()})
  else:a,b=[E[j]for j in n['args']];v=mul(a,b)if n['op']=='mul'else add(a,b,-1 if n['op']=='sub'else 1)
  E[i]=v
 memo={}
 def old(i):
  if i in memo:return memo[i]
  n=q['scalar_DAG']['nodes'][i];typ=n[0]
  if typ=='z':v=z
  elif typ=='one':v=one
  elif typ in['a','s','m']:v=mul(old(n[1]),old(n[2]))if typ=='m'else add(old(n[1]),old(n[2]),-1 if typ=='s'else 1)
  else:
   p,d=cp['matrix25x32'][n[1]][n[2]]if typ=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if typ=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   assert len(d)==1 and d[0][0][1:]==[0,0];v=atom({(m[0]-d[0][0][0],m[1],m[2]):c/F(d[0][1])for m,c in dec(p).items()})
  memo[i]=v;return v
 rows=[2,4,5,8];cols=q['new22_target_original32_columns'];roots=binding['column18_roots'];anchors=binding['anchor_roots'];outputs=binding['rawcross3_outputs'];assert roots==[21,31,41,51]and anchors==[56,57,58,59]and outputs==[69,72,75]
 for i,row in enumerate(rows):assert E[roots[i]]==old(q['remaining15x22_scalar_roots'][row][cols.index(18)])and E[anchors[i]]==atom(dec(bot[i][0]))
 for i,rt in enumerate(outputs,1):assert E[rt]==add(mul(E[anchors[0]],E[roots[i]]),mul(E[anchors[i]],E[roots[0]]),-1)
 assert load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW')['all_passed']
 # Verify literal producer ASCII corresponds to those source values, including
 # the harmless transfer of minimum a exponent to Laurent offset.
 ip=ROOT/'reports'/BASE/'native-input.txt';t=iter(ip.read_text().split());assert sum(len(p)for p in polys.values())==3564;assert(next(t),next(t))==(str(sum(len(p)for p in polys.values())),'76')
 for i,n in enumerate(N):
  k=next(t)
  if i not in need:assert k=='x';continue
  upper,D=int(next(t)),next(t);assert upper==(n['bound']['support']if n['bound']else 0)and D==n['QQ_denominator']
  if k=='p':
   shift,u,cnt=int(next(t)),next(t),int(next(t));raw=polys[i];m=min((e[0]for e,c in raw),default=0);assert shift==n['a_shift']+m and F(u)==F(n['QQ_unit'])and cnt==len(raw)
   p={(int(next(t)),int(next(t)),int(next(t))):F(next(t))for _ in range(cnt)};assert p=={(e[0]-m,e[1],e[2]):F(c)for e,c in raw}
  else:assert [int(next(t)),int(next(t))]==n['args']and dict(a='add',s='sub',m='mul')[k]==n['op']
 order=[int(next(t))for _ in range(11)];assert order==outputs+roots+anchors
 try:next(t);raise AssertionError('trailing')
 except StopIteration:pass
 saved=[load(BASE,f'root{i}.json')['polynomial']for i in order];assert sum(len(dec(p))for p in saved)+37236<150000
 uses=[0]*76
 for i in need:
  for j in N[i]['args']:uses[j]+=1
 own=OUT/'independent-input.txt'
 with own.open('w')as f:
  f.write('76 11\n')
  for rt,p in zip(order,saved):
   f.write(str(rt)+'\n'+str(len(p))+'\n')
   for m,c in p:f.write(' '.join(map(str,m+[c]))+'\n')
  f.write(' '.join(map(str,uses))+'\n')
  for i,n in enumerate(N):
   if i not in need:f.write(f'{i} x\n');continue
   if n['op']=='atom':
    p=polys[i];f.write(f"{i} p {n['a_shift']} {n['QQ_unit']}\n{len(p)}\n")
    for m,c in p:f.write(' '.join(map(str,m+[c]))+'\n')
   else:f.write(f'{i} '+dict(add='a',sub='s',mul='m')[n['op']]+' '+' '.join(map(str,n['args']))+'\n')
 cpp=Path(__file__).with_suffix('.cpp');exe=OUT/'independent-worker';cmd='source /workspace/cloud-setup/activate.sh && c++ -O2 -std=c++17 '+str(cpp)+' -o '+str(exe)+' $(pkg-config --cflags --libs flint)';s=subprocess.run(['bash','-lc',cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=max(1,55-(time.monotonic()-st)));(OUT/'compile.log').write_text(s.stdout);assert s.returncode==0
 cmd='source /workspace/cloud-setup/activate.sh && exec '+str(exe)+' '+str(own)+' '+str(OUT/'native.json');s=subprocess.run(['bash','-lc',cmd],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=max(1,58-(time.monotonic()-st)));(OUT/'native.log').write_text(s.stdout);assert s.returncode==0
 native=json.loads((OUT/'native.json').read_text());assert native['all_passed']and native['all3_raw_cross_QQ_products_zero']
 for p in[Path(__file__),cpp,ROOT/'code/s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp',ROOT/f'research/notes/{ID}.md']:pins[str(p.relative_to(ROOT))]=sha(p)
 report=dict(id=ID,all_passed=True,frozen=True,source_ops62=True,all11_saved_QQ_polynomials_match=True,three_actual_rawcross18_products_zero=True,no_kappa18_or_othercolumn=True,fixed_raw_pair_rank_at_most1_on_columns=[22,11,13,25,15,19,14,18],no_full17pure_or_rank_or_global57_claim=True,input_sha256=pins,evidence_sha256={p.name:sha(p)for p in OUT.iterdir()if p.is_file()and p.name not in['report.json','independent-worker']},seconds=time.monotonic()-st)
 write(OUT/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
