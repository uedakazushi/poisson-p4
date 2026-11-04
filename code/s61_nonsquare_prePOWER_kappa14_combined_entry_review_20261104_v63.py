#!/usr/bin/env python3
"""Independent fixed source73/raw proportionality arithmetic, no other column."""
import collections,hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63';GATE='S61_NONSQUARE_PREPOWER_KAPPA14_PREFLIGHT_20261104_v63';ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dec(p):return{tuple(m):F(c)for m,c in p if F(c)}
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 r=load(PRIMARY);assert pins[f'reports/{PRIMARY}/report.json']=='8f4212a18217214e9727596e9df606d31f3cc17ef68b54c5d413dc242d83f7d8'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 op=load(GATE,'operator-and-cost.json');profile=load(GATE,'profiles.json')[2];N=op['nodes'];need=set(profile['needed_nodes']);assert len(N)==76 and len(need)==73
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json')
 entry='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';assert load(entry+'_REVIEW')['all_passed'];assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW')['all_passed']
 bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6'];polys={}
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args'])
  if n['op']!='atom':continue
  x=n['reference']
  if 'constant'in x:p=[[[0,0,0],x['constant']]]
  elif 'named_FULLU_unit'in x:p=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
  elif x['path']=='accepted_Q2L_source_atom':
   a=x['definition'];assert q['scalar_DAG']['nodes'][x['old_DAG_root']]==a
   num,den=cp['matrix25x32'][a[1]][a[2]]if a[0]=='checkpoint'else pd['Lframe_dual6x2'][a[1]][1]if a[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   assert len(den)==1 and den[0][0][1:]==[0,0]and n['a_shift']==-den[0][0][0]and F(n['QQ_unit'])==F(den[0][1]);p=num
  else:
   p=load(*x['path'].split('/')[1:])[x['key']]
   for j in x.get('index',[]):p=p[j]
  polys[i]=p
 # Independent exact expression fingerprints for the four ACTUAL current15 rows.
 memo={}
 def key(x):
  if x not in memo:memo[x]=len(memo)
  return memo[x]
 def atom(p):return key(('poly',tuple(sorted((m,str(c))for m,c in p.items()if c))))
 z=atom({});one=atom({(0,0,0):F(1)});minus=atom({(0,0,0):F(-1)})
 def mul(a,b):return z if z in [a,b]else b if a==one else a if b==one else key(('mul',*sorted((a,b))))
 def add(a,b,s=1):return a if b==z else b if a==z and s==1 else mul(minus,b)if a==z else z if s<0 and a==b else key(('sub'if s<0 else'add',a,b))
 E=[]
 for i,n in enumerate(N):
  if n['op']=='atom':v=atom({(m[0]+n['a_shift'],m[1],m[2]):F(c)/F(n['QQ_unit'])for m,c in polys[i]if F(c)})
  else:
   a,b=[E[j]for j in n['args']];v=mul(a,b)if n['op']=='mul'else add(a,b,-1 if n['op']=='sub'else 1)
  E.append(v)
 oldmemo={}
 def old(i):
  if i in oldmemo:return oldmemo[i]
  n=q['scalar_DAG']['nodes'][i];k=n[0]
  if k=='z':v=z
  elif k=='one':v=one
  elif k in ['a','s','m']:v=mul(old(n[1]),old(n[2]))if k=='m'else add(old(n[1]),old(n[2]),-1 if k=='s'else 1)
  else:
   a,b=cp['matrix25x32'][n[1]][n[2]]if k=='checkpoint'else pd['Lframe_dual6x2'][n[1]][1]if k=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][n[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0];v=atom({(m[0]-b[0][0][0],m[1],m[2]):F(c)/F(b[0][1])for m,c in a if F(c)})
  oldmemo[i]=v;return v
 cols=q['new22_target_original32_columns'];rows=[2,4,5,8];M14=op['column14_roots'];cofs=op['kappa4_roots'];cross=op['raw_cross3_outputs'];out=op['kappa14_output'];anchors=[56,57,58,59]
 assert rows==op['current15_source_rows']and op['common_B_clearing_power']==0
 for row,root in zip(rows,M14):assert E[root]==old(q['remaining15x22_scalar_roots'][row][cols.index(14)])
 for root,oldroot in zip(cofs,[179,198,217,237]):assert dec(polys[root])==dec(load(entry,f'root{oldroot}.json')['polynomial'])
 for i,root in enumerate(anchors):assert E[root]==atom(dec(bot[i][0]))
 a=z
 for k,m in zip(cofs,M14):a=add(a,mul(E[k],E[m]))
 assert E[out]==a
 for i,root in enumerate(cross,1):assert E[root]==add(mul(E[anchors[0]],E[M14[i]]),mul(E[anchors[i]],E[M14[0]]),-1)
 # Every ASCII payload and denominator bound matches the frozen mathematical source.
 ip=ROOT/'reports'/PRIMARY/'native-input.txt';pins[str(ip.relative_to(ROOT))]=sha(ip);t=iter(ip.read_text().split());assert(next(t),next(t))==('8484','76');records={}
 for i,n in enumerate(N):
  k=next(t)
  if i not in need:assert k=='x';continue
  upper,D=int(next(t)),int(next(t));assert upper==n['bound']['support']and D==int(n['QQ_denominator'])
  if k=='p':
   s,u=int(next(t)),next(t);count=int(next(t));payload=[tuple(next(t)for _ in range(4))for j in range(count)];amin=min(m[0]for m,c in polys[i]);assert s==n['a_shift']+amin and F(u)==F(n['QQ_unit'])
   assert {(int(a),int(b),int(c)):F(v)for a,b,c,v in payload if F(v)}=={(m[0]-amin,m[1],m[2]):F(c)for m,c in polys[i]if F(c)};records[i]=(k,s,u,payload)
  else:
   a,b=int(next(t)),int(next(t));assert n['args']==[a,b]and n['op']==dict(a='add',s='sub',m='mul')[k];records[i]=(k,a,b)
 roots=[int(next(t))for _ in range(16)];assert roots==[out]+cross+cofs+M14+anchors
 try:next(t);raise AssertionError('trailing')
 except StopIteration:pass
 saved=[]
 for root in roots:
  d=load(PRIMARY,f'root{root}.json');assert d['root']==root;saved.append(d['polynomial'])
 actualterms=sum(len(dec(p))for p in saved);assert 75284+actualterms<150000
 save('preflight.json',{'source73_and16_saved_values_bound':True,'extra_held_saved_nonzero_terms':actualterms,'conservative_live_upper':75284+actualterms,'no_other_column_or_minor':True})
 use=collections.Counter(j for i in need for j in N[i]['args']);lines=['76 16']
 for root,p in zip(roots,saved):lines.extend([str(root),str(len(dec(p)))]);lines.extend(f'{m[0]} {m[1]} {m[2]} {c}'for m,c in dec(p).items())
 lines.append(' '.join(str(use[i])for i in range(76)))
 for i in range(76):
  if i not in need:lines.append(f'{i} x');continue
  rec=records[i]
  if rec[0]=='p':lines.extend([f'{i} p {rec[1]} {rec[2]}',str(len(rec[3]))]);lines+=[' '.join(x)for x in rec[3]]
  else:lines.append(f'{i} {rec[0]} {rec[1]} {rec[2]}')
 inp=OUT/'fixed-input.txt';inp.write_text('\n'.join(lines)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=OUT/'independent-checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 run=subprocess.run([str(exe),str(inp),str(OUT/'native-check.json')],capture_output=True,text=True,timeout=max(1,56-(time.monotonic()-start)));(OUT/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0;check=json.loads((OUT/'native-check.json').read_text());assert check['all_passed']
 for p in [cpp,ROOT/'code/s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp']:pins[str(p.relative_to(ROOT))]=sha(p)
 save('independent-proportionality-proof.json',{'all73_actual_source_expression_payloads_bound':True,'all16_saved_polynomials_exact_QQ_match':True,'kappa14_and_all3_distinct_raw_cross_polynomials_zero':True,'anchor_M2_22_nonzero_and_old6column_rank1_peer_bound':True,'seven_fixed_raw_pure_columns':[22,11,13,25,15,19,14],'whole_raw4x7_rank_at_most1_by_fractionfield_then_polynomial_domain_descent':True,'no_pointwise_anchor_inverse_needed':True,'kappa14_zero_alone_NOT_used_as_raw_proportionality':True,'all_other_columns_tails_and_full_rank_open':True})
 rep={'id':ID,'all_passed':True,'frozen':True,'kappa14_QQ_zero':True,'all3_raw_cross_QQ_zero':True,'raw4x7_rank_at_most1_on_exact_seven_columns':True,'different_monomial_order':'DEGREVLEX','no_other_column_rank_global57_or_quartic239_claim':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),'evidence_sha256':{n:sha(OUT/n)for n in ['preflight.json','native-check.json','independent-proportionality-proof.json','fixed-input.txt']},'seconds':time.monotonic()-start,'native_seconds':check['seconds'],'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024;save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':main()
