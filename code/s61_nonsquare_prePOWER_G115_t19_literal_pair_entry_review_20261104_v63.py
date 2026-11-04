#!/usr/bin/env python3
"""Independent fixed-source/Laplace/ordinary QQ audit of ONE normalized D6."""
import collections,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
PRIMARY='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_ENTRY_20261104_v63'
GATE='S61_NONSQUARE_PREPOWER_G115_T19_LITERAL_PAIR_PREFLIGHT_20261104_v63'
ID=PRIMARY+'_REVIEW';OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb')as f:return hashlib.file_digest(f,'sha256').hexdigest()
def save(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+'\n')
def dec(p):return {tuple(m):F(c)for m,c in p if F(c)}
def sign(p):return (-1)**sum(p[i]>p[j]for i in range(len(p))for j in range(i+1,len(p)))
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;k=str(p.relative_to(ROOT));pins[k]=sha(p)
  if k not in cache:cache[k]=json.loads(p.read_text())
  return cache[k]
 r=load(PRIMARY);assert pins[f'reports/{PRIMARY}/report.json']=='2a727a0a37b0d89ef04b18f95ae3aebf2c82146f14c073789f27008690680e42'
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 op=load(GATE,'operator-and-cost.json');cost=load(GATE,'cost.json');N=op['nodes'];need=set(op['needed_fused_nodes']);assert len(N)==1357 and len(need)==1325
 assert load('S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW')['all_passed']
 assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW')['all_passed']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json')
 hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');bot=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6']
 fx=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');fixed={v['old_DAG_root']:v['Laurent_polynomial']for v in fx['node_values']}
 B={(1,0,1):F(50),(0,2,1):F(18),(0,1,1):F(-9),(0,0,1):F(-2),(0,1,0):F(60),(0,0,0):F(-40)}
 polys={}
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args'])
  if n['op']!='atom':continue
  ref=n['reference']
  if 'constant'in ref:p=[[[0,0,0],ref['constant']]]
  elif 'named_FULLU_unit'in ref:assert ref['named_FULLU_unit']=='B';p=[[list(m),str(c)]for m,c in B.items()]
  elif ref['path']=='accepted_Q2L_source_atom':
   d=ref['definition'];assert q['scalar_DAG']['nodes'][ref['old_DAG_root']]==d
   a,b=cp['matrix25x32'][d[1]][d[2]]if d[0]=='checkpoint'else pd['Lframe_dual6x2'][d[1]][1]if d[0]=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0] and n['a_shift']==-b[0][0][0] and F(n['QQ_unit'])==F(b[0][1]) and ref['known_B_denominator_power']==0;p=a
  else:
   x=load(*ref['path'].split('/')[1:])
   if 'old_DAG_root'in ref:p=next(v['Laurent_polynomial']for v in x['node_values']if v['old_DAG_root']==ref['old_DAG_root'])
   elif 'h9_target'in ref:p=x['source_h9'][x['original_target_labels'].index(ref['h9_target'])]
   else:
    p=x[ref['key']]
    for j in ref.get('index',[]):p=p[j]
  polys[i]=p
 # Source-expression fingerprints, with only accepted cached-node boundaries.
 memo={};descs=[]
 def intern(x):
  if x not in memo:memo[x]=len(descs);descs.append(x)
  return memo[x]
 def atom(p):return intern(('p',tuple(sorted((m,str(c))for m,c in p.items()if c))))
 z=atom({});one=atom({(0,0,0):F(1)});minus=atom({(0,0,0):F(-1)})
 def mul(a,b):
  if z in [a,b]:return z
  if a==one:return b
  if b==one:return a
  return intern(('m',*sorted((a,b))))
 def add(a,b,s=1):
  if b==z:return a
  if a==z:return b if s==1 else mul(minus,b)
  if s<0 and a==b:return z
  return intern(('s'if s<0 else'a',a,b))
 E=[]
 for i,n in enumerate(N):
  if n['op']=='atom':v=atom({(m[0]+n['a_shift'],m[1],m[2]):F(c)/F(n['QQ_unit'])for m,c in polys[i]if F(c)})
  elif n['op']=='known_a_shift':v=intern(('shift',E[n['args'][0]],n['a_exponent_shift']))
  else:
   a,b=[E[j]for j in n['args']];v=mul(a,b)if n['op']=='mul'else add(a,b,-1 if n['op']=='sub'else 1)
  E.append(v)
 oldcache={};cols=q['new22_target_original32_columns'];labels=[22,11,13,25,15,19]
 # Literal current15 row2 is a PREFIX-combination, NOT bare actual135=20.
 stopped={q['remaining15x22_scalar_roots'][2][cols.index(c)]:p for c,p in zip(labels,bot[0])}
 def old(i):
  if i in oldcache:return oldcache[i]
  d=q['scalar_DAG']['nodes'][i];k=d[0]
  if i in stopped:v=atom(dec(stopped[i]))
  elif i in fixed:v=atom(dec(fixed[i]))
  elif k=='z':v=z
  elif k=='one':v=one
  elif k in ['a','s','m']:
   a,b=old(d[1]),old(d[2]);v=mul(a,b)if k=='m'else add(a,b,-1 if k=='s'else 1)
  else:
   a,b=cp['matrix25x32'][d[1]][d[2]]if k=='checkpoint'else pd['Lframe_dual6x2'][d[1]][1]if k=='accepted_pair_second'else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0]
   v=atom({(m[0]-b[0][0][0],m[1],m[2]):F(c)/F(b[0][1])for m,c in a if F(c)})
  oldcache[i]=v;return v
 M=op['matrix6x6'];assert op['row_a_factors']==[0,1,1,0,0,0]and op['row_B_clearings']==[0]*6 and op['original32_columns']==labels
 for rr,current in enumerate([1,3,6],1):
  for j,c in enumerate(labels):
   root=M[rr][j]
   if rr in [1,2]:assert N[root]['op']=='known_a_shift' and N[root]['a_exponent_shift']==-1;root=N[root]['args'][0]
   assert E[root]==old(q['remaining15x22_scalar_roots'][current][cols.index(c)])
 for j in range(6):assert E[M[4][j]]==atom(dec(bot[0][j]))
 g=atom(dec(hc['raw_compact14x5'][6][1]));h=atom(dec(hc['raw_compact14x5'][10][1]))
 c6=F(hc['known_row_units_rho14'][6]['constant']);c10=F(hc['known_row_units_rho14'][10]['constant'])
 for j,c in enumerate(labels):
  expected=add(mul(mul(g,atom({(0,0,0):1/c10})),old(q['remaining15x22_scalar_roots'][11][cols.index(c)])),mul(mul(h,atom({(0,0,0):1/c6})),old(q['remaining15x22_scalar_roots'][7][cols.index(c)])),-1)
  assert E[M[5][j]]==expected
 # Free matrix-entry polynomial checks of ALL15 fixed Laplace minors and signs.
 def padd(a,b,s=F(1)):
  o=dict(a)
  for m,c in b.items():
   v=o.get(m,F(0))+s*c
   if v:o[m]=v
   else:o.pop(m,None)
  return o
 def pmul(a,b):
  o={}
  for m,c in a.items():
   for n,d in b.items():o=padd(o,{tuple(sorted(m+n)):c*d})
  assert len(o)<1000;return o
 matrix=[[{(6*i+j,):F(1)}for j in range(6)]for i in range(6)]
 def det(A):
  out={}
  for p in itertools.permutations(range(len(A))):
   t={():F(sign(p))}
   for i,j in enumerate(p):t=pmul(t,A[i][j])
   out=padd(out,t)
  return out
 cuts={root:matrix[i][j]for i,row in enumerate(M)for j,root in enumerate(row)}
 # The one old P4 minor is an accepted exact coefficient identity, with two a-row factors.
 cuts[232]=det([row[:4]for row in matrix[:4]]);assert N[232]['a_shift']==-2
 values={}
 def formal(i):
  if i in cuts:return cuts[i]
  if i in values:return values[i]
  n=N[i]
  if n['op']=='atom':assert 'constant'in n['reference'];v={():F(n['reference']['constant'])}if F(n['reference']['constant'])else{}
  else:
   a,b=map(formal,n['args']);v=pmul(a,b)if n['op']=='mul'else padd(a,b,-1 if n['op']=='sub'else 1)
  values[i]=v;return v
 events=op['laplace15'];assert len(events)==15
 for chosen,event in zip(itertools.combinations(range(6),4),events):
  assert event['top4_columns']==list(chosen)and event['sign']==(-1)**sum(chosen)
  rem=[i for i in range(6)if i not in chosen]
  assert formal(event['top4_root'])==det([[row[j]for j in chosen]for row in matrix[:4]])
  assert formal(event['bottom2_root'])==det([[row[j]for j in rem]for row in matrix[4:]])
 # Bind every producer ASCII payload, division and event-scale without using its code.
 p=ROOT/'reports'/PRIMARY/'native-input.txt';pins[str(p.relative_to(ROOT))]=sha(p);t=iter(p.read_text().split());assert (next(t),next(t))==('7174','1357');records={}
 for i,n in enumerate(N):
  kind=next(t)
  if i not in need:assert kind=='x';continue
  upper,D=int(next(t)),int(next(t));assert upper==n['bound']['support']and D==int(n['QQ_denominator'])
  if kind=='p':
   shift,u=int(next(t)),next(t);count=int(next(t));payload=[tuple(next(t)for _ in range(4))for j in range(count)];amin=min(m[0]for m,c in polys[i]);assert shift==n['a_shift']+amin and F(u)==F(n['QQ_unit'])
   assert {(int(a),int(b),int(c)):F(v)for a,b,c,v in payload if F(v)}=={(m[0]-amin,m[1],m[2]):F(c)for m,c in polys[i]if F(c)};records[i]=(kind,upper,D,shift,u,payload)
  elif kind=='k':
   child,shift=int(next(t)),int(next(t));assert n['op']=='known_a_shift'and n['args']==[child]and n['a_exponent_shift']==shift;records[i]=(kind,upper,D,child,shift)
  else:
   a,b=int(next(t)),int(next(t));assert n['op']==dict(a='add',s='sub',m='mul')[kind]and n['args']==[a,b];records[i]=(kind,upper,D,a,b)
 tail=[int(next(t))for _ in range(9)];assert tail[:6]==[x for p in N[op['output']]['bound']['box']for x in p]and tail[6:8]==N[op['output']]['bound']['degree']and tail[8]==int(N[op['output']]['QQ_denominator']);assert int(next(t))==15
 schedule=[]
 for j in range(15):
  a=[int(next(t))for _ in range(6)];at,k,u,v,s,dscale=a;ev=events[k]
  assert (u,v,s)==(ev['top4_root'],ev['bottom2_root'],ev['sign'])and at==max(u,v)
  assert tail[8]==dscale*int(N[u]['QQ_denominator'])*int(N[v]['QQ_denominator']);schedule.append(a[:5])
 try:next(t);raise AssertionError('trailing')
 except StopIteration:pass
 bound=cost['total_whole_worker_coefficient_slots_upper']+cost['dense_output_lattice_slots'];assert bound==147444 and bound<150000
 save('preflight.json',{'independent_ordinary_QQ_product_extra_bound':38931,'total_live_upper':bound,'saved9047_loaded_only_after_source_release':True,'only_fixed1325_source15_Laplace_terms':True,'same60s2GB32MiB_guard':True})
 native=load(PRIMARY,'native.json');saved=native['exact_normalized_D6'];assert len(saved)==9047 and not native['normalized_D6_is_zero']
 use=collections.Counter(j for i in need for j in N[i]['args'])
 for e in events:use[e['top4_root']]+=1;use[e['bottom2_root']]+=1
 lines=['1357',' '.join(str(use[i])for i in range(1357))]
 for i in range(1357):
  if i not in need:lines.append(f'{i} x');continue
  rec=records[i]
  if rec[0]=='p':
   _,upper,D,shift,u,payload=rec;lines.extend([f'{i} p {upper} {D} {shift} {u}',str(len(payload))]);lines+=[' '.join(p)for p in payload]
  else:lines.append(f'{i} '+' '.join(map(str,rec)))
 lines.append('15');lines+=[' '.join(map(str,e))for e in schedule];lines.append(str(len(saved)));lines+=[f'{m[0]} {m[1]} {m[2]} {c}'for m,c in saved]
 inp=OUT/'fixed-input.txt';inp.write_text('\n'.join(lines)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=OUT/'independent-checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split()
 subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 run=subprocess.run([str(exe),str(inp),str(OUT/'native-check.json')],capture_output=True,text=True,timeout=max(1,56-(time.monotonic()-started)));(OUT/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0
 checked=json.loads((OUT/'native-check.json').read_text());assert checked['all_passed']
 for p in [cpp,ROOT/'code/s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp']:pins[str(p.relative_to(ROOT))]=sha(p)
 save('independent-source-and-Laplace-proof.json',{'all1325_source_nodes_atom_payloads_known_unit_divisions_and_a_offsets_bound':True,'actual_source_rows_current15':[0,1,3,6,2,11,7],'G20_G115_G85_are_current15_2_11_7_not_bare135_generators':True,'Ksmall_literal_gG115_over_c10_minus_hG85_over_c6':True,'all15_free_matrix_top4_bottom2_determinants_signs_and_integer_QQ_denominator_scales':True,'raw_D6_is_a_squared_times_normalized_by_actual_two_row_factors':True,'normalized_D6_all9047_QQ_coefficients_match_independently':True,'no_120_parameter_permutation_or_D5_expansion':True,'nonzero_entry_only_no_global_rank_claim':True})
 rep={'id':ID,'all_passed':True,'frozen':True,'all1325_actual_source_operations_and15_QQ_product_reproduction':True,'normalized_D6_nonzero_terms':9047,'raw_D6_equals_a_squared_times_normalized':True,'QQ_identity_scope_not_only_modp':True,'only_one_fixed_G115_target19_no_rank_global57_or_quartic239_claim':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{ID}.md'),'evidence_sha256':{n:sha(OUT/n)for n in ['preflight.json','native-check.json','independent-source-and-Laplace-proof.json','fixed-input.txt']},'seconds':time.monotonic()-started,'native_seconds':checked['seconds'],'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 assert sum(p.stat().st_size for p in OUT.iterdir())<32*1024*1024;save('report.json',rep);signal.alarm(0);print(json.dumps({k:v for k,v in rep.items()if k!='input_sha256'},indent=2))
if __name__=='__main__':main()
