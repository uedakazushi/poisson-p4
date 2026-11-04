#!/usr/bin/env python3
"""Independent fixed240 source/sign/B-clearing and four-product QQ verifier."""
import collections,hashlib,itertools,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];P='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_ENTRY_20261104_v63';S=P+'_REVIEW';O=R/'reports'/S;G='S61_NONSQUARE_PREPOWER_KAPPA_EXTRA2_PREFLIGHT_20261104_v63'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def save(n,d):(O/n).write_text(json.dumps(d,indent=2)+'\n')
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));O.mkdir(exist_ok=True);pins={}
 def load(st,name='report.json'):
  p=R/'reports'/st/name;pins[str(p.relative_to(R))]=sha(p);return json.loads(p.read_text())
 primary=load(P);assert pins[f'reports/{P}/report.json']=='dd9e28bc3907b2af492bb504e4a10dea593edbc8436cc0fcdaf61ecc212ae810' and primary['all_passed']
 for p,h in primary['input_sha256'].items():assert sha(R/p)==h
 op=load(G,'operator-and-cost.json');cert=load(G,'certificate.json');needed=set(cert['needed_nodes']);N=op['nodes'];assert len(N)==240 and len(needed)==238
 assert load('S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW')['all_passed']
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');bottom=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6']
 BP={(1,0,1):F(50),(0,2,1):F(18),(0,1,1):F(-9),(0,0,1):F(-2),(0,1,0):F(60),(0,0,0):F(-40)}
 def dec(p):return {tuple(m):F(c) for m,c in p if F(c)}
 rawpolys={};denpower={}
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args'])
  if n['op']!='atom':continue
  ref=n['reference']
  if 'constant' in ref:p=[[[0,0,0],ref['constant']]]
  elif 'named_FULLU_unit' in ref:assert ref['named_FULLU_unit']=='B';p=[[list(m),str(c)] for m,c in BP.items()]
  elif ref['path']=='accepted_Q2L_source_atom':
   d=ref['definition'];assert q['scalar_DAG']['nodes'][ref['old_DAG_root']]==d
   a,b=cp['matrix25x32'][d[1]][d[2]] if d[0]=='checkpoint' else pd['Lframe_dual6x2'][d[1]][1] if d[0]=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   dd=dec(b)
   if len(b)==1:assert b[0][0][1:]==[0,0] and n['a_shift']==-b[0][0][0] and F(n['QQ_unit'])==F(b[0][1]);bp=0
   else:assert dd=={m:F(n['QQ_unit'])*c for m,c in BP.items()} and n['a_shift']==0;bp=1
   assert ref['known_B_denominator_power']==bp;denpower[ref['old_DAG_root']]=bp;p=a
  elif 'raw_compact14x5'==ref['key']:
   r,c=ref['index'];assert r in [1,3,4,7] and c in [0,4];p=hc['raw_compact14x5'][r][c];assert F(n['QQ_unit'])==1/F(hc['known_row_units_rho14'][r]['constant'])
  else:assert ref['key']=='bottom4x6' and ref['index'][1]==0;p=bottom[ref['index'][0]][0];assert n['a_shift']==0 and F(n['QQ_unit'])==1
  rawpolys[i]=p
 # Symbolic scalar-DAG fingerprints: all four rows and B denominator alignments retained.
 memo={};desc=[]
 def key(x):
  if x not in memo:memo[x]=len(desc);desc.append(x)
  return memo[x]
 def atom(poly):return key(('p',tuple(sorted((m,str(c)) for m,c in poly.items() if c))))
 z=atom({});one=atom({(0,0,0):F(1)});minus=atom({(0,0,0):F(-1)})
 def mul(a,b):
  if z in [a,b]:return z
  if a==one:return b
  if b==one:return a
  return key(('m',*sorted((a,b))))
 def add(a,b,s=1):
  if b==z:return a
  if a==z:return b if s==1 else mul(minus,b)
  if s<0 and a==b:return z
  return key(('s' if s<0 else 'a',a,b))
 E=[]
 for i,n in enumerate(N):
  if n['op']=='atom':v=atom({(m[0]+n['a_shift'],m[1],m[2]):F(c)/F(n['QQ_unit']) for m,c in rawpolys[i] if F(c)})
  else:
   a,b=[E[j] for j in n['args']];v=mul(a,b) if n['op']=='mul' else add(a,b,-1 if n['op']=='sub' else 1)
  E.append(v)
 B=atom(BP);oldcache={}
 def timesB(x,n):
  for _ in range(n):x=mul(x,B)
  return x
 def old(i):
  if i in oldcache:return oldcache[i]
  d=q['scalar_DAG']['nodes'][i];k=d[0]
  if k=='z':v=(z,0)
  elif k=='one':v=(one,0)
  elif k in ['a','s','m']:
   a,b=old(d[1]);c,e=old(d[2]);p=max(b,e)
   v=(mul(a,c),b+e) if k=='m' else (add(timesB(a,p-b),timesB(c,p-e),-1 if k=='s' else 1),p)
  else:
   a,b=cp['matrix25x32'][d[1]][d[2]] if k=='checkpoint' else pd['Lframe_dual6x2'][d[1]][1] if k=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   if len(b)==1:shift=-b[0][0][0];unit=F(b[0][1]);bp=0
   else:unit=dec(b)[next(iter(BP))]/BP[next(iter(BP))];assert dec(b)=={m:unit*c for m,c in BP.items()};shift=0;bp=1
   v=(atom({(m[0]+shift,m[1],m[2]):F(c)/unit for m,c in a if F(c)}),bp)
  oldcache[i]=v;return v
 rows=[2,4,5,8];compact=[1,3,4,7];cols=q['new22_target_original32_columns'];power=op['common_B_clearing_power'];assert power==1
 for r,root in zip(rows,op['cleared_column_roots']):
  x,p=old(q['remaining15x22_scalar_roots'][r][cols.index(2)]);assert E[root]==timesB(x,power-p)
 A=[]
 for i,r in enumerate(compact):A.append([atom({tuple(m):F(c)*F(hc['known_row_units_rho14'][r]['constant']) for m,c in hc['raw_compact14x5'][r][j]}) for j in [0,4]]+[atom(dec(bottom[i][0]))])
 for i,root in enumerate(op['cofactor_roots']):
  rr=[r for r in range(4) if r!=i];v=z
  for pp in itertools.permutations(range(3)):
   t=one
   for j in range(3):t=mul(t,A[rr[j]][pp[j]])
   v=add(v,t,(-1)**(i+sum(pp[j]>pp[k] for j in range(3) for k in range(j+1,3))))
  assert E[root]==v
 out=z
 for c,x in zip(op['cofactor_roots'],op['cleared_column_roots']):out=add(out,mul(E[c],E[x]))
 assert E[op['output']]==out and op['source_current15_rows']==rows and op['anchor_original32']==22 and op['source_original32_target']==2
 # Exact native input coefficient payloads and240 source refs, independent of worker's generator.
 p=R/'reports'/P/'native-input.txt';pins[str(p.relative_to(R))]=sha(p);tokens=iter(p.read_text().split());assert (next(tokens),next(tokens))==('10456','240');records={}
 for i,n in enumerate(N):
  kind=next(tokens)
  if i not in needed:assert kind=='x';continue
  assert kind!='x';upper,den=int(next(tokens)),int(next(tokens));assert upper==n['bound']['support'] and den==int(n['QQ_denominator'])
  if kind=='p':
   shift,unit=int(next(tokens)),next(tokens);count=int(next(tokens));poly=[tuple(next(tokens) for _ in range(4)) for j in range(count)];rp=rawpolys[i];amin=min((m[0] for m,c in rp),default=0)
   assert shift==n['a_shift']+amin and F(unit)==F(n['QQ_unit'])
   assert {(int(a),int(b),int(c)):F(v) for a,b,c,v in poly if F(v)}=={(m[0]-amin,m[1],m[2]):F(c) for m,c in rp if F(c)};records[i]=(kind,shift,unit,poly)
  else:
   a,b=int(next(tokens)),int(next(tokens));assert n['op']=={'m':'mul','a':'add','s':'sub'}[kind] and n['args']==[a,b];records[i]=(kind,a,b)
 roots=[int(next(tokens)) for _ in range(9)];assert roots==[239,179,198,217,237,66,94,122,150]
 try:next(tokens);raise AssertionError('trailing')
 except StopIteration:pass
 saved=[]
 for i in roots:
  d=load(P,f'root{i}.json');assert d['root']==i;saved.append(d['polynomial'])
 terms=sum(map(len,saved));assert terms==18036 and len(saved[0])==7208
 bound=cert['last_use_peak_upper']+terms+15000;assert bound==126699 and bound<150000
 save('preflight.json',{'fixed238_source_and_all4_signed_cofactor_frames_bound':True,'held9_saved_terms':terms,'no_all240_resident_state':True,'conservative_core_upper':bound,'four_product_upper':112416,'hard_seconds':60,'no_large_D8_expansion':True})
 uses=collections.Counter(j for i in needed for j in N[i]['args']);text=['240 9']
 for i,poly in zip(roots,saved):
  text.extend([str(i),str(len(poly))]);text.extend(f'{m[0]} {m[1]} {m[2]} {c}' for m,c in poly)
 text.append(' '.join(str(uses[i]) for i in range(240)))
 for i in range(240):
  if i not in needed:text.append(f'{i} x');continue
  rec=records[i]
  if rec[0]=='p':
   text.extend([f'{i} p {rec[1]} {rec[2]}',str(len(rec[3]))]);text.extend(' '.join(row) for row in rec[3])
  else:text.append(f'{i} {rec[0]} {rec[1]} {rec[2]}')
 inp=O/'fixed-input.txt';inp.write_text('\n'.join(text)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=O/'independent-fixed-checker';flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 run=subprocess.run([str(exe),str(inp),str(O/'native-check.json')],capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-start)));(O/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0;checked=json.loads((O/'native-check.json').read_text());assert checked['all_passed']
 save('independent-source-proof.json',{'all67_literal_atoms_payloads_and_B_denominators_bound':True,'source240_238_needed_full4row_composition':True,'all4_signed3cofactors_have_actual_anchor_M22_no_extra_a':True,'cleared_column_is_B_times_actual_projected_M2':True,'all9_actual_saved_polynomials_QQ_verified':True,'four_product_order':[2,0,3,1],'T_nonzero_7208terms':True,'known6_annihilation':'Universal signed-cofactor kappa*A=0; accepted bottom4x6 rank<=1 and nonzero M22 anchor give kappa*all6=0 in fractionfield and then in integral source ring, with no pointwise anchor inversion.','no_other17cols_global57_or_actual_rawH7_lift_claim':True})
 pins[str(cpp.relative_to(R))]=sha(cpp);utility=R/'code/s61_nonsquare_prePOWER_raw8_zero_review_20261104_v63.cpp';pins[str(utility.relative_to(R))]=sha(utility)
 report={'id':S,'all_passed':True,'all240_source_and_four_signed_cofactor_B_alignment_bound':True,'independent_QQ_DEGREVLEX_source238_and_fourproducts':True,'T_nonzero_terms':7208,'all9_saved_polynomials_verified':True,'projected_source_only_rawH7_lift_gap_retained':True,'no_D8_expansion_ranksearch_othercolumn_point_prime_or_GB':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{n:sha(O/n) for n in ['preflight.json','native-check.json','independent-source-proof.json','fixed-input.txt']},'seconds':time.monotonic()-start,'native_seconds':checked['seconds'],'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};assert sum(p.stat().st_size for p in O.iterdir())<33554432;save('report.json',report);print(json.dumps(report))
if __name__=='__main__':main()
