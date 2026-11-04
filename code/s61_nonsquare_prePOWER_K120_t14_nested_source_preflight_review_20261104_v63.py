#!/usr/bin/env python3
"""Independent source/sign/integer-lattice gate. No parameter products or worker."""
import hashlib,itertools,json,math,resource,signal,time
from functools import lru_cache
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
BASE='S61_NONSQUARE_PREPOWER_K120_T14_NESTED_SOURCE_PREFLIGHT_20261104_v63'
ID=BASE+'_REVIEW'; OUT=ROOT/'reports'/ID
BP=[[[1,0,1],'50'],[[0,2,1],'18'],[[0,1,1],'-9'],[[0,0,1],'-2'],[[0,1,0],'60'],[[0,0,0],'-40']]
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dec(p):return {tuple(m):F(c) for m,c in p if F(c)}
def parity(p):return sum(a>b for i,a in enumerate(p) for b in p[i+1:])%2
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
@lru_cache(None)
def slots(box,degree):
 # Count integer exponent triples independently by total-degree slicing.
 return sum(max(0,min(box[2][1],s-a-b)-max(box[2][0],s-a-b)+1)
  for s in range(degree[0],degree[1]+1)
  for a in range(box[0][0],box[0][1]+1) for b in range(box[1][0],box[1][1]+1))
def main():
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);pins={};cache={}
 def load(st,name='report.json'):
  p=ROOT/'reports'/st/name;key=str(p.relative_to(ROOT));pins[key]=sha(p)
  if key not in cache:cache[key]=json.loads(p.read_text())
  return cache[key]
 def bind(r):
  for k,h in r.get('input_sha256',{}).items():assert sha(ROOT/k)==h,(k,'input hash')
  for k,h in r.get('evidence_sha256',{}).items():assert sha(ROOT/'reports'/r['id']/k)==h,(k,'evidence hash')
 r=load(BASE);assert pins[f'reports/{BASE}/report.json']=='10ff652bd664c2227399e7fb809917fb6e8b868d05294e8f5f29b712bcaff42c';bind(r)
 wr=load('S61_NONSQUARE_PREPOWER_K120_WEDGE_SIGNS_20261104_v63');wp=load(wr['id'],'proof.json');bind(wr)
 oldst='S61_NONSQUARE_PREPOWER_K120_T14_CANCELLED_SOURCE_PREFLIGHT_20261104_v63';orr=load(oldst);bind(orr);assert not orr['resource_gate_passes']
 old=load(oldst,'operator-and-cost.json');plan=load(BASE,'nested-plan.json');cost=load(BASE,'cost.json')
 assert plan['external_scalar_nodes']==dict(path=f'reports/{oldst}/operator-and-cost.json',sha256=pins[f'reports/{oldst}/operator-and-cost.json'],node_count=7323)
 N=old['nodes']+plan['appended_scalar_nodes'];assert len(N)==11523
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');hc=load('S61_NONSQUARE_CUBIC_GENUINE14_H5_PREFLIGHT_20261104_v63','certificate.json');h9=load('S61_NONSQUARE_PURE_KERNEL12_SOURCE_PREFLIGHT_20261104_v63','prePOWER-h9-source.json')
 fix=load('S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63','native.json');fixed={v['old_DAG_root']:v['Laurent_polynomial'] for v in fix['node_values']}
 bottom=load('S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63','native.json')['bottom4x6']
 for st in ['S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63_REVIEW','S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW','S61_NONSQUARE_PREPOWER_FULL57_SOURCE_LIFT_20261104_v63_REVIEW']:
  assert load(st)['all_passed']
 # Actual atom payloads, units, and metadata are reconstructed from pinned sources.
 polys={};stats=[]
 for i,n in enumerate(N):
  assert n['id']==i and all(j<i for j in n['args'])
  if n['op']=='atom':
   x=n['reference']
   if 'constant' in x:p=[[[0,0,0],x['constant']]]
   elif 'named_FULLU_unit' in x:assert x['named_FULLU_unit']=='B';p=BP
   elif 'derived_existing_coefficient_difference' in x:
    assert x['source']=='N11,1-N9,1';v=dec(hc['raw_compact14x5'][11][1])
    for m,c in dec(hc['raw_compact14x5'][9][1]).items():v[m]=v.get(m,F(0))-c
    p=[[list(m),str(c)] for m,c in v.items() if c];assert dec(p)==dec(x['polynomial'])
   elif x['path']=='accepted_Q2L_source_atom':
    a=x['definition'];assert a==q['scalar_DAG']['nodes'][x['old_DAG_root']]
    p,den=cp['matrix25x32'][a[1]][a[2]] if a[0]=='checkpoint' else pd['Lframe_dual6x2'][a[1]][1] if a[0]=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
    assert len(den)==1 and den[0][0][1:]==[0,0] and n['a_shift']==-den[0][0][0] and F(n['QQ_unit'])==F(den[0][1]);assert x['known_B_denominator_power']==0
   else:
    path=x['path'];assert path.startswith('reports/');v=load(*path.split('/')[1:])
    if 'h9_target' in x:p=v['source_h9'][v['original_target_labels'].index(x['h9_target'])]
    elif 'old_DAG_root' in x:p=fixed[x['old_DAG_root']]
    else:
     p=v[x['key']]
     for j in x.get('index',[]):p=p[j]
   p={(m[0]+n['a_shift'],m[1],m[2]):c/F(n['QQ_unit']) for m,c in dec(p).items()};polys[i]=p
   D=math.lcm(*(c.denominator for c in p.values())) if p else 1;L=sum(abs(c.numerator)*(D//c.denominator) for c in p.values())
   b=None if not p else dict(support=len(p),box=[[min(m[k] for m in p),max(m[k] for m in p)] for k in range(3)],degree=[min(map(sum,p)),max(map(sum,p))])
  else:
   args=n['args'];A,D1,L1=stats[args[0]]
   if n['op']=='known_a_shift':
    s=n['a_exponent_shift'];b=dict(support=A['support'],box=[[v+s for v in A['box'][0]]]+A['box'][1:],degree=[v+s for v in A['degree']]);D,L=D1,L1
   else:
    B,D2,L2=stats[args[1]];ism=n['op']=='mul';assert ism or n['op'] in ['add','sub'];D=D1*D2 if ism else math.lcm(D1,D2);L=L1*L2 if ism else L1*(D//D1)+L2*(D//D2)
    box=[[A['box'][k][0]+B['box'][k][0],A['box'][k][1]+B['box'][k][1]] if ism else [min(A['box'][k][0],B['box'][k][0]),max(A['box'][k][1],B['box'][k][1])] for k in range(3)]
    dg=[A['degree'][k]+B['degree'][k] for k in range(2)] if ism else [min(A['degree'][0],B['degree'][0]),max(A['degree'][1],B['degree'][1])]
    b=dict(support=min(A['support']*B['support'] if ism else A['support']+B['support'],slots(tuple(map(tuple,box)),tuple(dg)),L),box=box,degree=dg)
   assert A is not None
  assert b==n['bound'] and D==int(n['QQ_denominator']) and L==int(n['integer_l1_upper']),(i,'bound')
  stats.append((b,D,L))
 # Expression provenance, with fixed previously checked source values as leaves.
 ids={}
 def key(k):
  if k not in ids:ids[k]=len(ids)
  return ids[k]
 def atom(p):return key(('p',tuple(sorted((m,str(c)) for m,c in p.items() if c))))
 zero=atom({});one=atom({(0,0,0):F(1)});minus=atom({(0,0,0):F(-1)})
 def mul(a,b):return zero if zero in [a,b] else b if a==one else a if b==one else key(('m',*sorted((a,b))))
 def add(a,b,s=1):return a if b==zero else b if a==zero and s==1 else mul(minus,b) if a==zero else zero if s==-1 and a==b else key(('s' if s<0 else 'a',a,b))
 E=[]
 for i,n in enumerate(N):
  if n['op']=='atom':v=atom(polys[i])
  elif n['op']=='known_a_shift':v=key(('shift',E[n['args'][0]],n['a_exponent_shift']))
  else:a,b=[E[j] for j in n['args']];v=mul(a,b) if n['op']=='mul' else add(a,b,-1 if n['op']=='sub' else 1)
  E.append(v)
 memo={}
 def oldexpr(i):
  if i in memo:return memo[i]
  a=q['scalar_DAG']['nodes'][i];op=a[0]
  if i in fixed:v=atom(dec(fixed[i]))
  elif op=='z':v=zero
  elif op=='one':v=one
  elif op in ['a','s','m']:v=mul(oldexpr(a[1]),oldexpr(a[2])) if op=='m' else add(oldexpr(a[1]),oldexpr(a[2]),-1 if op=='s' else 1)
  else:
   p,d=cp['matrix25x32'][a[1]][a[2]] if op=='checkpoint' else pd['Lframe_dual6x2'][a[1]][1] if op=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][a[1]]
   assert len(d)==1 and d[0][0][1:]==[0,0];v=atom({(m[0]-d[0][0][0],m[1],m[2]):c/F(d[0][1]) for m,c in dec(p).items()})
  memo[i]=v;return v
 cols=q['new22_target_original32_columns'];labels=[22,11,13,25,15,19,14];top=old['top5_matrix'];bot=old['bottomABC_matrix'];assert old['source_pure_columns_original32']==labels
 assert old['known_a_top5_row_factors']==[0,1,1,0,0] and old['known_B_column_clearings']==[0]*7
 for ri,row in enumerate([0,1,3,6]):
  for j,label in enumerate(labels):
   v=atom(dec(h9['source_h9'][h9['original_target_labels'].index(label)])) if row==0 else oldexpr(q['remaining15x22_scalar_roots'][row][cols.index(label)])
   if old['known_a_top5_row_factors'][ri]:v=key(('shift',v,-1))
   assert E[top[ri][j]]==v,(ri,j,'source top')
 for j,label in enumerate(labels):
  p=load('S61_NONSQUARE_PREPOWER_KAPPA14_COMBINED_ENTRY_20261104_v63','root21.json')['polynomial'] if label==14 else bottom[0][labels.index(label)]
  assert E[top[4][j]]==atom(dec(p))
 c6=F(hc['known_row_units_rho14'][6]['constant']);c10=F(hc['known_row_units_rho14'][10]['constant']);c11=F(hc['known_row_units_rho14'][11]['constant']);assert hc['known_row_units_rho14'][11]==hc['known_row_units_rho14'][9]
 for k,(rows,c) in enumerate([([11],c10),([12,10],c11),([7],c6)]):
  for j,label in enumerate(labels):
   v=oldexpr(q['remaining15x22_scalar_roots'][rows[0]][cols.index(label)])
   if len(rows)>1:v=add(v,oldexpr(q['remaining15x22_scalar_roots'][rows[1]][cols.index(label)]),-1)
   assert E[bot[k][j]]==mul(v,atom({(0,0,0):1/c})),(k,j,'bottom')
 g,hp,h=old['weights_g_h120_h'];assert polys[g]==dec(hc['raw_compact14x5'][6][1]) and polys[h]==dec(hc['raw_compact14x5'][10][1])
 # Independent signed formal expansion: literal K-vectors, full5 and full7 permutations.
 wedge={}
 for i,ci in [(0,(1,'g')),(2,(-1,'h'))]:
  for j,cj in [(1,(1,'g')),(2,(-1,'hp'))]:
   if i==j:continue
   ij=tuple(sorted((i,j)));s=ci[0]*cj[0]*(1 if i<j else -1);weights=tuple(sorted((ci[1],cj[1])));wedge[(ij,weights)]=wedge.get((ij,weights),0)+s
 assert wedge=={((0,1),('g','g')):1,((0,2),('g','hp')):-1,((1,2),('g','h')):1}
 full5={p:(-1 if parity(p) else 1) for p in itertools.permutations(range(5))};seen5={}
 for ab in itertools.combinations(range(5),2):
  rem=[j for j in range(5) if j not in ab]
  for xx in [ab,ab[::-1]]:
   for rr in itertools.permutations(rem):
    p=(xx[0],*rr,xx[1]);sgn=(-1 if (sum(ab)+(xx!=ab)+parity(tuple(rem.index(x) for x in rr)))%2 else 1);seen5[p]=sgn
 assert seen5==full5
 for p in itertools.permutations(range(7)):
  rem=sorted(p[-2:]);ch=[j for j in range(7) if j not in rem];left=tuple(ch.index(x) for x in p[:5]);assert (11+sum(rem)+parity(left)+(p[5]>p[6]))%2==parity(p)
 def det2(U,V,j,k):return add(mul(E[U[j]],E[V[k]]),mul(E[U[k]],E[V[j]]),-1)
 def det3(M):
  v=zero
  for p in itertools.permutations(range(3)):
   t=one
   for i in range(3):t=mul(t,E[M[i][p[i]]])
   v=add(v,t,-1 if parity(p) else 1)
  return v
 @lru_cache(None)
 def closure(roots):
  s=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i not in s:s.add(i);todo.extend(N[i]['args'])
  return tuple(sorted(s))
 def life(roots):
  need=closure(tuple(roots));uses={i:0 for i in need}
  for i in need:
   for j in N[i]['args']:uses[j]+=1
  for i in roots:uses[i]+=1
  live=peak=0
  for i in need:
   live+=stats[i][0]['support'] if stats[i][0] else 0;sc=stats[i][0]['support'] if N[i]['op'] in ['mul','add','sub'] else 0;peak=max(peak,live+sc)
   for j in N[i]['args']:
    uses[j]-=1
    if uses[j]==0:live-=stats[j][0]['support'] if stats[j][0] else 0
  return list(need),peak
 allatoms=set();allneed=set();maxphase=0;innerpairs=0;outerpairs=0;globalD=int(plan['global_output_integer_denominator']);phase_records=[]
 assert len(plan['outer21'])==21
 for idx,(v,ev) in enumerate(zip(plan['outer21'],old['laplace21'])):
  assert v['index']==idx and v['top5_columns']==ev['top5_columns'] and v['bottom2_columns']==ev['bottom2_columns'] and v['weighted2_root']==ev['weighted_bottom2_root']
  ch=v['top5_columns'];j,k=v['bottom2_columns'];assert v['outer_sign']==(-1 if (11+j+k)%2 else 1)
  A,B,C=bot;w=add(add(mul(E[g],det2(A,B,j,k)),mul(E[hp],det2(A,C,j,k)),-1),mul(E[h],det2(C,B,j,k)),-1);assert w==E[v['weighted2_root']]
  need,peak=life([v['weighted2_root']]);assert need==v['weighted2_source_nodes'] and peak==v['weighted2_peak'];allneed.update(need);allatoms.update(i for i in need if N[i]['op']=='atom')
  topD=int(v['top5_integer_denominator']);maxordinary=0
  assert len(v['inner10'])==10
  for t,(a,b) in zip(v['inner10'],itertools.combinations(range(5),2)):
   ca,cb=ch[a],ch[b];rem=[x for x in ch if x not in [ca,cb]];assert t['source_columns']==[ca,cb] and t['literal3_columns']==rem and t['sign']==(-1 if (a+b)%2 else 1)
   assert E[t['left2minor']]==det2(top[0],top[4],ca,cb);assert E[t['right3minor']]==det3([[top[i][x] for x in rem] for i in [1,2,3]])
   nd,pk=life([t['left2minor'],t['right3minor']]);assert nd==t['source_nodes'] and pk==t['source_peak_with_scratch'];allneed.update(nd);allatoms.update(i for i in nd if N[i]['op']=='atom');maxordinary=max(maxordinary,pk)
   d=stats[t['left2minor']][1]*stats[t['right3minor']][1];assert topD%d==0 and int(t['integer_QQ_scale'])==topD//d
   innerpairs+=stats[t['left2minor']][0]['support']*stats[t['right3minor']][0]['support']
  assert v['top5_bound']==stats[ev['top5_root']][0];tb=v['top5_bound'];dense=slots(tuple(map(tuple,tb['box'])),tuple(tb['degree']));ws=stats[v['weighted2_root']][0]['support'];pk=max(peak,dense+ws+maxordinary)
  assert dense==v['top5_dense_slots'] and ws==v['retained_weighted2_slots'] and maxordinary==v['ordinary_inner_term_peak'] and pk==v['phase_peak_excludingouteratomswork'];maxphase=max(maxphase,pk)
  d=topD*stats[v['weighted2_root']][1];assert globalD%d==0 and int(v['outer_integer_QQ_scale'])==globalD//d;outerpairs+=dense*ws
  phase_records.append({'index':idx,'source_peak':peak,'top5_slots':dense,'W_slots':ws,'inner_peak':maxordinary,'phase_peak':pk})
 assert sorted(allatoms)==plan['all_atom_nodes']
 literal=sum(stats[i][0]['support'] if stats[i][0] else 0 for i in allatoms);ob=stats[old['output']][0];outslots=slots(tuple(map(tuple,ob['box'])),tuple(ob['degree']));total=outslots+literal+15000+maxphase
 oldD=stats[old['output']][1];assert globalD%oldD==0;L=stats[old['output']][2]*(globalD//oldD)
 inputbytes=0
 for i in allatoms:
  _,d,n=stats[i];u=F(N[i].get('QQ_unit','1'));inputbytes+=(stats[i][0]['support'] if stats[i][0] else 0)*(len(str(n))+len(str(d))+len(str(abs(u.numerator)))+len(str(u.denominator))+100)
 outputbytes=outslots*(len(str(L))+45)+len(str(globalD))+1000
 recipebytes=(ROOT/'reports'/BASE/'nested-plan.json').stat().st_size+(ROOT/'reports'/BASE/'formal-sign-proof.json').stat().st_size
 artifact=recipebytes+inputbytes+outputbytes+3_000_000
 assert (outslots,literal,maxphase,total,artifact)==(78338,11990,39549,144877,23622579)
 assert total<=150000 and artifact<32*1024*1024 and L.bit_length()==547 and globalD.bit_length()==232
 assert innerpairs+outerpairs==2_601_321_498
 assert plan['raw_D7_over_g_equals']=='a^2*output' and plan['no_numeric_g_division'] and plan['all_three_multilinear_terms_preserved']
 evidence=dict(source_atoms_checked=len(polys),all11523_metadata_nodes=True,source_top35_bottom21=True,wedge_computed_bilinearly=True,top120_signs=True,outer5040_signs=True,all210_inner_integer_scales=True,all21_outer_integer_scales=True,phase_records=phase_records,source_current15_not_bare135=True,raw7_proportionality_peer_bound=True,raw_normalization='raw(D7/g)=a^2*output',total_live_upper=total,artifact_byte_bound=artifact,convolution_pairs_upper=innerpairs+outerpairs,runtime_not_guaranteed=True,scalar_worker_not_run=True)
 write(OUT/'independent-source-sign-resource-proof.json',evidence)
 note=ROOT/f'research/notes/{ID}.md';pins[str(note.relative_to(ROOT))]=sha(note);pins[str(Path(__file__).relative_to(ROOT))]=sha(Path(__file__))
 report=dict(id=ID,all_passed=True,frozen=True,source_only=True,input_sha256=pins,evidence_sha256={'independent-source-sign-resource-proof.json':sha(OUT/'independent-source-sign-resource-proof.json')},seconds=time.monotonic()-start,scope='One fixed candidate source/resource gate only; no scalar output, nonzero, rank, or completion-in-60s claim.')
 write(OUT/'report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items() if k not in ['input_sha256','evidence_sha256']},indent=2))
if __name__=='__main__':main()
