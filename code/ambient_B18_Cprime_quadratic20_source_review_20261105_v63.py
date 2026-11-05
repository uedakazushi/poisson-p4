#!/usr/bin/env python3
"""Independent cached SOURCE checks for the complete Cprime quadratic map.
This file does not solve a canonical gauge or evaluate any new source value.
"""
import ast,collections,gc,hashlib,itertools,json,re,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
ID='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_20261105_v63_REVIEW'
OUT=ROOT/'reports'/ID
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return {str(k):enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [enc(x)for x in v]
 return v
def kernel20(d):
 cols=[{int(i):F(v)for i,v in c}for c in d['all44_actual_linear_columns27_sparse']]
 cp,rr=d['fixed_pivot_columns'],d['fixed_pivot_rows'];iv=[[F(v)for v in row]for row in d['fixed_minor_inverse']]
 assert len(cp)==len(rr)==len(iv)==24
 free=[i for i in range(44)if i not in cp];assert len(free)==20;K=[]
 for f in free:
  t=[sum((iv[i][j]*cols[f].get(rr[j],F(0))for j in range(24)),F(0))for i in range(24)];k={f:F(1)}
  for i,v in enumerate(t):
   if v:k[cp[i]]=-v
  assert all(sum((v*cols[j].get(i,F(0))for j,v in k.items()),F(0))==0 for i in range(27));K.append(k)
 assert [[k.get(j,F(0))for j in free]for k in K]==[[F(i==j)for j in range(20)]for i in range(20)]
 assert sum(map(len,K))==49
 return free,K
def ring_projection(d):
 b27=[()]+[(i,)for i in range(9)]+[m for m in itertools.combinations_with_replacement(range(9),2)if m[1]>=7]
 assert len(b27)==27;b26=[m for m in b27 if m!=(8,8)];ix={m:i for i,m in enumerate(b26)};mp=[]
 for z,r,e in d['basis33']:
  m=tuple(([z-1]if z else[])+[7]*r+[8]*e);mp.append({ix[m]:F(1)}if m in ix else{})
 for i in range(33):
  for j in range(33):
   got={};want={}
   for k,c in d['table_columns_for_each_left_basis'][i][j]:
    for t,v in mp[k].items():got[t]=got.get(t,F(0))+F(c)*v
   got={t:v for t,v in got.items()if v}
   for a,c in mp[i].items():
    for b,e in mp[j].items():
     m=tuple(sorted(b26[a]+b26[b]))
     if m in ix:want[ix[m]]=want.get(ix[m],F(0))+c*e
   assert got==want
 return b27,b26,mp

def minimal_full9(path):
 """Decode only the required unique top-level blocks, then prune unused QQ."""
 text=path.read_text();decoder=json.JSONDecoder()
 def block(key):
  hits=list(re.finditer(r'"'+re.escape(key)+r'"\s*:',text));assert len(hits)==1
  pos=hits[0].end()
  while text[pos].isspace():pos+=1
  return decoder.raw_decode(text,pos)[0]
 both=block('both57_source2');mons=block('source2_monomials')
 names=['unknown57','equations57','F9_nodes','all57_source2jets','all57_equations_all55_coefficients_verified']
 result={f:{k:both[f][k]for k in names}for f in ['CA','CB']}
 del both,text;gc.collect()
 return mons,result

def source_DAG_profile(source,both):
 ids=[set(i for i,v in both[f]['all57_source2jets'])for f in ['CA','CB']]
 assert all(len(x)==57 for x in ids)and len(ids[0]&ids[1])==9 and len(ids[0]|ids[1])==105
 assert ids[0]&ids[1]==set(range(2,11))
 for f in ['CA','CB']:
  b=both[f];assert b['unknown57']==[i for i,v in b['all57_source2jets']]
  assert b['equations57'][-9:]==b['F9_nodes']and len(b['equations57'])==57
  assert b['all57_equations_all55_coefficients_verified']and all(len(v)==55 for i,v in b['all57_source2jets'])
 nodes=source['DAG_scalar_nodes'];roots=[i for f in source['families'].values()for pair in f['complete_normalized_physical_w_coefficients'][0]for i in pair]+[f['actual_tracefree_nu']for f in source['families'].values()]+[i for f in ['CA','CB']for i in both[f]['equations57']]+[7790,28539]
 bind=ids[0]|ids[1];seen=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in seen:continue
  seen.add(i)
  if i in bind:continue
  n=nodes[i];assert n[0]in ['c','a','m','n','i']
  if n[0]!='c':todo.extend(n[1:])
 inverse={i:(j,F(c))for i,j,c in source['all_unit_inverse_nodes_and_origin_values']}
 for i in seen-bind:
  if nodes[i][0]=='i':assert i in inverse and inverse[i][0]==nodes[i][1]and inverse[i][1]
 uses=collections.Counter(roots)
 for i in seen-bind:
  if nodes[i][0]!='c':uses.update(nodes[i][1:])
 live=set();peak=0
 for i in sorted(seen):
  live.add(i);peak=max(peak,len(live)*27)
  if i not in bind and nodes[i][0]!='c':
   for j in nodes[i][1:]:
    uses[j]-=1
    if not uses[j]:live.remove(j)
 assert len(seen)==20383 and peak==13689
 return dict(common105=sorted(bind),shared9=sorted(ids[0]&ids[1]),retained_roots=sorted(set(roots)),all114_residual_labels={f:both[f]['equations57']for f in ['CA','CB']},F9_target_subtraction_required={f:both[f]['F9_nodes']for f in ['CA','CB']},ancestors=len(seen),source_QQ_peak=peak,all_source_inverses_have_recorded_nonzero_origin=True)

def numeric_leaves(v):
 if isinstance(v,bool)or v is None:return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return int(re.fullmatch(r'-?\d+(?:/\d+)?',v)is not None)
 if isinstance(v,(list,tuple)):return sum(map(numeric_leaves,v))
 if isinstance(v,dict):return sum(map(numeric_leaves,v.values()))
 return 0

def verify_subset(receipt,subset):
 p=ROOT/receipt['source_path'];raw=p.read_bytes();assert sha(p)==receipt['source_sha256']
 text=raw.decode();assert text.isascii();decoder=json.JSONDecoder();selected={}
 for key,rec in receipt['raw_decode_fields_only'].items():
  hits=list(re.finditer(r'"'+re.escape(key)+r'"\s*:',text));assert len(hits)==1;start=hits[0].end()
  while text[start].isspace():start+=1
  value,end=decoder.raw_decode(text,start);assert (start,end)==(rec['start_byte'],rec['end_byte'])
  assert hashlib.sha256(raw[start:end]).hexdigest()==rec['source_slice_sha256'];selected[key]=value
 assert numeric_leaves(selected['both57_source2'])==receipt['largest_parsed_block_numeric_leaves']==29310
 names=['unknown57','equations57','F9_nodes','f63','determinant','all57_source2jets','all57_equations_all55_coefficients_verified']
 wanted={'source2_monomials':selected['source2_monomials'],'both57_source2':{f:{k:selected['both57_source2'][f][k]for k in names}for f in ['CA','CB']}}
 assert subset==wanted and numeric_leaves(subset)==receipt['retained_numeric_leaves_including_labels']==6733
 assert sum(len(v)for f in subset['both57_source2'].values()for i,v in f['all57_source2jets'])==6270
 return dict(all_subset_fields_literal=True,extraction_peak_numeric_leaves=29310,retained_numeric_leaves=6733,retained_rational_tensor_entries=6270,disjoint_input_parse_and_arithmetic_phases=True,whole_source_not_decoded=True)

def main():
 st=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));OUT.mkdir(exist_ok=True);assert not(OUT/'report.json').exists();pins={}
 def load(stem,name='report.json',expected=None):
  p=ROOT/'reports'/stem/name;h=sha(p);pins[str(p.relative_to(ROOT))]=h
  if expected:assert h==expected
  return json.loads(p.read_text())
 def checkpins(d):
  for p,h in d.get('input_sha256',{}).items():assert sha(ROOT/p)==h;pins[p]=h
 def save(n,v):(OUT/n).write_text(json.dumps(enc(v),separators=(',',':'))+'\n')
 primary='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE_PREFLIGHT_20261105_v63'
 r=load(primary,expected='5c3bf7163906cfbd94c0af76b1c5d489ce936e986527c677235415235273f6e3');checkpins(r)
 for n,h in r['evidence_sha256'].items():assert sha(ROOT/'reports'/primary/n)==h
 p=load(primary,'complete-four-by20-and-ell-square-cache-source-gate.json')
 for s,h in [('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63_REVIEW','4108f86bd3e122b92fd0488ec970fac0a38fae494194b4454617b906d3a91450'),('AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW','92d5c943b12da089eb504862043a6e811ad17bef9d07e2b693e4eaf07a539b30'),('AMBIENT_B18_CPRIME_PRODUCT_ANN_SOURCE_20261105_v63_REVIEW','4283817886115e59a158afa331b415a83161056a6e6c1dd9444e4c39a47d89c7')]:assert load(s,expected=h)['all_passed']
 sub='AMBIENT_B18_CPRIME_QUADRATIC20_SOURCE2_SUBSET_20261105_v63';sr=load(sub,expected='bb05969e852de56ff58fa3d808033898141ee3ca4d7b015e8c60b754447b4caf');ss=load(sub,'minimal-two57-source2-and-labels.json',expected='9c3cd8d7a91b9f079bf36b348daabaa5348b91d0ce5f4674cd9d9cbfbcaeaf1a');custody=verify_subset(sr,ss)
 v1=load('AMBIENT_B18_CPRIME_QUADRATIC20_PREPARED_WORKER_20261105_v63',expected='6bd6f687d3c618fc57d764b0ebb64dedd0a67f73a1b20d5dba5a312c0629e7d3');checkpins(v1)
 v2=load('AMBIENT_B18_CPRIME_QUADRATIC20_PREPARED_WORKER_V2_20261105_v63',expected='5b2086b79e1c11729187fc09b4825b4cbbddc026a1a8ead187319cb77b5e7302');checkpins(v2)
 old=ROOT/'code/ambient_B18_Cprime_quadratic20_values_20261105_v63.py';new=ROOT/v2['corrected_worker'];assert sha(old)==v1['prepared_worker_sha256']==v2['old_prepared_worker_unchanged_sha256']=='1c3d658491d481a03d9156120008055d270afc1f8511283831a6879eafeb9b4e';assert sha(new)==v2['corrected_worker_sha256']=='8c0c791a47419c59a520a06255e7cf1a48a4db8cb5aca3d734b0f6a2563d45f8'
 for name,manifest in [('repository_preserved_prepared_source',v1),('repository_preserved_corrected_worker',v2)]:
  q=ROOT/manifest[name];assert sha(q)==manifest['snapshot_sha256'];pins[str(q.relative_to(ROOT))]=sha(q)
 oldtext=old.read_text();newtext=new.read_text();assert newtext==oldtext.replace('+[7790];needed=', '+[7790,28539];needed=').replace(' pairs={f:[[vals[x],vals[y]]',' assert overlap==164\n pairs={f:[[vals[x],vals[y]]')
 tree=ast.parse(newtext);fieldmap=next(ast.literal_eval(n.value)for n in ast.walk(tree)if isinstance(n,ast.Assign)and any(isinstance(t,ast.Name)and t.id=='keys'for t in n.targets));assert 'source2'not in fieldmap and set(fieldmap)=={'source','interface','proper','second5','J30','firstpreflight','source33','Vfirst','gauge33','proper33','adapter2','radical'}
 assert 'cubic_stream'not in newtext and 'actual_complete_rank' in newtext and '--source-peer-report' in newtext
 for literal in ['signal.alarm(60)','peak<=118472','2048','32<<20','del vals,own,YB,J,old105',"del D['source33'],D['adapter2']",'assert overlap==164','assert res==[b.Z]*57','assert all(v==a54.Z for v in constraints())']:
  assert literal in newtext,literal
 first=load('AMBIENT_B18_CPRIME_LINEAR_INITIAL_VALUES_20261105_v63','phase1-actual-fullnine-first-V-columns.json');free,K=kernel20(first)
 assert free==p['full9_V1_kernel_free_source_labels']and K==[{int(i):F(c)for i,c in col}for col in p['all20_constant_kernel_columns44_sparse']]
 table=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63','exact-existing-raw33-multiplication-table.json');b27,b26,maps=ring_projection(table)
 assert b27==list(map(tuple,p['base27_monomials']))and b26==list(map(tuple,p['base26_monomials']))and maps==[{int(i):F(c)for i,c in col}for col in p['actual33_to26_columns']]
 assert ss['source2_monomials']==[[]]+[[i]for i in range(9)]+[list(m)for m in itertools.combinations_with_replacement(range(9),2)]
 source=load('AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63','actual-source-coefficients.json');profile=source_DAG_profile(source,ss['both57_source2']);assert profile['ancestors']==p['source_DAG_ancestor_count']and profile['source_QQ_peak']==p['budget']['stream_actual_physical_DAG_peak']
 old33=load('AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63','actual-source105-ell2-raw33-values.json');oldroots=set(i for i,v in old33['all_requested164_root_values']);assert len(oldroots)==164 and oldroots.issubset(profile['retained_roots']);assert {f:ss['both57_source2'][f]['f63']for f in ['CA','CB']}=={'CA':7790,'CB':28539}
 q4=load('AMBIENT_B18_CPRIME_PRODUCT_ANN_QUOTIENT_SOURCE_PREFLIGHT_20261105_v63');assert p['actual_NF4_rows_by135_sparse']==q4['literal_fixed_pi4_rows_by135_sparse'];mons2=list(itertools.combinations_with_replacement(range(9),2));support={mons2[i%45]for row in p['actual_NF4_rows_by135_sparse']for i,c in row};assert support=={(7,7),(7,8),(8,8)}
 J=load('AMBIENT_B18_NATIVE128_FIRST_20261104_v63','phase2-fixed30-Jacobian.json');a0=[[F(c)for c in row]for row in J['A0']];a1=[[F(c)for c in row]for row in J['A1_eta_delta']];ji=[[F(c)for c in row]for row in J['J30_inverse']];jj=[a0[i]+[F(0)]*15 for i in range(15)]+[a1[i]+a0[i]for i in range(15)];assert F(J['det_J30'])==F(4,9)
 assert all(sum((jj[i][k]*ji[k][j]for k in range(30)),F(0))==F(i==j)and sum((ji[i][k]*jj[k][j]for k in range(30)),F(0))==F(i==j)for i in range(30)for j in range(30))
 interface=load('AMBIENT_B18_PROPER_PHYSICAL_INTERFACE_FIRSTJET_20261104_v63','actual-source-physical-firstjet-interface.json');proper=load('AMBIENT_B18_PROPER_POLYNOMIAL_MODEL_PREFLIGHT_20261104_v63','literal-proper-source14-preflight.json');assert len(interface['actual_group15_positions'])==15 and len(proper['CA_fiber10_fixed_coefficient_indices'])==10 and len(proper['CA_free_fiber9_coefficient_indices'])==9 and interface['retained_source9_order']==['a','d','v2','v3','v4','v5','v6','v7','v8']
 budget=p['budget'];assert sum(budget.values())==r['whole_saved_live_QQ_upper']==v1['whole_QQ_bound']==v2['whole_QQ_bound_unchanged']==118472<150000
 assert budget['cached_full9_source2_used']==2*57*55 and budget['common105_base27']==105*27 and budget['physical150_pairs_and_nus']==602*27 and budget['gauge30_base27']==30*27 and budget['actual14_kappa9_C108']==23*108 and budget['finite_GL_source_work']==400*108
 assert p['hard_seconds']==60 and p['hard_total_bytes']==32<<20 and p['derived_bit_guard']==2048
 assert p['optional_cubic_source_count_bound']==44*45+4*9+3+27==2046
 save('independent-complete4by20-source-frame-and-immutable-worker-v2.json',dict(free20=free,K20=K,all540_V1_kernel_and400_tail_identity_entries=True,all1089_A33_B26_products=True,source114=profile,selective_input_custody=custody,unit_J30_left_right_entries=1800,J30_determinant='4/9',all164_diagnostic_roots_retained=True,old_prepared_driver_preserved=True,narrow_v2_change='Add retainedCBf63=28539 and require overlap164. Ancestor20383/peak13689 unchanged.',whole_live_QQ_upper=118472,separate_extraction_phase_peak=29310,finite_GL_C108_envelope=43200,whole_complete_quadratic_source_theorem='ForV0=0, everydegree2initial is V2*a0+V1*a1 withV1*a0=0. WholeconstantkernelK20 yields the exactNF4*V2*K20 map. No all-order kernel/freeness claim.',mandatory_value_gates=['Both cached own57s compose in actual full9 YA/YB with ALLnine adapter2 curvature; sharedphysical9 must agree.','Check all114 equations, subtracting full YA/YB targets in the last9 per family; check105source and all164B26 root overlaps.','Only ell² gauge30 is new derived data. SameJ30 inverse, all30 constraints exactly0 and fullB26 gauge/proper overlap required.','Keepnu carry and fullproper rational functions beforetrace. Plane6 projection is an algebra quotient only after the actual fullB27 proper values are established.','All44 origins and plane firsts match prioractualV1; applyliteralNF4 and wholeK20. Fixed4x20 rank requires savednonzero minor plusall20span checks.','Keep118472live/2048bit/60s/2GB/32MiB guards and preservedphase onfailure. No source2implicit solve, source3, oroldr/g quotient.'],optional_cubic_NOT_in_this_worker=True,actual_gauge_V2_or_rank_NOT_computed=True))
 for q in[Path(__file__),ROOT/f'research/notes/{ID}.md',new,old]:pins[str(q.relative_to(ROOT))]=sha(q)
 raw=sum(q.stat().st_size for q in OUT.rglob('*')if q.is_file());assert raw<32<<20 and time.monotonic()-st<60
 save('report.json',dict(id=ID,all_passed=True,frozen=True,source_frame_budget_and_prepared_V2_gate_passes=True,complete_actual_quadratic_shape=[4,20],whole_live_QQ_upper=118472,new_derived_gauge_coefficients_allowed=30,new_physical_source_or_jet_coefficients=0,actual_value_or_rank_job_run=False,optional_cubic_separately_gated=True,prepared_worker_sha256=sha(new),input_sha256=pins,evidence_sha256={q.name:sha(q)for q in OUT.iterdir()if q.is_file()and q.name!='report.json'},actual_RAW=raw,seconds=time.monotonic()-st,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
 print(json.dumps(dict(all_passed=True,whole_QQ_upper=118472,actual_value_job_run=False,seconds=time.monotonic()-st),indent=2))
if __name__=='__main__':main()
