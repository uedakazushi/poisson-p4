#!/usr/bin/env python3
"""SOURCE: faithful B63/plane10 cache composition and positive bounds."""
import hashlib,json,math,signal,time
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_20261105_v63'
PATHS={
 'B63':'reports/AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63/finite63-socle-safe-output-theorem-and-actual-lower-input-overlaps.json',
 'B63_peer':'reports/AMBIENT_B18_C_CONNECTING_B63_SAFE_OUTPUT_SOURCE_20261105_v63_REVIEW/report.json',
 'head29':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63/coherent-head29-high-operator-and-three-selected-plane4-source-gate.json',
 'head29_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_SOURCE_20261105_v63_REVIEW/report.json',
 'head30':'reports/AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63/full-coherent-head30-and-safe63-cached-callbacks-source.json',
 'head30_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63_REVIEW/report.json',
 'head30_zero':'reports/AMBIENT_B18_C_CONNECTING_HEAD30_SOURCE_20261105_v63/inactive-beta3-native-carry-literal-zero-companion.json',
 'original':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63/compiled-original164-DAG-and-provided-pure-ell-inputs.json',
 'original_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_SOURCE_20261105_v63_REVIEW/report.json',
 'C2':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-physical12-and-C2-right-exact-quotient.json',
 'oldroots':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',
 'oldroots_peer':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json',
 'known_lambda':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63/phase1-single-root-existing-cache-readback.json',
 'known_lambda_peer':'reports/AMBIENT_B18_C_CONNECTING_HEAD23_CACHED_VALUES_20261105_v63_REVIEW/report.json',
 'ell':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63/all164-cached-pure-ell-readbacks-and-physical12-second-values.json',
 'ell_peer':'reports/AMBIENT_B18_C3_PURE_ELL_READBACK_VALUES_20261105_v63_REVIEW/report.json',
 'availability':'reports/AMBIENT_B18_C_CONNECTING_HEAD29_CACHE_AVAILABILITY_20261105_v63/literal-raw33-plane9-map-and-one-selected-fourth-head29-gate.json',
}
H=0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def ceilabs(q):return(abs(q.numerator)+q.denominator-1)//q.denominator
def ht(v):return max(abs(v[0].numerator).bit_length(),v[0].denominator.bit_length())if isinstance(v[0],Q)else(v[0]*v[1]).bit_length()+1
def track(v):
 global H
 H=max(H,ht(v));assert H<=384 and 2*H+1<=1153 and 2*H+1<=2048,'SOURCE_POSITIVE_HEIGHT_STOP';return v
def bd(q):return track((ceilabs(q),q.denominator))
def plus(a,b):
 if isinstance(a[0],Q)and isinstance(b[0],Q):return track((a[0]+b[0],1))
 if isinstance(a[0],Q):a=bd(a[0])
 if isinstance(b[0],Q):b=bd(b[0])
 return track((a[0]+b[0],math.lcm(a[1],b[1])))
def times(a,b):
 if isinstance(a[0],Q)and isinstance(b[0],Q):return track((a[0]*b[0],1))
 if isinstance(a[0],Q):a=bd(a[0])
 if isinstance(b[0],Q):b=bd(b[0])
 if not a[0]or not b[0]:return(0,1)
 return track((a[0]*b[0],a[1]*b[1]))
def main():
 start=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 peers={k:json.loads((ROOT/p).read_text())for k,p in PATHS.items()if k.endswith('_peer')};assert all(p['all_passed']for p in peers.values())
 original=S.selective(ROOT/PATHS['original'],['operator_records_originalID_op_args'],receipts)
 omap={r[0]:r for r in original['operator_records_originalID_op_args']}
 C2=S.selective(ROOT/PATHS['C2'],['physical12_coefficients'],receipts)
 phys={r['physical_row']:r for r in C2['physical12_coefficients']};oddrows=[18,19,28,29,60,63,67,72]
 groots=[20144,40922];oroots=[phys[r]['roots']['bB_node']for r in oddrows];roots=groots+oroots
 def cl(rr):
  wanted=set();todo=list(rr)
  while todo:
   i=todo.pop()
   if i in wanted:continue
   wanted.add(i);r=omap[i]
   if r[1]not in['c','input']:todo.extend(r[2:])
  return wanted
 gneed=cl(groots);allneed=cl(roots);records=[omap[i]for i in sorted(allneed)];precision=[63 if r[0]in gneed else 10 for r in records]
 # Shared original records are counted once; lists are borrowed views plus
 # index metadata.  Nothing derived is a requested coefficient response.
 phase0=S.count(original)+2*len(omap)+len(gneed)+len(allneed)+len(records)+len(precision)+S.count(C2)+S.count(peers)+12000
 assert phase0<=150000,('SOURCE_ORIGINAL_CUSTODY_STOP',phase0)
 del original,omap,gneed,allneed
 b=S.selective(ROOT/PATHS['B63'],['basis63_monomials','provided105_lower_extension63','socle32_indices'],receipts)
 mon63=[tuple(m)for m in b['basis63_monomials']];mon10=[m for m in mon63 if all(j>=7 for j in m)and len(m)<=3]
 mons={63:mon63,10:mon10};indices={p:{m:i for i,m in enumerate(mm)}for p,mm in mons.items()};idx10=[indices[63][m]for m in mon10]
 pp={}
 for p,mm in mons.items():
  table=[]
  for i,a in enumerate(mm):
   for j,c in enumerate(mm):
    k=indices[p].get(tuple(sorted(a+c)))
    if k is not None:table.append([i,j,k])
  pp[p]=table
 def zero(p):return[(Q(0),1)]+[(0,1)]*(p-1)
 def restrict(v,p):return v if len(v)==p else[v[j]for j in idx10]
 def conv(a,c,p):
  v=zero(p)
  for i,j,k in pp[p]:v[k]=plus(v[k],times(a[i],c[j]))
  return v
 given={}
 for node,arr in b['provided105_lower_extension63']:
  v=zero(63)
  for i,s in enumerate(arr):
   q=Q(s);v[i]=(q,1)if i==0 else bd(q)
  given[node]=v
 assert len(given)==105
 uses=Counter(roots)
 for r in records:
  if r[1]not in['c','input']:uses.update(r[2:])
 cache={};outputs={};max_live=0;maxQQ=0;constants=[];units=[];majorhash=hashlib.sha256();retire=[]
 for r,p in zip(records,precision):
  i,op=r[:2]
  if op=='input':v=restrict(given[i],p)
  elif op=='c':v=zero(p);v[0]=(Q(r[2],r[3]),1);constants.append(r)
  elif op=='a':v=[plus(a,c)for a,c in zip(restrict(cache[r[2]],p),restrict(cache[r[3]],p))]
  elif op=='n':v=restrict(cache[r[2]],p).copy();v[0]=(-v[0][0],1)
  elif op=='m':v=conv(restrict(cache[r[2]],p),restrict(cache[r[3]],p),p)
  elif op=='i':
   a=restrict(cache[r[2]],p);assert a[0][0]==1;v=zero(p);v[0]=(Q(1),1)
   for k in range(1,p):
    for x,y,t in pp[p]:
     if t==k and x>0 and y<k:v[k]=plus(v[k],times(a[x],v[y]))
   conv(a,v,p);conv(v,a,p);units.append([i,r[2],'1'])
  else:raise AssertionError(op)
  cache[i]=v;max_live=max(max_live,len(cache));maxQQ=max(maxQQ,sum(map(len,cache.values())))
  majorhash.update((json.dumps([i,p,[[str(a),str(c)]for a,c in v]],separators=(',',':'))+'\n').encode())
  if i in roots:outputs[i]=v.copy();uses[i]-=roots.count(i)
  gone=[]
  if op not in['c','input']:
   for j in r[2:]:
    uses[j]-=1
    if not uses[j]:del cache[j];gone.append(j)
  if not uses[i]:del cache[i];gone.append(i)
  retire.append(gone)
 assert not cache and set(outputs)==set(roots)
 # Positive-only source bound for subtraction; no output coefficient dot.
 delta=[plus(a,c)for a,c in zip(outputs[groots[0]],outputs[groots[1]])]
 root_H=H;root_pre=2*H+1
 old=S.selective(ROOT/PATHS['oldroots'],['all_requested164_root_values'],receipts);old=dict(old['all_requested164_root_values'])
 ell=S.selective(ROOT/PATHS['ell'],['all164_actual_pure_ell_root_values3'],receipts);ell=dict(ell['all164_actual_pure_ell_root_values3'])
 requested=[]
 for row,node in zip(oddrows,oroots):
  vals=dict(old[node]);checks=[]
  for m,rawix in[((),0),((7,),8),((8,),9),((7,7),24),((7,8),25),((7,7,7),29),((7,7,8),30)]:
   checks.append([list(m),vals.get(rawix,'0')])
  # Pure ell² is absent as an independent raw33 slot. Read the accepted
  # literal pureell164 cache rather than assign it by a support guess.
  checks.append([[8,8],ell[node][2]])
  requested.append({'original_node':node,'physical_row':row,'trusted_precision':10,'root_lower_checks':checks})
 known=S.selective(ROOT/PATHS['known_lambda'],['requested_bB63_r_ell2'],receipts)
 heads={}
 for tag in['head29','head30']:
  heads[tag]=S.selective(ROOT/PATHS[tag],['known_cached_inputs','four_genuine_and_canonical_high_polynomials_UNPAIRED'],receipts)
 leafbounds={}
 for h in heads.values():
  for label,q,locator in h['known_cached_inputs']:
   label=tuple(label)
   if label in leafbounds:assert leafbounds[label]==bd(Q(q))
   leafbounds[label]=bd(Q(q))
 for label,mon in[(('NEW_delta63_plane4','r4'),(7,7,7,7)),(('NEW_delta63_plane4','r3ell'),(7,7,7,8)),(('NEW_delta63_plane4','r2ell2'),(7,7,8,8))]:leafbounds[label]=delta[indices[63][mon]]
 leafbounds[('NEW_CACHED_bB63_ell3',)]=outputs[phys[63]['roots']['bB_node']][indices[10][(8,8,8)]]
 for m,j in indices[63].items():leafbounds[('SAFE_delta63',)+m]=delta[j]
 for row,node in zip(oddrows,oroots):leafbounds[('TRUSTED_plane3','bB',row,1,2)]=outputs[node][indices[10][(7,8,8)]]
 post=[]
 for tag,h in heads.items():
  bounds=[]
  for poly in h['four_genuine_and_canonical_high_polynomials_UNPAIRED']:
   v=(0,1)
   for labels,weight in poly:
    term=(Q(weight),1)
    for label in labels:term=times(term,leafbounds[tuple(label)])
    v=plus(v,term)
   bounds.append(v)
  post.append([tag,[[str(a),str(c)]for a,c in bounds]])
 post_H=H;post_pre=2*H+1
 # All phases use disjoint parsed custody; compilation is serialized before
 # root/post bound work.  Runtime phase1 releases records/input vectors
 # before loading the two symbolic postoperators.
 assert sum(map(len,retire))==len(records)
 runtime1={'records_precision_and_precompiled_retirement':S.count(records)+len(precision)+S.count(retire),'runtime_dependency_Counter':0,
  'given105_vectors_and_IDs':105*63+105,'mixed_precision_lastuse_QQ_and_keys':maxQQ+2*max_live,
  'root_outputs_and_product_tables':sum(map(len,outputs.values()))+S.count(pp),
  'all_lowercheck_values_and_metadata':S.count(requested)+S.count(known),'controllers_copies_guard_scratch':10000}
 runtime2={'post_polys_cached_leaves_locators':S.count(heads),'root_outputs_and_safe_difference':sum(map(len,outputs.values()))+63,
  'leafbound_and_working_term_views':2*len(leafbounds)+5000,'serialization_fixed_kernel_basis_controller':10000}
 future=max(sum(runtime1.values()),sum(runtime2.values()));assert future<=75000,('FUTURE_VALUE_STATIC_STOP',runtime1,runtime2)
 phase1=S.count(records)+len(precision)+S.count(retire)+S.count(b)+S.count(given)+2*len(records)+2*maxQQ+S.count(pp)+S.count(outputs)+S.count(peers)+S.count(C2)+2*len(phys)+12000
 # The SOURCE compiler keeps parsed records/input bounds while computing
 # postmajorants. Count every such still-owned array, including the old
 # overlap tables. Runtime will instead release phase1 before post.
 phase2=S.count(records)+len(precision)+S.count(retire)+S.count(b)+S.count(given)+2*len(records)+S.count(pp)+S.count(C2)+2*len(phys)+S.count(old)+S.count(ell)+S.count(known)+S.count(peers)+S.count(heads)+S.count(leafbounds)+S.count(outputs)+S.count(post)+S.count(requested)+12000
 source_whole=max(phase0,phase1,phase2);assert source_whole<150000,('SOURCE_WHOLE_STOP',phase0,phase1,phase2)
 compiled={'original_records_precision_and_retirement':[[r,p,gone]for r,p,gone in zip(records,precision,retire)],'basis63_monomials':[list(m)for m in mon63],
  'plane10_monomials':[list(m)for m in mon10],'provided105_lower_extension63':b['provided105_lower_extension63'],
  'delta_safe_root_pair':groots,'odd_trusted_root_requests':requested,'unit_inverse_centers':units,
  'known_bB63_r_ell2':known['requested_bB63_r_ell2'],'known_delta63_plane4_r4_r3ell':['-1236125/1728','-91325/288'],
  'output_precision_warning':'Only difference20144-40922 is safe on entireB63. Each bB output is used only in plane10. No wholeB63 implicit105-equation claim.',
  'original_provided_input_provenance':'B63SOURCE all630 input overlaps and lower31 exact; originalQ20 105/114 equations and jointplane3 are inherited by immutable pins, not re-solved.',
  'head29_head30_paths':[PATHS['head29'],PATHS['head30']]}
 profile={'scope':'SOURCE/cache-only simultaneous fullHEAD29/30 fourcanonicalhighs. No new physical/IFT/gauge/jet, no requested root value or head response evaluated.',
  'records':len(records),'precision63_records':precision.count(63),'precision10_records':precision.count(10),
  'lastuse_peak_nodes':max_live,'lastuse_peak_QQ_slots':maxQQ,'all_node_positive_majorant_sha256':majorhash.hexdigest(),
  'literal_n_d_constants':constants,'root_stored_bound':root_H,'root_before_operation_bound':root_pre,
  'post_partial_stored_bound':post_H,'post_before_operation_bound':post_pre,'post_output_U_D_bounds':post,
  'future_runtime1_numeric_components':runtime1,'future_runtime2_numeric_components':runtime2,'future_static_owned_upper':future,
  'SOURCE_phase_numeric_upper':[phase0,phase1,phase2],'whole_SOURCE_owned_upper':source_whole,
  'static_cap':75000,'hard_value_owned':100000,'global_owned':150000,'committed_QQ_bits':384,'temporary_QQ_bits':1153,'hard_preop_bits':2048,'seconds':60,'memory_bytes':2147483648,
  'fresh_RAW_cap':2000000,'global_RAW_cap':33554432,'all_new_output_values_UNCOMPUTED':True,
  'further_worker_gate_OPEN':'Need immutable guarded worker/code/custody, all phase lifetime/RAW/operational prechecks, fresh Mark SOURCE/code/wholepeer and ROOT explicitONEGO before any requested coefficient composition.',
  'selective_source_receipts':receipts}
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False)
 cp=out/'compiled-ten-roots-mixed63-plane10-and-actual-provided-lower-inputs.json';cp.write_text(json.dumps(compiled,separators=(',',':'))+'\n')
 ppout=out/'positive-before-actual-twohead-cache-and-whole-budget.json';ppout.write_text(json.dumps(profile,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_TEN_ROOTS_MIXED63_PLANE10_BOTH_HEADS_POSITIVE_NO_VALUES',
  'actual_root_or_edge_values_evaluated':False,'future_static_owned_upper':future,'whole_SOURCE_owned_upper':source_whole,
  'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),
  'evidence_sha256':{p.name:sha(p)for p in[cp,ppout]},'seconds':time.monotonic()-start}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
