#!/usr/bin/env python3
"""Independent source extraction binding, QQ190+90 check; no rank search."""
import hashlib,json,resource,signal,subprocess,time
from fractions import Fraction as F
from pathlib import Path
R=Path(__file__).resolve().parents[1];S='S61_NONSQUARE_PREPOWER_RAW8_ZERO_20261104_v63_REVIEW';O=R/'reports'/S
B='S61_NONSQUARE_PREPOWER_RAW8_BOTTOM24_SOURCE_20261104_v63';W='S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_ENTRY_20261104_v63';G='S61_NONSQUARE_PREPOWER_NORMALIZED_RAW8_FUSED_PREFLIGHT_20261104_v63'
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
 parent=load(B);assert pins[f'reports/{B}/report.json']=='fbbcd518f798255d93764d872d9acebeec7ad4f75e5ed87a50ac7997c7a16573';worker=load(W);assert pins[f'reports/{W}/report.json']=='5d166d7004ee01a93ff160fc2e521550bf2a861d836f5b35603f1eeb4bf18db2' and worker['all_passed'] and worker['normalized_minor_is_zero']
 assert load('S61_NONSQUARE_PREPOWER_NESTED_RAW8_SOURCE_20261104_v63_REVIEW')['all_passed']
 for p,h in parent['input_sha256'].items():assert sha(R/p)==h
 op=load(G,'operator-and-cost.json');N=op['scalar_nodes'];binding=load(B,'source-binding.json');values=load(B,'native.json');assert values['all_passed']
 roots=[i for row in op['matrix8x8_after_row_a_normalization'][4:] for i in row[:6]];needed=set();todo=roots[:]
 while todo:
  i=todo.pop()
  if i in needed:continue
  needed.add(i);todo.extend(N[i]['args'])
 order=sorted(needed);assert len(order)==190 and binding['minimal_closure190']==order and binding['source_roots24']==roots
 assert binding['source_current15_rows']==[2,4,5,8] and binding['original32_columns']==[22,11,13,25,15,19] and op['row_known_a_factors'][4:]==[0,0,0,0]
 p=R/'reports'/W/'native-input.txt';pins[str(p.relative_to(R))]=sha(p);assert sha(p)==binding['saved_native_input_sha256'];t=iter(p.read_text().split());assert (next(t),next(t))==('8048','2094');raw={}
 for i in range(2094):
  kind=next(t)
  if kind=='x':continue
  upper,den=next(t),next(t)
  if kind=='p':
   shift,unit=next(t),next(t);count=int(next(t));poly=[tuple(next(t) for j in range(4)) for k in range(count)];record=(kind,shift,unit,tuple(poly))
  else:record=(kind,int(next(t)),next(t) if kind=='k' else int(next(t)))
  raw[i]=record
 # Independently bind all46 atom coefficient payloads to source JSON, not only worker bytes.
 q=load('S61_NONSQUARE_CUBIC_Q2_L_UNIT_COMPLETION_20261104_v63','certificate.json');cp=load('S61_NONSQUARE_CUBIC_ZW25_UNITS_20261104_v63','checkpoint.json');pd=load('S61_NONSQUARE_CUBIC_ZW25_PAIR_DUAL_PREFLIGHT_20261104_v63','certificate.json');source_cache={}
 def source_atom(n):
  ref=n['reference']
  if 'constant' in ref:return [[[0,0,0],ref['constant']]]
  if ref['path']=='accepted_Q2L_source_atom':
   d=ref['definition'];assert q['scalar_DAG']['nodes'][ref['old_DAG_root']]==d
   a,b=cp['matrix25x32'][d[1]][d[2]] if d[0]=='checkpoint' else pd['Lframe_dual6x2'][d[1]][1] if d[0]=='accepted_pair_second' else q['scalar_DAG']['new_unit_atoms_num_den'][d[1]]
   assert len(b)==1 and b[0][0][1:]==[0,0] and F(n['QQ_unit'])==F(b[0][1]) and n['a_shift']==-b[0][0][0];return a
  path=ref['path']
  if path not in source_cache:
   p=R/path;pins[path]=sha(p);source_cache[path]=json.loads(p.read_text())
  d=source_cache[path]
  if 'old_DAG_root' in ref:return next(x['Laurent_polynomial'] for x in d['node_values'] if x['old_DAG_root']==ref['old_DAG_root'])
  if 'h9_target' in ref:return d['source_h9'][d['original_target_labels'].index(ref['h9_target'])]
  a=d[ref['key']]
  for i in ref.get('index',[]):a=a[i]
  return a
 # Own byte-token extraction comparison, and normalized2094 source op binding.
 p=R/'reports'/B/'input.txt';pins[str(p.relative_to(R))]=sha(p);tt=iter(p.read_text().split());assert next(tt)=='190';local={old:j for j,old in enumerate(order)};atomcount=0
 for i in order:
  assert int(next(tt))==i;kind=next(tt);rec=raw[i];assert rec[0]==kind
  if kind=='p':
   shift,unit=next(tt),next(tt);count=int(next(tt));poly=tuple(tuple(next(tt) for j in range(4)) for k in range(count));assert (kind,shift,unit,poly)==rec;atomcount+=1
   assert N[i]['op']=='atom' and unit==N[i]['QQ_unit']
   ap=source_atom(N[i]);amin=min((m[0] for m,c in ap),default=0);assert int(shift)==N[i]['a_shift']+amin
   wanted={(m[0]-amin,m[1],m[2]):F(c) for m,c in ap if F(c)};got={(int(a),int(b),int(c)):F(v) for a,b,c,v in poly if F(v)};assert got==wanted
  elif kind=='k':assert (int(next(tt)),next(tt))==(local[rec[1]],rec[2]);assert N[i]['op']=='known_a_shift' and N[i]['args']==[rec[1]] and int(rec[2])==N[i]['a_exponent_shift']
  else:
   assert (int(next(tt)),int(next(tt)))==(local[rec[1]],local[rec[2]]);assert N[i]['op']=={'a':'add','s':'sub','m':'mul'}[kind] and N[i]['args']==[rec[1],rec[2]]
 assert atomcount==46 and [int(next(tt)) for i in range(24)]==[local[i] for i in roots]
 try:next(tt);raise AssertionError('trailing tokens')
 except StopIteration:pass
 # New fixed input retains exact190 records, then appends all24 actual saved polynomial values.
 text=p.read_text();out=[text.rstrip()]
 for row in values['bottom4x6']:
  assert len(row)==6
  for poly in row:
   out.append(str(len(poly)));out.extend(f'{m[0]} {m[1]} {m[2]} {c}' for m,c in poly)
 nativein=O/'fixed-input.txt';nativein.write_text('\n'.join(out)+'\n');cpp=Path(__file__).with_suffix('.cpp');exe=O/'independent-fixed-checker'
 save('phase.json',{'stage':'fixed190_and24_source_labels_bound_before_QQ_products','input_sha256':pins,'no1898_worker_rerun':True})
 flags=subprocess.check_output(['pkg-config','--cflags','--libs','flint'],text=True).split();subprocess.run(['g++','-O2','-std=c++17',str(cpp),'-o',str(exe),*flags],check=True,timeout=20)
 run=subprocess.run([str(exe),str(nativein),str(O/'native-check.json')],capture_output=True,text=True,timeout=max(1,55-(time.monotonic()-start)));(O/'native.log').write_text(run.stdout+run.stderr);assert run.returncode==0
 checked=json.loads((O/'native-check.json').read_text());assert checked['all_passed'] and checked['all90_QQ_two_minors_zero']
 proof={'actual_bottom4x6_source_rows':[2,4,5,8],'actual_pure_target_columns':[22,11,13,25,15,19],'all190_extracted_source_operations_bound_to_original2094':True,'all24_saved_values_independently_QQ_reproduced':True,'all90_two_minors_literal_QQ_zero':True,'bottom_pure_block_rank_at_most1_every_geometric_point':True,'bottom_heavy2_block_rank_at_most2':True,'top4_heavy2_literal_zero':True,'actual8x8_rank_at_most7_so_exact_determinant_zero':True,'normalized_D8_and_raw_D8_both_zero_known_unit_transport':True,'no_new11pivot_chart_exists_for_this_selected_minor':True,'no_claim_for_other17pure_columns_or_other_candidate_or_full57':True}
 save('independent-raw8-zero-proof.json',proof)
 pins[str(cpp.relative_to(R))]=sha(cpp)
 report={'id':S,'all_passed':True,'independent_QQ_fixed190_plus90_products':True,'exact_selected_raw8_minor_zero':True,'source_bottom_pure6_rank_le1':True,'no_11pivot_promotion':True,'no_producer_import_rerun_full1898_sourceworker_ranksearch_point_prime_or_GB':True,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(R/'research/notes'/f'{S}.md'),'evidence_sha256':{n:sha(O/n) for n in ['native-check.json','independent-raw8-zero-proof.json','fixed-input.txt']},'seconds':time.monotonic()-start,'native_seconds':checked['seconds'],'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};assert sum(p.stat().st_size for p in O.iterdir())<33554432;save('report.json',report);save('phase.json',{'stage':'FINAL','all_passed':True});print(json.dumps(report))
if __name__=='__main__':main()
