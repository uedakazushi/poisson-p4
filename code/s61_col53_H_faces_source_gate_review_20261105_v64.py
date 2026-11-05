#!/usr/bin/env python3
"""Independent FULL2 H-face support/plans/integer-code/whole-budget audit."""
import ast,hashlib,itertools,json,math,pathlib,re,resource,signal,time
from fractions import Fraction as Q
from collections import Counter
from ambient_B18_C_genuine_projection7_source_review_20261105_v63 import count,sha
from s61_col53_vertex_source_values_review_20261105_v64 import decoder,literal_record
ROOT=pathlib.Path(__file__).resolve().parents[1]
SRC='S61_COL53_H_FACES_SOURCE_20261105_v64';GATE='S61_COL53_H_FACES_VALUE_GATE_20261105_v64'
ID='S61_COL53_H_FACES_SOURCE_CODE_20261105_v64_REVIEW'
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def legacy(x):
 if isinstance(x,bool):return 0
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',x) else 0
 if isinstance(x,(list,tuple)):return sum(map(legacy,x))
 if isinstance(x,dict):return sum(map(legacy,x.values()))
 return 0
def main():
 started=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2<<30,2<<30))
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pins={};RAW=0
 def pin(p,h=None,parse=True):
  nonlocal RAW
  p=pathlib.Path(p);k=str(p.relative_to(ROOT))
  if k not in pins:pins[k]=sha(p);RAW+=p.stat().st_size
  if h:assert pins[k]==h,(k,pins[k],h)
  return json.loads(p.read_text()) if parse else None
 sr=pin(ROOT/'reports'/SRC/'report.json','5e66838fb845def468a9a1c883fd08a5836f21fdfd3a88de8f4d62109013d718')
 gr=pin(ROOT/'reports'/GATE/'report.json','fe261a87001601168e3f52eae66fe4bd5e1e80ddda810836de7b656736ac09fe')
 for report in [sr,gr]:
  for p,h in report['input_sha256'].items():pin(ROOT/p,h,False)
 def doc(st,n,r):return pin(ROOT/'reports'/st/n,r['evidence_sha256'][n])
 D=doc(SRC,'compact-face-graph.json',sr);P=doc(SRC,'positive-face-profiles-and-VALUE-lifetimes.json',sr)
 F=doc(SRC,'two-complete-positive-H-faces.json',sr);B=doc(GATE,'whole-integer-face-VALUE-budget.json',gr)
 proofmeta=doc(GATE,'actual-code-lifetime-and-codec-source-proof.json',gr);helpdoc=doc(GATE,'worker-help-receipt.json',gr)
 nodes=D['nodes'];N=len(nodes);root=D['root'];assert N==1203 and root==1202
 # Exact previous original-source SSA, independently accepted, with only
 # in-place vertex->H compaction and the NEW integer metadata added.
 vp=ROOT/'reports/S61_COL53_VERTEX_SOURCE_GATE_V2_20261105_v64';vr=pin(vp/'report.json','6ad0a98af6dd6539bf2d9cdf16ec8e431ee70df66dbafef540c6aac32c73a1b8')
 VG=pin(vp/'graph.json',vr['evidence_sha256']['graph.json'])
 accepted=pin(ROOT/'reports/S61_COL53_VERTEX_SOURCE_VALUES_V2_20261105_v64_REVIEW/report.json','d7e8134bbed2dc1483bef0933d9c892223222583619f72583d564ea2dbb66831');assert accepted['all_passed']
 compaction_peak=count([sr,gr,D,P,F,B,proofmeta,helpdoc,vr,VG,accepted,pins])+15000;assert compaction_peak<150000
 maxq=0;maxM=0;weights_count=0
 for k,(old,n) in enumerate(zip(VG['nodes'],nodes)):
  assert {x:y for x,y in n.items() if x not in ['H_extrema','integer_weights','actual_zero_by_source_q','leaf_integer_ratio']}=={x:y for x,y in old.items() if x!='leads'}
  assert n['H_extrema']==[None if old['leads'][i] is None else old['leads'][i][0] for i in [0,1]]
  q=Q(n['source_profile']['known_QQ_content']);M=int(n['source_profile']['primitive_integer_abs_upper'])
  assert n['actual_zero_by_source_q']==(not q) and all(i<k for i in n['args'])
  maxq=max(maxq,abs(q.numerator).bit_length(),q.denominator.bit_length());maxM=max(maxM,M.bit_length())
  if not q:assert n['integer_weights']==[];continue
  args=n['args'];kind=n['op'];qq=[Q(nodes[i]['source_profile']['known_QQ_content']) for i in args]
  if kind=='input':w=[]
  elif kind in ['copy','immutable_SOURCE_view']:w=[qq[0]*Q(n['scalar'])/q]
  elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:w=[qq[0]*qq[1]*Q(n['scalar'])/q]
  else:w=[x*(Q(n.get('scalarB','1')) if j else 1)/q for j,x in enumerate(qq)]
  assert all(x.denominator==1 for x in w) and list(map(str,(x.numerator for x in w)))==n['integer_weights'];weights_count+=len(w)
 assert maxq==92 and maxM==401 and weights_count==876
 del VG
 # Guard independent reader before full support arithmetic; all current
 # structural metadata plus bounded live5tuples/leafdecode reserve.
 own_base=count([sr,gr,D,P,F,B,proofmeta,helpdoc,vr,accepted,pins])+18000
 assert own_base+5*6000+20000<150000
 binding_terms=max(n['literal_profile']['terms_upper'] for n in nodes if n['op']=='input' and 'binding' in n['source'])
 binary_terms=max(n['literal_profile']['terms_upper'] for n in nodes if n['op']=='input' and 'file' in n['source'])
 assert max(19*binding_terms+1024,11*binary_terms+128)<=20000
 save(out/'before-independent-full-positive-face-supports.json',{'SOURCE_only_no_signed_faces':True,'metadata_and_controller':own_base,'live_support_reserve_5x6000':30000,'literal_decode_reserve':20000,'whole_bound':own_base+50000,'hard_numeric':150000})
 def input_support(n):
  s=n['source'];p=n['literal_profile']
  if 'binding' in s:
   b=s['binding'];path=ROOT/OLD/b['scope'];pin(path,vr['input_sha256'][str(path.relative_to(ROOT))],False)
   rec=literal_record(path,b['label']);lit=rec['original_literal'];content=Q(lit['content']);offset=lit['shift'];terms=[(int(t[0]),tuple(t[1:])) for t in lit['terms']]
   scalar=Q(rec['QQ_transform_scalar']);shift=rec['transform_monomial_exponent5']
  else:
   path=ROOT/s['file'];pin(path,vr['input_sha256'][s['file']],False)
   with path.open('rb') as f:
    if 'offset' in s:f.seek(s['offset']);raw=f.read(s['bytes']);assert hashlib.sha256(raw).hexdigest()==s['payload_sha256']
    else:raw=f.read();assert raw[:9]==b'Q5CELL1\n\x01' and pins[s['file']]==s['sha256'];raw=raw[9:]
   content,offset,terms=decoder(raw);scalar=Q(1);shift=[0]*5
   if 'offset' not in s:
    low=[offset[j]+min(e[j] for c,e in terms) for j in range(5)];shift=[p['exponent_box5'][j][0]-low[j] for j in range(5)]
    scalar=Q(p['known_QQ_content'])/content
    assert (scalar,shift)==((Q(353894400000),[1,0,0,0,1]) if 'col56-ninth-cleared' in s['file'] else (Q(1),[0]*5))
  q=Q(n['source_profile']['known_QQ_content'])
  if q:
   ratio=content*scalar/q;assert ratio.denominator==1 and str(ratio.numerator)==n['leaf_integer_ratio']
   assert max((abs(c)*abs(ratio.numerator) for c,e in terms),default=0)<=int(n['source_profile']['primitive_integer_abs_upper'])
  exponents=[tuple(offset[j]+e[j]+shift[j] for j in range(5)) for c,e in terms]
  assert len(exponents)==len(set(exponents));return exponents
 roots=[];source_pairs=0;source_peak=own_base;allcounts=[];required_modes=[];forecast=[];peakvalue=B['VALUE_metadata_base'];maxterms=maxmaps=maxmap=0;product_bounds=0
 assert [legacy(D),legacy(P),legacy(F)]==[28362,6681,4904]
 assert B['VALUE_metadata_parts']=={'runtime_and_python_code_reserve':25000,'loaded_compact_graph':28362,'loaded_SOURCE_profiles':6681,'loaded_both_root_faces':4904,'indices_Mlimits_required_active_counts_uses_cache_bookkeeping':14436,'root_allowed_5tuple_and_index_duplicate':4596,'ROOT_document_numeric_upper':2000,'gate_budget_new_receipts_numeric_reserve':4000,'failure_numeric_reserve':8192}
 assert sum(B['VALUE_metadata_parts'].values())==98171==B['VALUE_metadata_base']
 for mode,H in enumerate([-3,21]):
  plan=D['mode_plans'][mode];needed=set();active={};todo=[root]
  while todo:
   k=todo.pop()
   if k in needed:continue
   needed.add(k);n=nodes[k];h=n['H_extrema'][mode];args=n['args'];kind=n['op']
   if h is None or kind=='input':child=[]
   elif kind in ['copy','immutable_SOURCE_view']:
    child=args;assert h==nodes[child[0]]['H_extrema'][mode]+n['shift'][0]
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
    child=args;assert h==sum(nodes[i]['H_extrema'][mode] for i in child)+n['shift'][0]
   else:
    assert n['shift']==[0]*5;child=[i for i in args if nodes[i]['H_extrema'][mode]==h]
   active[k]=child;todo.extend(child)
  assert sorted(needed)==plan['SOURCE_needed'] and active=={int(k):v for k,v in plan['active_tied_args'].items()}
  uses=Counter(i for args in active.values() for i in args);uses[root]+=1;live={};counts={}
  for k in sorted(needed):
   n=nodes[k];h=n['H_extrema'][mode];kind=n['op'];args=active[k];support=set();resident=sum(map(len,live.values()))
   if h is None:pass
   elif kind=='input':
    support={e for e in input_support(n) if e[0]==h}
   elif kind in ['copy','immutable_SOURCE_view']:
    support={tuple(e[j]+n['shift'][j] for j in range(5)) for e in live[args[0]]}
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
    source_pairs+=len(live[args[0]])*len(live[args[1]])
    for e in live[args[0]]:
     for f in live[args[1]]:support.add(tuple(e[j]+f[j]+n['shift'][j] for j in range(5)))
   else:
    for child in args:support.update(live[child])
   assert all(len(e)==5 and e[0]==h for e in support)
   source_peak=max(source_peak,own_base+5*(resident+len(support)+sum(map(len,roots)))+6*len(support)+20000);assert source_peak<150000
   counts[k]=len(support);assert counts[k]==P['profiles'][mode]['face_terms_by_node'][str(k)]
   digest=hashlib.sha256((json.dumps([list(e) for e in sorted(support)],separators=(',',':'))+'\n').encode()).hexdigest()
   assert digest==P['profiles'][mode]['face_hash_by_node'][str(k)]
   live[k]=support
   for child in args:
    uses[child]-=1
    if not uses[child]:del live[child]
   if not uses[k]:del live[k]
  rootface=live[root];assert rootface=={tuple(e) for e in F['faces'][mode]['full5_exponents']} and len(rootface)==[766,214][mode]
  roots.append(rootface);allcounts.append(counts)
  # VALUEclosure: ONLY knownsourceqzero additionally pruned, never a
  # vanished vertex coefficient, actual support or unknown cancellation.
  req=set();todo=[root]
  while todo:
   k=todo.pop()
   if k in req:continue
   req.add(k)
   if not nodes[k]['actual_zero_by_source_q']:todo.extend(active[k])
  assert sorted(req)==plan['VALUE_needed'];required_modes.append(req)
  uses=Counter(child for k in req if not nodes[k]['actual_zero_by_source_q'] for child in active[k]);uses[root]+=1;live={};rows=[];modepeak=98171
  for k in sorted(req):
   n=nodes[k];sz=0 if n['actual_zero_by_source_q'] else counts[k];old=sum(live.values());decode=16152 if n['op']=='input' and sz else 0
   whole=98171+7*(old+sz)+decode+2048;modepeak=max(modepeak,whole);peakvalue=max(peakvalue,whole);maxterms=max(maxterms,old+sz);maxmaps=max(maxmaps,len(live)+1);maxmap=max(maxmap,sz)
   rows.append({'node':k,'existing_live_terms':old,'new_output_before_release':sz,'leaf_decode_numeric_upper':decode,'whole_VALUE_numeric_upper':whole})
   if not n['actual_zero_by_source_q']:
    M=int(n['source_profile']['primitive_integer_abs_upper']);children=active[k];ww=list(map(int,n['integer_weights']));kind=n['op']
    if kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
     a,b=children;assert min(counts[a],counts[b])*int(nodes[a]['source_profile']['primitive_integer_abs_upper'])*int(nodes[b]['source_profile']['primitive_integer_abs_upper'])*abs(ww[0])<=M;product_bounds+=1
    elif kind in ['copy','immutable_SOURCE_view']:assert abs(ww[0])*int(nodes[children[0]]['source_profile']['primitive_integer_abs_upper'])<=M
    elif kind!='input':assert sum(abs(w)*int(nodes[c]['source_profile']['primitive_integer_abs_upper']) for c,w in zip(n['args'],ww))<=M
   live[k]=sz
   if not n['actual_zero_by_source_q']:
    for child in active[k]:
     uses[child]-=1
     if not uses[child]:live.pop(child)
   if not uses[k]:live.pop(k,None)
  codec=98171+10*counts[root]+4096;peakvalue=max(peakvalue,codec)
  forecast.append({'mode':mode,'H':H,'schedule':rows,'max_normal_owned':modepeak,'codec_owned':codec})
 assert forecast==proofmeta['reconstructed_forecasts'] and peakvalue==126262 and (maxterms,maxmaps,maxmap)==(2912,26,766)
 assert product_bounds==proofmeta['positive_coefficient_majorant_product_checks']==388
 # The final FULLfive-coordinate snapshot proof is also checked directly.
 endpoint=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json';end=pin(endpoint,vr['input_sha256'][str(endpoint.relative_to(ROOT))])
 snap=endpoint.parent/end['all_live_snapshot']['file'];pin(snap,vr['input_sha256'][str(snap.relative_to(ROOT))],False);observed=[set(),set()]
 with snap.open('rb') as f:
  assert f.read(8)==b'Q5SUPP1\n';header=json.loads(f.read(int.from_bytes(f.read(4),'big')))
  for rec in header['maps']:
   for _ in range(rec['terms']):
    b=f.read(5);assert len(b)==5
    if rec['id']==end['accumulator']:
     word=int.from_bytes(b,'big');e=tuple(((word>>(8*j))&255)-64 for j in range(5))
     if e[0] in [-3,21]:observed[int(e[0]==21)].add(e)
  assert not f.read(1)
 assert observed==roots
 # Literal code, sourcepins, failure/codec and unchanged hard limits.
 worker=ROOT/'code/s61_col53_H_faces_values_20261105_v64.py';pin(worker,'d7e9a87cda9384a31f4899e887d2a3816e4530fa815d94786708b6abf0934357',False)
 code=worker.read_text();functions={n.name:ast.get_source_segment(code,n) for n in ast.walk(ast.parse(code)) if isinstance(n,ast.FunctionDef)}
 assert '7*(sum(map(len,cache.values()))+len(work)' in functions['own'] and "transient+extra" in functions['own']
 assert functions['commit'].index('own(newterm=True)')<functions['commit'].index('work[e]=v')
 assert "actual_zero_by_source_q" in functions['main'] and 'S61_COL53_VERTEX_VALUES' not in code
 assert 'cache[k]=work;work={};own()' in functions['main'] and "all(e[0]==target" in functions['main']
 assert 'META-=reserve' in functions['checkpoint_maps'] and 'failure-cache-index.jsonl' in functions['checkpoint_maps']
 assert "if work and 'node'in context and context['node']not in cache" in functions['checkpoint_maps']
 assert "signal.alarm(55)" in functions['main'] and 'time.monotonic()-start<50' in functions['main']
 assert 'len(e)==5' in functions['commit'] and 'abs(v).bit_length()<=640' in functions['commit']
 assert 'abs(v).bit_length()<=1921' in functions['rawtemp']
 codecpath=ROOT/'code/s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64.py';pin(codecpath,gr['input_sha256'][str(codecpath.relative_to(ROOT))],False)
 cc=codecpath.read_text();cf={n.name:ast.get_source_segment(cc,n) for n in ast.walk(ast.parse(cc)) if isinstance(n,ast.FunctionDef)}
 assert "MAGIC=b'Q5CELL1\\n'" in cc
 assert 'keys=descending_alias_keys(values)' in cf['write_cell'] and 'tmp.rename(path)' in cf['write_cell']
 assert "for j in range(5)" in cf['write_cell'] and 'values.get(e)==c*G' in cf['validate_cell']
 assert 'gc==1 and not f.read(1)' in cf['validate_cell'] and "not values" in cf['validate_cell']
 assert B['normal_commit_bits']==640 and B['normal_temp_bits']==1921 and B['stored_VALUE_bits_unchanged']==2048
 assert maxq+2*640==B['codec_naive_numerator_temp_bits']==1372<1921
 assert 3*766+4096+512==6906<=8192 and B['failure_numeric_reserve']==8192
 assert B['face_binary_RAW_upper']==[512+128*766,512+128*214] and B['generated_RAW_upper']==2<<20
 assert B['fresh_RAW_upper']==17661238<32<<20 and B['whole_VALUE_numeric_upper']==126262<150000
 assert helpdoc['exit_code']==0 and helpdoc['no_VALUE_execution']
 assert not (ROOT/'reports/S61_COL53_H_FACES_VALUES_20261105_v64').exists()
 assert all(sha(ROOT/p)==h for p,h in pins.items())
 proof={'all1203_originalSSA_and_H_compaction_checked':True,'all_integer_q_rescaling_identities':876,
  'all_SOURCE_modeplans_and_positive_full5_supports_regenerated':True,'root_positive_face_terms':[766,214],
  'all_mode_node_support_hashes_and_full_snapshot_faces_match':True,
  'all_VALUE_qzero_only_closures_and_lastuse_lifetimes':True,'integer_coefficient_majorant_products_checked':388,
  'face_projection_proof':'ForfiniteLaurent support, productHextreme=hA+hB andtheWHOLEfacepolynomialistheproductoftheWHOLEextremefaces (plusliteralshift/scalar). Atsumextreme onlyH-tiedchildrencontribute. Sourceq0 permitsactualzero pruning; vanishedlexvertexcoefficients doNOTpermit anypruning. ExtremesremainSOURCEextremes whenactualfacescancel.',
  'full5_coordinate_codec_and_zero_nonzero_payload_semantics_checked':True,
  'worker_after_operation_checks_safe_by_prior_integer_majorants':'Forintegercopy/product/sum,sourceM andexactintegerweightsboundeveryactualterm andpartialsum beforeexecution. Standardtemporarybound404 andcodecnaivebound1372 areundercommitted1921/hard2048. No newunknowninverse orpoint isintroduced.',
  'prospective_each_map_charge7_and_move_without_duplicate':True,
  'normal_VALUE_META':98171,'normal_whole_VALUE_upper':126262,'max_complete_live_and_work_terms':2912,'max_maps':26,
  'failure_reserve_and_atomic_prefix_preservation':True,'failure_reserve':8192,'codec_extra_upper':6906,
  'fresh_RAW_bound':17661238,'generated_RAW_cap':2<<20,'math_signal_hard_seconds':[50,55,60],
  'ROOT_report_cap_bytes_numeric':[65536,2000],
  'own_SOURCE_support_peak':source_peak,'own_compaction_peak':compaction_peak,
  'no_signed_face_coefficients_evaluated':True,'parent_conditional_ONE_twofaces_readiness':True,
  'scope':'TwofixedFULLHfacesofrow6/col53clearednumeratoronly. Bothzeroisnotentryzero; anynonzeroimpliesentrynonzeroonly. No wholeentry/rank/fullinverse/exhaustion,parameterpoint,primeornewtargetjob.'}
 ep=out/'independent-full-two-Hface-support-mode-code-and-resource-proof.json';save(ep,proof)
 report={'id':ID,'all_passed':True,'frozen':True,'SOURCE_CODE_only':True,'input_sha256':pins,'code_sha256':sha(__file__),
  'evidence_sha256':{ep.name:sha(ep)},'parent_conditional_ONE_twofaces_readiness':True,'actual_signed_face_values_evaluated':False,'seconds':time.monotonic()-started}
 save(out/'report.json',report);signal.alarm(0);print(json.dumps({'all_passed':True,'sha256':sha(out/'report.json'),'own_source_peak':source_peak,'seconds':report['seconds']}))
if __name__=='__main__':main()
