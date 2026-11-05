#!/usr/bin/env python3
"""Prepared SOURCE continuation ONLY from immutable LIVE07, remaining3 components."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,re,resource,signal,time
from s61_sparse_generic5_source_view_profile_20261105_v64 import Engine,count_numeric,unpack
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
PLAN='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
SCALAR='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_SCALAR_SOURCE_V2_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.profile and a.source_gate_file and a.source_gate_sha and a.root_review_file and a.root_review_sha
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 inputs={};termsaved=[];phases=[];eng=None;budget={'generated_RAW_upper':16384};counter=188
 def load(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def write(name,v,cap=65536):
  b=enc(v);assert len(b)<=cap;assert sum(p.stat().st_size for p in out.iterdir()if p.is_file())+len(b)<=budget['generated_RAW_upper']
  p=out/name;assert not p.exists();p.write_bytes(b);return p
 def checkpoint(e,n):
  if getattr(e,'inside_checkpoint',False):raise RuntimeError('nested SOURCE checkpoint guard')
  e.inside_checkpoint=True
  shot=None
  if sum(p.stat().st_size for p in out.iterdir()if p.is_file())+budget.get('failure_view_snapshot_RAW_upper',0)+262144<=budget['generated_RAW_upper']:
   try:shot=e.failure_view_snapshot(out/'failure-live-views.q5view')
   except Exception:shot={'saved':False,'last_atomic_boundary_and_receipts_retained':True}
  write('failure.json',{'diagnostic_live_view_snapshot':shot,'all_passed':False,'source_only':True,'stop':str(n),'context':e.context,'completed_terms':termsaved,'completed_phases':phases,'resume_old_boundary07_preserved':True,'no_automatic_retry':True},262144)
 def timeout(s,f):
  if eng is not None:checkpoint(eng,'SOURCE math55s; last atomic snapshot retained')
  raise TimeoutError('SOURCE55s; preserve prefix within hard60')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(55)
 def fresh(tag):
  nonlocal counter
  counter+=1;return 'new_'+str(counter)+'_'+tag
 try:
  gp=ROOT/a.source_gate_file;assert sha(gp)==a.source_gate_sha;g=load(gp);assert g['all_passed']
  bp=gp.parent/'resume-views-fused-budget.json';assert sha(bp)==g['evidence_sha256'][bp.name];budget=load(bp);assert sha(Path(__file__))==budget['worker_sha256']
  helper=ROOT/'code/s61_sparse_generic5_source_view_profile_20261105_v64.py';assert sha(helper)==budget['helper_sha256']
  rp=ROOT/a.root_review_file;assert sha(rp)==a.root_review_sha and rp.stat().st_size<=65536
  rootdoc=load(rp);assert rootdoc['all_passed']and count_numeric(rootdoc)<=2000
  for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]=h
  assert budget['fresh_RAW_upper']<=32*1024*1024
  phase=load(ROOT/'reports'/PREV/'phase-07.json');accumulator=phase['accumulator'];assert accumulator=='new_188_sum'
  plan=load(ROOT/'reports'/PLAN/'exact-11-row-coefficient-operator-and-callback-owners.json')
  scalar=load(ROOT/'reports'/SCALAR/'old-live07-and-remaining93-scalar-H-majorants.json')
  assert scalar['max_positive_M_bits_upper']==404 and scalar['max_QQ_field_bits_upper']==80
  planned_majorants={(x['component'],x['term']):x for x in scalar['scalar_and_H_majorants93']}
  old=load(ROOT/'reports'/OLD/'staged-literal-col56-program-and-positive-profiles.json');registry={o['output']:o['binding']for o in old['literal_ops']if o['op']=='input'};del old
  base=25000+count_numeric(phase)+count_numeric(plan)+count_numeric(scalar)+count_numeric(registry)+4096
  assert base<=budget['SOURCE_metadata_base_upper'];eng=Engine(metadata_base=base,checkpoint=checkpoint)
  snap=ROOT/'reports'/PREV/phase['complete_live_snapshot']['file'];assert sha(snap)==phase['complete_live_snapshot']['sha256'];inputs[str(snap.relative_to(ROOT))]=sha(snap)
  with snap.open('rb')as f:
   assert f.read(8)==b'Q5SUPP1\n';n=int.from_bytes(f.read(4),'big');assert n<=64*4096+256
   eng.own(8192);header=json.loads(f.read(n));assert len(header['maps'])==22
   for rec in header['maps']:
    k=rec['id'];q=Q(rec['q']);M=int(rec['M']);eng.fields(q,M);S=set();eng.live[k]={'S':S,'q':q,'M':M,'off':0}
    for _ in range(rec['terms']):
     b=f.read(5);assert len(b)==5;v=int.from_bytes(b,'big');assert v<2**40 and v not in S;eng.own(1+8192);S.add(v)
    p=phase['nodes'][k];assert q==Q(p['known_QQ_content'])and M==int(p['primitive_integer_abs_upper'])and len(S)==p['terms_upper']
    eng.live.pop(k);eng.save(k,S,q,M,{'op':'resume_existing_support_snapshot','snapshot_sha256':sha(snap)})
    for field in['terms_upper','known_QQ_content','primitive_integer_abs_upper','exponent_box5','H_slice_sizes']:assert eng.nodes[k][field]==p[field]
   assert not f.read(1)
  eng.boundary_hash(accumulator,phase['canonical_accumulator_hash'])
  funcs={}
  alltokens={x for row in plan['raw_coefficient_side11']for factors,c in row for x,power in factors if x not in['H','s','q','t']}
  for token in alltokens:
   name=token if token in ['UF10','C10_Ni56','C9_C10_alpha_i']else'leaf_'+token
   if name in eng.live:funcs[token]=name
  assert len(funcs)==21 and set(eng.live)==set(funcs.values())|{accumulator}
  write('resume-live07-proof.json',{'all_passed':True,'source_only':True,'snapshot':str(snap.relative_to(ROOT)),'total_keys':44555,'maps':22,'old_accumulator_hash_exact':True,'counter':188,'no_old8_source_arithmetic_replayed':True})
  eng.nodes={k:v for k,v in eng.nodes.items()if k in eng.live};eng.ops=[];eng.meta=sum(count_numeric(v)+12 for v in eng.nodes.values())
  def function(token):
   if token in funcs:return funcs[token]
   name=token
   if token.endswith('i')and token[:-1]in['U','V','C','F','N','M']:name=token[:-1]+'6'
   elif token.startswith('Ai_'):name='A6_'+token[3:]
   binding=registry[name];eng.own(8500+2048);eng.base+=8500
   try:
    z=load(ROOT/'reports'/OLD/binding['scope']);assert count_numeric(z)<=8500
    z=z['common_operands'][binding['label']]if 'common_operands'in z else z['row_operands'][binding['label']]
    p=z['original_literal'];assert len({tuple(t[1:])for t in p['terms']})==len(p['terms'])
    node='leaf_'+token;eng.leaf(node,p,source={'old_source_node':name,'binding':binding},scalar=Q(z['QQ_transform_scalar']),shift=tuple(z['transform_monomial_exponent5']));funcs[token]=node;return node
   finally:eng.base-=8500
  # Three exact old input literals; their payload hashes remain tied to E.
  eng.own(12000+2048)
  E=load(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']
  E={(x['row'],x['column']):x for x in E if x['row']in[11,44,52]and x['column']==53}
  for component in[8,9,10]:
   row=plan['raw_row_ids11'][component];r=E[row,53];ref=r['source'];p=ROOT/ref['file']
   if str(p.relative_to(ROOT))not in inputs:inputs[str(p.relative_to(ROOT))]=sha(p)
   with p.open('rb')as f:f.seek(ref['offset']);b=f.read(ref['bytes'])
   assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];eng.own(8192);eng.leaf('raw53',decode(b),source=ref)
   for ti,(factors,coefficient)in enumerate(plan['raw_coefficient_side11'][component]):
    shift=[0]*5;names=[]
    for token,power in factors:
     if token in['H','s','q','t']:shift[{'H':0,'s':2,'q':3,'t':4}[token]]+=power
     else:assert power>=0;names.extend([function(token)]*power)
    names.sort(key=lambda k:len(eng.live[k]['S']))
    if names:
     node=eng.copy(fresh('term'),names[0])
     for x in names[1:]:
      nxt=eng.mul(fresh('term_product'),node,x);eng.drop(node);node=nxt
     nxt=eng.copy(fresh('term_scaled'),node,Q(coefficient),tuple(shift));eng.drop(node);node=nxt
    else:node=eng.leaf(fresh('constant'),{'content':str(Q(coefficient)),'shift':shift,'terms':[['1',0,0,0,0,0]]},source={'QQ_monomial':True})
    dot=fresh('dot');nxt=fresh('sum');eng.fused_union_product(nxt,accumulator,node,'raw53',dot);eng.drop(node);accumulator=nxt
    m=planned_majorants[component,ti];v=eng.nodes[accumulator]
    assert Q(v['known_QQ_content'])==Q(m['q_exact'])and int(v['primitive_integer_abs_upper'])<=int(m['positive_M_upper'])
    assert v['exponent_box5'][0][0]>=m['H_BOX_upper'][0]and v['exponent_box5'][0][1]<=m['H_BOX_upper'][1]
    newids={o['output']for o in eng.ops};record={'source_only':True,'component':component,'term':ti,'nodes':{k:eng.nodes[k]for k in newids},'literal_ops':eng.ops,'accumulator':accumulator,'scalar_H_majorants_verified':True}
    f=write(f'term-{component}-{ti:02d}.json',record,budget['term_receipt_RAW_upper']);termsaved.append({'component':component,'term':ti,'file':f.name,'sha256':sha(f)})
    eng.nodes={k:v for k,v in eng.nodes.items()if k in eng.live};eng.ops=[];eng.meta=sum(count_numeric(v)+12 for v in eng.nodes.values());eng.own()
    assert time.monotonic()-start<52
   eng.drop('raw53');assert set(eng.live)==set(funcs.values())|{accumulator}
   h=eng.boundary_hash(accumulator);assert len(eng.live)<=39
   generated=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
   assert generated+budget['one_boundary_snapshot_RAW_upper']+budget['phase_receipt_RAW_upper']+262144<=budget['generated_RAW_upper']
   shot=eng.snapshot(out/f'live-{component:02d}.q5supports');assert shot['bytes']<=budget['one_boundary_snapshot_RAW_upper']
   p=write(f'phase-{component:02d}.json',{'source_only':True,'accumulator':accumulator,'profile':eng.nodes[accumulator],'boundary_hash':h,'all_live_snapshot':shot},budget['phase_receipt_RAW_upper'])
   phases.append({'component':component,'file':p.name,'sha256':sha(p),'snapshot':shot});gc.collect()
  final=eng.nodes[accumulator];raw=sum((ROOT/p).stat().st_size for p in inputs)+sum(p.stat().st_size for p in out.iterdir()if p.is_file())+65536
  assert raw<=budget['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'source_only':True,'old8_not_replayed':True,'resume44555keys_verified':True,'remaining3_components_complete':True,'93_terms_complete':len(termsaved)==93,'immutable_views_no_copy_of_support_sets':True,'fused_SOURCE_union_no_dot_set':True,'SOURCE_owned_peak':eng.peak,'literal_pair_majorant':eng.pair_count,'full_row6_col53_positive_support':final,'phases':phases,'all_three_new_snapshots_roundtrip':True,'actual_READ_output_RAW_upper':raw,'all167x14_87plus66_RHS_oneP290_F9_F10_complements_kept':True,'no_VALUE_or_all14_row_claim':True,'input_sha256':inputs,'seconds':time.monotonic()-start};write('report.json',report,262144);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','phases','full_row6_col53_positive_support')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():write('failure.json',{'all_passed':False,'source_only':True,'error':repr(exc),'completed_terms':termsaved,'completed_phases':phases,'no_automatic_retry':True},262144)
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile',action='store_true');p.add_argument('--source-gate-file');p.add_argument('--source-gate-sha');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');main(p.parse_args())
