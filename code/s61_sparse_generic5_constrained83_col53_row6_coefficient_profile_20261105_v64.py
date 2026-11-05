#!/usr/bin/env python3
"""Prepared changed-association SOURCE profile, fixed first row6 and col53 only."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,gc,hashlib,json,resource,signal,time
from s61_sparse_generic5_source_packed5_profile_20261105_v64 import Engine,count_numeric
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
PLAN='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def main(a):
 assert a.profile
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={};phasefiles=[];counter=0;eng=None
 def checkpoint(e,n):
  d={'all_passed':False,'source_only':True,'forecast':n,'limit':150000,'context':e.context,'completed_phase_files':phasefiles,'no_new_matrix_coefficients':True,'no_automatic_retry':True}
  (out/'failure.json').write_bytes(enc(d));(out/'partial-profile.json').write_bytes(enc({'nodes':e.nodes,'ops':e.ops,'input_sha256':inputs}))
 def timeout(s,f):
  if eng is not None:checkpoint(eng,'hard60s')
  raise TimeoutError('SOURCE60s; no automatic retry')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
 def load(p):
  p=Path(p);inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def write(name,v):
  b=enc(v);total=sum(p.stat().st_size for p in out.iterdir()if p.is_file())
  assert total+len(b)<=8*1024*1024
  p=out/name;assert not p.exists();p.write_bytes(b);return p
 def fresh(tag):
  nonlocal counter
  counter+=1;return 'new_'+str(counter)+'_'+tag
 try:
  proof=load(ROOT/'reports'/PLAN/'exact-11-row-coefficient-operator-and-callback-owners.json')
  assert proof['all_passed']and proof['stage8_vector_identity_exact']
  old=load(ROOT/'reports'/OLD/'staged-literal-col56-program-and-positive-profiles.json')
  assert old['all_passed']
  byname={o['output']:o for o in old['literal_ops']};all_input={k:o for k,o in byname.items()if o['op']=='input'}
  target=['row6_N8_sum_m','row58_N8_sum_m'];needed=set()
  def visit(k):
   if k in needed:return
   needed.add(k);o=byname[k]
   for key in('A','B'):
    if key in o:visit(o[key])
  for k in target:visit(k)
  subset=[o for o in old['literal_ops']if o['output']in needed]
  expected={k:old['positive_support_and_height_profiles'][k]for k in needed}
  c8='S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64'
  a8=load(ROOT/'reports'/c8/'actual-cleared-col82-and-commonQ3RE-representation.json')
  b8={r['row']:r for r in a8['actual16_col82_numerator_values']}
  c9='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64'
  a9=load(ROOT/'reports'/c9/'actual15-col56-and-once-shared-row35.json')
  n9={r['row']:r for r in a9['actual15_ninth_col56_numerators']}
  Edata=load(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/actual-fixed24-entry-source-and-exponent-plan.json')
  E={(r['row'],r['column']):r for r in Edata['actual_entry_sources576']};del Edata
  allowedrows={6 if r=='i'else r for r in proof['raw_row_ids11']}
  E={k:v for k,v in E.items()if k[0]in allowedrows and k[1]==53}
  parse_numeric=count_numeric(old)+count_numeric(proof)+count_numeric(subset)+count_numeric(expected)+count_numeric(all_input)+count_numeric(a8)+count_numeric(a9)+count_numeric(E)+4096+20000
  assert parse_numeric<150000
  del old,byname;gc.collect()
  controls=count_numeric(proof)+count_numeric(subset)+count_numeric(expected)+count_numeric(all_input)+count_numeric(a8)+count_numeric(a9)+count_numeric(E)+4096
  # SOURCE only: these are metadata and existing one-leaf parse allowances,
  # not packed VALUE coefficient counts.
  base=25000+controls;assert base<=50000
  planned=set(ROOT/p for p in inputs)
  for o in all_input.values():planned.add(ROOT/'reports'/OLD/o['binding']['scope'])
  for row in(6,58):
   planned.add(ROOT/'reports'/c8/b8[row]['file']);planned.add(ROOT/'reports'/c9/n9[row]['file'])
  for v in E.values():
   if v['terms']:planned.add(ROOT/v['source']['file'])
  for p in[Path(__file__),ROOT/'code/s61_sparse_generic5_source_packed5_profile_20261105_v64.py',ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py',ROOT/f'research/notes/{ID}.md']:
   inputs[str(p.relative_to(ROOT))]=sha(p);planned.add(p)
  planned_RAW=sum(p.stat().st_size for p in planned);assert planned_RAW<24*1024*1024
  write('SOURCE-plan-budget-before-support.json',{'source_only':True,'cold_parse_numeric':parse_numeric,'metadata_base':base,'leaf_document_transient':8500,'cached_binary_decode_transient':20000,'SOURCE_owned_hard':150000,'planned_ALL_READ_RAW':planned_RAW,'all_new_versions_RAW_upper':8*1024*1024,'RAW_total_hard':32*1024*1024,'VALUE_no_GO':True,'fixed_row6_col53':True})
  eng=Engine(metadata_base=base,checkpoint=checkpoint)
  uses={k:0 for k in needed}
  for o in subset:
   for key in('A','B'):
    if key in o:uses[o[key]]+=1
  for k in target:uses[k]+=1
  def oldoperand(o):
   b=o['binding'];p=ROOT/'reports'/OLD/b['scope'];z=load(p);assert count_numeric(z)<=8500
   return z['common_operands'][b['label']]if 'common_operands'in z else z['row_operands'][b['label']]
  def oldleaf(name,k):
   eng.own(8500+2048);eng.base+=8500
   try:
    z=oldoperand(all_input[k]);p=z['original_literal']
    return eng.leaf(name,p,source={'old_program_node':k,'binding':all_input[k]['binding']},scalar=Q(z['QQ_transform_scalar']),shift=tuple(z['transform_monomial_exponent5']))
   finally:eng.base-=8500
  for o in subset:
   k=o['output'];name='O_'+k
   if o['op']=='input':oldleaf(name,k)
   elif o['op']=='copy':eng.copy(name,'O_'+o['A'])
   elif o['op']=='multiply':eng.mul(name,'O_'+o['A'],'O_'+o['B'],Q(o['scalar']),tuple(o['shift5']))
   else:assert o['op']=='add';eng.add(name,'O_'+o['A'],'O_'+o['B'],Q(o.get('scalarB','1')))
   assert eng.nodes[name]['support_sha256']==expected[k]['support_sha256']
   assert eng.nodes[name]['known_QQ_content']==expected[k]['known_QQ_content']
   for key in('A','B'):
    if key in o:
     uses[o[key]]-=1
     if not uses[o[key]]:eng.drop('O_'+o[key])
  # The callback values are not computed; only their exact source supports.
  def cache_leaf(name,stem,rec,scalar=Q(1),shift=(0,0,0,0,0)):
   p=ROOT/'reports'/stem/rec['file'];assert sha(p)==rec['sha256'];inputs[str(p.relative_to(ROOT))]=rec['sha256'];eng.own(20000+2048);eng.base+=20000
   try:
    b=p.read_bytes();poly=decode(b[9:]);assert count_numeric(poly)<=20000
    return eng.leaf(name,poly,source={'file':str(p.relative_to(ROOT)),'sha256':rec['sha256']},scalar=scalar,shift=shift)
   finally:eng.base-=20000
  cache_leaf('B58',c8,b8[58]);cache_leaf('B6',c8,b8[6])
  A=eng.mul('alphaA','B58','O_'+target[0]);B=eng.mul('alphaB','B6','O_'+target[1],Q(-1));eng.add('alpha',A,B,consume=True);eng.drop('B58','B6',*['O_'+k for k in target])
  # Existing cached readbacks, not newly evaluated parameter coefficients.
  cache_leaf('UF10',c9,n9[58],Q(353894400000),(1,0,0,0,1))
  cache_leaf('C10_Ni56',c9,n9[6],Q(353894400000),(1,0,0,0,1))
  eng.copy('C9_C10_alpha_i','alpha',Q(62208000000*353894400000),(1,0,1,1,3));eng.drop('alpha')
  # Resolve generic function tokens to already frozen leaf source nodes.
  funcs={k:k for k in ['UF10','C10_Ni56','C9_C10_alpha_i']}
  def function(token):
   if token in funcs:return funcs[token]
   k=token
   if token[-1:]=='i'and token[:-1]in['U','V','C','F','N','M']:k=token[:-1]+'6'
   elif token.startswith('Ai_'):k='A6_'+token[3:]
   name='leaf_'+token
   assert k in all_input,(token,k)
   oldleaf(name,k);funcs[token]=name;return name
  # Read only the selected raw source entry per one dot-product component.
  actual_source_data=count_numeric(proof)+count_numeric(subset)+count_numeric(expected)+count_numeric(all_input)+count_numeric(a8)+count_numeric(a9)+count_numeric(E)+4096
  assert 25000+actual_source_data<=base+20000
  source_READ_max=sum((ROOT/p).stat().st_size for p in inputs)
  rawfiles={ROOT/E[(6 if row=='i'else row),53]['source']['file']for row in proof['raw_row_ids11']if E[(6 if row=='i'else row),53]['terms']}
  planned_RAW=source_READ_max+sum(p.stat().st_size for p in rawfiles if str(p.relative_to(ROOT))not in inputs)
  assert planned_RAW<24*1024*1024
  write('SOURCE-preflight.json',{'source_only':True,'fixed_row':6,'fixed_column':53,'parse_numeric_upper':parse_numeric,'engine_metadata_base':base,'dynamic_source_own_limit':150000,'planned_READ_RAW':planned_RAW,'all_outputs_RAW_limit':8*1024*1024,'total_RAW_limit':32*1024*1024,'no_coefficient_worker_GO':True})
  accumulator=None
  for position,(row,terms)in enumerate(zip(proof['raw_row_ids11'],proof['raw_coefficient_side11'])):
   row=6 if row=='i'else row;ent=E[row,53]
   if ent['terms']:
    ref=ent['source'];p=ROOT/ref['file']
    if str(p.relative_to(ROOT))not in inputs:inputs[str(p.relative_to(ROOT))]=sha(p)
    with p.open('rb')as f:f.seek(ref['offset']);b=f.read(ref['bytes'])
    assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];eng.leaf('raw53',decode(b),source=ref)
   else:eng.leaf('raw53',{'content':'0','shift':[0]*5,'terms':[]},source={'known_source_zero':True})
   for factors,coefficient in terms:
    shift=[0]*5;names=[]
    for token,power in factors:
     if token in['H','s','q','t']:shift[{'H':0,'s':2,'q':3,'t':4}[token]]+=power
     else:
      assert power>=0;names.extend([function(token)]*power)
    names.sort(key=lambda k:len(eng.live[k]['S']))
    if names:
     node=eng.copy(fresh('term'),names[0])
     for x in names[1:]:
      nxt=eng.mul(fresh('term_product'),node,x);eng.drop(node);node=nxt
     nxt=eng.copy(fresh('term_scaled'),node,Q(coefficient),tuple(shift));eng.drop(node);node=nxt
    else:node=eng.leaf(fresh('constant'),{'content':str(Q(coefficient)),'shift':shift,'terms':[['1',0,0,0,0,0]]},source={'known_QQ_Laurent_monomial':True})
    dot=eng.mul(fresh('dot'),node,'raw53');eng.drop(node)
    if accumulator is None:accumulator=dot
    else:
     nxt=eng.add(fresh('sum'),accumulator,dot,consume=True);accumulator=nxt
   eng.drop('raw53')
   pf=write(f'phase-{position:02d}.json',{'source_only':True,'nodes':eng.nodes,'ops':eng.ops,'accumulator':accumulator});phasefiles.append({'file':pf.name,'sha256':sha(pf)})
   keep=set(eng.live);eng.nodes={k:v for k,v in eng.nodes.items()if k in keep};eng.ops=[];eng.meta=sum(count_numeric(v)+12 for v in eng.nodes.values());gc.collect()
   assert time.monotonic()-start<55
  final=eng.nodes[accumulator];report={'id':ID,'all_passed':True,'source_only':True,'fixed_row':6,'fixed_col53':True,'raw_components':11,'output_support_upper':final['terms_upper'],'output_profile':final,'SOURCE_owned_peak':eng.peak,'H_slice_hash_alias_peak':eng.hash_H_alias_peak,'literal_pair_majorant':eng.pair_count,'phase_files':phasefiles,'future_VALUE_profile_7_numeric_output_alone':7*final['terms_upper'],'no_VALUE_GO_or_all14_row_claim':True,'all167x14_87plus66_RHS_P290once_F9_F10_complements_kept':True,'input_sha256':inputs,'seconds':time.monotonic()-start}
  write('report.json',report);signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','phase_files')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  if not(out/'failure.json').exists():write('failure.json',{'all_passed':False,'source_only':True,'error':repr(exc),'completed_phase_files':phasefiles,'no_new_matrix_coefficients':True,'no_automatic_retry':True})
  raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--profile',action='store_true');main(p.parse_args())
