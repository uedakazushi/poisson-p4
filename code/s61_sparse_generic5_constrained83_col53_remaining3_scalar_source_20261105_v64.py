#!/usr/bin/env python3
"""Read old LIVE07 losslessly and bound remaining3 scalar/H-BOX calculations."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,hashlib,json,resource,signal,time
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_SCALAR_SOURCE_20261105_v64'
PREV='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64'
OLD='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
PLAN='S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_COEFFICIENT_SIDE_SOURCE_20261105_v64'
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def unpack(k):return tuple((k//256**j)%256-64 for j in range(5))
def main(a):
 assert a.freeze
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('SOURCE60s')));signal.alarm(60)
 out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={}
 def read(p):
  inputs[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 try:
  previous=ROOT/'reports'/PREV
  fail=read(previous/'failure.json');assert fail['forecast']==155866
  phase=read(previous/'phase-07.json');snap=previous/phase['complete_live_snapshot']['file']
  assert sha(snap)==phase['complete_live_snapshot']['sha256']=='b4cecb1b7ff1fd78cd02268e20e53a794cb05b4e643bca25aa2d48dda2e97f61'
  inputs[str(snap.relative_to(ROOT))]=sha(snap);stats={};decoded=0;max_reader_set=0;snapshot_checks=[]
  with snap.open('rb')as f:
   assert f.read(8)==b'Q5SUPP1\n';n=int.from_bytes(f.read(4),'big');assert n<=64*4096+256
   header=json.loads(f.read(n));assert len(header['maps'])==22
   for rec in header['maps']:
    k=rec['id'];q=Q(rec['q']);M=int(rec['M']);p=phase['nodes'][k]
    assert q==Q(p['known_QQ_content'])and M==int(p['primitive_integer_abs_upper'])and rec['terms']==p['terms_upper']
    S=set();low=[191]*5;high=[-64]*5;hs={}
    for _ in range(rec['terms']):
     b=f.read(5);assert len(b)==5;e=int.from_bytes(b,'big');assert e<2**40 and e not in S;S.add(e);v=unpack(e)
     for j in range(5):low[j]=min(low[j],v[j]);high[j]=max(high[j],v[j])
     hs[v[0]]=hs.get(v[0],0)+1
    assert [[low[j],high[j]]for j in range(5)]==p['exponent_box5']
    assert sorted(hs.items())==[tuple(v)for v in p['H_slice_sizes']]
    if k==phase['accumulator']:
     hh=hashlib.sha256();hh.update(b'[');first=True
     for h in sorted(hs):
      bucket=sorted(unpack(e)for e in S if unpack(e)[0]==h)
      for v in bucket:
       if not first:hh.update(b',')
       first=False;hh.update(enc(list(v))[:-1])
     hh.update(b']\n');assert hh.hexdigest()==phase['canonical_accumulator_hash']
    decoded+=len(S);max_reader_set=max(max_reader_set,len(S))
    stats[k]={'q':q,'M':M,'n':len(S),'lo':low[0],'hi':high[0]}
    snapshot_checks.append({'id':k,'terms':len(S),'q_M_BOX_Hcounts_exact':True});del S
   assert not f.read(1)
  assert decoded==44555 and max_reader_set==35473 and phase['accumulator']=='new_188_sum'
  plan=read(ROOT/'reports'/PLAN/'exact-11-row-coefficient-operator-and-callback-owners.json')
  old=read(ROOT/'reports'/OLD/'staged-literal-col56-program-and-positive-profiles.json')
  registry={o['output']:o['binding']for o in old['literal_ops']if o['op']=='input'};del old
  funcs={token:stats[name]for token,name in [('UF10','UF10'),('C10_Ni56','C10_Ni56'),('C9_C10_alpha_i','C9_C10_alpha_i')]}
  for token in {x for row in plan['raw_coefficient_side11']for fac,c in row for x,n in fac if x not in['H','s','q','t']}:
   if 'leaf_'+token in stats:funcs[token]=stats['leaf_'+token]
  assert len(funcs)==21
  owners={}
  def leaf(token):
   if token in funcs:return funcs[token]
   name=token
   if token.endswith('i')and token[:-1]in['U','V','C','F','N','M']:name=token[:-1]+'6'
   elif token.startswith('Ai_'):name='A6_'+token[3:]
   b=registry[name];p=ROOT/'reports'/OLD/b['scope'];z=read(p);z=z['common_operands'][b['label']]if 'common_operands'in z else z['row_operands'][b['label']]
   poly=z['original_literal'];shift=poly['shift'][0]+z['transform_monomial_exponent5'][0];h=[shift+t[1]for t in poly['terms']]
   v={'q':Q(poly['content'])*Q(z['QQ_transform_scalar']),'M':max([abs(int(t[0]))for t in poly['terms']]+[0]),'n':len(poly['terms']),'lo':min(h)if h else 0,'hi':max(h)if h else 0};funcs[token]=v;owners[token]={'old_source_node':name,'binding':b,'existing_leaf_only':True};return v
  common=read(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_TENTH_COL53_STAGED_SOURCE_20261105_v64/common-source-operands.json')['common_operands']
  def rawstat(label):
   z=common[label];p=z['original_literal'];h=[p['shift'][0]+t[1]+z['transform_monomial_exponent5'][0]for t in p['terms']]
   return{'q':Q(p['content'])*Q(z['QQ_transform_scalar']),'M':max([abs(int(t[0]))for t in p['terms']]+[0]),'n':len(p['terms']),'lo':min(h)if h else 0,'hi':max(h)if h else 0}
  raw=[rawstat('A11_53'),rawstat('S'),rawstat('A52_53')]
  records=[];state=stats[phase['accumulator']].copy();maxbits=0;maxq=0;maxspan=0;maxraw=0
  def track(v):
   nonlocal maxbits,maxq,maxspan
   maxbits=max(maxbits,v['M'].bit_length());maxq=max(maxq,abs(v['q'].numerator).bit_length(),v['q'].denominator.bit_length());maxspan=max(maxspan,v['hi']-v['lo']+1)
   assert -64<=v['lo']<=v['hi']<=191;return v
  def product(A,B,materialized=True):
   n=A['n']*B['n']
   return track({'q':A['q']*B['q'],'M':min(A['n'],B['n'])*A['M']*B['M'],'n':min(125000,n)if materialized else n,'lo':A['lo']+B['lo'],'hi':A['hi']+B['hi']})
  for component,(terms,r)in enumerate(zip(plan['raw_coefficient_side11'][8:],raw),8):
   for ti,(factors,coefficient)in enumerate(terms):
    shift=0;v=[]
    for token,power in factors:
     if token in['H','s','q','t']:
      if token=='H':shift+=power
     else:
      assert power>=0;v.extend([leaf(token)]*power)
    v.sort(key=lambda x:x['n'])
    term=v[0].copy()if v else{'q':Q(1),'M':1,'n':1,'lo':0,'hi':0}
    for b in v[1:]:term=product(term,b)
    term['q']*=Q(coefficient);term['lo']+=shift;term['hi']+=shift;track(term)
    dot=product(term,r,False);d=state['q'].denominator*dot['q'].denominator//gcd(state['q'].denominator,dot['q'].denominator)
    aa=int(state['q']*d);bb=int(dot['q']*d);g=gcd(abs(aa),abs(bb));wa=aa//g if g else 0;wb=bb//g if g else 0
    state=track({'q':Q(g,d)if g else Q(0),'M':abs(wa)*state['M']+abs(wb)*dot['M'],'n':min(125000,state['n']+dot['n']),'lo':min(state['lo'],dot['lo']),'hi':max(state['hi'],dot['hi'])})
    records.append({'component':component,'term':ti,'q_exact':str(state['q']),'positive_M_upper':str(state['M']),'H_BOX_upper':[state['lo'],state['hi']],'materialized_node_card_bound_conditional_on_SAME_SOURCE_guard':125000,'dot_not_materialized_pair_card_upper':dot['n']})
  assert len(records)==93 and len(funcs)==38
  assert maxbits<=1024 and maxq<=256 and maxspan<=64
  ev={'all_passed':True,'source_only_no_profile_or_VALUE':True,'snapshot_all22_maps44555_keys_checked':snapshot_checks,'old_accumulator_hash_exact':True,'resume_phase':8,'resume_counter':188,'resume_function_registry_size':21,'remaining_terms':93,'coefficient_side_raw_components_remaining':[8,9,10],'new_function_owners':owners,'scalar_and_H_majorants93':records,'max_positive_M_bits_upper':maxbits,'max_QQ_field_bits_upper':maxq,'max_H_BOX_width_upper':maxspan,'conditional_materialized_set_cap':125000,'qM_node_storage_guard_proposed':[1024,256],'SOURCE_temp16384_unchanged':True,'views_and_fused_union_source_ONLY_not_VALUE_GO':True}
  ep=out/'old-live07-and-remaining93-scalar-H-majorants.json';ep.write_bytes(enc(ev))
  for p,h in read(previous/'partial-profile.json')['input_sha256'].items():assert sha(ROOT/p)==h;inputs[p]=h
  for p in[Path(__file__),ROOT/f'research/notes/{ID}.md']:inputs[str(p.relative_to(ROOT))]=sha(p)
  raw=sum((ROOT/p).stat().st_size for p in inputs)+ep.stat().st_size+65536;assert raw<32*1024*1024
  source_numeric_upper=25000+max_reader_set+6*2970+8192+12000;assert source_numeric_upper<150000
  report={'id':ID,'all_passed':True,'source_only':True,'snapshot22maps44555_keys_verified':True,'accumulator_hash_verified':True,'remaining_terms':93,'max_positive_M_bits_upper':maxbits,'max_QQ_field_bits_upper':maxq,'max_H_BOX_width_upper':maxspan,'SOURCE_reader_numeric_upper':source_numeric_upper,'no_profile_VALUE_or_GO':True,'READ_output_RAW_upper':raw,'input_sha256':inputs,'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'source_only':True,'error':repr(exc),'no_profile_or_VALUE_run':True,'no_automatic_retry':True}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true');main(p.parse_args())
