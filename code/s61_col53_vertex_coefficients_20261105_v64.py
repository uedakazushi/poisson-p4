"""SOURCE then bounded exact coefficients at ten exposed positive-support vertices."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,hashlib,json,signal,resource,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'code'))
from s61_sparse_generic5_source_view_profile_20261105_v64 import count_numeric
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
OLD='reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64'
ORDERS=[([ (i+j)%5 for j in range(5)],s)for i in range(5)for s in [1,-1]]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main(a):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));o=ROOT/'reports'/('S61_COL53_VERTEX_SOURCE_GATE_20261105_v64'if a.prepare else'S61_COL53_VERTEX_VALUES_20261105_v64');o.mkdir(exist_ok=False);pins={}
 def read(p):
  p=Path(p);pins[str(p.relative_to(ROOT))]=sha(p);return json.loads(p.read_text())
 def put(n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 def polynomial(src,profile):
  if 'binding'in src:
   b=src['binding'];doc=read(ROOT/OLD/b['scope']);v=(doc['common_operands']if 'common_operands'in doc else doc['row_operands'])[b['label']];return v['original_literal'],Q(v['QQ_transform_scalar']),v['transform_monomial_exponent5']
  if 'file'in src:
   p=ROOT/src['file'];pins[str(p.relative_to(ROOT))]=sha(p)
   if 'offset'in src:
    with p.open('rb')as f:f.seek(src['offset']);b=f.read(src['bytes'])
    assert hashlib.sha256(b).hexdigest()==src['payload_sha256']
   else:
    assert sha(p)==src['sha256'];b=p.read_bytes()
   return decode(b),Q(1),[0]*5
  assert src.get('QQ_monomial') and profile['terms_upper']==1
  return {'content':profile['known_QQ_content'],'shift':[x[0]for x in profile['exponent_box5']],'terms':[['1',0,0,0,0,0]]},Q(1),[0]*5
 def vertices(p,shift):
  ex=[tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5))for t in p['terms']];assert len(ex)==len(set(ex))
  return [min(ex,key=lambda m:tuple(sign*m[j]for j in order))if ex else None for order,sign in ORDERS]
 if a.prepare:
  endpoint=read(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/report.json');assert endpoint['all_passed'];graph=[];binding={};records=[]
  d0=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_ROW6_COEFFICIENT_PROFILE_V2_20261105_v64'
  records +=[d0/f'phase-{i:02d}.json'for i in range(8)]
  for st,components in [('S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_REMAINING3_VIEW_PROFILE_20261105_v64',[8]),('S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64',[9,10])]:
   for component in components:records+=sorted((ROOT/'reports'/st).glob(f'term-{component}-*.json'))
  for rp in records:
   rec=read(rp);ops=rec.get('ops',rec.get('literal_ops'));assert ops is not None
   for op in ops:
    name=op['output'];kind=op['op'];prof=rec['nodes'][name];node={'name':name,'op':kind,'source_profile':{k:prof[k]for k in ['known_QQ_content','primitive_integer_abs_upper']},'args':[],'scalar':op.get('scalar','1'),'shift':op.get('shift5',[0]*5)}
    if kind=='input':
     node['source']=op['source'];p,sc,shift=polynomial(op['source'],prof);leads=vertices(p,shift)
     if op['source'].get('QQ_monomial'):node['constant_profile']=prof
    elif kind in ['copy','immutable_SOURCE_view']:
     k=binding[op['A']];node['args']=[k];leads=[None if m is None else tuple(x+y for x,y in zip(m,node['shift']))for m in graph[k]['leads']]
    elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:
     i,j=binding[op['A']],binding[op['B']];node['args']=[i,j];leads=[None if x is None or y is None else tuple(x[t]+y[t]+node['shift'][t]for t in range(5))for x,y in zip(graph[i]['leads'],graph[j]['leads'])]
    elif kind in ['add','SOURCE_fused_union_product']:
     if kind=='add':args=[binding[op['A']],binding[op['B']]];node['scalarB']=op.get('scalarB','1')
     else:args=([binding[op['C']]]if op['C']is not None else[])+[binding[op['literal_dot']]];node['scalarB']='1'
     node['args']=args;leads=[]
     for t,(order,sign)in enumerate(ORDERS):
      v=[graph[k]['leads'][t]for k in args if graph[k]['leads'][t]is not None];leads.append(min(v,key=lambda m:tuple(sign*m[j]for j in order))if v else None)
    else:raise ValueError(kind)
    node['leads']=leads;binding[name]=len(graph);graph.append(node)
    assert 18000+count_numeric(graph)+10000<150000
  root=binding[endpoint['phases'][-1]['snapshot']['file'] if False else read(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json')['accumulator']]
  # Verify endpoint vertices against ALL exact positive-support snapshot keys.
  phase=read(ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64/phase-10.json');snap=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_COL53_LAST2_VIEW_PROFILE_20261105_v64'/phase['all_live_snapshot']['file'];pins[str(snap.relative_to(ROOT))]=sha(snap)
  from s61_sparse_generic5_source_view_profile_20261105_v64 import unpack
  best=[None]*10
  with snap.open('rb')as f:
   assert f.read(8)==b'Q5SUPP1\n';head=json.loads(f.read(int.from_bytes(f.read(4),'big')))
   for rec in head['maps']:
    for _ in range(rec['terms']):
     raw=f.read(5);assert len(raw)==5
     if rec['id']==phase['accumulator']:
      m=unpack(int.from_bytes(raw,'big'))
      for i,(order,s)in enumerate(ORDERS):
       if best[i]is None or tuple(s*m[j]for j in order)<tuple(s*best[i][j]for j in order):best[i]=m
   assert not f.read(1)
  assert [tuple(x)for x in graph[root]['leads']]==best
  required=set();stack=[root]
  while stack:
   k=stack.pop()
   if k not in required:required.add(k);stack.extend(graph[k]['args'])
  remap={k:i for i,k in enumerate(sorted(required))};graph=[{**graph[k],'args':[remap[j]for j in graph[k]['args']]}for k in sorted(required)];root=remap[root]
  bits=max((int(n['source_profile']['primitive_integer_abs_upper']).bit_length()+abs(Q(n['source_profile']['known_QQ_content']).numerator).bit_length())for n in graph);assert bits<=640
  own=count_numeric(graph)+len(graph)*20+40000;assert own<150000
  put('graph.json',dict(nodes=graph,root=root,orders=ORDERS,root_vertices=best))
  raw=sum((ROOT/n).stat().st_size for n in pins);assert raw+2000000<33554432
  put('report.json',dict(all_passed=True,SOURCE_ONLY_no_new_signed_entry_coefficient=True,ten_exposed_vertices_verified_against_full_snapshot=True,SSA_nodes=len(graph),raw53_each_redefinition_distinguished_by_SSA=True,owned_upper=own,hard_owned150000=150000,coefficient_prior_bits=bits,committed640_temp1921_within_VALUE2048=True,RAW_upper=raw+2000000,proof='Every finite Laurent support has a minimum for each signed cyclic lex order. This order is translation invariant, so product extremal exponents add. At a sum extremum only children with that same extremum contribute. Recursively computing that coefficient remains correct when cancellation makes it zero; it never proves whole polynomial zero. Full original literal sources and common monomial/scalar transforms retained.',input_sha256=pins,code_sha256=sha(Path(__file__)),evidence_sha256={'graph.json':sha(o/'graph.json')}));print(sha(o/'report.json'));return
 assert a.gate_sha;gp=ROOT/'reports/S61_COL53_VERTEX_SOURCE_GATE_20261105_v64/report.json';assert sha(gp)==a.gate_sha;g=read(gp);assert g['all_passed'] and sha(Path(__file__))==g['code_sha256']
 for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h
 ep=gp.parent/'graph.json';assert sha(ep)==g['evidence_sha256']['graph.json'];D=read(ep);graph=D['nodes'];values=[];put('phase-before-values.json',dict(ONE_bounded10coefficients_GO=True,gate_sha=a.gate_sha,owned=g['owned_upper'],prior=g['coefficient_prior_bits']))
 def checked(q):assert max(abs(q.numerator).bit_length(),q.denominator.bit_length())<=640;return q
 try:
  for node in graph:
   kind=node['op'];leads=node['leads'];args=node['args'];row=[]
   if kind=='input':
    p,sc,shift=polynomial(node['source'],node.get('constant_profile',{}));coeff={tuple(p['shift'][j]+int(t[j+1])+shift[j]for j in range(5)):int(t[0])for t in p['terms']}
    row=[checked(Q(p['content'])*sc*coeff.get(tuple(m),0))if m is not None else Q(0)for m in leads]
   elif kind in ['copy','immutable_SOURCE_view']:row=[checked(x*Q(node['scalar']))for x in values[args[0]]]
   elif kind in ['multiply','unmaterialized_multiply_for_SOURCE_union']:row=[checked(x*y*Q(node['scalar']))for x,y in zip(values[args[0]],values[args[1]])]
   else:
    for i,m in enumerate(leads):
     x=Q(0)
     for j,k in enumerate(args):
      if graph[k]['leads'][i]==m:x=checked(x+values[k][i]*(Q(node['scalarB'])if j else Q(1)))
     row.append(x)
   values.append(row);assert time.monotonic()-start<45
  result=[dict(monomial=m,coefficient=str(v),nonzero=bool(v))for m,v in zip(D['root_vertices'],values[D['root']])];put('ten-actual-vertex-coefficients.json',result)
  put('report.json',dict(all_passed=True,all10_exact_entry_coefficients_evaluated=True,nonzero_count=sum(v['nonzero']for v in result),entry_nonzero_proved=any(v['nonzero']for v in result),zero_coefficients_do_not_imply_zero_entry=True,scope='row6/col53 numerator only; full167x14/87+66/oneP290/F9F10/allcomplements retained, no full inverse or exhaustion claim',input_sha256={**pins,**g['input_sha256']},code_sha256=sha(Path(__file__)),evidence_sha256={p.name:sha(p)for p in o.iterdir()},seconds=time.monotonic()-start));print(sha(o/'report.json'))
 except Exception as e:put('failure-prefix.json',dict(error=repr(e),completed_values=[[str(x)for x in row]for row in values],next_node=len(values),no_auto_retry=True));raise
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--gate-sha');main(p.parse_args())
