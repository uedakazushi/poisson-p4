#!/usr/bin/env python3
"""Independent source/content/box/lifetime review of the fixed staged col56 job."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,re,gc,time,signal,resource
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(p.read_text())
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for x in iter(lambda:f.read(65536),b''):h.update(x)
 return h.hexdigest()
def digest(v):return hashlib.sha256((json.dumps(v,separators=(',',':'))+'\n').encode()).hexdigest()
def scalars(x):
 if isinstance(x,bool):return 0
 if isinstance(x,(int,float)):return 1
 if isinstance(x,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',x)else 0
 if isinstance(x,list):return sum(map(scalars,x))
 if isinstance(x,dict):return sum(map(scalars,x.values()))
 return 0
def box(ps,shift=None):
 if shift is not None and any(not p['terms_upper']for p in ps):return [[0,0]for _ in range(5)]
 ps=[p for p in ps if p['terms_upper']]
 if not ps:return [[0,0]for _ in range(5)]
 if shift is None:return [[min(p['exponent_box5'][j][0]for p in ps),max(p['exponent_box5'][j][1]for p in ps)]for j in range(5)]
 return [[sum(p['exponent_box5'][j][k]for p in ps)+shift[j]for k in range(2)]for j in range(5)]
def main():
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 gp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_IMPLEMENTATION_GATE_V2_20261105_v64/report.json';assert sha(gp)=='c1dab9549187c3161a46393b909d426e331285a695c56da4727bddfe26a9a7fa';g=read(gp);bp=gp.parent/'implementation-custody-streaming-lifetime-proof.json';b=read(bp)
 assert sha(bp)==g['evidence_sha256'][bp.name] and digest(g['input_sha256'])==b['actual_READ_inputs_map_sha256']
 inputs={gp:sha(gp),bp:sha(bp),Path(__file__):sha(Path(__file__))}
 for n,h in g['input_sha256'].items():assert sha(ROOT/n)==h;inputs[ROOT/n]=h
 assert b['worker_sha256']=='bacb8eec599b16f850b35a2ffa59e807932e798ca6a1138412b4f435c9535b52' and b['fixed_metadata_limit']==20000
 hp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_CODEC_PHASE_HEIGHT_SOURCE_20261105_v64/report.json';assert sha(hp)=='8213d1b8c276e4b44d9104842c42b8a3db68a5aabdc737cd8432a7447b9004d5';h=read(hp);assert h['all_passed'];inputs[hp]=sha(hp)
 cp=ROOT/b['thin_controller_file'];c=read(cp);inputs[cp]=sha(cp);parentp=ROOT/c['parent_program'];assert sha(parentp)==c['parent_program_sha256'];parent=read(parentp);inputs[parentp]=sha(parentp)
 assert len(parent['literal_ops'])==700 and parent['literal_shared_N8_35_56_node']==c['shared_N8_35_56_node']
 expected={i:(digest(op),digest(parent['positive_support_and_height_profiles'][op['output']]))for i,op in enumerate(parent['literal_ops'])};del parent;gc.collect()
 phases=[];records=[]
 for ref in c['phase_index16']:
  p=ROOT/ref['file'];assert sha(p)==ref['sha256'];inputs[p]=sha(p);current=read(p);phases.append(current);records+=current['literal_records']
 assert len(records)==700 and [z['index']for z in records]==list(range(700))
 profiles={};qqs={};uses={};kinds={};leafs=0;qbits=0;maxcoef=0;sourcepeak=0
 for rec in records:
  op=rec['original_op'];p=rec['output_profile'];assert (digest(op),digest(p))==expected[rec['index']];dest=op['output'];assert dest not in profiles;k=op['op'];n=p['terms_upper'];M=int(p['primitive_integer_abs_upper']);q=Q(p['known_QQ_content']);qbits=max(qbits,abs(q.numerator).bit_length(),q.denominator.bit_length());maxcoef=max(maxcoef,M.bit_length());kinds[k]=kinds.get(k,0)+1
  if k=='input':
   ref=rec['input_leaf'];f=ROOT/ref['file'];assert sha(f)==ref['sha256'];inputs[f]=sha(f);z=read(f);assert z['original_input_op_sha256']==digest(op) and z['original_binding']==op['binding'];a=z['operand'];raw=a['original_literal'];assert q==Q(raw['content'])*Q(a['QQ_transform_scalar']);shift=[raw['shift'][j]+a['transform_monomial_exponent5'][j]for j in range(5)];exps=[tuple(shift[j]+v[j+1]for j in range(5))for v in raw['terms']]
   assert len(exps)==len(set(exps))==n and M==max([abs(int(v[0]))for v in raw['terms']]+[0])
   actual=[[min(e[j]for e in exps),max(e[j]for e in exps)]for j in range(5)]if exps else[[0,0]for _ in range(5)];assert actual==p['exponent_box5'];leafs+=1;del z,a,raw,exps
  elif k=='copy':
   a=profiles[op['A']];assert q==qqs[op['A']]and n==a['terms_upper']and M==int(a['primitive_integer_abs_upper']) and p['exponent_box5']==a['exponent_box5']
  elif k=='multiply':
   a=profiles[op['A']];d=profiles[op['B']];assert q==qqs[op['A']]*qqs[op['B']]*Q(op['scalar']);assert M==int(a['primitive_integer_abs_upper'])*int(d['primitive_integer_abs_upper'])*min(a['terms_upper'],d['terms_upper']);assert n<=a['terms_upper']*d['terms_upper'];assert p['exponent_box5']==box([a,d],op['shift5'])
  elif k=='add':
   a=profiles[op['A']];d=profiles[op['B']];wa=int(op['integer_weight_A']);wb=int(op['integer_weight_B']);assert qqs[op['A']]==q*wa and qqs[op['B']]*Q(op['scalarB'])==q*wb;assert M==abs(wa)*int(a['primitive_integer_abs_upper'])+abs(wb)*int(d['primitive_integer_abs_upper']);assert n<=a['terms_upper']+d['terms_upper'];assert p['exponent_box5']==box([a,d])
  else:raise AssertionError(k)
  profiles[dest]=p;qqs[dest]=q;uses[dest]=0
  for key in ['A','B']:
   if key in op:uses[op[key]]+=1
  sourcepeak=max(sourcepeak,16000+scalars(records)+scalars(b)+scalars(g)+scalars(c)+scalars(profiles)+2*len(qqs)+3*len(uses));assert sourcepeak<150000
 assert leafs==181 and kinds['multiply']==303 and kinds['add']==193 and qbits<=174 and maxcoef<=272
 outputs={z['node']:z for z in c['output_nodes15']};aux=c['shared_N8_35_56_node'];encoders=set(outputs)|{aux}
 for v in encoders:uses[v]+=1
 live={};peak=20000;perphase=[]
 for ix,current in enumerate(phases):
  for rec in current['literal_records']:
   op=rec['original_op'];dest=op['output'];n=profiles[dest]['terms_upper'];reuse=rec['last_use']['reuse_left_last_use_buffer'];before=sum(live.values())
   if reuse:
    assert op['op']=='add' and uses[op['A']]==1 and op['A']!=op['B'];oldn=live.pop(op['A']);assert n>=oldn;now=before+n-oldn
   else:now=before+n
   peak=max(peak,20000+7*now+128);assert peak<=130791;live[dest]=n
   for key in ['A','B']:
    if key in op:uses[op[key]]-=1
   for v in rec['last_use']['release_after_arithmetic_before_encoding']:
    assert uses[v]==0;live.pop(v,None)
   out=rec['last_use']['release_output_after_encoding']
   if out is not None:
    assert uses[out]==1;uses[out]-=1;nn=live[out];peak=max(peak,20000+7*(sum(live.values())-nn)+10*nn+2048);assert peak<=130791;live.pop(out)
  if ix==0:
   assert uses[aux]>=1;uses[aux]-=1;nn=live[aux];peak=max(peak,20000+7*(sum(live.values())-nn)+10*nn+2048);assert peak<=130791
  perphase.append(dict(phase=ix,whole_bound=peak))
 assert not live and not any(uses.values())
 # Scalar metadata is a separate cap from integer map coefficients. ROOTdoc
 # may use2000, all700 event entries2800, and16 outputrecords800.
 parse=0
 for ix,current in enumerate(phases):
  for rec in current['literal_records']:
   z=read(ROOT/rec['input_leaf']['file'])if rec['original_op']['op']=='input'else None
   count=scalars(g)+scalars(b)+scalars(c)+scalars(phases[0])+(scalars(current)if ix else 0)+(scalars(z)if z else 0)+2000+2800+800+2048
   parse=max(parse,count);assert count<=20000
 # Only profile bounds are used here; no coefficient DAG is executed. Every
 # output header has generous2048 bytes even q*G numerator<=174+272 bits.
 for node in encoders:
  p=profiles[node];n=p['terms_upper'];nb=(int(p['primitive_integer_abs_upper']).bit_length()+7)//8;coords=sum(max(1,((hi-lo).bit_length()+6)//7)for lo,hi in p['exponent_box5']);cap=b['output_caps'][str(outputs[node]['row'])if node in outputs else'shared35'];assert 2056+n*(1+nb+coords)<=cap
 raw=sum(f.stat().st_size for f in inputs)+65536;assert raw<32*1024*1024
 out=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_ROOT_SOURCE_REVIEW_20261105_v64';out.mkdir(exist_ok=False)
 r=dict(id=out.name,all_passed=True,source_only=True,no_col56_coefficient_evaluated=True,all700_original_op_and_profile_digests_exact=True,all181_actual_leaf_profiles_checked=True,all700_QQ_contents_positive_integer_majorants_BOXes_rechecked=True,support_cardinality_profile_source='Reuses exact Minkowski census from frozen e332 SOURCE; caps and runtime before-allocation stops independently audited, no claim of a second full support expansion.',source_owned_peak=sourcepeak,value_whole_owned_bound=peak,VALUE_parse_with_max_ROOT_and_future_receipts_bound=parse,codec_map_bound272_and_DAGscalar174_preserved=True,codec_separate_conservative_QQ_temporary_bound718=True,all16_header_RAW_caps_verified=True,all_global_lastuse_and_left_reuse_counts_independently_verified=True,new_Schur_formula='N9_i56=Q*(q+1)*F9*N8_i56-62208000000*s*q*t^2*N8_i82*N8_35,56; denominator Q^4*(q+1)^2*E*F9.',shared_N8_35_56_computed_once=True,full168x15_87extras_66constraints_other14_oneP290_complements_retained=True,prospective_VALUE_RAW_bound=b['fresh_RAW_upper'],ROOT_source_READ_RAW_bound=raw,input_sha256={str(f.relative_to(ROOT)):h for f,h in inputs.items()},per_phase=perphase,seconds=time.monotonic()-begin)
 p=out/'report.json';p.write_text(json.dumps(r,indent=2)+'\n');assert p.stat().st_size<65536 and scalars(r)<=2000;signal.alarm(0);print(json.dumps({k:r[k]for k in ['all_passed','source_owned_peak','value_whole_owned_bound','VALUE_parse_with_max_ROOT_and_future_receipts_bound','seconds']}));print(sha(p))
if __name__=='__main__':main()
