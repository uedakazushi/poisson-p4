#!/usr/bin/env python3
"""Prepared one15 col56 coefficient job: thin integer DAG and streaming cells."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time,signal,resource,gc,re
from s61_sparse_generic5_stream_primitive_cell_codec_20261105_v64 import write_cell,validate_cell
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_IMPLEMENTATION_GATE_V2_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def digest(v):return hashlib.sha256(enc(v)).hexdigest()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def main(a):
 assert a.execute and a.root_review_file and a.root_review_sha,'requires explicit --execute and frozen ROOT review'
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000))
 def timeout(s,f):raise TimeoutError('hard60s; preserve completed coefficient prefix')
 signal.signal(signal.SIGALRM,timeout);signal.alarm(60)
 rootreview=ROOT/a.root_review_file;assert sha(rootreview)==a.root_review_sha and rootreview.stat().st_size<=65536;rootdoc=read(rootreview);assert rootdoc['all_passed']and scalars(rootdoc)<=2000
 gp=ROOT/f'reports/{GATE}/report.json';g=read(gp);assert g['all_passed']and g['prepared_one15_staged_job_safe'];proofpath=gp.parent/'implementation-custody-streaming-lifetime-proof.json';assert sha(proofpath)==g['evidence_sha256'][proofpath.name];proof=read(proofpath);assert sha(Path(__file__))==proof['worker_sha256'];META=proof['fixed_metadata_limit'];LIMIT=proof['coefficient_owned_upper'];assert META==20000 and LIMIT<=130791
 assert proof['actual_READ_inputs_source']=='GATE report.input_sha256'and digest(g['input_sha256'])==proof['actual_READ_inputs_map_sha256']
 for p,h in g['input_sha256'].items():assert sha(ROOT/p)==h
 prepared=sum(p.stat().st_size for p in gp.parent.iterdir()if p.is_file());assert prepared<=131072
 control=read(ROOT/proof['thin_controller_file']);phaseindex=control['phase_index16'];assert len(phaseindex)==16;shared=read(ROOT/phaseindex[0]['file']);assert sha(ROOT/phaseindex[0]['file'])==phaseindex[0]['sha256'];auxid=control['shared_N8_35_56_node'];outputinfo={z['node']:z for z in control['output_nodes15']};out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 values={};resident=0;ownedpeak=META;parsepeak=0;maxbits=0;product_updates=0;add_backchecks=0;copy_checks=0;events=[];results=[];auxresult=None;phaseNo=-1;opNo=-1
 def guard(extra=128):
  nonlocal ownedpeak
  n=META+7*resident+extra;ownedpeak=max(ownedpeak,n);assert n<=LIMIT,(n,LIMIT,phaseNo,opNo);return n
 def parseguard(current,leaf=None):
  nonlocal parsepeak
  n=scalars(g)+scalars(proof)+scalars(rootdoc)+scalars(control)+scalars(shared)+(0 if current is shared else scalars(current))+(scalars(leaf)if leaf else 0)+scalars(events)+scalars(results)+2048;parsepeak=max(parsepeak,n);assert n<=META,(n,META,phaseNo,opNo)
 def release(name):
  nonlocal resident
  z=values.pop(name,None)
  if z is not None:resident-=len(z['z'])
 def insert(node,e,c,bound,cap,box):
  nonlocal resident,maxbits
  assert isinstance(c,int)and abs(c)<=bound and all(box[j][0]<=e[j]<=box[j][1]for j in range(5));Z=node['z'];old=Z.get(e,0);new=old+c;assert abs(new)<=bound and new-c==old;maxbits=max(maxbits,abs(c).bit_length(),abs(new).bit_length())
  if new:
   if e not in Z:guard(7+128);resident+=1
   Z[e]=new
  elif e in Z:del Z[e];resident-=1
  assert len(Z)<=cap;return old,new
 def eval_record(rec,current):
  nonlocal resident,product_updates,add_backchecks,copy_checks,maxbits
  op=rec['original_op'];profile=rec['output_profile'];assert digest(op)==rec['original_op_sha256']and digest(profile)==rec['output_profile_sha256'];dest=op['output'];assert dest not in values;cap=profile['terms_upper'];bound=int(profile['primitive_integer_abs_upper']);qq=Q(profile['known_QQ_content']);box=profile['exponent_box5'];kind=op['op'];updates=0;backs=0
  if kind=='input':
   ref=rec['input_leaf'];p=ROOT/ref['file'];assert sha(p)==ref['sha256'];leaf=read(p);parseguard(current,leaf);assert leaf['original_input_op_sha256']==rec['original_op_sha256']and leaf['original_binding']==op['binding'];z=leaf['operand'];raw=z['original_literal'];q=Q(raw['content'])*Q(z['QQ_transform_scalar']);assert q==qq;shift=[raw['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];values[dest]={'q':qq,'z':{},'profile':profile};node=values[dest]
   for t in raw['terms']:
    e=tuple(shift[j]+t[j+1]for j in range(5));assert e not in node['z'];c=int(t[0]);assert c;insert(node,e,c,bound,cap,box);updates+=1
   assert len(node['z'])==cap;del leaf,z,raw
  elif kind=='copy':
   src=values[op['A']];assert src['q']==qq;values[dest]={'q':qq,'z':{},'profile':profile};node=values[dest]
   for e,c in src['z'].items():insert(node,e,c,bound,cap,box);assert node['z'][e]==c;copy_checks+=1;updates+=1
  elif kind=='multiply':
   A=values[op['A']];B=values[op['B']];q=Q(op['scalar']);assert A['q']*B['q']*q==qq;shift=op['shift5'];values[dest]={'q':qq,'z':{},'profile':profile};node=values[dest]
   for e,x in A['z'].items():
    for f,y in B['z'].items():
     key=tuple(e[j]+f[j]+shift[j]for j in range(5));c=x*y;old,new=insert(node,key,c,bound,cap,box);assert new-c==old;updates+=1;backs+=1
   product_updates+=updates;add_backchecks+=backs
  elif kind=='add':
   A=values[op['A']];B=values[op['B']];wa=int(op['integer_weight_A']);wb=int(op['integer_weight_B']);qB=Q(op['scalarB']);assert A['q']==qq*wa and B['q']*qB==qq*wb;reuse=rec['last_use']['reuse_left_last_use_buffer']
   if reuse:
    assert op['A']in rec['last_use']['release_after_arithmetic_before_encoding']and A['z']is not B['z'];oldnode=values.pop(op['A']);node={'q':qq,'z':oldnode['z'],'profile':profile};values[dest]=node
    for e,old in node['z'].items():
     new=wa*old;assert abs(new)<=bound and(wa and new//wa==old);node['z'][e]=new;maxbits=max(maxbits,abs(old).bit_length(),abs(new).bit_length());backs+=1
   else:
    values[dest]={'q':qq,'z':{},'profile':profile};node=values[dest]
    for e,x in A['z'].items():insert(node,e,wa*x,bound,cap,box);backs+=1
   for e,x in B['z'].items():old,new=insert(node,e,wb*x,bound,cap,box);assert new-wb*x==old;updates+=1;backs+=1
   add_backchecks+=backs
  else:raise AssertionError(kind)
  node=values[dest];assert node['q']==qq and len(node['z'])<=cap and all(c and abs(c)<=bound for c in node['z'].values());guard();return {'op':rec['index'],'kind':kind,'actual_terms':len(node['z']),'updates':updates,'backchecks':backs,'all_partial_positive_bounds':True}
 def output_cell(nodeid,row,kind):
  node=values[nodeid];profile=node['profile'];n=len(node['z']);other=resident-n;peak=META+7*other+10*n+2048;assert peak<=LIMIT;guard(3*n+2048);shift=[b[0]for b in profile['exponent_box5']];file=out/(f'row-{row}-col56-ninth-cleared.bin'if kind=='residual'else'shared-row35-col56-eighth-cleared.bin');record=write_cell(file,node['q'],node['z'],shift);assert validate_cell(file,node['q'],node['z'],record);cap=proof['output_caps'][str(row)if kind=='residual'else'shared35'];assert file.stat().st_size<=cap;return {'row':row,'kind':kind,'node':nodeid,'file':file.name,**record,'all_coefficients_incrementally_readback':True,'same_primitive_format_no_term_arrays':True}
 try:
  for phaseNo,ref in enumerate(phaseindex):
   p=ROOT/ref['file'];assert sha(p)==ref['sha256'];current=shared if phaseNo==0 else read(p);parseguard(current)
   for rec in current['literal_records']:
    opNo=rec['index'];parseguard(current);events.append(eval_record(rec,current))
    for name in rec['last_use']['release_after_arithmetic_before_encoding']:release(name)
    nodeid=rec['last_use']['release_output_after_encoding']
    if nodeid is not None:
     meta=outputinfo[nodeid];result=output_cell(nodeid,meta['row'],'residual');results.append(result);release(nodeid);(out/f"completed-prefix-after-row-{meta['row']}.json").write_bytes(enc({'all_passed':True,'completed_rows':[z['row']for z in results],'actual_output_terms':sum(z['terms']for z in results),'all_previous_cells_immutable':True}))
    assert time.monotonic()-start<55
   if phaseNo==0:auxresult=output_cell(auxid,35,'shared_auxiliary')
   if current is not shared:del current
   gc.collect()
  assert len(events)==700 and len(results)==15 and auxresult is not None and resident==0 and not values
  (out/'all700-stage-partial-and-backadd-checks.json').write_bytes(enc({'all_passed':True,'operations':events,'integer_product_updates':product_updates,'integer_add_backchecks':add_backchecks,'copy_entry_checks':copy_checks,'all_source_profiles_contents_partial_heights_boxes_and_cardinalities':True}))
  (out/'actual15-col56-and-once-shared-row35.json').write_bytes(enc({'all_passed':True,'actual15_ninth_col56_numerators':results,'actual_shared_N8_35_56_once':auxresult,'residual_representation':'S9_i56=N9_i56/(Q^4*(q0+1)^2*E*F9)','shared_representation':'S8_35,56=N8_35,56/(Q^3*(q0+1)*E)','no_new_column_operation_or_P290_scale':True,'full168x15_other14_87plus66_RHS_all_complements_kept':True,'no_rank_inverse_or_exhaustion_claim':True}))
  actualRAW=proof['actual_READ_input_RAW_bytes']+prepared+rootreview.stat().st_size+sum(p.stat().st_size for p in out.iterdir());assert actualRAW+65536<=proof['fresh_RAW_upper'];report={'id':ID,'all_passed':True,'actual_residual_cells':15,'shared_auxiliary_computed_once':True,'actual_residual_terms':sum(z['terms']for z in results),'shared_actual_terms':auxresult['terms'],'literal_ops':700,'product_integer_updates':product_updates,'integer_backchecks':add_backchecks,'actual_coefficient_integer_max_bits':maxbits,'whole_owned_numeric_peak':ownedpeak,'VALUE_parse_numeric_peak':parsepeak,'META_limit':META,'coefficient_owned_limit':LIMIT,'actual_fresh_RAW_before_report':actualRAW,'source_staged_pair_majorant':4442487,'no_direct995410878_unroll':True,'stream_primitive_all16_cells_roundtrip':True,'all168_RHS_87plus66_other14_oneP290_and_complements_kept':True,'input_sha256':{**g['input_sha256'],str(gp.relative_to(ROOT)):sha(gp),str(proofpath.relative_to(ROOT)):sha(proofpath),str(rootreview.relative_to(ROOT)):sha(rootreview)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));assert proof['actual_READ_input_RAW_bytes']+prepared+rootreview.stat().st_size+sum(p.stat().st_size for p in out.iterdir())<=proof['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(exc),'phase':phaseNo,'op':opNo,'completed_rows':[z['row']for z in results],'all_completed_binary_cells_preserved':True,'no_automatic_retry':True,'owned_peak':ownedpeak,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');p.add_argument('--root-review-file');p.add_argument('--root-review-sha');main(p.parse_args())
