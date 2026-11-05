#!/usr/bin/env python3
"""ONE16 row-streamed eighth-col82 QQ Laurent numerators over Q^3(q+1)E."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
from itertools import product
import argparse,hashlib,json,time,resource,signal,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1]
ID='S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64'
GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_IMPLEMENTATION_GATE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def vhash(p):
 h=hashlib.sha256();h.update((str(len(p['terms']))+'\n').encode());q=Q(p['content'])
 for t in p['terms']:h.update((str(q*int(t[0]))+' '+' '.join(str(p['shift'][j]+t[j+1])for j in range(5))+'\n').encode())
 return h.hexdigest()
def add(A,B,sgn=1):
 R=dict(A)
 for e,z in B.items():
  R[e]=R.get(e,0)+sgn*z
  if not R[e]:del R[e]
 return R
def compile_leaf(z):
 p=z['original_literal'];s=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)]
 return {'content':Q(p['content'])*Q(z['QQ_transform_scalar']),'terms':[(int(t[0]),tuple(s[j]+t[j+1]for j in range(5)))for t in p['terms']]}
def main(a):
 assert a.execute,'use --execute only after immutable gate'
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 gp=ROOT/f'reports/{GATE}/report.json';g=read(gp);assert g['all_passed']and g['ONE16_column82_job_safe'];ap=gp.parent/'implementation-custody-storage-bit-addendum.json'
 for name,h in g['evidence_sha256'].items():assert sha(gp.parent/name)==h
 b=read(ap);assert sha(Path(__file__))==b['worker_sha256']
 for p,h in b['actual_READ_input_sha256'].items():assert sha(ROOT/p)==h
 prepared=sum(p.stat().st_size for p in gp.parent.iterdir()if p.is_file());assert prepared<=131072
 index=read(ROOT/b['row_index_file']);assert len(index['row_records16'])==16
 data=read(ROOT/index['common_operand_file']);common={k:compile_leaf(z)for k,z in data['common_operands'].items()};del data;gc.collect()
 records=[];events=[];peak=40000;total_tuples=0;max_integer_bits=0;backadds=0;parse_peak=0
 try:
  for ref in index['row_records16']:
   row=ref['row'];rp=ROOT/ref['row_operand_file'];pp=ROOT/ref['plan_file'];assert sha(rp)==ref['row_operand_sha256']and sha(pp)==ref['plan_sha256']
   raw=read(rp);plan=read(pp);assert raw['row']==plan['row']==row
   leaves={**common,**{k:compile_leaf(z)for k,z in raw['row_operands'].items()}}
   fixed=numeric(raw)+numeric(plan)+numeric(b)+numeric(g)+numeric(index)+numeric(records)+numeric(events)+7*sum(len(z['terms'])for z in leaves.values())+2048
   parse_peak=max(parse_peak,fixed);assert fixed<40000;del raw;gc.collect()
   D=int(plan['integer_common_D']);M=int(plan['integer_numerator_abs_upper']);allowed={tuple(e)for e in plan['exact_source_support_sorted']};acc={};row_peak=40000;row_tuples=0;row_back=0
   for ip,p in enumerate(plan['products']):
    factors=[leaves[k]for k in p['factor_labels']];q=Q(p['full_parameter_monomial_prefactor']['QQ']);shift=p['full_parameter_monomial_prefactor']['exponent5']
    for z in factors:q*=z['content']
    assert q==Q(p['source_product_content'])and D%q.denominator==0;weight=q.numerator*(D//q.denominator);assert weight==int(p['exact_scalar_integer_weight_after_D']);block={};tuples=0
    for combo in product(*(z['terms']for z in factors)):
     c=weight;e=list(shift)
     for z,f in combo:
      c*=z;assert abs(c)<=M
      for j in range(5):e[j]+=f[j]
     e=tuple(e);assert e in allowed;block[e]=block.get(e,0)+c;assert abs(block[e])<=M;tuples+=1;max_integer_bits=max(max_integer_bits,abs(c).bit_length(),abs(block[e]).bit_length())
    block={e:z for e,z in block.items()if z};assert tuples==p['literal_tuple_count']and len(block)<=p['exact_product_support_cardinality'];before=dict(acc);acc=add(acc,block);assert set(acc)<=allowed and all(abs(z)<=M for z in acc.values());back=add(acc,block,-1);assert back==before
    extra=7*(len(acc)+len(block)+len(before)+len(back))+128;row_peak=max(row_peak,40000+extra);assert row_peak<=plan['future_owned_extra_including40kmeta']<=106064;row_back+=1;backadds+=1;row_tuples+=tuples
    events.append({'row':row,'product_index':ip,'literal_tuple_count':tuples,'all_partial_integer_coefficients_abs_le_M':True,'source_product_backadd_exact':True});del block,before,back
   assert row_tuples==ref['source_tuple_majorant'];total_tuples+=row_tuples;assert len(acc)<=plan['exact_union_support_cardinality']and set(acc)<=allowed;G=0
   for z in acc.values():G=gcd(G,abs(z));assert abs(z)<=M
   n=len(acc)
   if n:
    base=plan['expected_shift_lower5'];order=sorted(acc,key=lambda e:(sum(e),e),reverse=True);value={'content':str(Q(G,D)),'shift':base,'terms':[[str(acc[e]//G),*[e[j]-base[j]for j in range(5)]]for e in order]};expected=vhash(value);row_peak=max(row_peak,40000+20*n+128);assert row_peak<=plan['future_owned_extra_including40kmeta'];del order,acc,allowed;gc.collect();payload=encode(value);raw=b'Q5CELL1\n\x01'+payload;back=decode(payload);assert encode(back)==payload and vhash(back)==expected;del payload,back
   else:value={'content':'0','shift':[0]*5,'terms':[]};expected=vhash(value);raw=b'Q5CELL1\n\x00';del acc,allowed
   file=out/f'row-{row}-col82-cleared.bin';assert not file.exists();tmp=file.with_suffix('.bin.partial');tmp.write_bytes(raw);tmp.rename(file);assert file.stat().st_size<=plan['uncompressed_binary_RAW_upper'];blob=file.read_bytes();assert blob==raw
   if n:
    p=decode(blob[9:]);assert vhash(p)==expected and encode(p)==blob[9:];q=Q(p['content'])
    for t in p['terms']:z=q*int(t[0])*D;assert z.denominator==1 and abs(z.numerator)<=M
    del p
   records.append({'row':row,'column':82,'file':file.name,'sha256':sha(file),'actual_QQ_numerator_value_sha256':expected,'terms':n,'value_denominator':'(5q0+1)^3*(q0+1)*(32s0+15q0+35)','row_source_backadds':row_back,'literal_primitive_tuples':row_tuples,'all_partial_source_height_checks':True,'all_support_and_binary_roundtrip_checks':True});peak=max(peak,row_peak);assert peak<=106064;del value,raw,blob,leaves,plan;gc.collect();assert time.monotonic()-start<55
  assert total_tuples==2728788 and len(records)==16 and not list(out.glob('*.partial'))
  (out/'all16-source-product-backadds.json').write_bytes(enc({'all_passed':True,'products':events,'all_partial_integer_coefficients_bounded':True}))
  api={'all_passed':True,'actual16_col82_numerator_values':records,'actual_value_representation':'S8_i82 = numerator_file / ((5q0+1)^3*(q0+1)*(32s0+15q0+35)) on the existing eighth Q*(q0+1)*E chart; representation ONLY.','common_denominator_polynomial':'(5q0+1)^3*(q0+1)*(32s0+15q0+35)','no_Q_cancellation_assumed':True,'new_Q_column_scale_NOT_applied':True,'P290_never_multiplied_again':True,'all_other15_columns_still_typed_no_updates':True,'full169_top103_plus66_RHS_87extras_retained':True,'selected16_square_NOTsamecokernel_as_full169':True,'fixed24_det_and_exhaustion_unresolved':True}
  (out/'actual-cleared-col82-and-commonQ3RE-representation.json').write_bytes(enc(api));base=b['actual_READ_input_RAW_bytes']+prepared;rawtotal=base+sum(p.stat().st_size for p in out.iterdir());assert rawtotal<=b['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'actual_cleared_QQ_cells':16,'actual_output_terms':sum(z['terms']for z in records),'literal_primitive_tuples':total_tuples,'source_product_backadds':backadds,'all16_binary_roundtrips_valuehashes_and_height_checks':True,'integer_partial_max_bits':max_integer_bits,'whole_coefficient_owned_upper':106064,'checked_coefficient_owned_peak':peak,'checked_streamed_parse_numeric_peak':parse_peak,'source_full_parse_numeric_upper_reference':81715,'actual_fresh_RAW_before_small_report':rawtotal,'commonQ3RE_representation_NOT_new_column_operation':True,'all15_other_columns_RHS_and_P290_scope_kept':True,'input_sha256':{**b['actual_READ_input_sha256'],str(gp.relative_to(ROOT)):sha(gp),str(ap.relative_to(ROOT)):sha(ap)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));assert base+sum(p.stat().st_size for p in out.iterdir())<=b['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as e:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(e),'all_completed_cells_preserved':True,'no_automatic_retry':True,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true');main(p.parse_args())
