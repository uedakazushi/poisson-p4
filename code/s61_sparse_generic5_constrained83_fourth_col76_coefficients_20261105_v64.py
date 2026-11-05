#!/usr/bin/env python3
"""ONE20 source-bound clearedLaurent numerators, integerconvolution only."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
from itertools import product
import argparse,hashlib,json,time,resource,signal,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_COEFFICIENTS_20261105_v64';GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_FOURTH_COL76_IMPLEMENTATION_GATE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def metadata(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(metadata,v))
 if isinstance(v,dict):return sum(map(metadata,v.values()))
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

def main(a):
 assert a.execute,'use --execute only after immutable source/storage gate'
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);gp=ROOT/f'reports/{GATE}/report.json';g=read(gp);assert g['all_passed']and g['ONE20_cleared_column_job_safe'];ap=gp.parent/'implementation-custody-storage-bit-addendum.json'
 for name,h in g['evidence_sha256'].items():assert sha(gp.parent/name)==h
 b=read(ap);assert sha(Path(__file__))==b['worker_sha256']
 for p,h in b['actual_READ_input_sha256'].items():assert sha(ROOT/p)==h
 prepared=sum(p.stat().st_size for p in gp.parent.iterdir()if p.is_file());assert prepared<=131072
 data=read(ROOT/b['operand_bundle']);plans=read(ROOT/b['source_plans'])['all20_source_plans'];assert len(plans)==20;leaves={}
 for k,z in data['source_operands'].items():
  p=z['original_literal'];shift=[p['shift'][j]+z['transform_monomial_exponent5'][j]for j in range(5)];leaves[k]={'content':Q(p['content'])*Q(z['QQ_transform_scalar']),'terms':[(int(t[0]),tuple(shift[j]+t[j+1]for j in range(5)))for t in p['terms']]}
 assert metadata(data)+metadata(plans)+metadata(b)+sum(7*len(z['terms'])for z in leaves.values())<25000
 records=[];events=[];peak=25000;total_tuples=0;max_integer_bits=0;backadds=0
 try:
  for plan in plans:
   row=plan['row'];D=int(plan['integer_common_D']);M=int(plan['integer_numerator_abs_upper']);allowed={tuple(e)for e in plan['exact_source_support_sorted']};acc={};row_peak=25000;row_tuples=0;row_back=0
   for ip,p in enumerate(plan['products']):
    factors=[leaves[k]for k in p['factor_labels']];q=Q(1)if len(p['factor_labels'])==2 else(Q(405,8)if p['factor_labels'][0].startswith('U')else Q(-405,8))
    for z in factors:q*=z['content']
    assert q==Q(p['source_product_content'])and D%q.denominator==0;weight=q.numerator*(D//q.denominator);assert weight==int(p['exact_scalar_integer_weight_after_D']);block={};tuples=0
    for combo in product(*(z['terms']for z in factors)):
     c=weight;e=list(p['prefactor_monomial_exponent5'])
     for z,f in combo:
      c*=z;assert abs(c)<=M
      for j in range(5):e[j]+=f[j]
     e=tuple(e);assert e in allowed;block[e]=block.get(e,0)+c;assert abs(block[e])<=M;tuples+=1;max_integer_bits=max(max_integer_bits,abs(c).bit_length(),abs(block[e]).bit_length())
    block={e:z for e,z in block.items()if z};assert tuples==p['literal_tuple_count']and len(block)<=p['exact_product_support_cardinality'];before=dict(acc);acc=add(acc,block);assert set(acc)<=allowed and all(abs(z)<=M for z in acc.values());back=add(acc,block,-1);assert back==before;extra=7*(len(acc)+len(block)+len(before)+len(back))+128;row_peak=max(row_peak,25000+extra);assert row_peak<=plan['future_owned_extra_including25kmeta']<=28968;row_back+=1;backadds+=1;row_tuples+=tuples
    events.append({'row':row,'product_index':ip,'literal_tuple_count':tuples,'integer_D':str(D),'all_partial_integer_coefficients_abs_le_M':True,'source_product_backadd_exact':True});del block,before,back
   assert row_tuples==sum(p['literal_tuple_count']for p in plan['products']);total_tuples+=row_tuples;assert len(acc)<=plan['exact_union_support_cardinality']and set(acc)<=allowed;G=0
   for z in acc.values():G=gcd(G,abs(z));assert abs(z)<=M
   n=len(acc)
   if n:
    base=plan['expected_shift_lower5'];order=sorted(acc,key=lambda e:(sum(e),e),reverse=True);value={'content':str(Q(G,D)),'shift':base,'terms':[[str(acc[e]//G),*[e[j]-base[j]for j in range(5)]]for e in order]};expected=vhash(value);row_peak=max(row_peak,25000+20*n+128);assert row_peak<=plan['future_owned_extra_including25kmeta'];del order,acc,allowed;gc.collect();payload=encode(value);raw=b'Q5CELL1\n\x01'+payload;back=decode(payload);assert encode(back)==payload and vhash(back)==expected;del payload,back
   else:value={'content':'0','shift':[0]*5,'terms':[]};expected=vhash(value);raw=b'Q5CELL1\n\x00'
   file=out/f'row-{row}-col76-cleared.bin';assert not file.exists();tmp=file.with_suffix('.bin.partial');tmp.write_bytes(raw);tmp.rename(file);assert file.stat().st_size<=plan['uncompressed_binary_RAW_upper'];blob=file.read_bytes();assert blob==raw
   if n:
    p=decode(blob[9:]);assert vhash(p)==expected and encode(p)==blob[9:];q=Q(p['content'])
    for t in p['terms']:z=q*int(t[0])*D;assert z.denominator==1 and abs(z.numerator)<=M
    del p
   records.append({'row':row,'column':76,'file':file.name,'sha256':sha(file),'actual_QQ_numerator_value_sha256':expected,'terms':n,'value_denominator':'5q0+1','row_source_backadds':row_back,'literal_primitive_tuples':row_tuples,'all_partial_source_height_checks':True,'all_support_and_binary_roundtrip_checks':True});peak=max(peak,row_peak);assert peak<=28968;del value,raw,blob;gc.collect();assert time.monotonic()-start<55
  assert total_tuples==1068 and len(records)==20 and not list(out.glob('*.partial'))
  (out/'all20-source-product-backadds.json').write_bytes(enc({'all_passed':True,'products':events,'all_partial_integer_coefficients_bounded':True}))
  api={'all_passed':True,'actual20_col76_numerator_values':records,'actual_value_representation':'S4_i76 = numerator_file / (5q0+1) on the EXISTINGaddedQchart. This isNOTa new columnoperation.','common_denominator_polynomial':'5q0+1','new_Q_column_scale_NOT_applied':True,'P290_never_multiplied_again':True,'all_other19_columns_still_typed_no_updates':True,'full173_top107_plus66_RHS_87extras_retained':True,'selected20_square_NOTsamecokernel_as_full173':True,'fixed24_det_and_exhaustion_unresolved':True}
  (out/'actual-cleared-col76-and-commonQ-representation.json').write_bytes(enc(api));base=b['actual_READ_input_RAW_bytes']+prepared;rawtotal=base+sum(p.stat().st_size for p in out.iterdir());assert rawtotal<=b['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'actual_cleared_QQ_cells':20,'actual_output_terms':sum(z['terms']for z in records),'literal_primitive_tuples':total_tuples,'source_product_backadds':backadds,'all20_binary_roundtrips_valuehashes_and_height_checks':True,'integer_partial_max_bits':max_integer_bits,'whole_coefficient_owned_upper':28968,'checked_coefficient_owned_peak':peak,'source_full_parse_numeric_upper_reference':87029,'actual_fresh_RAW_before_small_report':rawtotal,'commonQ_representation_NOT_new_column_operation':True,'all19_other_columns_RHS_and_P290_scope_kept':True,'input_sha256':{**b['actual_READ_input_sha256'],str(gp.relative_to(ROOT)):sha(gp),str(ap.relative_to(ROOT)):sha(ap)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));assert base+sum(p.stat().st_size for p in out.iterdir())<=b['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as e:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(e),'all_completed_cells_preserved':True,'no_automatic_retry':True,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--execute',action='store_true',help='run ONLY after immutable implementation/source/storage gate');main(p.parse_args())
