#!/usr/bin/env python3
"""ONE21 shifted literal cached differences ONLY, exactQQ/binary/source bounds."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import hashlib,json,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_COEFFICIENTS_20261105_v64';GATE='S61_SPARSE_GENERIC5_CONSTRAINED83_SCHUR3_COL48_IMPLEMENTATION_GATE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def fingerprint(p):
 h=hashlib.sha256();h.update((str(len(p['terms']))+'\n').encode());q=Q(p['content'])
 for t in p['terms']:
  c=q*int(t[0]);e=[p['shift'][j]+t[j+1]for j in range(5)];h.update((str(c)+' '+' '.join(map(str,e))+'\n').encode())
 return h.hexdigest()
def metadata(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(metadata,v))
 if isinstance(v,dict):return sum(map(metadata,v.values()))
 return 0
def actual_poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def add(A,B):
 R=dict(A)
 for e,q in B.items():
  R[e]=R.get(e,Q(0))+q
  if not R[e]:del R[e]
 return R

def main():
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);gatep=ROOT/f'reports/{GATE}/report.json';g=read(gatep);assert g['all_passed']and g['ONE21_difference_job_safe'];ap=gatep.parent/'implementation-storage-bit-addendum.json';assert sha(ap)==g['evidence_sha256'][ap.name];a=read(ap)
 for p,h in a['actual_READ_input_sha256'].items():assert sha(ROOT/p)==h
 assert sha(Path(__file__))==a['worker_sha256'];prepared=gatep.stat().st_size+ap.stat().st_size;assert prepared<=65536
 frame=read(ROOT/a['cell_plans']);plans=frame['all21_exact_difference_plans'];assert len(plans)==21;assert metadata(plans)+metadata(a)<25000;records=[];peak=25000;zeros=[];monomials=[]
 try:
  for plan in plans:
   assert plan['column']==48;row=plan['row'];A={};B={};counts={};parsepeaks=[]
   for op in plan['operand_refs']:
    assert op['column']in(48,50);ref=op['source'];assert ref['kind']=='unchanged20_column_literal';p=ROOT/ref['file']
    with p.open('rb')as f:f.seek(ref['offset']);payload=f.read(ref['bytes'])
    assert hashlib.sha256(payload).hexdigest()==ref['payload_sha256'];value=decode(payload);assert encode(value)==payload and len(value['terms'])==op['terms']
    sign=1 if op['column']==48 else-1;transform=Q(1)if op['column']==48 else Q(256,3);assert Q(value['content'])*sign*transform==Q(op['QQ_content_with_sign']);V=actual_poly(value)
    if op['column']==50:V={tuple(e[j]+[1,0,0,0,-1][j]for j in range(5)):q*transform for e,q in V.items()}
    counts[op['column']]=len(V)
    # metadata25k contains all descriptor/exponent/index structures; this
    # extra term bound charges QQ coefficients and live value representations.
    parsepeaks.append(20*len(V)+7*(len(A)+len(B))+32)
    if op['column']==48:A=V
    else:B=V
    del value,payload,V
   na,nb=len(A),len(B);R=dict(A)
   for e,q in B.items():
    R[e]=R.get(e,Q(0))-q
    if not R[e]:del R[e]
   nr=len(R);assert nr<=plan['support_union_cardinality_upper'];back=add(R,B);assert back==A
   backpeak=7*(na+nb+nr+len(back))+32;del back
   D=int(plan['integer_common_denominator']);M=int(plan['integer_numerator_abs_upper']);source_support=set(A)|set(B);assert set(R)<=source_support
   for e,q in R.items():
    assert all(plan['actual_exponent_lower5'][j]<=e[j]<=plan['actual_exponent_upper5'][j]for j in range(5));v=q*D;assert v.denominator==1 and abs(v.numerator)<=M
   if not R:
    value={'content':'0','shift':[0]*5,'terms':[]};raw=b'Q5CELL1\n\x00';expected=fingerprint(value);localpeak=max(parsepeaks+[backpeak,128])
   else:
    base=plan['actual_exponent_lower5'];ints={e:(q*D).numerator for e,q in R.items()};G=0
    for v in ints.values():G=gcd(G,abs(v))
    assert G;ordered=sorted(R,key=lambda e:(sum(e),e),reverse=True);value={'content':str(Q(G,D)),'shift':base,'terms':[[str(ints[e]//G),*[e[j]-base[j]for j in range(5)]]for e in ordered]};expected=fingerprint(value)
    convertpeak=13*nr+32;del ints,ordered,R,A,B,source_support;gc.collect()
    payload=encode(value);raw=b'Q5CELL1\n\x01'+payload;decoded=decode(payload);assert encode(decoded)==payload and fingerprint(decoded)==expected;localpeak=max(parsepeaks+[backpeak,convertpeak,14*nr+32]);del decoded,payload
   peak=max(peak,25000+localpeak);assert 25000+localpeak<=plan['whole_owned_upper']<=30000
   file=out/f'row-{row}-col48.bin';assert not file.exists();tmp=file.with_suffix('.bin.partial');tmp.write_bytes(raw);tmp.rename(file);assert file.stat().st_size<=plan['binary_RAW_upper'];blob=file.read_bytes();assert blob==raw
   if blob[8]:again=decode(blob[9:]);assert fingerprint(again)==expected and encode(again)==blob[9:];del again
   else:assert len(blob)==9 and nr==0
   if nr==0:zeros.append(row)
   if nr==1:
    t=value['terms'][0];monomials.append({'row':row,'column':48,'QQ_coefficient':str(Q(value['content'])*int(t[0])),'parameter_exponent5':[value['shift'][j]+t[j+1]for j in range(5)],'unit_chart':'Nonzero displayed parameter monomial; all variables are K5 base parameters, not unknownnormal coordinates.'})
   records.append({'row':row,'column':48,'file':file.name,'sha256':sha(file),'actual_QQ_value_sha256':expected,'terms':nr,'backadd_result_plus_lambda_old50_equals_old48':True,'all_scalar_height_and_support_checks':True,'binary_roundtrip':True,'new_polynomial_products':0});del value,raw,blob;gc.collect()
   assert time.monotonic()-begin<55
  assert len(records)==21 and not list(out.glob('*.partial'))
  api={'all_passed':True,'actual21_residual_col48_entries':records,'zero_rows':zeros,'source_monomial_units':monomials,'all_zero':len(zeros)==21,'all21_backadds_and_binary_roundtrips':True,'layer':'ONLYselected21 post-unit3 residualcol48. Response96 remainsinput; oldpost-unit2col81cache is now ONLYinput to further typedcallbacks, not automaticallycurrentcoefficients.','full174_transport':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_THIRD_UNIT_SOURCE_20261105_v64/report.json','all87_other_top_RHS_plus66_constraints_kept':True,'chosen21_square_not_same_cokernel_as_full174':True,'other_residual20_columns_still_typed_no_updates':True,'if_all_zero':'If all21differences vanish, post-unit2col48=lambda*post-unit2col50 on all22rows, including pivotrow83. Since both originalpivotrows have0 in48/50, the actualfixed24columns48/50 aredependent. This only provesTHISfixed24detZERO, notfulltoprank.','if_nonzero':'Only recorded monomial/source structural conclusions; no automatic further field elimination.'}
  (out/'actual-residual-col48-and-structural-frontier.json').write_bytes(enc(api));readbase=a['actual_READ_input_RAW_bytes']+prepared;rawtotal=readbase+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert rawtotal<=a['fresh_RAW_upper']
  report={'id':ID,'all_passed':True,'actual_QQ_cells':21,'actual_output_terms':sum(z['terms']for z in records),'all21_backadds':True,'all21_binary_roundtrips_and_valuehashes':True,'all_zero':len(zeros)==21,'zero_cells':len(zeros),'monomial_unit_cells':len(monomials),'whole_owned_upper':30000,'checked_owned_peak':peak,'actual_fresh_RAW_before_small_report':rawtotal,'full174_87plus66_RHS_scope_retained':True,'no_P_double_or_other_matrix_updates':True,'input_sha256':{**a['actual_READ_input_sha256'],str(gatep.relative_to(ROOT)):sha(gatep),str(ap.relative_to(ROOT)):sha(ap)},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));assert readbase+sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<=a['fresh_RAW_upper'];signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as e:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(e),'all_completed_cells_preserved':True,'no_automatic_retry':True,'seconds':time.monotonic()-begin}));raise
if __name__=='__main__':main()
