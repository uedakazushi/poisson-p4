#!/usr/bin/env python3
"""Distinct staged QQ Laurent VALUE validation, no integer-DAG producer calls."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,resource,signal,re,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_QQ_VALUE_REVIEW_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def scalars(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,str):return 3 if re.fullmatch(r'[+-]?[0-9]+(?:/[0-9]+)?',v)else 0
 if isinstance(v,list):return sum(map(scalars,v))
 if isinstance(v,dict):return sum(map(scalars,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def literal(z):return{tuple(e[j]+z['transform_monomial_exponent5'][j]for j in range(5)):q*Q(z['QQ_transform_scalar'])for e,q in poly(z['original_literal']).items()}
def stats(A):
 d=1
 for q in A.values():d=d*q.denominator//gcd(d,q.denominator)
 return d,max([abs((q*d).numerator)for q in A.values()]+[0])
def vhash(p):
 h=hashlib.sha256();h.update((str(len(p['terms']))+'\n').encode());q=Q(p['content'])
 for t in p['terms']:h.update((str(q*int(t[0]))+' '+' '.join(str(p['shift'][j]+t[j+1])for j in range(5))+'\n').encode())
 return h.hexdigest()
def main(a):
 assert a.check
 start=time.monotonic();resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));signal.signal(signal.SIGALRM,lambda s,f:(_ for _ in()).throw(TimeoutError('hard60s')));signal.alarm(60);out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
 pp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_COEFFICIENTS_20261105_v64/report.json';assert sha(pp)=='20524cd69e98bc7f5334a07544035fc04ea4fc1486ac64ce885670c0ce70b173';pr=read(pp);assert pr['all_passed'];rootpeer=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_ROOT_SOURCE_REVIEW_20261105_v64/report.json';assert sha(rootpeer)=='14597dfdab8c51626900f47a2b7c1dedb9f245917858b24871e5ea683a051868';rootdoc=read(rootpeer);assert rootdoc['all_passed'];inputs={ROOT/p:h for p,h in pr['input_sha256'].items()};inputs.update({pp:sha(pp),rootpeer:sha(rootpeer),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')})
 for name,h in pr['evidence_sha256'].items():inputs[pp.parent/name]=h
 for p,h in inputs.items():assert sha(p)==h
 api=read(pp.parent/'actual15-col56-and-once-shared-row35.json');recordsExpected=api['actual15_ninth_col56_numerators'];auxExpected=api['actual_shared_N8_35_56_once'];assert len(recordsExpected)==15;src=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_NINTH_COL56_STAGED_SOURCE_V4_20261105_v64';oldAPIpath=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/actual-cleared-col82-and-commonQ3RE-representation.json';assert oldAPIpath in inputs;oldAPI={z['row']:z for z in read(oldAPIpath)['actual16_col82_numerator_values']}
 P={};META=30000;peak=META;parsepeak=0;bits=0;events=[];results=[];products=0;backchecks=0
 def owned(extra=0):
  nonlocal peak
  n=META+8*sum(len(A)for A in P.values())+extra;peak=max(peak,n);assert n<=150000,(n,150000);return n
 def put(k,A):P[k]=A;owned()
 def drop(*names):
  for k in names:P.pop(k,None)
 def addinto(dest,block):
  nonlocal backchecks,bits
  A=P[dest];B=P[block];dA,mA=stats(A);dB,mB=stats(B);D=dA*dB//gcd(dA,dB);M=mA*(D//dA)+mB*(D//dB);bits=max(bits,D.bit_length(),M.bit_length());owned(8*len(A)+128);R=dict(A)
  for e,q in B.items():
   if e not in R:owned(8*(len(R)+1)+128)
   R[e]=R.get(e,Q(0))+q;z=R[e];assert D%z.denominator==0 and abs((z*D).numerator)<=M
   if not z:del R[e]
  owned(16*len(R)+128);back=dict(R)
  for e,q in B.items():
   if e not in back:owned(8*(len(R)+len(back)+1)+128)
   back[e]=back.get(e,Q(0))-q
   if not back[e]:del back[e]
  assert back==A;P[dest]=R;backchecks+=1;events.append({'kind':'QQ_stage_backadd','destination':dest,'terms':len(R),'exact':True});del back;owned()
 def pmul(dest,Aname,Bname,q=Q(1),shift=(0,0,0,0,0)):
  nonlocal products,bits
  assert dest not in(Aname,Bname);drop(dest);A=P[Aname];B=P[Bname];S=set()
  for e in A:
   for f in B:
    g=tuple(e[j]+f[j]+shift[j]for j in range(5))
    if g not in S:owned(13*(len(S)+1)+128);S.add(g)
  dA,mA=stats(A);dB,mB=stats(B);D=dA*dB*q.denominator;M=min(len(A),len(B))*mA*mB*abs(q.numerator);bits=max(bits,D.bit_length(),M.bit_length());R={}
  for e,x in A.items():
   for f,y in B.items():
    g=tuple(e[j]+f[j]+shift[j]for j in range(5));z=q*x*y;R[g]=R.get(g,Q(0))+z;assert D%R[g].denominator==0 and abs((R[g]*D).numerator)<=M
  R={e:z for e,z in R.items()if z};assert set(R)<=S;owned(5*len(S)+8*len(R)+128);del S;put(dest,R);products+=1
 def madd(dest,A,B,q,shift):pmul('_rv_block',A,B,q,shift);addinto(dest,'_rv_block');drop('_rv_block')
 def m2(dest,k):
  put(dest,dict(P[f'A{k}_56']));madd(dest,f'A{k}_49','S',Q(324000),(0,0,-1,0,-1));madd(dest,f'A{k}_49','Delta',Q(36000),(0,0,-1,0,-1));madd(dest,f'A{k}_51','Delta',Q(-675,8),(-1,0,0,-1,0))
 def loadrow(row):
  nonlocal parsepeak
  path=src/'rows'/f'row-{row}-source-operands.json';assert path in inputs and sha(path)==inputs[path];raw=read(path);assert raw['row']==row;count=scalars(raw)+scalars(pr)+scalars(rootdoc)+scalars(api)+scalars(events)+scalars(results)+6000;parsepeak=max(parsepeak,count);assert count<META
  assert all(not k.startswith('_rv_')for k in raw['row_operands'])
  for k,z in raw['row_operands'].items():put(k,literal(z))
  del raw
 def stage8(row):
  loadrow(row);m2('_rv_BI',row);pmul('_rv_H4','Q','_rv_BI');drop('_rv_BI');madd('_rv_H4','U'+str(row),'B31',Q(-1),(0,0,0,0,0))
  # Different association: combine U/V curvature first, then multiply B83 once.
  put('_rv_T',{tuple(e[j]+(1,0,0,0,-2)[j]for j in range(5)):q*Q(-1,40)for e,q in P['U'+str(row)].items()});put('_rv_TV',{tuple(e[j]+(0,0,0,0,-2)[j]for j in range(5)):q*Q(-405,8)for e,q in P['V'+str(row)].items()});addinto('_rv_T','_rv_TV');drop('_rv_TV');pmul('_rv_QT','Q','_rv_T');drop('_rv_T');madd('_rv_H4','_rv_QT','B83',Q(1),(0,0,0,0,0));drop('_rv_QT');madd('_rv_H4','C'+str(row),'G',Q(-6075,2),(-1,0,0,0,0));pmul('_rv_F6','Q','_rv_H4');drop('_rv_H4');madd('_rv_F6','F'+str(row),'B104',Q(-765450000),(-1,0,0,-1,0));pmul('_rv_N7','P7','_rv_F6');drop('_rv_F6');madd('_rv_N7','N'+str(row),'B51',Q(3645,64),(-2,0,1,0,1));pmul('_rv_N8','E','_rv_N7');drop('_rv_N7');madd('_rv_N8','M'+str(row),'A52_56',Q(-273375,32768),(-3,0,0,0,2));return '_rv_N8'
 def verify(node,record,kind):
  file=pp.parent/record['file'];assert sha(file)==record['sha256'];blob=file.read_bytes();assert blob[:8]==b'Q5CELL1\n';old=decode(blob[9:])if blob[8]else{'content':'0','shift':[0]*5,'terms':[]};expected=poly(old);owned(6*len(old['terms'])+8*len(expected)+128);assert P[node]==expected and len(expected)==record['terms'];assert(not blob[8]and blob==b'Q5CELL1\n\x00')or b'Q5CELL1\n\x01'+encode(old)==blob;assert vhash(old)==record['value_sha256'];results.append({'row':record['row'],'kind':kind,'terms':len(expected),'all_QQ_coefficients_formula_equal':True,'old_codec_full_roundtrip_wire_equal':True,'valuehash':record['value_sha256']});return len(expected)
 try:
  cpath=src/'common-source-operands.json';assert cpath in inputs;common=read(cpath)['common_operands'];assert all(not k.startswith('_rv_')for k in common);parsepeak=max(parsepeak,scalars(common)+scalars(pr)+scalars(rootdoc)+scalars(api)+6000);assert parsepeak<META
  for k,z in common.items():put(k,literal(z))
  del common
  for k in(31,83,14,35,104):m2('B'+str(k),k)
  put('G',dict(P['B14']));put('_rv_sixB35',{e:6*q for e,q in P['B35'].items()});addinto('G','_rv_sixB35');drop('_rv_sixB35');put('B51',dict(P['A51_56']));put('_rv_sum',dict(P['A11_56']));addinto('_rv_sum','S');pmul('_rv_sumQ','Q','_rv_sum',Q(320,27),(1,0,-1,0,-1));addinto('B51','_rv_sumQ');drop('_rv_sum','_rv_sumQ');pmul('_rv_K','P7','F9')
  node=stage8(35);verify(node,auxExpected,'shared_N8_35_56');put('_rv_shared',P.pop(node));kept={'Q','P7','E','S','Delta','A52_56','B31','B83','B104','G','B51','_rv_K','_rv_shared'}
  for k in list(P):
   if k not in kept:drop(k)
  for record in recordsExpected:
   row=record['row'];node=stage8(row)
   for k in list(P):
    if k not in kept and k!=node:drop(k)
   pmul('_rv_N9','_rv_K',node);drop(node);v=oldAPI[row];oldfile=oldAPIpath.parent/v['file'];assert oldfile in inputs and sha(oldfile)==v['sha256'];b=oldfile.read_bytes();old=decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]};put('_rv_Z82',poly(old));del b,old;madd('_rv_N9','_rv_Z82','_rv_shared',Q(-62208000000),(0,0,1,1,2));drop('_rv_Z82');verify('_rv_N9',record,'ninth_residual');drop('_rv_N9');(out/f'checked-through-row-{row}.json').write_bytes(enc({'all_passed':True,'verified_residual_rows':[z['row']for z in results if z['kind']=='ninth_residual'],'shared_verified':True}));gc.collect();assert time.monotonic()-start<55
  assert len(results)==16 and sum(z['terms']for z in results if z['kind']=='ninth_residual')==28886 and results[0]['terms']==113
  (out/'all16-distinct-staged-QQ-formulas-and-backadds.json').write_bytes(enc({'all_passed':True,'results':results,'staged_backadd_records':events,'distinct_association':'combine U/V curvature before B83 product; direct QQ coefficients and separate M2/shared row, no integer known-content DAG evaluator','same_agent_alternate_arithmetic_validation_not_extra_external_peer':True,'independent_source_code_peer_sha256':sha(rootpeer)}));base=sum(p.stat().st_size for p in inputs);assert base+sum(p.stat().st_size for p in out.iterdir())+131072<32*1024*1024;report={'id':ID,'all_passed':True,'all15_residual_and_shared113_QQ_coefficients_equal':True,'residual_terms':28886,'shared_terms':113,'all16_binary_roundtrips_and_valuehashes':True,'distinct_QQ_products':products,'exact_stage_backadds':backchecks,'owned_numeric_peak':peak,'metadata_parse_peak':parsepeak,'prearithmetic_scaled_integer_bit_upper':bits,'same_agent_distinct_arithmetic':True,'producer_not_imported_or_reexecuted':True,'all168_RHS_87plus66_other14_oneP290_complements_kept':True,'no_matrix_rank_or_inverse_claim':True,'actual_READ_RAW_bytes':base,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
 except Exception as exc:
  (out/'failure.json').write_bytes(enc({'all_passed':False,'error':repr(exc),'verified_rows':results,'no_automatic_retry':True,'owned_peak':peak,'seconds':time.monotonic()-start}));raise
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');main(p.parse_args())
