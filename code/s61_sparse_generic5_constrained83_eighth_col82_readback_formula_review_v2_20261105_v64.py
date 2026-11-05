#!/usr/bin/env python3
"""Alternate staged QQ Laurent verification of all16 completed col82 cells."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import argparse,json,hashlib,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import encode,decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_READBACK_FORMULA_REVIEW_V2_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def literal(z):return{tuple(e[j]+z['transform_monomial_exponent5'][j]for j in range(5)):q*Q(z['QQ_transform_scalar'])for e,q in poly(z['original_literal']).items()}
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def stats(A):
 d=1
 for q in A.values():d=d*q.denominator//gcd(d,q.denominator)
 return d,max([abs((q*d).numerator)for q in A.values()]+[0])
def valuehash(p):
 h=hashlib.sha256();h.update((str(len(p['terms']))+'\n').encode());q=Q(p['content'])
 for t in p['terms']:h.update((str(q*int(t[0]))+' '+' '.join(str(p['shift'][j]+t[j+1])for j in range(5))+'\n').encode())
 return h.hexdigest()
def main(a):
 assert a.check
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);pp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64/report.json';assert sha(pp)=='f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb';pr=read(pp);assert pr['all_passed'];rp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_ROOT_SOURCE_REVIEW_20261105_v64/report.json';assert sha(rp)=='0dc56acb9944d845dc537302f29cda8a5d57608d873fe815496f4a513fb9a64a';assert read(rp)['all_passed']
 inputs={ROOT/p:h for p,h in pr['input_sha256'].items()};inputs.update({pp:sha(pp),rp:sha(rp),Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')})
 for name,h in pr['evidence_sha256'].items():inputs[pp.parent/name]=h
 for p,h in inputs.items():assert sha(p)==h
 src=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_UNIT_COL82_SOURCE_20261105_v64';idx=read(src/'row-stream-input-and-plan-index16.json');api=read(pp.parent/'actual-cleared-col82-and-commonQ3RE-representation.json');vals={z['row']:z for z in api['actual16_col82_numerator_values']};P={};events=[];records=[];peak=40000;bits=0;muls=0;backs=0
 def owned(extra=0):
  nonlocal peak
  n=40000+7*sum(len(A)for A in P.values())+extra;peak=max(peak,n);assert n<=150000;return n
 def drop(*ks):
  for k in ks:P.pop(k,None)
 def put(k,A):P[k]=A;owned()
 def addmaps(A,B,sign=Q(1)):
  dA,mA=stats(A);dB,mB=stats(B);d=dA*dB//gcd(dA,dB);M=mA*(d//dA)+abs(sign.numerator)*mB*(d//dB)*sign.denominator
  # All callers use integer sign, so this is also an unreduced temporary bound.
  assert sign.denominator==1;R=dict(A)
  for e,q in B.items():
   R[e]=R.get(e,Q(0))+sign*q
   if not R[e]:del R[e]
  assert all(d%q.denominator==0 and abs((q*d).numerator)<=M for q in R.values());return R
 def pmul(dest,Aname,Bname,q=Q(1),shift=(0,0,0,0,0)):
  nonlocal muls,bits
  assert dest not in(Aname,Bname);drop(dest);A=P[Aname];B=P[Bname];S=set()
  for e in A:
   for f in B:
    S.add(tuple(e[j]+f[j]+shift[j]for j in range(5)));owned(12*len(S)+64)
  dA,mA=stats(A);dB,mB=stats(B);D=dA*dB*q.denominator;M=min(len(A),len(B))*mA*mB*abs(q.numerator);bits=max(bits,D.bit_length(),M.bit_length());R={}
  for e,x in A.items():
   for f,y in B.items():
    g=tuple(e[j]+f[j]+shift[j]for j in range(5));R[g]=R.get(g,Q(0))+q*x*y;z=R[g];assert D%z.denominator==0 and abs((z*D).numerator)<=M
  R={e:z for e,z in R.items()if z};assert set(R)<=S;owned(5*len(S)+7*len(R)+64);del S;put(dest,R);muls+=1
 def addinto(dest,block,sign=Q(1)):
  nonlocal backs
  old=P[dest];B=P[block];R=addmaps(old,B,sign);owned(14*len(R)+7*len(old)+64);back=addmaps(R,B,-sign);assert back==old;P[dest]=R;backs+=1;events.append({'destination':dest,'block':block,'sign':str(sign),'exact_backadd':True});del back;owned()
 def madd(dest,A,B,q,shift):
  pmul('_block',A,B,q,shift);addinto(dest,'_block');drop('_block')
 def m2(dest,k):
  put(dest,dict(P[f'A{k}_82']));madd(dest,f'A{k}_49','S',Q(324000),(0,0,-1,0,-1));madd(dest,f'A{k}_49','Delta',Q(36000),(0,0,-1,0,-1));madd(dest,f'A{k}_51','Delta',Q(-675,8),(-1,0,0,-1,0))
 commonraw=read(src/'common-source-operands.json')['common_operands'];assert all(not k.startswith('_stage_')for k in commonraw);P.update({k:literal(z)for k,z in commonraw.items()});assert numeric(commonraw)+7*sum(len(z)for z in P.values())<40000;del commonraw
 for k in(31,83,14,35,104):m2('B'+str(k),k)
 put('G',dict(P['B14']));put('_sixB35',{e:6*q for e,q in P['B35'].items()});addinto('G','_sixB35');drop('_sixB35')
 put('B51',dict(P['A51_82']));put('_sum',dict(P['A11_82']));addinto('_sum','S');pmul('_qsum','Q','_sum',Q(320,27),(1,0,-1,0,-1));addinto('B51','_qsum');drop('_sum','_qsum')
 fixedkeys=set(P)
 for ref in idx['row_records16']:
  row=ref['row'];raw=read(ROOT/ref['row_operand_file']);plan=read(ROOT/ref['plan_file']);assert sha(ROOT/ref['row_operand_file'])==ref['row_operand_sha256']and sha(ROOT/ref['plan_file'])==ref['plan_sha256'];assert numeric(raw)+numeric(plan)+numeric(events)+numeric(records)+7*sum(len(literal(z))for z in raw['row_operands'].values())+8000<40000
  assert all(not k.startswith('_stage_')for k in raw['row_operands'])
  for k,z in raw['row_operands'].items():put(k,literal(z))
  del raw;m2('_stage_BI',row);pmul('_stage_H4','Q','_stage_BI');madd('_stage_H4','U'+str(row),'B31',Q(-1),(0,0,0,0,0));pmul('_qu','Q','U'+str(row));madd('_stage_H4','_qu','B83',Q(-1,40),(1,0,0,0,-2));drop('_qu');pmul('_qv','Q','V'+str(row));madd('_stage_H4','_qv','B83',Q(-405,8),(0,0,0,0,-2));drop('_qv','_stage_BI')
  put('_stage_H5',P.pop('_stage_H4'));madd('_stage_H5','C'+str(row),'G',Q(-6075,2),(-1,0,0,0,0));pmul('_stage_F6','Q','_stage_H5');drop('_stage_H5');madd('_stage_F6','F'+str(row),'B104',Q(-765450000),(-1,0,0,-1,0));pmul('_stage_N7','P7','_stage_F6');drop('_stage_F6');madd('_stage_N7','N'+str(row),'B51',Q(3645,64),(-2,0,1,0,1));pmul('_stage_N8','E','_stage_N7');drop('_stage_N7');madd('_stage_N8','M'+str(row),'A52_82',Q(-273375,32768),(-3,0,0,0,2))
  v=vals[row];file=pp.parent/v['file'];b=file.read_bytes();assert sha(file)==v['sha256']and b[:8]==b'Q5CELL1\n';literal_out=decode(b[9:])if b[8]else{'content':'0','shift':[0]*5,'terms':[]};assert not b[8]or encode(literal_out)==b[9:];assert valuehash(literal_out)==v['actual_QQ_numerator_value_sha256'];expected=poly(literal_out);owned(7*len(expected)+6*len(literal_out['terms'])+64);assert P['_stage_N8']==expected and len(expected)==v['terms'];D=int(plan['integer_common_D']);M=int(plan['integer_numerator_abs_upper']);assert all(D%q.denominator==0 and abs((q*D).numerator)<=M for q in expected.values());records.append({'row':row,'terms':len(expected),'value_sha256':v['actual_QQ_numerator_value_sha256'],'binary_roundtrip_and_canonical_hash':True,'alternate_staged_QQ_Schur_formula_equal':True,'final_source_height_bound':True});del plan,expected,literal_out,b
  for k in list(P):
   if k not in fixedkeys:del P[k]
  gc.collect();assert time.monotonic()-start<55
 assert len(records)==16 and sum(z['terms']for z in records)==4160
 (out/'alternate-staged-QQ-formula-and-backadds.json').write_bytes(enc({'all_passed':True,'records16':records,'all_staged_backadds':events,'QQ_products_checked':muls,'no_cartesian_tuple_producer_import_or_replay':True,'source_peer':'0dc56acb9944d845dc537302f29cda8a5d57608d873fe815496f4a513fb9a64a','same_agent_alternate_arithmetic_validation_not_new_external_peer':True}));base=sum(p.stat().st_size for p in inputs);assert base+sum(p.stat().st_size for p in out.iterdir())+100000<32*1024*1024
 report={'id':ID,'all_passed':True,'actual_cells16':16,'all4160_QQ_coefficients_staged_formula_equal':True,'all16_binary_roundtrips_and_valuehashes':True,'staged_QQ_products':muls,'exact_stage_backadds':backs,'owned_numeric_peak':peak,'derived_temporary_integer_bits_upper':bits,'root_independent_source_peer_hash':'0dc56acb9944d845dc537302f29cda8a5d57608d873fe815496f4a513fb9a64a','same_agent_alternate_arithmetic_validation':True,'all169_RHS_87plus66_other15_oneP290_chart_scope_kept':True,'actual_READ_RAW_bytes':base,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');main(p.parse_args())
