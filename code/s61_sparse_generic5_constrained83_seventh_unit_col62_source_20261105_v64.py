#!/usr/bin/env python3
"""SOURCE-only row-stream supports/heights for seventeenth next-col62 layer."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd,prod
import argparse,json,hashlib,time,signal,resource,gc
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_UNIT_COL62_SOURCE_20261105_v64'
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
def primitive(A):
 if not A:return{'content':'0','shift':[0]*5,'terms':[]}
 D=1
 for q in A.values():D=D*q.denominator//gcd(D,q.denominator)
 Z={e:(q*D).numerator for e,q in A.items()};G=0
 for z in Z.values():G=gcd(G,abs(z))
 lo=[min(e[j]for e in A)for j in range(5)]
 return{'content':str(Q(G,D)),'shift':lo,'terms':[[str(Z[e]//G),*[e[j]-lo[j]for j in range(5)]]for e in sorted(A,key=lambda e:(sum(e),e),reverse=True)]}
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def add(A,B,k=Q(1)):
 R=dict(A)
 for e,q in B.items():
  R[e]=R.get(e,Q(0))+k*q
  if not R[e]:del R[e]
 return R
def mm(A,B):return A[0]*B[0],tuple(x+y for x,y in zip(A[1],B[1]))
def main(a):
 assert a.build
 begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);(out/'rows').mkdir();inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def ev(p,r,name):
  f=p.parent/name;assert sha(f)==r['evidence_sha256'][name];inputs[f]=sha(f);return read(f)
 p,s7=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL74_SEVENTH_UNIT_SOURCE_20261105_v64','23dc07db98639cc1a418192064e922d6c467ea8003b11b4b4dc55113dc1c5319');op=ev(p,s7,'actual-row51-seventh-unit-and-full170x17-operator.json');cost=ev(p,s7,'remaining17-literal-source-cost-comparison.json');assert cost['selected_minimum_in_this_frozen_source_model']==62;rows=op['selected17_row_ids'];j=62
 p,s6=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_UNIT_COL74_SOURCE_20261105_v64','19a897e17f89aa55e68a0487d5455425f397b169e07b1cb8c6766fed54f8feb4');oldcommon=ev(p,s6,'common-source-operands.json')['common_operands'];oldindex=ev(p,s6,'row-stream-input-and-plan-index18.json');unit=ev(p,s6,'old-unit2-monomial-source-prefactors.json');oldrefs={z['row']:z for z in oldindex['row_records18']}
 ap,ar=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_COL74_IMPLEMENTATION_GATE_V2_20261105_v64','a3e8e6b84de0400f2980ac5bd541f2a6ab28270031f73aee17af79e77a9a02d1');ac=ev(ap,ar,'implementation-custody-storage-bit-addendum.json');assert ac['source_custody_counts']=={'raw':67,'explicit':2,'delta':1,'reuse':54,'F':18};del ac
 cp,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SIXTH_COL74_COEFFICIENTS_20261105_v64','96c24538f7f6dfce75b78a25cf755895116023b9fd46861d8a04410373d467cd');api=ev(cp,cr,'actual-cleared-col74-and-commonQ2-representation.json');Fsix={z['row']:z for z in api['actual18_col74_numerator_values']}
 ep,er=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');E={(z['row'],z['column']):z for z in ev(ep,er,'actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 checked=set()
 def raw(k,c):
  z=E[k,c]
  if not z['terms']:return{'content':'0','shift':[0]*5,'terms':[]},{'kind':'sourcezero','row':k,'column':c}
  ref=z['source'];f=ROOT/ref['file'];assert ref['kind']=='unchanged20_column_literal'
  if f not in checked:assert sha(f)==er['input_sha256'][str(f.relative_to(ROOT))];inputs[f]=sha(f);checked.add(f)
  with f.open('rb')as h:h.seek(ref['offset']);b=h.read(ref['bytes'])
  assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];v=decode(b);assert len(v['terms'])==z['terms'];return v,ref
 inner=[(Q(z['QQ']),tuple(z['exponent5']))for z in unit['all_four_inner_prefactors']];common={};cpfs={}
 def profile(z):
  p=z['original_literal'];s=[p['shift'][h]+z['transform_monomial_exponent5'][h]for h in range(5)];ts=[(int(t[0]),tuple(s[h]+t[h+1]for h in range(5)))for t in p['terms']];return{'n':len(ts),'q':Q(p['content'])*Q(z['QQ_transform_scalar']),'S':{e for q,e in ts},'M':max([abs(q)for q,e in ts]+[0])}
 def commonleaf(k,p,src):
  z={'original_literal':p,'source':src,'QQ_transform_scalar':'1','transform_monomial_exponent5':[0]*5};common[k]=z;cpfs[k]=profile(z)
 for k in('Q','Q2'):common[k]=oldcommon[k];cpfs[k]=profile(common[k])
 commonleaf('P7',{'content':'1','shift':[0]*5,'terms':[['5',0,0,0,2,0],['6',0,0,0,1,0],['1',0,0,0,0,0]]},{'kind':'explicit_Q_times_qplus1'})
 z,ref=raw(44,j);commonleaf('S',z,ref);d=add(poly(raw(11,j)[0]),poly(z),Q(-8));commonleaf('Delta',primitive(d),{'kind':'fresh_raw_A11_minus8A44','column':j})
 for k in(31,83,14,35,104):
  for c in(j,49,51):
   z,ref=raw(k,c);commonleaf(f'A{k}_{c}',z,ref)
 for k in(51,11):
  z,ref=raw(k,j);commonleaf(f'A{k}_{j}',z,ref)
 commonpath=out/'common-source-operands.json';commonpath.write_bytes(enc({'all_passed':True,'parameter_order':['H','L','s0','q0','t0'],'common_operands':common}));index=[];alltuples=0;caps=0;sourcepeak=0;maxowned=0
 fixed=25000+numeric(E)+numeric(op)+numeric(cost)+numeric(oldcommon)+numeric(oldindex)+numeric(unit)+numeric(common)+numeric(api)+4096;assert fixed<150000
 for i in rows:
  oldref=oldrefs[i];rp=ROOT/oldref['row_operand_file'];assert sha(rp)==oldref['row_operand_sha256'];inputs[rp]=sha(rp);oldraw=read(rp)['row_operands'];row={};pf=dict(cpfs)
  for c in(j,49,51):
   k=f'A{i}_{c}'
   if k not in common:
    z,ref=raw(i,c);row[k]={'original_literal':z,'source':ref,'QQ_transform_scalar':'1','transform_monomial_exponent5':[0]*5}
  for prefix in('U','V','C','F'):row[prefix+str(i)]=oldraw[prefix+str(i)]
  v=Fsix[i];f=cp.parent/v['file'];assert sha(f)==v['sha256'];inputs[f]=sha(f);blob=f.read_bytes();z=decode(blob[9:])if blob[8]else{'content':'0','shift':[0]*5,'terms':[]};assert len(z['terms'])==v['terms'];row['N'+str(i)]={'original_literal':z,'source':{'kind':'actual_sixth_col74_numerator','row':i,'path':str(f.relative_to(ROOT)),'sha256':v['sha256'],'value_representation':'S6_i74=N_i/Q^2'},'QQ_transform_scalar':'1','transform_monomial_exponent5':[0]*5}
  pf.update({k:profile(z)for k,z in row.items()});del oldraw,z,blob
  outer=[(i,['P7','Q2'],(Q(1),(0,0,0,0,0))),(31,['P7','Q','U'+str(i)],(Q(-1),(0,0,0,0,0))),(83,['P7','Q2','U'+str(i)],(Q(-1,40),(1,0,0,0,-2))),(83,['P7','Q2','V'+str(i)],(Q(-405,8),(0,0,0,0,-2))),(14,['P7','Q','C'+str(i)],(Q(-6075,2),(-1,0,0,0,0))),(35,['P7','Q','C'+str(i)],(Q(-18225),(-1,0,0,0,0))),(104,['P7','F'+str(i)],(Q(-765450000),(-1,0,0,-1,0)))]
  active=[]
  for k,ol,opre in outer:
   inn=[([f'A{k}_{j}'],inner[0]),([f'A{k}_49','S'],inner[1]),([f'A{k}_49','Delta'],inner[2]),([f'A{k}_51','Delta'],inner[3])]
   for il,ipre in inn:
    ls=ol+il
    if all(pf[l]['n']for l in ls):active.append((k,ls,mm(opre,ipre)))
  for ls,pre in [(['N'+str(i),f'A51_{j}'],(Q(3645,64),(-2,0,1,0,1))),(['N'+str(i),'Q',f'A11_{j}'],(Q(675),(-1,0,0,0,0))),(['N'+str(i),'Q','S'],(Q(675),(-1,0,0,0,0)))]:
   if all(pf[l]['n']for l in ls):active.append((51,ls,pre))
  dc=1;aa=[]
  for k,ls,pre in active:
   q=pre[0]
   for l in ls:q*=pf[l]['q']
   dc=dc*q.denominator//gcd(dc,q.denominator);aa.append((k,ls,pre,q))
  union=set();M=0;products=[];rowtuples=0;peaksets=0
  for k,ls,pre,q in aa:
   # Forecast support-only scratch before each distinct-factor Minkowski sum.
   forecast=fixed+numeric(row)+numeric(products)+numeric(index)+7*sum(z['n']for z in pf.values())+6*len(union)+4096;assert forecast<150000
   S={(0,0,0,0,0)};stages=[]
   for l in ls:
    predicted=len(S)*pf[l]['n'];assert forecast+6*(len(S)+predicted)<150000
    T={tuple(e[h]+f[h]for h in range(5))for e in S for f in pf[l]['S']};stages.append({'factor':l,'previous_support':len(S),'factor_terms':pf[l]['n'],'exact_next_support':len(T)});peaksets=max(peaksets,len(S)+len(T)+len(union));S=T
   S={tuple(e[h]+pre[1][h]for h in range(5))for e in S};union|=S;nt=prod(pf[l]['n']for l in ls);h=(nt//max(pf[l]['n']for l in ls))*prod(pf[l]['M']for l in ls);w=q.numerator*(dc//q.denominator);M+=abs(w)*h;rowtuples+=nt;products.append({'source_M2_row':k,'factor_labels':ls,'full_parameter_monomial_prefactor':{'QQ':str(pre[0]),'exponent5':list(pre[1])},'source_product_content':str(q),'exact_scalar_integer_weight_after_D':str(w),'literal_tuple_count':nt,'per_output_assignment_upper':nt//max(pf[l]['n']for l in ls),'primitive_product_integer_abs_upper':str(h),'support_stages':stages,'exact_product_support_cardinality':len(S)})
  U=len(union);lower=[min(e[h]for e in union)for h in range(5)]if union else[0]*5;upper=[max(e[h]for e in union)for h in range(5)]if union else[0]*5;eb=sum(max(1,((upper[h]-lower[h]).bit_length()+6)//7)for h in range(5));rawcap=2048+(dc.bit_length()+7)//8+U*((M.bit_length()+7)//8+eb+16);nleaves=sum(z['n']for z in pf.values());owned=25000+16*(nleaves+U+max([p['exact_product_support_cardinality']for p in products]+[0]))+2048;assert owned<150000
  whole=fixed+numeric(row)+numeric(products)+numeric(index)+7*nleaves+6*peaksets+4096;sourcepeak=max(sourcepeak,whole);assert whole<150000
  plan={'row':i,'column':j,'value_representation':'S7_i62=Laurent_numerator_over_Qcubed_times_qplus1','integer_common_D':str(dc),'integer_numerator_abs_upper':str(M),'products':products,'exact_source_support_sorted':[list(e)for e in sorted(union,key=lambda e:(sum(e),e),reverse=True)],'exact_union_support_cardinality':U,'expected_shift_lower5':lower,'upper5':upper,'uncompressed_binary_RAW_upper':rawcap,'future_owned_extra_including25kmeta':owned,'future_common_and_row_leaf_terms':nleaves,'no_Q_cancellation_or_new_columnoperation':True}
  rf=out/'rows'/f'row-{i}-source-operands.json';pp=out/'rows'/f'row-{i}-support-height-plan.json';rf.write_bytes(enc({'all_passed':True,'row':i,'row_operands':row}));pp.write_bytes(enc(plan));index.append({'row':i,'row_operand_file':str(rf.relative_to(ROOT)),'row_operand_sha256':sha(rf),'plan_file':str(pp.relative_to(ROOT)),'plan_sha256':sha(pp),'source_tuple_majorant':rowtuples,'binary_RAW_upper':rawcap,'future_owned_upper':owned});alltuples+=rowtuples;caps+=rawcap;maxowned=max(maxowned,owned);del row,pf,aa,active,outer,products,plan,union;gc.collect();assert time.monotonic()-begin<55
  (out/'completed-source-prefix.json').write_bytes(enc({'all_passed':True,'source_only':True,'completed_rows':[z['row']for z in index],'no_valueworker':True}))
 expected=next(z['new17_fixed_source_tuple_majorant']for z in cost['all17_fixed_source_tuple_majorants']if z['column']==j);assert alltuples==expected==894424
 (out/'row-stream-input-and-plan-index17.json').write_bytes(enc({'all_passed':True,'source_only':True,'common_operand_file':str(commonpath.relative_to(ROOT)),'common_operand_sha256':sha(commonpath),'row_records17':index,'all17_source_tuple_majorant':alltuples,'VALUE_must_stream_one_row':True,'full170_top104_plus66_87extra_RHS_scope':True,'other16_columns_still_typed':True,'oneP290_Qcubed_times_qplus1_representation_ONLY':True,'future_immutable_implementation_and_ROOTGO_REQUIRED':True}));base=sum(p.stat().st_size for p in inputs);sourcebytes=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());future=base+sourcebytes+caps+2*1024*1024;assert future<32*1024*1024
 report={'id':ID,'all_passed':True,'source_only':True,'next_fixed_column':j,'source_tuple_majorant':alltuples,'future_coeff_owned_upper':maxowned,'future_coeff_RAW_upper':future,'whole_source_numeric_upper':sourcepeak,'old14MBmatrix_actualREAD_correctly_charged':True,'no_new_column_values_or_worker_authorized':True,'future_immutable_implementation_and_ROOTGO_REQUIRED':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{str(p.relative_to(out)):sha(p)for p in out.rglob('*')if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
