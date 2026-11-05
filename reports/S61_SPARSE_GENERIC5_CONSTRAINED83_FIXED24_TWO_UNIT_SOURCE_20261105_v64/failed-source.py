#!/usr/bin/env python3
"""Exact two source units, typed22x22 Schur, no polynomial coefficient update."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TWO_UNIT_SOURCE_20261105_v64'
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
ZERO=(0,0,0,0,0)
def mono(q,e):return{tuple(e):Q(q)}if q else{}
def mm(A,B):
 R={}
 for e,a in A.items():
  for f,b in B.items():
   g=tuple(x+y for x,y in zip(e,f));R[g]=R.get(g,Q(0))+a*b
 return{e:q for e,q in R.items()if q}
def mn(A):return{e:-q for e,q in A.items()}
def ma(A,B):
 R=dict(A)
 for e,q in B.items():R[e]=R.get(e,Q(0))+q
 return{e:q for e,q in R.items()if q}
def mi(A):assert len(A)==1;e,q=next(iter(A.items()));return mono(1/q,[-x for x in e])
def serial(A):return[[str(q),*e]for e,q in sorted(A.items())]
def bit_profile(blob):
 at=0
 def u():
  nonlocal at
  v=s=0
  while True:
   b=blob[at];at+=1;v|=(b&127)<<s
   if not b&128:return v
   s+=7
 def z(read_value=False):
  nonlocal at
  tag=u();n=tag//2;assert at+n<=len(blob)
  val=int.from_bytes(blob[at:at+n],'big')if read_value else 0;at+=n
  return(-val if tag&1 else val,n)
 num,_=z(True);den,_=z(True);q=Q(num,den);shift=[]
 for _ in range(5):
  v=u();shift.append(v//2 if not v&1 else-(v//2)-1)
 n=u();maxbytes=0
 for _ in range(n):
  _,nb=z();maxbytes=max(maxbytes,nb)
  for _ in range(5):u()
 assert at==len(blob);return n,q,maxbytes*8
begin=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False)
inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
planrp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64/report.json';assert sha(planrp)=='a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a';rr=read(planrp);inputs[planrp]=sha(planrp);planp=planrp.parent/'actual-fixed24-entry-source-and-exponent-plan.json';assert sha(planp)==rr['evidence_sha256'][planp.name];inputs[planp]=sha(planp);v=read(planp);rows=v['rows24_fixed_before_any_rank_or_matching_test'];cols=v['columns24'];entries={(z['row'],z['column']):z for z in v['actual_entry_sources576']};del v
crp=ROOT/'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64/report.json';assert sha(crp)=='0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371';cr=read(crp);inputs[crp]=sha(crp);cp=crp.parent/'whole-cached-rows11-44-exact-comparison.json';assert sha(cp)==cr['evidence_sha256'][cp.name];inputs[cp]=sha(cp);comparison=read(cp);delta={z['column']:z['difference_terms']for z in comparison['all24_complete_actual_row_comparisons']};del comparison
peak=numeric(entries)+numeric(delta)+25000;assert peak<150000
for z in entries.values():
 if z['terms']:
  p=ROOT/z['source']['file'];h=rr['input_sha256'][str(p.relative_to(ROOT))];assert sha(p)==h;inputs[p]=h
def get(z):
 if not z['terms']:return b''
 ref=z['source'];p=ROOT/ref['file']
 if ref['kind']=='unchanged20_column_literal':
  with p.open('rb')as f:f.seek(ref['offset']);b=f.read(ref['bytes'])
  assert hashlib.sha256(b).hexdigest()==ref['payload_sha256'];return b
 b=p.read_bytes();assert sha(p)==ref['sha256']and b[:9]==b'Q5CELL1\n\x01';return b[9:]
def mon(z):
 if not z['terms']:return{}
 assert z['terms']==1;p=decode(get(z));assert len(p['terms'])==1;t=p['terms'][0];return mono(Q(p['content'])*int(t[0]),[p['shift'][j]+t[j+1]for j in range(5)])
a=mon(entries[44,49]);bb=mon(entries[44,51]);assert len(delta[51])==1;d=mono(delta[51][0][0],delta[51][0][1:]);assert a==mono('-1/324000',[0,0,1,0,1])and d==mono('8/675',[1,0,0,1,0]);ia,id=mi(a),mi(d)
Uinv=[[ia,mn(mm(mm(ia,bb),id))],[{},id]];U=[[a,bb],[{},d]]
for i in range(2):
 for j in range(2):
  val={}
  for k in range(2):val=ma(val,mm(U[i][k],Uinv[k][j]))
  assert val==mono(int(i==j),ZERO)
  val={}
  for k in range(2):val=ma(val,mm(Uinv[i][k],U[k][j]))
  assert val==mono(int(i==j),ZERO)
restrows=[i for i in rows if i not in(11,44)];restcols=[j for j in cols if j not in(49,51)];assert len(restrows)==len(restcols)==22
coeff={};unitprofiles=[]
for row in restrows:
 c49=mon(entries[row,49]);c51=mon(entries[row,51]);alpha=mn(mm(c49,ia));beta=ma(mm(mm(mm(c49,ia),bb),id),mn(mm(c51,id)));coeff[row]=(alpha,beta);unitprofiles.append({'row':row,'minus_Ai49_over_a':serial(alpha),'Ai49_b_over_ad_minus_Ai51_over_d':serial(beta),'both_scalar_QQ_only':all(e==ZERO for e in list(alpha)+list(beta))})
# Operand coefficient heights, never any residual polynomial is constructed.
profiles={}
for key,z in entries.items():
 if not z['terms']:profiles[key]={'terms':0,'content':'0','primitive_bits':0,'lower5':None,'upper5':None};continue
 n,q,bits=bit_profile(get(z));assert n==z['terms'];profiles[key]={'terms':n,'content':str(q),'primitive_bits':bits,'lower5':z['lower5'],'upper5':z['upper5']}
dstats={}
for col,terms in delta.items():
 if not terms:dstats[col]={'terms':0,'content':'0','primitive_bits':0,'lower5':None,'upper5':None};continue
 D=1
 for t in terms:q=Q(t[0]);D=D*q.denominator//gcd(D,q.denominator)
 M=max(abs(Q(t[0]).numerator*(D//Q(t[0]).denominator))for t in terms);dstats[col]={'terms':len(terms),'content':str(Q(1,D)),'primitive_bits':M.bit_length(),'lower5':[min(t[j+1]for t in terms)for j in range(5)],'upper5':[max(t[j+1]for t in terms)for j in range(5)]}
# Delta's actual coefficients stay only in immutablecp, not copied to output.
del delta;gc_peak=numeric(entries)+numeric(profiles)+numeric(dstats)+numeric(unitprofiles)+25000;peak=max(peak,gc_peak);assert peak<150000
roots=[];budgets=[];possible_mono=[]
for row in restrows:
 alpha,beta=coeff[row]
 for col in restcols:
  operands=[('actual',profiles[row,col],mono(1,ZERO)),('row44',profiles[44,col],alpha),('row11_minus8row44',dstats[col],beta)];active=[]
  for kind,p,scale in operands:
   for e,q in scale.items():
    if p['terms']:active.append((kind,p,e,q*Q(p['content'])))
  D=1
  for _,p,e,q in active:D=D*q.denominator//gcd(D,q.denominator)
  M=sum(abs(q.numerator)*(D//q.denominator)*((1<<p['primitive_bits'])-1)for _,p,e,q in active);Ubound=sum(p['terms']for _,p,e,q in active)
  if active:
   lower=[min(p['lower5'][j]+e[j]for _,p,e,q in active)for j in range(5)];upper=[max(p['upper5'][j]+e[j]for _,p,e,q in active)for j in range(5)];expbytes=sum(max(1,(upper[j]-lower[j]).bit_length()+6)//7+1 for j in range(5))
  else:lower=upper=[0]*5;expbytes=5
  inputn=profiles[row,col]['terms']+profiles[44,col]['terms']+dstats[col]['terms'];monon=max([p['terms']for _,p,e,q in active]+[0]);owned=25000+inputn+4*Ubound+4*monon+512
  binarycap=2048+((D.bit_length()+7)//8)+Ubound*((M.bit_length()+7)//8+expbytes+16)
  root={'row':row,'column':col,'op':'actual_source_entry_plus_two_shared_row_callbacks','actual_entry':entries[row,col]['source'],'row44_entry':entries[44,col]['source'],'row44_coefficient':serial(alpha),'delta_row11_minus8row44_column':col,'delta_source':str(cp.relative_to(ROOT)),'delta_coefficient':serial(beta)};roots.append(root)
  budgets.append({'row':row,'column':col,'active_scaled_operand_count':len(active),'exact_support_union_cardinality_upper':Ubound,'integer_common_denominator':str(D),'integer_numerator_abs_upper':str(M),'output_lower5':lower,'output_upper5':upper,'three_operand_input_terms':inputn,'future_owned_with_parser_and_backadd_upper':owned,'uncompressed_binary_output_upper':binarycap,'coefficient_product_pairs':0,'scalar_monomial_times_polynomial_only':True})
  if len(active)==1 and active[0][1]['terms']==1:possible_mono.append({'row':row,'column':col,'reason':'Single nonzero cached monomial operand times a nonzero source monomial; no cancellation or new coefficient calculation needed.'})
peak=max(peak,numeric(entries)+numeric(profiles)+numeric(dstats)+numeric(roots)+numeric(budgets)+25000);assert peak<150000
operator={'all_passed':True,'source_only':True,'field':'QQ(H,L,s0,q0,t0)','fixed24_rows':rows,'fixed24_columns':cols,'unimodular_row_difference':'row11 <- row11-8*row44, determinant1','pivot_row_order':[44,11],'pivot_columns':[49,51],'actual_triangular_U2':[[serial(x)for x in row]for row in U],'actual_parameter_monomial_U2_inverse':[[serial(x)for x in row]for row in Uinv],'all8_left_right_inverse_identities_checked':True,'residual_rows22':restrows,'residual_columns22':restcols,'source_coefficients22':unitprofiles,'all484_typed_residual_entries':roots,'residual_formula':'Mij=Aij-(Ai49/a)*A44j+(Ai49*b/(a*d)-Ai51/d)*(A11j-8*A44j). No residual coefficient was evaluated.','substitution':'x51=h_delta/d-sum_j(delta_j/d)*xj; x49=(h44-b*h_delta/d)/a+sum_j(-A44j/a+b*delta_j/(a*d))*xj.','RHS':'h_delta=h11-8*h44; residual h_i=h_i-(Ai49/a)*h44+(Ai49*b/(a*d)-Ai51/d)*h_delta. Retain all unselected87top and66constraint RHS equations.','all_actual_entries_and20unchanged_column_refs_preserved':True,'all177_rows_and66_zero_constraint_consistency_preserved':True,'chart':'Old base monomial units*P290*detA*detB and H*s0*q0*t0 nonzero; record any added coordinate divisor as complement. No normal-coordinate inverse.','no_fixed24_basis_actual_rank24_or_inverse_claim':True,'no_new_row_or_weight_search':True,'no_new_polynomial_update':True,'further_exact_source_monomial_callbacks':possible_mono}
(out/'actual-unit2-and-typed22x22-Schur.json').write_bytes(enc(operator));(out/'future-cell-support-height-ownership-RAW-upper.json').write_bytes(enc({'all_passed':True,'source_only':True,'all484_cell_forecasts':budgets,'not_a_frozen_implementation_or_arithmetic_GO':True,'future_worker_obligations':'Stream one cell+two shared operand callbacks; common-base exact monomial scaling, add/backadd/valuehash, canonical primitive binary roundtrip; source scalar and coefficient-byte bounds must be included. Full input/reference controller must be released before native.','future_all484_max_owned_upper':max(z['future_owned_with_parser_and_backadd_upper']for z in budgets),'future_all484_binary_output_upper':sum(z['uncompressed_binary_output_upper']for z in budgets)}))
ledger=[{'path':str(p.relative_to(ROOT)),'sha256':h,'bytes':p.stat().st_size}for p,h in inputs.items()];raw=sum(z['bytes']for z in ledger)+sum(p.stat().st_size for p in out.rglob('*')if p.is_file());assert raw<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'typed_residual_shape':[22,22],'exact_two_parameter_monomial_units':True,'no_polynomial_coefficient_update':True,'further_single_operand_monomial_callbacks':len(possible_mono),'max_future_cell_owned_upper':max(z['future_owned_with_parser_and_backadd_upper']for z in budgets),'future_binary_output_upper_all484':sum(z['uncompressed_binary_output_upper']for z in budgets),'source_owned_peak':peak,'actual_source_RAW_before_small_report':raw,'input_sha256':{z['path']:z['sha256']for z in ledger},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
