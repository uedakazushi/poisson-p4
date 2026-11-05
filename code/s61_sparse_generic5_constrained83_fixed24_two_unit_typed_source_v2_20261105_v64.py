#!/usr/bin/env python3
"""Exact two source units, typed22x22 Schur, no polynomial coefficient update."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
import json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TWO_UNIT_TYPED_SOURCE_V2_20261105_v64'
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

# GENERAL Ai49/Ai51 are polynomial callbacks, not assumedmonomials.
assert mon(entries[44,81])==bb
assert len(delta[81])==1 and mono(delta[81][0][0],delta[81][0][1:])==d
nodes=[]
def node(op,*args):nodes.append([op,*args]);return len(nodes)-1
roots=[];coefrows=[]
for row in restrows:
 a49=node('actual_poly_ref',entries[row,49]['source'],entries[row,49]['terms']);a51=node('actual_poly_ref',entries[row,51]['source'],entries[row,51]['terms'])
 alpha=node('negative_monomial_scale',serial(ia),a49)
 beta=node('subtract',node('monomial_scale',serial(mm(mm(ia,bb),id)),a49),node('monomial_scale',serial(id),a51))
 coefrows.append({'row':row,'alpha_node':alpha,'beta_node':beta,'Ai49_terms':entries[row,49]['terms'],'Ai51_terms':entries[row,51]['terms']})
 for col in restcols:
  if col==81:root=node('subtract_cached_actual_entries',entries[row,81]['source'],entries[row,51]['source'])
  else:root=node('three_term_callback_sum',entries[row,col]['source'],alpha,entries[44,col]['source'],beta,col,str(cp.relative_to(ROOT)))
  roots.append({'row':row,'column':col,'typed_root':root})
# SOURCE forecasts for the cheap EXACT post-Schur column81 only.
plans=[]
for row in restrows:
 pp=[]
 for col,sgn in [(81,1),(51,-1)]:
  z=entries[row,col]
  if z['terms']:
   n,q,bits=bit_profile(get(z));assert n==z['terms'];pp.append((col,n,sgn*q,bits,z))
 D=1
 for col,n,q,bits,z in pp:D=D*q.denominator//gcd(D,q.denominator)
 M=sum(abs(q.numerator)*(D//q.denominator)*((1<<bits)-1)for col,n,q,bits,z in pp);Ubound=sum(n for _,n,_,_,_ in pp)
 if pp:lower=[min(z['lower5'][j]for _,_,_,_,z in pp)for j in range(5)];upper=[max(z['upper5'][j]for _,_,_,_,z in pp)for j in range(5)]
 else:lower=upper=[0]*5
 eb=sum(max(1,((upper[j]-lower[j]).bit_length()+6)//7)for j in range(5));cap=2048+(D.bit_length()+7)//8+Ubound*((M.bit_length()+7)//8+eb+16)
 owner=25000+Ubound+4*Ubound+4*max([n for _,n,_,_,_ in pp]+[0])+512
 plans.append({'row':row,'column':81,'formula':'actual(row,81)-actual(row,51)','operand_refs':[{'column':col,'terms':n,'QQ_content_with_sign':str(q),'primitive_integer_bit_upper':bits,'source':z['source']}for col,n,q,bits,z in pp],'support_union_cardinality_upper':Ubound,'integer_common_denominator':str(D),'integer_numerator_abs_upper':str(M),'actual_exponent_lower5':lower,'actual_exponent_upper5':upper,'binary_RAW_upper':cap,'whole_owned_upper':owner,'no_multiplication_of_polynomials':True})
# Full22 remains factored; no unauthorizedpolynomial update/search.
operator={'all_passed':True,'source_only':True,'fixed24_rows':rows,'fixed24_columns':cols,'unimodular_row_difference':'row11 <- row11-8row44','pivot_rows':[44,11],'pivot_columns':[49,51],'actual_triangular_U2':[[serial(x)for x in row]for row in U],'actual_U2_inverse_parameter_monomials':[[serial(x)for x in row]for row in Uinv],'all8_left_right_inverse_identities_checked':True,'residual_rows22':restrows,'residual_columns22':restcols,'typed_nodes':nodes,'all484_actual_residual_roots':roots,'general_Ai49_Ai51_polynomial_coefficient_callbacks':coefrows,'residual_formula':'Mij=Aij-(Ai49/a)A44j+(Ai49*b/(ad)-Ai51/d)*delta_j','column81_exact_simplification':'(A44,81,delta81)=(b,d)=U2*e2, so remainingM_i81=A_i81-A_i51 EXACTLY. No coefficient update evaluated.','column81_source_vector_identity_checked':True,'substitution':'x51=h_delta/d-sum(delta_j/d)xj; x49=(h44-b*h_delta/d)/a+sum(-A44j/a+b*delta_j/(ad))xj','RHS':'h_delta=h11-8h44; h_remaining_i=h_i-(Ai49/a)h44+(Ai49*b/(ad)-Ai51/d)h_delta','all24_actual_entries_retained_before_two_unit_operation':True,'all177_rows_other20_sourcecols_and66_zero_RHS_constraints_retained':True,'chart':'Old parameter-monomial units*P290*detA*detB and H*s0*q0*t0; coordinate complement obligations explicit. No unknown-normal inverse.','no_new_polynomial_coefficient_update_or_inverse_claim':True,'no_fixed24_rank24_basis_claim':True,'previous_SOURCE_only_assumption_failure_preserved':'reports/S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_TWO_UNIT_SOURCE_20261105_v64/failure.json'}
(out/'actual-unit2-and-typed22x22-general-callback-Schur.json').write_bytes(enc(operator));(out/'column81-source-support-bit-and-storage-plans.json').write_bytes(enc({'all_passed':True,'source_only':True,'all22_exact_difference_plans':plans,'no_coefficients_updated':True,'future_implementation_and_storage_addendum_required':True}))
ledger=[{'path':str(p.relative_to(ROOT)),'sha256':h,'bytes':p.stat().st_size}for p,h in inputs.items()];sourcebytes=sum(z['bytes']for z in ledger);outputbytes=sum(p.stat().st_size for p in out.rglob('*')if p.is_file());future=sourcebytes+outputbytes+sum(p['binary_RAW_upper']for p in plans)+2*1024*1024;ownedfuture=max(p['whole_owned_upper']for p in plans);peak=max(peak,numeric(entries)+numeric(delta)+numeric(nodes)+numeric(plans)+numeric(roots)+25000);assert peak<150000 and future<32*1024*1024
report={'id':ID,'all_passed':True,'source_only':True,'typed_residual_shape':[22,22],'all_actual8_unit_inverse_identities':True,'column81_exact_cached_difference_simplification':True,'no_polynomial_coefficient_update':True,'max_general_Ai49_terms':max(entries[row,49]['terms']for row in restrows),'max_general_Ai51_terms':max(entries[row,51]['terms']for row in restrows),'source_owned_peak':peak,'future_column81_owned_upper':ownedfuture,'future_column81_RAW_upper':future,'implementation_addendum_required_before_coefficients':True,'input_sha256':{z['path']:z['sha256']for z in ledger},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()if p.is_file()},'seconds':time.monotonic()-begin};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
