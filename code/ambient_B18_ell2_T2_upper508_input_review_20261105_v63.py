#!/usr/bin/env python3
"""Independent dense nonnegative-polynomial proof of ONE508 zero envelopes.

No signed upper response and no physical/honest function jet is evaluated.
"""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_UPPER508_INPUT_GATE_20261105_v63';STEM=A+'_REVIEW'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(v):
 if isinstance(v,dict):return sum(count(k)+count(x) for k,x in v.items())
 if isinstance(v,(list,tuple)):return sum(map(count,v))
 if isinstance(v,str):
  try:F(v);return 1
  except(ValueError,ZeroDivisionError):return 0
 return int(isinstance(v,(int,float,F)))
def fields(p,names):
 txt=p.read_text();dec=json.JSONDecoder();i=txt.index('{')+1;z={}
 while True:
  while txt[i].isspace() or txt[i]==',':i+=1
  if txt[i]=='}':break
  k,i=dec.raw_decode(txt,i)
  while txt[i].isspace() or txt[i]==':':i+=1
  v,i=dec.raw_decode(txt,i)
  if k in names:z[k]=v
  del v
 assert set(z)==set(names);return z
def main(out):
 t0=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};raw=0;maximum=0;operations=0
 def read(p,names=None):
  nonlocal raw
  q=ROOT/p;pins[p]=sha(q);raw+=2*q.stat().st_size
  return json.loads(q.read_text()) if names is None else fields(q,names)
 def dump(n,v):(out/n).write_text(json.dumps(v,separators=(',',':'))+'\n')
 rp=f'reports/{A}/report.json';r=read(rp);assert r['all_passed'] and pins[rp]=='806eeeea7454edbf4c34ee95ab45fe1ce985dcef7aca976a4b85f84376c2419d'
 ep='ONE508-positive-full-cache-input-and-honest-extension-certificate.json';d=read(f'reports/{A}/{ep}');assert pins[f'reports/{A}/{ep}']==r['evidence_sha256'][ep]
 sp='reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63/ONE508-source-closure-B21-and-exact-carry-zeros.json';s=read(sp)
 pr='reports/AMBIENT_B18_ELL2_T2_UPPER508_SOURCE_20261105_v63_REVIEW/report.json';peer=read(pr);assert peer['all_passed'] and pins[pr]=='8da6bd72d7ca9cf1f62b80d3a273d04b406bbf6302a8780cdce82a498a28c95f'
 for k in ['ordinary_basis_B21','Gamma_two','D_two_by10_by56','literal_C0_program38','full_typed_program86','Gamma_old595_source_roots','D_old595_source_roots']:assert d[k]==s[k],k
 assert s['literal_original768_source']==[[595,'1']] and s['a1_all27times9_exact_empty'] and s['native_c0_exact_empty'] and s['same617_source_recipe']['subtract_actual_basis563609']==['0','0']
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json'
 ff=read(fp,['first640'])['first640']['all79_atom_values_full10'];firsts={i:dict(v) for i,v in ff};del ff
 pp='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json';pv=read(pp,['all_requested164_root_values'])['all_requested164_root_values'];phys={i:dict(v) for i,v in pv};del pv
 sel=s['raw33_to_B21_selector21'];catalog={};physical=[];formal=[]
 for i,kind,v in d['physical164_actual_B21_and_only79_truefirst_extension_catalog39']:
  f=[firsts[i].get(j,'0') for j in range(10)];assert f[0]=='0'
  if i in phys:
   assert kind=='original_physical164_actual_B21' and list(map(F,v))==[F(phys[i].get(j,'0')) for j in sel];physical.append(i)
  else:
   assert kind=='whole_functional_higher0_extension' and list(map(F,v))==list(map(F,f))+[F(0)]*11;formal.append(i)
  assert list(map(F,v[:10]))==list(map(F,f));catalog[i]=v
 assert physical==d['physical_function_ids23'] and formal==d['arbitrary_higher_extension_ids16'] and len(physical)==23 and len(formal)==16
 assert sorted(catalog)==s['function_leaf39_ids']
 j1=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD_J1_FORMAL_20261105_v63/617-formal-direct-old-third-and-second-symbol-sensitivities.json',['all617_formal_directOLD_J0_Jr_Jell'])['all617_formal_directOLD_J0_Jr_Jell'];j=next(v['two_formal_J0_Jr_Jell_symbol_rows'] for v in j1 if v['free_case619']==508)
 o2=read('reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_FORMAL_20261105_v63/actual617-formal-combined-OLD2-second-symbol-coefficients.json',['all617_combined_OLD2_independent_symbol_coefficients'])['all617_combined_OLD2_independent_symbol_coefficients'];o=next(v for i,v in o2 if i==508)
 assert j==d['original_case508_J1_certificate'] and o==d['original_case508_OLD2_certificate']
 for row in j:
  for i,co in row:
   if i in formal:assert i in [299650,505938] and co[0]=='0'
 for row in o:
  for pair,i,q in row:
   if i in formal:assert i in [299650,505938] and pair in [8,9]
 for i in [299650,505938]:assert all(F(q)==0 for q in catalog[i][:10])
 order=read('reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW/report.json');assert order['all_passed']
 qpath='reports/AMBIENT_B18_ELL2_T2_UPPER617_FULL_COORDINATE_SOURCE_20261105_v63/known-Psi-point-correction-and-full-b1-source-annihilation.json';Q=read(qpath,['Q_two6'])['Q_two6'];assert Q==d['fixed_Q_embedding_two6']
 owner0=count([r,d,s,peer,firsts,phys,catalog,j1,o2,Q,pins])+20000;assert owner0<150000 and raw+1_000_000<32*1024*1024
 dump('phase-before-positive-replay.json',{'phase':'literal_cache_all39firsts23physical_and_whole_extension_scope_PASS','source90_temporary271':True,'owned_bound':owner0,'no_signed_function_response':True})
 # Dense21 exact nonnegative rationals differ from the producer's sparse(D,U)
 # implementation; all operations are positive upper bounds only.
 def ck(q):
  nonlocal maximum
  assert q>=0;maximum=max(maximum,q.numerator.bit_length(),q.denominator.bit_length());assert maximum<=90;return q
 def add(a,b):return ck(a+b)
 def mul(a,b):
  nonlocal operations
  operations+=1;return ck(a*b)
 b=[tuple(v) for v in d['ordinary_basis_B21']];bi={v:i for i,v in enumerate(b)}
 table=[[bi.get(tuple(x+y for x,y in zip(a,c))) for c in b] for a in b]
 def plus(a,c):return [add(x,y) for x,y in zip(a,c)]
 def scale(a,q):return [mul(x,q) for x in a]
 def product(a,c):
  z=[F(0)]*21
  for i,x in enumerate(a):
   if not x:continue
   for k,y in enumerate(c):
    if y and table[i][k] is not None:q=table[i][k];z[q]=add(z[q],mul(x,y))
  return z
 vals={i:[ck(abs(F(q))) for q in v] for i,v in catalog.items()};C={};T={}
 for i,n in d['literal_C0_program38']:
  if n[0]=='c':z=[ck(abs(F(n[1],n[2])))]+[F(0)]*20
  elif n[0]=='n':z=list(vals[n[1]] if n[1] in vals else C[n[1]])
  elif n[0]=='a':z=plus(vals[n[1]] if n[1] in vals else C[n[1]],vals[n[2]] if n[2] in vals else C[n[2]])
  else:assert n[0]=='m';z=product(vals[n[1]] if n[1] in vals else C[n[1]],vals[n[2]] if n[2] in vals else C[n[2]])
  C[i]=z
 for i,n in d['full_typed_program86']:
  if n[0]=='actual_C0_scalar_function':z=list(vals[n[1]] if n[1] in vals else C[n[1]])
  elif n[0]=='sum':
   z=[F(0)]*21
   for j,a,c in n[1]:z=plus(z,scale(T[j],ck(abs(F(a,c)))))
  else:
   assert n[0]=='product';z=[ck(abs(F(n[1],n[2])))]+[F(0)]*20
   for j in n[3]:z=product(z,T[j])
  T[i]=z
 for key,own in [('positive_C0_B21',C),('positive_typed_B21',T)]:
  for i,cols in d[key]:
   want={k:F(u,den) for k,den,u in cols};assert own[i]==[want.get(k,F(0)) for k in range(21)]
 Gm=[dict(v) for v in d['Gamma_two']];D=d['D_two_by10_by56'];gp=[];dp=[]
 for t in range(2):
  a=F(0);u=F(0)
  for row,(kind,node) in d['Gamma_old595_source_roots']:
   v=T[node] if kind=='typed' else vals[node] if node in vals else C[node];a=add(a,mul(ck(abs(F(Gm[t][row]))),v[20]))
  for row,node in d['D_old595_source_roots']:
   for pair in range(10):u=add(u,mul(ck(abs(F(D[t][pair][row]))),T[node][10+pair]))
  gp.append(a);dp.append(u)
 assert gp==dp==[F(0),F(0)]
 total=[add(x,y) for x,y in zip(gp,dp)];six=[]
 for k in range(6):
  q=F(0)
  for t in range(2):q=add(q,mul(ck(abs(F(dict(Q[t]).get(k,'0')))),total[t]))
  six.append(q)
 for key,own in [('positive_Gamma_direct_cubic_two',gp),('positive_D_full_firstcross_second_two',dp),('positive_actual_upper_two',total),('positive_full_six',six)]:assert [F(u,den) for den,u in d[key]]==own==[F(0)]*len(own)
 for p0,ha in r['input_sha256'].items():
  if p0 not in pins:raw+=(ROOT/p0).stat().st_size;pins[p0]=sha(ROOT/p0)
  assert pins[p0]==ha,p0
 owner=count([r,d,s,peer,firsts,phys,catalog,j1,o2,Q,pins,vals,C,T,table,gp,dp,total,six])+15000
 assert owner<150000 and raw+1_000_000<32*1024*1024
 ep2='independent-dense-positive124-ONE508-zero-envelope-proof.json';e={'all39_truefull9firsts_and23_actual_physical_B21_inputs_verified':True,
  'whole_functional_extension_of16_only_higher_slots_certified':True,'actual_conductor_seconds_only_rr_rell_order3zero_thirdcoefficient0':True,
  'all38_C0_86typed_positive_polynomials_match':True,'both_directGamma_cubic_and_D_firstcross_envelopes_exact0':True,
  'all_two_upper_six_embedding_positive_envelopes0':True,'signed_zero_only_value_job_unnecessary':True,
  'actual_ONE508_upper_class_zero_follows_from_nonnegative_bound':True,'only_case508_no_other605_HF_or_length_claim':True,
  'no_signed_function_response_or_upper_evaluation':True,'own_positive_bits':maximum,'positive_multiplications':operations,'owned_upper':owner}
 dump(ep2,e);rep={'id':STEM,'all_passed':True,'positive_zero_proof_only':True,'input_sha256':pins,'evidence_sha256':{ep2:sha(out/ep2)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'whole_numeric_upper':owner,'whole_RAW_upper':raw+1_000_000,'seconds':time.monotonic()-t0}
 dump('report.json',rep);print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'owned':owner,'positive_bits':maximum,'seconds':rep['seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
