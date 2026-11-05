#!/usr/bin/env python3
"""P7/C28 source and positive bounds only. No signed G8 third response."""
import argparse,hashlib,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_SOURCE_20261105_v63'
FILES={
 'canonical_source':'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/report.json',
 'canonical_input':'reports/AMBIENT_B18_ELL2_NATIVE_CANONICAL33_INPUT_GATE_V2_20261105_v63/canonical33-to-B17-C68-G8-source-input-and-response-budget.json',
 'canonical_source_peer':'reports/AMBIENT_B18_NATIVE_CANONICAL68_ROOT_SOURCE_PEER_20261105_v63/report.json',
 'canonical_values_peer':'reports/AMBIENT_B18_NATIVE_CANONICAL68_DIRECT_REVIEW_20261105_v63/report.json',
 'MS_source_values':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/report.json',
 'MS_source_peer':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63_REVIEW/report.json',
 'MS_kappa':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase3-actual-kappa9-and-normalized-nilframe.json',
 'MS_monic':'reports/AMBIENT_B18_CPRIME_ELL2_NILIMAGE_VALUES_20261105_v63/phase1-source150-readback-fingerprints.json',
 'table':'reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_PREFLIGHT_20261105_v63/exact-existing-raw33-multiplication-table.json',
 'plane_peer':'reports/AMBIENT_B18_ELL2_T2_UPPER617_BOUNDARY_PLANE_SOURCE_20261105_v63_REVIEW/report.json'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def count(x):
 if isinstance(x,dict):return sum(count(k)+count(v) for k,v in x.items())
 if isinstance(x,(list,tuple)):return sum(map(count,x))
 if isinstance(x,str):
  try:F(x);return 1
  except (ValueError,ZeroDivisionError):return 0
 return int(isinstance(x,(int,float,F)))
def main(out):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));out.mkdir(parents=True,exist_ok=False);pins={};D={};raw=0;maxbits=0;opcount=0;cursor={}
 def dump(n,x):(out/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
 for k,p in FILES.items():
  q=ROOT/p;pins[p]=sha(q);raw+=q.stat().st_size;D[k]=json.loads(q.read_text())
 for k in ['canonical_source','canonical_source_peer','canonical_values_peer','MS_source_values','MS_source_peer','plane_peer']:assert D[k]['all_passed']
 assert pins[FILES['canonical_source']]=='48ecdc989f97b4d59a383d87feddaa59b671c86a6391b2a13c9177a8137820cd'
 assert pins[FILES['canonical_input']]==D['canonical_source']['evidence_sha256'][Path(FILES['canonical_input']).name]
 assert pins[FILES['MS_source_values']]=='cffac50257d126e9f0b4cae908228ec7d1b8ea2f547d5765cf331bdfd6633782'
 assert pins[FILES['MS_kappa']]==D['MS_source_values']['evidence_sha256'][Path(FILES['MS_kappa']).name]
 assert pins[FILES['MS_monic']]==D['MS_source_values']['evidence_sha256'][Path(FILES['MS_monic']).name]
 assert count(D)+40000<150000 and raw+4_000_000<32*1024*1024
 dump('phase.json',{'phase':'literal_inputs_and_source_owner_PASS_before_ring_or_positive_products','no_signed_G8_third_response':True})
 ci=D['canonical_input'];tab=D['table'];basis=[(0,0),(1,0),(0,1),(2,0),(1,1),(3,0),(2,1)];index={v:i for i,v in enumerate(basis)}
 sel=[0,8,9,24,25,29,30];proj={b:i for i,b in enumerate(sel)}
 assert [tab['basis33'][i] for i in sel]==[[0,a,e] for a,e in basis]
 columns=[[] if j not in proj else [[proj[j],'1']] for j in range(33)]
 bt={(i,j):index[(a+c,e+f)] for i,(a,e) in enumerate(basis) for j,(c,f) in enumerate(basis) if (a+c,e+f) in index}
 for i in range(33):
  for j in range(33):
   actual={proj[k]:v for k,v in tab['table_columns_for_each_left_basis'][i][j] if k in proj}
   expected={bt[proj[i],proj[j]]:'1'} if i in proj and j in proj and (proj[i],proj[j]) in bt else {}
   assert actual==expected,(i,j,actual,expected)
 kappa=D['MS_kappa']['kappa9_C132'];nu33=D['MS_monic']['nu_A_B33']
 assert kappa==ci['actual_kappa9_C132_cached_UNEVALUATED'] and nu33==ci['actual_nu_A_B33_cached_UNEVALUATED']
 K=[[q[4*b+m] for b in sel for m in range(4)] for q in kappa]
 nu=[[q[b] for b in sel] for q in nu33]
 assert all(q[0]=='0' for q in K) and all(q[:3]==['0','0','0'] for q in nu)
 # Weighted source proof: base and trace weights1; ν starts at base degree2.
 # Largest basis weight3+2=5. Hence C28 radical^6=0, κ^6=0 and d^6=0.
 # The old independently verified radical7 polynomial is therefore filtered
 # at totalκdegree5, which is exact here; no numerical G8 value is computed.
 oldprogram=ci['G8_formal_kappa9_polynomials_mod_radical7']
 program=[]
 for p in oldprogram:
  assert all(len(m)==9 for m,a,b in p)
  program.append([[m,a,b] for m,a,b in p if sum(m)<=5])
 assert len(program)==8 and sum(map(len,program))<1000
 payload={'physical_R0_variables':ci['physical_R0_variables'],'raw33_to_P7_columns':columns,'P7_raw33_selector_indices':sel,
  'P7_basis_r_ell':basis,'P7_unital_multiplication_entries':[[i,j,k] for (i,j),k in bt.items()],
  'all1089_raw33_P7_product_identities_checked':True,
  'canonical_branch':ci['canonical30_branch'],'free4_basis_mask_order':ci['free4_basis_mask_order'],
  'C28_index_rule':'4*P7_basis+mask; same tauA,tauB monics and ordered physicalR0',
  'actual_kappa9_C28_cached_inputs_UNEVALUATED_AS_G8':K,'actual_nu_two_P7_cached_inputs':nu,
  'source_G8_order':ci['source_G8_order'],'G8_formal_kappa9_polynomials_mod_radical6':program,
  'radical6_proof':'wt(r,ell,tauA,tauB)=1, nu plane order>=2, max basisweight5. No weight-lowering carry; radical6=0, P=5+d² inverse three terms.',
  'requested_response_inventory':{'G8':{'base_monomial':'rA^2*ell','base7_index':6,'all8_generators_and4_masks':32},'monic':{'base7_index':6,'two_scalars':2}},
  'no_pure_r3_or_full_G8_C28_response_requested':True,
  'no_honest_quotient_N3_or_individual_pivot_carry_descent_claim':True,
  'future_value_worker_prepared_or_authorized':False}
 # Positive-only C28 structure: retain the original ν operands as absolute
 # envelopes. No signed product/sum of κ or G8 response operands is evaluated.
 Z=(1,0)
 def ck(n):
  nonlocal maxbits
  maxbits=max(maxbits,abs(n).bit_length())
  if abs(n).bit_length()>384:
   dump('positive-guard-frontier.json',{'cursor':cursor,'uncommitted_positive_integer':str(n),'source_inputs':pins,'no_signed_response_evaluated':True});raise RuntimeError('positive384 STOP')
  return n
 def ex(x):
  q=F(x);return (ck(q.denominator),ck(abs(q.numerator))) if q else Z
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=ck(math.lcm(a[0],b[0]));return d,ck(ck(a[1]*(d//a[0]))+ck(b[1]*(d//b[0])))
 def times(a,b):
  nonlocal opcount
  opcount+=1
  return (ck(a[0]*b[0]),ck(a[1]*b[1])) if a[1] and b[1] else Z
 routes={}
 for x in range(28):
  bx,mx=divmod(x,4)
  for y in range(28):
   by,my=divmod(y,4);a,e=basis[bx];b,f=basis[by];collision=mx&my;outmask=mx^my;r=[]
   if collision==0:
    if (a+b,e+f) in index:r=[(4*index[a+b,e+f]+outmask,ex(1))]
   elif collision in [1,2]:
    for q,(u,v) in enumerate(basis):
     if F(nu[collision-1][q]) and (a+b+u,e+f+v) in index:r.append((4*index[a+b+u,e+f+v]+outmask,ex(nu[collision-1][q])))
   else:assert collision==3 # νAνB starts degree4, zero in P7.
   if r:routes[x,y]=r
 def prod(x,y):
  z={}
  for a,p in x.items():
   for b,q in y.items():
    for k,v in routes.get((a,b),[]):z[k]=plus(z.get(k,Z),times(times(p,q),v))
  return z
 kval=[{i:ex(q) for i,q in enumerate(v) if F(q)} for v in K];powers={(j,1):v for j,v in enumerate(kval)}
 for j,v in enumerate(kval):
  for n in range(2,6):cursor.update({'phase':'positive_cached_kappa_power','kappa':j,'power':n});powers[j,n]=prod(powers[j,n-1],v)
 bound=[]
 for g,p in enumerate(program):
  z={}
  for term,(m,a,b) in enumerate(p):
   cursor={'phase':'positive_G8_polynomial_bound_ONLY','G8':g,'term':term};v={0:ex(F(a,b))}
   for j,n in enumerate(m):
    if n:v=prod(v,powers[j,n])
   for k,q in v.items():z[k]=plus(z.get(k,Z),q)
  bound.append([[str(z.get(24+m,Z)[0]),str(z.get(24+m,Z)[1])] for m in range(4)])
 assert maxbits<=384 and 3*maxbits+1<=1153
 payload['positive_G8_requested32_D_U']=bound
 payload['positive_monic_requested2_D_U']=[[str(a),str(b)] for a,b in map(ex,[nu[0][6],nu[1][6]])]
 payload['positive_propagation_operations_upper']=opcount
 payload['positive_committed_bits']=maxbits;payload['positive_temporary_bits']=3*maxbits+1
 # Include all JSON/parsed input copies, at most45 κ powers, current terms,
 # coefficient outputs+codec, integer scratch, source metadata and guards.
 n=count(payload);future=2*(n+100)+count([K,nu,powers,routes])+2*34+6*28*2+8000
 assert future<60000
 payload['prospective_value_whole_numeric_bound']=future
 source=count([D,payload,kval,powers,routes,bound,program,pins])+25000
 assert source<150000 and raw+4_000_000<32*1024*1024
 ep='native-P7-C28-canonical-third-source-and-positive34-response-input.json';dump(ep,payload)
 rep={'id':STEM,'all_passed':True,'source_and_positive_only_no_signed_G8_third':True,
  'all1089_P7_products_and_same_canonical_C132_monic_branch_bound':True,'P7_dimension':7,'C28_dimension':28,
  'radical6_and_P_inverse_three_terms_from_literal_plane_monic_orders':True,
  'requested_G8_r2ell_masks32_and_monic2_only':True,'positive_committed_bits':maxbits,'positive_temporary_bits':3*maxbits+1,
  'source_owned_numeric_upper':source,'future_value_numeric_upper':future,'whole_RAW_read_plus_output_upper':raw+4_000_000,
  'future_value_worker_prepared_or_authorized':False,'input_sha256':pins,'evidence_sha256':{ep:sha(out/ep)},
  'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 dump('report.json',rep);dump('phase.json',{'phase':'FINAL_SOURCE_ONLY','report_sha256':sha(out/'report.json')});print(json.dumps({k:rep[k] for k in ['all_passed','positive_committed_bits','positive_temporary_bits','source_owned_numeric_upper','future_value_numeric_upper','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
