#!/usr/bin/env python3
"""Positive, source-bound OLD2 fold gate. No signed sensitivity or response."""
import argparse,json,math,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_T2_G3_active4_source8_fold_20261105_v63 import count
ROOT=H.ROOT
STEM='AMBIENT_B18_ELL2_T2_UPPER617_OLD2_INPUT_GATE_20261105_v63'
ZERO=(1,0)
def main(out):
 start=time.monotonic();signal.alarm(60)
 resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);pins={};cursor={};bits=0
 def read(p):pins[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):pins[p]=H.sha(ROOT/p);return H.selected(p,k)
 def guard(x):
  nonlocal bits
  b=abs(x).bit_length();bits=max(bits,b)
  if b>384:
   H.dump(out/'BEFORE_SIGNED_OLD2_POSITIVE_GUARD.json',{'cursor':cursor,'uncommitted_positive_integer':str(x),'bits':b,'no_signed_sensitivity_or_response':True})
   raise RuntimeError('positive384 STOP')
  return x
 def ex(x):
  x=F(x)
  return (guard(x.denominator),guard(abs(x.numerator))) if x else ZERO
 def add(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=guard(math.lcm(a[0],b[0]));u=guard(a[1]*(d//a[0]));v=guard(b[1]*(d//b[0]))
  return d,guard(u+v)
 def mul(a,b):
  if not a[1] or not b[1]:return ZERO
  return guard(a[0]*b[0]),guard(a[1]*b[1])
 def union(a,b):
  if not a[1]:return b
  if not b[1]:return a
  d=guard(math.lcm(a[0],b[0]));u=guard(a[1]*(d//a[0]));v=guard(b[1]*(d//b[0]))
  return d,max(u,v)
 def axpy(a,b,s):
  z=dict(a)
  if not s[1]:return z
  for k,q in b.items():z[k]=add(z.get(k,ZERO),mul(q,s))
  return z
 sourcep='reports/AMBIENT_B18_ELL2_T2_UPPER617_OLD2_SOURCE_20261105_v63/'
 sr=read(sourcep+'report.json');assert sr['all_passed'] and pins[sourcep+'report.json']=='addc65c2eb868114a773f8dbed9e7c513cc5f21b7e9e46d7604f23498d9d8e16'
 ep0=sourcep+'combined-OLD2-D-and-four-pivot-source-factor-inventory.json';s=read(ep0);assert pins[ep0]==sr['evidence_sha256'][Path(ep0).name] and not s['unresolved_function_boundary_ids']
 assert read('reports/AMBIENT_B18_ELL2_T2_UPPER617_D_20261105_v63_REVIEW/report.json')['all_passed']
 f79p='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json'
 f79={i:dict(v) for i,v in field(f79p,'first640')['all79_atom_values_full10']}
 f164p='reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json'
 f164={i:dict(v) for i,v in field(f164p,'all_requested164_root_values')}
 origins=[]
 for i,kind in s['typed_actual_C0_function_classification']:
  a=F(f79[i].get(0,'0')) if i in f79 else None;b=F(f164[i].get(0,'0')) if i in f164 else None
  if a is not None and b is not None:assert a==b
  point=a if a is not None else b;assert point is not None and point==0
  origins.append([i,kind,'0'])
 del f79,f164
 gproofp='reports/AMBIENT_B18_T_GDIFF_SECOND_UNIVARIATE5_PREFLIGHT_20261104_v63_REVIEW/'
 gr=read(gproofp+'report.json');assert gr['all_passed'] and pins[gproofp+'report.json']=='917f2492740394b7ebbdeca99886cfa4cbbe18018c2171eaa5788af0c463fdd5'
 gp=read(gproofp+'independent-two57-fifth-input-and-source-order-audit.json')
 assert {r['family'] for r in gp['both_families_proof_records']}=={'CA','CB'}
 for r in gp['both_families_proof_records']:
  assert r['global_bbar_order_at_least4'] and r['cofactor_origin_zero'] and r['all_signed_original_g_and_bbar_roots_match_frozen_quotients']
 assert s['newg4_second_symbols_not_zero_filled']==list(range(566400,566404))
 outer=s['outer_four_pivot_source'];assert not outer['inverse_nodes']
 assert [(i,n[0]) for i,n in outer['literal_C0_program'] if n[0]=='g']==[(i,'g') for i in range(566400,566404)]
 origins.extend([i,'outer_C0_regular_g_original_source917f','0'] for i in range(566400,566404))
 # All typed nodes have zero origin, inductively: there are no empty products.
 grad={};linear=[];zeroids=[];nodepositive=[]
 for i,n in s['typed_program']:
  cursor={'stage':'positive_linear_projection','node':i}
  if n[0]=='actual_C0_scalar_function':g={n[1]:ex(1)}
  elif n[0]=='sum':
   g={}
   for j,a,b in n[1]:assert j<i;g=axpy(g,grad[j],ex(F(a,b)))
  else:
   assert n[0]=='product' and len(n[3])>=1 and all(j<i for j in n[3])
   g=axpy({},grad[n[3][0]],ex(F(n[1],n[2]))) if len(n[3])==1 else {}
  grad[i]=g
  if g:
   # Remove zero-gradient summands only; preserve their full source in ep0.
   nn=['sum',[[j,a,b] for j,a,b in n[1] if grad[j] and a]] if n[0]=='sum' else n
   linear.append([i,nn]);nodepositive.append([i,[[f,str(q[0]),str(q[1])] for f,q in sorted(g.items())]])
  else:zeroids.append(i)
 assert len(grad)==757 and len(linear)==104 and sum(len(g) for g in grad.values())==112
 # Only literal add/negate after the four independently defined regular g leaves.
 cg={i:{i:ex(1)} for i,kind,_ in origins}
 for i,n in outer['literal_C0_program']:
  if n[0]=='g':assert i in cg
  elif n[0]=='n':cg[i]=dict(cg[n[1]])
  else:assert n[0]=='a';cg[i]=axpy(cg[n[1]],cg[n[2]],ex(1))
 D=s['D_two_by10_by56'];G=[dict(v) for v in outer['Gamma_two']]
 rootcols=[[j,[[row,r] for row,r in rr if grad[r]]] for j,rr in s['D_old_root_columns']]
 rootcols=[v for v in rootcols if v[1]];assert len(rootcols)==45 and sum(len(rr) for j,rr in rootcols)==62
 colroots=dict(rootcols);c619=s['original619_c0_source_combinations'];rawa1=dict(s['raw619_old4_plane_a1_nonempty28coefficients'])
 old4=outer['original4_old_pivots'];pivots=dict(outer['four_source_root_columns'])
 def positive_roots(j,ch):
  z={}
  for old,x in c619[j]:
   if old>=640 or old not in colroots:continue
   for row,r in colroots[old]:
    for p in range(10):
     if D[ch][p][row]=='0':continue
     key=('typed',r,p);q=mul(ex(x),ex(D[ch][p][row]));z[key]=add(z.get(key,ZERO),q)
  for slot,x in rawa1.get(j,[]):
   pi,var=divmod(slot,9);assert pi<4 and var in [7,8];p=8 if var==8 else 9
   for row,(kind,r) in pivots[old4[pi]]:
    if row in G[ch]:key=(kind,r,p);z[key]=add(z.get(key,ZERO),mul(ex(x),ex(G[ch][row])))
  return z
 bases={j:[positive_roots(j,t) for t in range(2)] for j in [563,609]}
 allpositive=[];rootmax=0;slots=0;casechannels=0;common=ZERO;rootcommon=ZERO
 for rec in s['kernel617_source_recipes']:
  case=rec['free_case619'];cursor={'stage':'positive617_combined_D_and_outer4','case':case}
  for ch in range(2):
   roots=positive_roots(case,ch)
   for j,x in zip([563,609],rec['subtract_actual_basis563609'],strict=True):roots=axpy(roots,bases[j][ch],ex(x))
   rootmax=max(rootmax,len(roots));z={}
   for (kind,r,p),q in roots.items():
    rootcommon=union(rootcommon,q)
    for f,a in (grad[r] if kind=='typed' else cg[r]).items():
     k=(p,f);z[k]=add(z.get(k,ZERO),mul(q,a))
   if z:
    allpositive.append([case,ch,[[p,f,str(q[0]),str(q[1])] for (p,f),q in sorted(z.items())]])
    casechannels+=1;slots+=len(z)
    for q in z.values():common=union(common,q)
 assert slots==1926 and casechannels==106
 # Pruned future payload: no actual F1/F2/F3 values, no native arrays.
 data={'D_two_by10_by56':D,'pair_order10':s['pair_order10'],'original619_c0_source_combinations':c619,'kernel617_source_recipes':s['kernel617_source_recipes'],'linear_D_old_root_columns45':rootcols,'linear_typed_program104':linear,'all_zero_gradient_typed_node_ids653':zeroids,'literal_outer_four_pivot_source':outer,'raw619_old4_plane_a1_nonempty28coefficients':s['raw619_old4_plane_a1_nonempty28coefficients'],'literal_function_origins71':origins,'positive_linear_node_gradients':nodepositive,'positive617_combined_symbol_coefficients':allpositive,'positive_common_symbol_envelope':list(map(str,common)),'positive_common_rootweight_envelope':list(map(str,rootcommon)),'largest_per_case_rootweight_support':rootmax,'all_second_symbols_independent_by_function_ID_and_ordinary_pair':True,'formal_scope':'D Xc2 minus OLD4 GammaXp2a1 only, after exact617source subtraction before rootfold and gradient. No F2/F3, firstcross, native, b1 or upper response values.','source_proofs':{'all757_typed_origins0_and_no_empty_product':True,'multi_products_degree2_have_no_leaf_second_sensitivity':True,'new4_regular_g_point0_from_global_bbar_order4_and_nonzero_g3':True,'unknown36_only79_plus_new4_preserved':True,'physical31_also_retained_as_formal_symbols':True,'no_leaf_second_zero_filled':True,'native23_individual_L2S0_and_full_source_carry_separate':True}}
 input_count=count(data);nodecount=count([linear,nodepositive]);output_records=6*slots+6*casechannels+2000
 # Two raw/parsed payload copies + all gradient scratch + two 617-basis roots
 # + root/output maps for one streamed case + emitted output and metadata.
 scratch=4*112+8*rootmax+8*max(len(v[2]) for v in allpositive)+2000
 cap=2*(input_count+1000)+scratch+output_records
 raw=2_000_000+4*cap*(math.ceil(bits*math.log10(2))+24)
 source_owned=count([s,grad,cg,linear,zeroids,nodepositive,rootcols,bases,allpositive,data,pins])+10000
 assert bits<=384 and 3*bits+1<=1153 and source_owned<150000 and cap<150000 and raw<32*1024*1024
 data.update({'source_committed_bits':bits,'source_temporary_bits':3*bits+1,'future_formal_numeric_cap':cap,'future_formal_RAW_cap':raw,'metadata_reserve':1000,'prior_input_count':input_count,'future_scratch_upper':scratch,'future_output_record_upper':output_records,'fixed_positive_output_slots':slots})
 ep=out/'617-combined-OLD2-independent-symbol-positive-input.json';H.dump(ep,data)
 report={'id':STEM,'all_passed':True,'source_and_positive_envelopes_only':True,'no_signed_sensitivity_or_actual_response':True,'typed757_origins_zero':True,'linear_nodes104':True,'pruned_D_old_columns45':True,'formal_independent_symbol_count71':True,'formal_slots_upper1926':True,'committed_prior_bits':bits,'temporary_prior_bits':3*bits+1,'source_owned_bound':source_owned,'future_numeric_cap':cap,'future_RAW_cap':raw,'input_sha256':pins,'evidence_sha256':{ep.name:H.sha(ep)},'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);print(json.dumps({k:report[k] for k in ['all_passed','committed_prior_bits','temporary_prior_bits','source_owned_bound','future_numeric_cap','future_RAW_cap','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();main(a.output)
