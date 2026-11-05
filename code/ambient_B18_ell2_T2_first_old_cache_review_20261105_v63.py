#!/usr/bin/env python3
"""Independent sparse-Rell4 substitution of saved source/copy data; no solves."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_first_X2_source_fold_review_20261105_v63 import jetmul
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale
ROOT=Path(__file__).resolve().parents[1];A='AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_20261105_v63';STEM=A+'_REVIEW'
G='AMBIENT_B18_ELL2_T2_FIRST_OLD_CACHE_INPUT_GATE_20261105_v63'
P={'primary':f'reports/{A}/report.json','values':f'reports/{A}/actual104-source-rell-values-and-two-by619-old-knownparts.json',
'gate_report':f'reports/{G}/report.json','gate':f'reports/{G}/actual104-old-C0-ancestor-cache-bindings-and-Rell4-gate.json',
'gate_peer':f'reports/{G}_REVIEW/report.json',
'formal':'reports/AMBIENT_B18_ELL2_T2_FIRST_OLD_AD_20261105_v63/actual-two-by619-old-formal-rell-second-sensitivities.json',
'formal_peer':'reports/AMBIENT_B18_ELL2_T2_FIRST_X2_SOURCE_FOLD_20261105_v63_REVIEW/report.json',
'common':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase1-actual105-B56-cache-composition.json',
'base':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63/phase2-all105-equations-and-primitive56-source-columns.json',
'common_peer':'reports/AMBIENT_B18_ELL2_V_CONORMAL_EXISTING_CACHE_VALUES_20261105_v63_REVIEW/report.json',
'FJ':'reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json',
'C0':'reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json',
'quartic':'reports/AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63/quotient-firsts-and-actual-D0-firstjets.json',
'atom':'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/phase1-actual-atom-firsts.json',
'mu_gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json'}
MON=[(0,0),(1,0),(0,1),(1,1)];ONE={(0,0):Q(1)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def decode(v):return {m:Q(c) for m,c in zip(MON,v) if Q(c)}
def encode(v):return [str(v.get(m,Q(0))) for m in MON]
def inverse(v):
 c=v[(0,0)];h={m:a for m,a in v.items() if m!=(0,0)};out=add(scale(ONE,1/c),scale(h,-1/c**2),scale(jetmul(h,h),1/c**3));assert jetmul(v,out)==ONE;return out

def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'INDEPENDENT_SAVED_SOURCE_AND_HONEST_SYMBOLIC_SUBSTITUTION','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k not in {'FJ','C0','quartic'}}
 for k in ['primary','gate_report','gate_peer','formal_peer','common_peer']:assert D[k]['all_passed']
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 E=D['values'];I=D['gate'];nodes=dict(I['source_node_closure']);basis=list(map(tuple,D['common']['basis56']));assert basis==list(map(tuple,D['base']['basis56']))
 targetindex={m:i for i,m in enumerate(MON)};images=[]
 for m in basis:images.append({(m[7],m[8]):Q(1)} if not any(m[:7]) and (m[7],m[8]) in targetindex else {})
 ix={m:i for i,m in enumerate(basis)}
 for i,a in enumerate(basis):
  for j,b in enumerate(basis):
   k=ix.get(tuple(x+y for x,y in zip(a,b)));assert jetmul(images[i],images[j])==(images[k] if k is not None else {})
 def project(v):return add(*(scale(images[i],c) for i,c in sp(v).items()))
 known={i:project(v) for i,v in D['base']['all164_original_source_roots']}
 for i,v in zip(D['common']['unknown105_order'],D['common']['all105_cached_values_B56']):
  p=project(v)
  if i in known:assert p==known[i]
  known[i]=p
 assert D['base']['all105_equations_exact'] and all(not v for v in D['base']['all105_equation_residuals'])
 fj=fields(ROOT/P['FJ'],['both57_source2','source2_monomials','all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']);cb=fj['both57_source2']['CB']
 def evaluate(roots,bind,opaque=None):
  opaque=opaque or {};need=set();todo=list(roots)
  while todo:
   i=todo.pop()
   if i in need:continue
   need.add(i)
   if i in bind or i in opaque:continue
   n=nodes[i]
   if n[0] in ('a','m'):todo.extend(n[1:3])
   elif n[0] in ('n','i'):todo.append(n[1])
   else:assert n[0]=='c',(i,n)
  vals={};grads={}
  for i in sorted(need):
   if i in bind:v=bind[i];g={}
   elif i in opaque:v=opaque[i];g={i:Q(1)}
   else:
    n=nodes[i];op=n[0]
    if op=='c':v={(0,0):Q(n[1],n[2])} if n[1] else {};g={}
    elif op=='a':v=add(vals[n[1]],vals[n[2]]);g=add(grads[n[1]],grads[n[2]])
    elif op=='n':v=scale(vals[n[1]],-1);g=scale(grads[n[1]],-1)
    elif op=='m':
     a,b=n[1:];v=jetmul(vals[a],vals[b]);g=add(scale(grads[a],vals[b].get((0,0),Q(0))),scale(grads[b],vals[a].get((0,0),Q(0))))
    else:assert op=='i';v=inverse(vals[n[1]]);g=scale(grads[n[1]],-1/vals[n[1]][0,0]**2)
   vals[i]=v;grads[i]=g
  return {i:vals[i] for i in roots},{i:grads[i] for i in roots}
 y,_=evaluate(cb['F9_nodes'][7:9],known);args=[{}]*7+[y[i] for i in cb['F9_nodes'][7:9]];bindings=dict(I['copied57_exact_source_bindings'])
 for i,p in cb['all57_source2jets']:
  v={}
  for word,c in zip(fj['source2_monomials'],p):
   if not Q(c):continue
   term=scale(ONE,Q(c))
   for j in word:term=jetmul(term,args[j])
   v=add(v,term)
  known[bindings[i]]=v
 res,_=evaluate(I['copied57_equations'],known);assert all(not v for v in res.values())
 assert [[i,encode(y[i])] for i in cb['F9_nodes'][7:9]]==E['actual_YB7_8_from_common105']
 opaque={}
 for field in ['all792_CA_h_firstjets','all784_CB_h_firstjets_in_YA']:
  for i,p in fj[field].items():opaque[int(i)]=decode([Q(p[0]),Q(p[8])+Q(p[9]),Q(p[9])/5,Q(0)])
 prefix={b['node']:sp(b['origin_and_ds']) for b in D['mu_gate']['all35_ordered_prefix_bindings']}
 for b in D['atom']['all35_ordered_CA_prefix_firsts']:
  q=sp(b['first9_bar']);opaque[b['node']]=decode([q.get(0,Q(0)),q.get(8,Q(0)),prefix[b['node']].get(1,Q(0))/5,Q(0)])
 quartic=fields(ROOT/P['quartic'],['original_regular_function_nodes_unmodified','all8_regular_function_first_gradients_in_FULL9_CA_F9']);regular=quartic['original_regular_function_nodes_unmodified']
 for i,p in zip(regular['old_CB']+regular['new_CA'],quartic['all8_regular_function_first_gradients_in_FULL9_CA_F9']):opaque[i]=decode([Q(0),Q(p[7])+Q(p[8]),Q(p[8])/5,Q(0)])
 for i in fields(ROOT/P['C0'],['actual_C0_basis8_scalar_nodes'])['actual_C0_basis8_scalar_nodes']:known[i]={}
 values,grads=evaluate(I['formal_original_C0_second_labels104'],known,opaque);assert all(not g for g in grads.values())
 assert [[i,encode(values[i]),enc(grads[i])] for i in sorted(values)]==E['all104_original_C0_values_and_opaque_gradients']
 actual=[]
 for row in D['formal']['all619_results']:
  parts=[]
  for p in row['two_formal_old_functionals']:
   v=Q(p['known_lower_first_response_only'])
   for i,c in sp(p['literal_original_C0_second_rell_sensitivity']).items():v+=c*values[i].get((1,1),Q(0))
   parts.append(str(v))
  actual.append(parts)
 assert actual==E['all619_two_old_knownparts_including_frozen_pivot_first_correction'];nz=[i for i,p in enumerate(actual) if any(Q(c) for c in p)];assert len(nz)==6
 ep=o/'independent104-Rell4-source-values-and1238-old-responses.json';dump(o,ep.name,{'inputs':inputs,'all3136_ring_products_and_all105_equations_homomorphic_from_accepted_actualB56':True,'actualYB7_8_and_all57_copied_equations_verified':True,'all104_originalC0_values_and_zero_opaque_second_gradients_verified':True,'all1238_old_responses_and_actualpivot_corrections_verified':True,'nonzero_old_response_indices':nz,'all_six_old_nonzero_vectors':[[i,actual[i]] for i in nz],'new_source_or_implicit_or_gauge_coefficients':0,'native_contributions_remain_external':True,'actual_T2_HF1_and_length_not_computed':True})
 report={'id':STEM,'all_passed':True,'all104_old_source_rell_values_and1238_partial_responses_accepted':True,'complete_native_or_T2_HF1_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t,'max_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps({'all_passed':True,'report_sha256':sha(o/'report.json'),'seconds':report['seconds'],'nonzero_old_rows':nz}))
if __name__=='__main__':main()
