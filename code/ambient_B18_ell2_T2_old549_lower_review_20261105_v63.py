#!/usr/bin/env python3
"""Independent548 functional postcheck and case508 symbolic order proof."""
import argparse,hashlib,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as OWN
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def sp(v):return {int(i):F(c) for i,c in v if F(c)}
def main(out,orderout):
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));out.mkdir(parents=True,exist_ok=False);orderout.mkdir(parents=True,exist_ok=False);hs={}
 def read(p):hs[p]=sha(ROOT/p);return json.loads((ROOT/p).read_text())
 def field(p,k):hs[p]=sha(ROOT/p);return OWN.selected(p,k)
 stem='AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_20261105_v63_REVIEW';base='reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_20261105_v63';primary=read(base+'/report.json');assert hs[base+'/report.json']=='5ac0e171508b0994c9210ca6392affe39fdf73448fe6223bef52a4fe9b729d76' and primary['all_passed'];saved=read(base+'/actual548-first-only-lower-prefix-and-witness-if-any.json');assert primary['evidence_sha256']['actual548-first-only-lower-prefix-and-witness-if-any.json']==hs[base+'/actual548-first-only-lower-prefix-and-witness-if-any.json']
 gp='reports/AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_INPUT_20261105_v63';g=read(gp+'/report.json');d=read(gp+'/pruned548-first-only-cache-source-and-positive-bound.json');assert read(gp+'_REVIEW/report.json')['all_passed'];assert read('reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63_REVIEW/report.json')['all_passed']
 for p,h in primary['input_sha256'].items():assert sha(ROOT/p)==h
 pairs=[tuple(v) for v in d['ordinary_pair_order8']];allowed={()}|{(i,) for i in range(9)}|set(pairs);first={i:sp(v) for i,v in d['actual_C0_leaf_first10']};ops=dict(d['typed_closure_nodes']);maximum=0
 def ck(x):
  nonlocal maximum
  maximum=max(maximum,abs(x.numerator).bit_length(),x.denominator.bit_length());assert maximum<=60;return x
 def add(a,b,c=F(1)):
  v=dict(a)
  for k,x in b.items():
   y=ck(v.get(k,F(0))+ck(c*x))
   if y:v[k]=y
   else:v.pop(k,None)
  return v
 def product(a,b,keep):
  v={}
  for u,x in a.items():
   for w,y in b.items():
    q=tuple(sorted(u+w))
    if keep(q):v=add(v,{q:ck(x*y)})
  return v
 jets={}
 for i,n in sorted(ops.items()):
  if n[0]=='actual_C0_scalar_function':v={(k-1,) if k else ():x for k,x in first[n[1]].items()}
  elif n[0]=='c':v={():F(n[1],n[2])} if n[1] else {}
  elif n[0]=='sum':
   v={}
   for j,a,b in n[1]:v=add(v,jets[j],F(a,b))
  else:
   assert n[0]=='product';v={():F(n[1],n[2])} if n[1] else {}
   for j in n[3]:v=product(v,jets[j],lambda q:q in allowed)
  jets[i]=v
 G=list(map(sp,d['Gamma4_source56_times8']));B={j:dict(v) for j,v in d['oldB_source_rows_pruned']};R={j:sp(v) for j,v in d['actual_pivot_r_first_pruned']};E={j:sp(v) for j,v in d['actual_old_pivot_ell_first_pruned']};L0=list(map(sp,d['native_L0_rows_pruned']));Le=list(map(sp,d['native_L1ell_rows_pruned']));S0={j:sp(v) for j,v in d['native_pivot_S0_columns']};Se={j:sp(v) for j,v in d['native_pivot_S1ell_columns']}
 for j,local in d['native_pivot_local128_ids']:
  q={}
  for matrix,col in [(L0,Se[j]),(Le,S0[j])]:
   for k,x in col.items():q=add(q,matrix[k],x)
  E[j]=q
 expected=[]
 for n,rec in enumerate(d['selected548_source_records']):
  src=rec['original_old640_source'];a1=sp(rec['actual_a1_27_times9']);assert all(k%9<7 for k in a1);v=[]
  for ch,row in enumerate(G):
   x=F(0)
   for k,c in row.items():
    target,p=divmod(k,8);u,w=pairs[p];root=B[src].get(target)
    if root is not None:x=ck(x+ck(c*jets[root].get((u,w),F(0))))
    if u<7:
     for j in d['needed_pivot_indices']:
      a=a1.get(9*j+u,F(0))
      if a:x=ck(x-ck(c*ck(a*(R[j].get(target,F(0)) if w==7 else E[j].get(target,F(0))))))
   D,U=map(int,d['all548_four_images_global_D_U'][ch]);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U;v.append(x)
  r={'case619':rec['case619'],'original_old640_source':src,'actual_lower4_values':list(map(str,v))};assert r==saved['computed_source_prefix'][n] and not any(v);expected.append(r)
 assert len(expected)==saved['computed_prefix_count']==548 and saved['all548_zero'] and saved['nonzero_witness'] is None
 ep=out/'independent-all2192-functional-identities.json';ep.write_text(json.dumps({'all548_literal_case_and_source_rows':[[r['case619'],r['original_old640_source']] for r in expected],'all548_times4_functional_values_match_and_zero':True,'ordinary_polynomial_convolution_different_from_producer_vector_order':True,'arbitrary_second_extension_applied_only_after_formal_weight_independence':True,'true_leaf_second_jets_not_claimed_zero':True,'all_firstcross_and_individual_LeS0_carries_retained':True,'no_new_physical_gauge_or_unknown_second_coefficients':True},separators=(',',':'))+'\n')
 rep={'id':stem,'all_passed':True,'all548_times4_exact_functional_identities':True,'lower548_allzero':True,'no_true_leaf_second_values_evaluated':True,'actual_maximum_bits':maximum,'input_sha256':dict(hs),'evidence_sha256':{ep.name:sha(ep)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{stem}.md'),'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(rep,separators=(',',':'))+'\n')
 # Separate symbolic truncated-free-polynomial proof for source595/case508.
 ostem='AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63_REVIEW';op='reports/AMBIENT_B18_ELL2_T2_OLD_CASE508_ORDER_20261105_v63';orderprimary=read(op+'/report.json');assert hs[op+'/report.json'].startswith('98463d9f') and orderprimary['all_passed'];order_ev=read(op+'/case508-original595-source-order-and-plane-cut-proof.json')
 for r in orderprimary['input_sha256'].values():assert sha(ROOT/r['path'])==r['sha256']
 c0='reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json';basis=field(c0,'actual_C0_basis8_scalar_nodes');div=field(c0,'all35_regular_coordinate_divisions');assert basis[7]==div[0]['restriction_r_s_node']==505938
 rel='reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json';gb=field(rel,'CB_g_restriction_node');oldf=field(rel,'CB_f63_node');instances=field(rel,'all16_pointed_unit57_inverse_instances');tail_instances=[r for r in instances if r['retained_Y_indices']==[7,8]];assert len(tail_instances)==1 and F(tail_instances[0]['unit_Jacobian']) and len(tail_instances[0]['all57_pointed_equations'])==57
 tail=next(r for r in field('reports/AMBIENT_B18_NORMALIZER_CHAIN_LIFT_640_NEGATIVE_BOUNDARY_20261104_v63/actual-UP-QB-640-typed-operator.json','accepted_actual_regular_function_identities') if r['kind']=='accepted_actual_CB_f63_tail');assert tail['gB']==gb==299650 and tail['function']==oldf==28539
 fp='reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json';mons=field(fp,'source2_monomials');families=field(fp,'both57_source2');cut=[i for i,m in enumerate(mons) if set(m)<=set([7,8])];assert len(cut)==6;assert families['CA']['f63']==div[0]['original_pure_CA_source_function'] and families['CB']['f63']==oldf
 for family in ['CA','CB']:
  jet=families[family]['pure_function2jets'][str(families[family]['f63'])];assert all(F(jet[i])==0 for i in cut)
 cc=field('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json','all_zero_bottom_and_Y0_source_combinations')[508];assert cc['original768_combination']==[[595,'1']] and not cc['actual_a1_27_times9']
 mug='reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json';allops=dict(field(mug,'old640_typed_closure_nodes'));cb=dict(field(mug,'old640_source_columns')[595]);gamma=list(map(sp,field('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json','Gamma4_on_source56_times8')));coord=set(basis[:7]);high={505938,299650};memo={}
 def truncated(i,plane):
  key=(i,plane)
  if key in memo:return memo[key]
  n=allops[i]
  if n[0]=='c':q={():F(n[1],n[2])} if n[1] else {}
  elif n[0]=='actual_C0_scalar_function':q={} if n[1] in high or plane and n[1] in coord else {(n[1],):F(1)}
  elif n[0]=='sum':
   q={}
   for j,a,b in n[1]:q=add(q,truncated(j,plane),F(a,b))
  else:
   assert n[0]=='product';q={():F(n[1],n[2])} if n[1] else {}
   for j in n[3]:q=product(q,truncated(j,plane),lambda w:len(w)<=2)
  memo[key]=q;return q
 receipts=[]
 for ch,row in enumerate(gamma):
  for k,c in row.items():
   target,p=divmod(k,8);root=cb.get(target)
   if root is None:continue
   plane=all(i in [7,8] for i in pairs[p]);assert not truncated(root,plane);receipts.append([ch,target,p,root,plane])
 ordinary=truncated(cb[29],False);assert cb[29]==281 and ordinary and all(any(i in coord for i in word) for word in ordinary);assert not truncated(281,True)
 oe=orderout/'independent-case508-truncated-free-polynomial-and-zero7-proof.json';oe.write_text(json.dumps({'case619':508,'source_old640':595,'a1_empty':True,'literal_gA_gB_distinct_ordinary_f63':[505938,299650,28539],'both_pointed_zero7_conductors_all6_source2_coefficients_zero':True,'all_supported_gamma_coefficients_zero_by_symbolic_order':receipts,'ordinary_root281_terms_contain_actual_z_coordinate':[[list(w),str(x)] for w,x in ordinary.items()],'plane_root281_truncated_free_polynomial_zero':True,'no_numerical_quadratic_image_or_newjet':True,'source_proof':'Opaqueotherleaforiginszeroarealreadycertified7d2918;assigngA/gBweight3fromliteralpointedzero7unit57source2. Truncatefreepolynomialsmodm3;coordinateplanequotientkillsz7. EveryGamma-supportedsecondcoefficientvanishes.'},separators=(',',':'))+'\n')
 orderreport={'id':ostem,'all_passed':True,'case508_all_four_lower_values_zero_by_source_order':True,'root281_zideal_mod_m3_checked_symbolically':True,'literal_conductor_order_bindings_verified':True,'no_new_second_or_physical_coefficients':True,'input_sha256':hs,'evidence_sha256':{oe.name:sha(oe)},'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{ostem}.md'),'combined_with_peered548_closes549_OLD_lower_only':True,'remaining70_native_mixed_and_all619_upper2_open':True,'seconds':time.monotonic()-t};(orderout/'report.json').write_text(json.dumps(orderreport,separators=(',',':'))+'\n');print(json.dumps({'548_all_passed':True,'all2192_values_zero':True,'508_order_passed':True,'548_report_sha256':sha(out/'report.json'),'508_report_sha256':sha(orderout/'report.json'),'seconds':time.monotonic()-t}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--order-output',required=True,type=Path);a=p.parse_args();main(a.output.resolve(),a.order_output.resolve())
