#!/usr/bin/env python3
"""First-only548 source/cache/positive gate, no actual relation image."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_OLD548_FIRST_ONLY_INPUT_20261105_v63'
PAIRS=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)]
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p,keys=None):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 formal='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63';fr=read(formal+'/report.json');assert fr['all_passed']
 fv=read(formal+'/actual-formal619-oldsecond-weights-and-symbol-support.json');nonzero=set(fv['nonzero_case619_indices'])
 C=read('reports/AMBIENT_B18_ELL2_T2_KNOWN_FIRST_PROFILE_20261105_v63/actual-known-bottom-first-rank-and-corrected-source-kernel.json',['all_zero_bottom_and_Y0_source_combinations'])['all_zero_bottom_and_Y0_source_combinations']
 cases=[]
 for i,c in enumerate(C):
  if i not in nonzero and all(j<640 for j,x in c['original768_combination']):
   assert len(c['original768_combination'])==1 and Q(c['original768_combination'][0][1])==1
   assert all(k%9<7 for k,x in c['actual_a1_27_times9'])
   cases.append({'case619':i,'original_old640_source':c['original768_combination'][0][0],'actual_a1_27_times9':c['actual_a1_27_times9']})
 assert len(cases)==548 and sum(len(c['actual_a1_27_times9']) for c in cases)==103
 gg=read('reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json');Gamma=list(map(sp,gg['Gamma4_on_source56_times8']))
 mug=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','native_Se_columns','native_S1ell_columns','native_S0_columns','L0_Res4_Ntranspose_columns']) if False else read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','native_S1ell_columns','native_S0_columns','L0_Res4_Ntranspose_columns'])
 B=list(map(dict,mug['old640_source_columns']));ops=dict(mug['old640_typed_closure_nodes']);usedrows=sorted({k//8 for row in Gamma for k in row});assert max(usedrows)<48
 prunedB={rec['original_old640_source']:{i:B[rec['original_old640_source']][i] for i in usedrows if i in B[rec['original_old640_source']]} for rec in cases}
 roots={k for v in prunedB.values() for k in v.values()};closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=ops[k]
  for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(j)
 for root in roots:visit(root)
 leaves=sorted({ops[k][1] for k in closure if ops[k][0]=='actual_C0_scalar_function'})
 fp='reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json';cache=read(fp)['first640'];firstall=dict(cache['all79_atom_values_full10']);assert set(leaves)<=set(firstall)
 first={i:sp(firstall[i]) for i in leaves};assert all(not v.get(0,Q(0)) for v in first.values())
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['original27_X_source_columns'])
 pivids=mini['original27_X_source_columns'];neededpiv=sorted({slot//9 for c in cases for slot,x in c['actual_a1_27_times9']})
 R={}
 for off,path,key in [(0,'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),(640,'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  raw=read(path,[key])[key]
  for p in neededpiv:
   j=pivids[p]
   if off<=j<off+len(raw):R[p]={row:sp(v)[8] for row,v in raw[j-off] if row in usedrows and sp(v).get(8)}
 Eold=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns']
 Eold={p:{row:x for row,x in sp(Eold[pivids[p]]).items() if row in usedrows} for p in neededpiv if pivids[p]<640}
 nativeids={p:pivids[p]-640 for p in neededpiv if pivids[p]>=640}
 nf=read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',['actual_L1ell_columns'])
 L0=[{i:x for i,x in sp(v).items() if i in usedrows} for v in mug['L0_Res4_Ntranspose_columns']]
 Le=[{i:x for i,x in sp(v).items() if i in usedrows} for v in nf['actual_L1ell_columns']]
 Se={p:sp(mug['native_S1ell_columns'][j]) for p,j in nativeids.items()};S0={p:sp(mug['native_S0_columns'][j]) for p,j in nativeids.items()}
 data={'selected548_source_records':cases,'ordinary_pair_order8':PAIRS,'Gamma4_source56_times8':gg['Gamma4_on_source56_times8'],
 'oldB_source_rows_pruned':[[j,[[i,k] for i,k in sorted(v.items())]] for j,v in sorted(prunedB.items())],'typed_closure_nodes':[[k,ops[k]] for k in sorted(closure)],
 'actual_C0_leaf_first10':[[i,[[k,str(x)] for k,x in sorted(v.items())]] for i,v in sorted(first.items())],
 'point_pivot27_original_source_ids':pivids,'needed_pivot_indices':neededpiv,'actual_pivot_r_first_pruned':[[p,[[i,str(x)] for i,x in sorted(v.items())]] for p,v in sorted(R.items())],
 'actual_old_pivot_ell_first_pruned':[[p,[[i,str(x)] for i,x in sorted(v.items())]] for p,v in sorted(Eold.items())],
 'native_pivot_local128_ids':[[p,j] for p,j in sorted(nativeids.items())],
 'native_L0_rows_pruned':[[[i,str(x)] for i,x in sorted(v.items())] for v in L0],'native_L1ell_rows_pruned':[[[i,str(x)] for i,x in sorted(v.items())] for v in Le],
 'native_pivot_S0_columns':[[p,[[i,str(x)] for i,x in sorted(v.items())]] for p,v in sorted(S0.items())],
 'native_pivot_S1ell_columns':[[p,[[i,str(x)] for i,x in sorted(v.items())]] for p,v in sorted(Se.items())],
 'zero_extension_proof':'EachselectedcasehasnonnativeoldsinglecolumnandALL4weightedhonestsecondcoefficientszero. Theonlypossiblehonestseconddependenceislinear. HenceANYsecondextensionoftheseEXACTfirstscomputesthetrueweightedimage;settingextension0doesNOTclaimtrueleafsecondszero. Firstcrossproductsremain.',
 'native_pivot_first_rule':'Xp_ell=L0*S1ell+L1ell*S0INDIVIDUALLY;neverdroptheLeS0carry.',
 'no_actual548_lower_values':True}
 base=count(data);assert base<60000
 maximum=0
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());return b
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator)))
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else (1,0)
 zero=(1,0);one=(1,1)
 def jm(a,b):
  v=[mul(a[0],b[0])]
  for k in range(1,10):v.append(plus(mul(a[0],b[k]),mul(b[0],a[k])))
  for j,(u,w) in enumerate(PAIRS):
   c=plus(mul(a[0],b[10+j]),mul(b[0],a[10+j]));c=plus(c,mul(a[u+1],b[w+1]))
   if u!=w:c=plus(c,mul(a[w+1],b[u+1]))
   v.append(c)
  return v
 jets={}
 for k in sorted(closure):
  n=ops[k];op=n[0]
  if op=='actual_C0_scalar_function':v=[ex(first[n[1]].get(i,Q(0))) for i in range(10)]+[zero]*8
  elif op=='c':v=[ex(Q(n[1],n[2]))]+[zero]*17
  elif op=='sum':
   v=[zero]*18
   for j,a,b in n[1]:v=[plus(x,mul(ex(Q(a,b)),y)) for x,y in zip(v,jets[j],strict=True)]
  else:
   assert op=='product';v=[ex(Q(n[1],n[2]))]+[zero]*17
   for j in n[3]:v=jm(v,jets[j])
  jets[k]=v
 ER={p:{i:ex(x) for i,x in v.items()} for p,v in Eold.items()}
 for p in nativeids:
  v={}
  for A,b in [(L0,Se[p]),(Le,S0[p])]:
   for col,c in b.items():
    for i,x in A[col].items():v[i]=plus(v.get(i,zero),mul(ex(x),ex(c)))
  ER[p]=v
 globalbound=[zero]*4
 for rec in cases:
  cb=B[rec['original_old640_source']];a1=sp(rec['actual_a1_27_times9'])
  for target,row in enumerate(Gamma):
   q=zero
   for flat,c in row.items():
    source,pair=divmod(flat,8);root=cb.get(source)
    if root is not None:q=plus(q,mul(ex(c),jets[root][10+pair]))
    u,v=PAIRS[pair]
    if u<7:
     for p in neededpiv:
      a=a1.get(9*p+u,Q(0))
      if a:
       firstbound=ex(R[p].get(source,Q(0))) if v==7 else ER[p].get(source,zero)
       q=plus(q,mul(ex(c),mul(firstbound,ex(a))))
   D,U=globalbound[target];dd,uu=q;nd=math.lcm(D,dd);globalbound[target]=mark((nd,max(U*(nd//D),uu*(nd//dd))))
 assert maximum<=384
 sourceowned=count([data,jets,ER,C,mug,nf,cache])+8000;future=2*(base+20)+18*len(closure)+8*548+6*48*2+8000
 raw=2_000_000+future*(2*math.ceil(maximum*math.log10(2))+40)
 ready=sourceowned<150000 and future<150000 and raw<32*1024*1024
 data.update({'positive_committed_bit_bound':maximum,'positive_temporary_bit_bound':3*maximum+1,'all548_four_images_global_D_U':[[str(D),str(U)] for D,U in globalbound],
 'future_first_only_numeric_bound':future,'future_first_only_RAW_bound':raw,'source_positive_owned_bound':sourceowned,'source_prior_ready':ready,
 'future_scope':'Witness-firstscanof548lower4imagesONLY.Noranksearch.NativeX2andallactualsecondresponsesexcluded.NovalueworkerGOimpliedbythissourcegate.'})
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'pruned548-first-only-cache-source-and-positive-bound.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':ready,'source_only_no_actual548_values':True,'selected_count':548,'old_typed_closure_nodes':len(closure),'cached_leaf_firsts':len(leaves),'positive_prior_bits':maximum,'future_numeric_bound':future,'future_RAW_bound':raw,'source_owned_bound':sourceowned,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':ready,'closure':len(closure),'leaves':len(leaves),'prior':maximum,'future':future,'RAW':raw,'sourceowned':sourceowned,'report_sha256':sha(out/'report.json')}))
if __name__=='__main__':main()
