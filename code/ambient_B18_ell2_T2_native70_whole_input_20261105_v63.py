#!/usr/bin/env python3
"""Whole native70 cache/source/positive gate, never a signed image scan."""
import json,math,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count
ROOT=Path(__file__).resolve().parents[1];STEM='AMBIENT_B18_ELL2_T2_NATIVE70_WHOLE_INPUT_GATE_20261105_v63'
PAIRS=[(7,7),(7,8),(0,7),(0,8),(4,7),(4,8),(5,7),(5,8)];ZERO=(1,0)
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));pins={};maximum=0
 def read(p,keys=None):pins[p]={'path':p,'sha256':sha(ROOT/p)};return fields(ROOT/p,keys) if keys else json.loads((ROOT/p).read_text())
 nb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_POSITIVE_SOURCE_20261105_v63';nr=read(nb+'/report.json');assert nr['all_passed'] and nr['positive_envelopes_only_no_native70_values']
 nd=read(nb+'/positive-C72-program-and-native70-partial-coefficient-envelopes.json');assert sha(ROOT/nb/'positive-C72-program-and-native70-partial-coefficient-envelopes.json')==nr['evidence_sha256']['positive-C72-program-and-native70-partial-coefficient-envelopes.json']
 frbase='reports/AMBIENT_B18_ELL2_T2_NATIVE70_C72_FRAME_SOURCE_20261105_v63';fr=read(frbase+'/report.json');assert fr['all_passed'];frame=read(frbase+'/actual-C72-cache-frame-and-nine-direction-Nfirst-source-bundle.json')
 kb='reports/AMBIENT_B18_ELL2_T2_NATIVE70_A18_SOURCE_20261105_v63';kr=read(kb+'/report.json');assert kr['all_passed'] and kr['all70_fullC4_point_kernel_verified'];kernel=read(kb+'/all70-C4-kernels-and-full-raw33-A18-ring-map-source.json');recs=kernel['all70_post27_native_point_kernel_records'];assert [v['case619'] for v in recs]==list(range(549,619))
 gp='reports/AMBIENT_B18_ELL2_T2_MIXED619_LOWER4_X2_SOURCE_20261105_v63/actual-four-Gamma-X2-source-rows-and-lazy619-recipe.json';gg=read(gp,['Gamma4_on_source56_times8','ordinary_second_pair_order8']);assert list(map(tuple,gg['ordinary_second_pair_order8']))==PAIRS;Gamma=list(map(sp,gg['Gamma4_on_source56_times8']));rows=sorted({k//8 for row in Gamma for k in row});assert max(rows)<48
 mug=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',['old640_source_columns','old640_typed_closure_nodes','native_S1ell_columns','native_S0_columns','L0_Res4_Ntranspose_columns'])
 B=list(map(dict,mug['old640_source_columns']));ops=dict(mug['old640_typed_closure_nodes']);oldids=sorted({i for r in recs for i,c in r['original768_constant_combination'] if i<640})
 prows={j:{i:B[j][i] for i in rows if i in B[j]} for j in oldids};closure=set()
 def visit(k):
  if k in closure:return
  closure.add(k);n=ops[k]
  for j in ([j for j,a,b in n[1]] if n[0]=='sum' else n[3] if n[0]=='product' else []):visit(j)
 for v in prows.values():
  for k in v.values():visit(k)
 leaves=sorted({ops[k][1] for k in closure if ops[k][0]=='actual_C0_scalar_function'})
 allfirst=dict(read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/phase1-independent-full9-first640-source-proof.json',['first640'])['first640']['all79_atom_values_full10']);assert set(leaves)<=set(allfirst)
 first={i:sp(allfirst[i]) for i in leaves};assert all(not v.get(0,Q(0)) for v in first.values())
 fs='reports/AMBIENT_B18_ELL2_T2_MIXED619_OLD_LINEAR_SENSITIVITY_20261105_v63';fsr=read(fs+'/report.json');assert fsr['all_passed'];sens=read(fs+'/actual-formal619-oldsecond-weights-and-symbol-support.json',['nonzero_original_C0_symbols'])['nonzero_original_C0_symbols']
 honest13=[7199,7454,7505,8045,27897,27948,28203,28254,28437,28692,28794,299650,505938];assert sorted(sens)==sorted(honest13)
 assert read(fs+'_REVIEW/report.json')['all_passed']
 raw=read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63/actual-source105-ell2-raw33-values.json',['all_requested164_root_values','all105_equations_exactly_zero','all6_unit_inverse_products_exact']);assert raw['all105_equations_exactly_zero'] and raw['all6_unit_inverse_products_exact']
 assert read('reports/AMBIENT_B18_ELL2_C2_FINITE_SOURCE_20261105_v63_REVIEW/report.json')['all_passed']
 vals=dict(raw['all_requested164_root_values']);copy=frame['raw33_to_A18_coeff_copy_source_indices18'];primitive7=[27897,27948,28203,28254,28437,28692,28794];seconds={i:[Q(0)]*8 for i in honest13};cacheoverlap=[]
 for i in primitive7:
  v=sp(vals[i]);a=[v.get(j,Q(0)) for j in copy];assert len(a)==18
  assert all(a[j]==first.get(i,sp(allfirst[i])).get(j,Q(0)) for j in range(10))
  seconds[i]=a[10:];cacheoverlap.append({'root':i,'actual18_from_raw33':list(map(str,a))})
 assert all(not first[i] for i in [299650,505938] if i in first)
 mini=read('reports/AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63/actual27-plus99-unit-complements-and-source-combinations.json',['original27_X_source_columns']);piv=mini['original27_X_source_columns'];assert len(piv)==27
 needed=sorted({slot//9 for r in recs for slot,c in r['actual_a1_27_times9']});Xp={p:[{} for _ in range(8)] for p in needed}
 for off,path,key in [(0,'reports/AMBIENT_B18_COMPARISON_X_FIRST640_LINEAR_20261104_v63/actual-first640-X-mod-m2.json','source_columns640_rows56_sparse_dual9'),(640,'reports/AMBIENT_B18_NATIVE128_FIRST_20261104_v63/actual-native128-X-mod-m2.json','source_columns128_rows48_sparse_dual9')]:
  arr=read(path,[key])[key]
  for p in needed:
   if off<=piv[p]<off+len(arr):
    for row,v in arr[piv[p]-off]:
     if row in rows:
      for idx,x in sp(v).items():
       if 1<=idx<=8:Xp[p][idx-1][row]=x
 Eraw=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/phase1-actual-first640-ell-from-existing-full9-cache.json',['actual_first640_ell_columns'])['actual_first640_ell_columns'];Eold={p:{i:x for i,x in sp(Eraw[piv[p]]).items() if i in rows} for p in needed if piv[p]<640}
 L0=[{i:x for i,x in sp(v).items() if i in rows} for v in mug['L0_Res4_Ntranspose_columns']]
 Le=[{i:x for i,x in sp(v).items() if i in rows} for v in read('reports/AMBIENT_B18_ELL2_T2_FIRST_NATIVE_CACHE_20261105_v63/phase3-actual-native-S2rell-L1ell-and619-coherent-parts.json',['actual_L1ell_columns'])['actual_L1ell_columns']]
 nativep={p:piv[p]-640 for p in needed if piv[p]>=640};S0={p:sp(mug['native_S0_columns'][j]) for p,j in nativep.items()};Se={p:sp(mug['native_S1ell_columns'][j]) for p,j in nativep.items()}
 def enc(v):return [[i,str(x)] for i,x in sorted(v.items())]
 data={'source70_records':recs,'ordinary_pair_order8':PAIRS,'Gamma4_source56_times8':gg['Gamma4_on_source56_times8'],'oldB_source_rows_pruned':[[j,sorted(v.items())] for j,v in sorted(prows.items())],'typed_closure_nodes':[[k,ops[k]] for k in sorted(closure)],'actual_C0_leaf_first10':[[i,enc(v)] for i,v in first.items()],
 'actual_13_honest_second8':[[i,list(map(str,seconds[i]))] for i in honest13],'actual7_primitive_raw33_A18_cache_overlap':cacheoverlap,
 'point_pivot27_original_source_ids':piv,'needed_pivot_indices':needed,'actual_pivot_old8_first_pruned':[[p,[enc(v) for v in Xp[p]]] for p in needed],'actual_old_pivot_ell_first_pruned':[[p,enc(v)] for p,v in Eold.items()],
 'native_pivot_local128_ids':sorted(nativep.items()),'native_L0_rows_pruned':list(map(enc,L0)),'native_L1ell_rows_pruned':list(map(enc,Le)),'native_pivot_S0_columns':[[p,enc(v)] for p,v in S0.items()],'native_pivot_S1ell_columns':[[p,enc(v)] for p,v in Se.items()],
 'A18_basis_monomials':frame['A18_basis_monomials'],'C72_multiplication_records':nd['C72_nonzero_multiplication_records'],'actual_kappa9_C72':frame['actual_kappa9_C72'],'actual_nu_A_B_A18':frame['actual_nu_A_B_A18'],'regular_G8_source_program':frame['regular_G8_source_program'],'actual_G8_point_and_first40':frame['actual_G8_point_and_first40'],
 'actual_needed_N_point_and_first40':[[i,j,frame['actual_N12_by12_point_and_first40'][i][j]] for i in range(8) for j in [0,4,5,7]],
 'old_second_extension_scope':'Only13linearhonestsecondlabelspropagatetothewholefourΓOLDfunctionalsbytheacceptedcompleteformal619fold. TheotherleafsecondsareanARBITRARYextension0,notassertedactual0. Allfirstcrossproductsremain. Sevenprimitive13labelsuseactualfinite105readbacks;fourliteralcoordinatesandtwoorder3conductorshaveactualsecond0.',
 'native_extension_scope':'SameactualN0/N1andnu;all70fullC4S0cN0. AnysecondscalarfunctionextensionofNhasidenticalWEIGHTEDnativeX2. ActualN2notclaimed0;individualXp_ell=L0Se+LeS0stillretained.',
 'zero6_second_source_binding':{'literal_coordinate_nodes':[7199,7454,7505,8045],'actual_zero7_conductor_gA_gB_order3_nodes':[505938,299650],'accepted_existing_root_ancestry':True},'no_actual70_images':True}
 base=count(data);assert base<60000
 def mark(b):
  nonlocal maximum
  maximum=max(maximum,b[0].bit_length(),b[1].bit_length());assert maximum<=384;return b
 def ex(x):x=Q(x);return mark((x.denominator,abs(x.numerator))) if x else ZERO
 def plus(a,b):
  if not a[1]:return b
  if not b[1]:return a
  D=math.lcm(a[0],b[0]);return mark((D,a[1]*(D//a[0])+b[1]*(D//b[0])))
 def mul(a,b):return mark((a[0]*b[0],a[1]*b[1])) if a[1] and b[1] else ZERO
 def jm(a,b):
  v=[mul(a[0],b[0])]+[plus(mul(a[0],b[k]),mul(b[0],a[k])) for k in range(1,10)]
  for j,(u,w) in enumerate(PAIRS):
   z=plus(mul(a[0],b[10+j]),mul(b[0],a[10+j]));z=plus(z,mul(a[u+1],b[w+1]))
   if u!=w:z=plus(z,mul(a[w+1],b[u+1]))
   v.append(z)
  return v
 jets={}
 for k in sorted(closure):
  n=ops[k];op=n[0]
  if op=='actual_C0_scalar_function':v=[ex(first[n[1]].get(i,Q(0))) for i in range(10)]+[ex(x) for x in seconds.get(n[1],[Q(0)]*8)]
  elif op=='c':v=[ex(Q(n[1],n[2]))]+[ZERO]*17
  elif op=='sum':
   v=[ZERO]*18
   for j,a,b in n[1]:v=[plus(x,mul(ex(Q(a,b)),y)) for x,y in zip(v,jets[j],strict=True)]
  else:
   assert op=='product';v=[ex(Q(n[1],n[2]))]+[ZERO]*17
   for child in n[3]:v=jm(v,jets[child])
  jets[k]=v
 E={p:{i:ex(x) for i,x in v.items()} for p,v in Eold.items()}
 for p in nativep:
  v={}
  # SUM, never a union/max: individual native pivot carry is not eligible
  # for the full70 S0-kernel simplification.
  for A,b in [(L0,Se[p]),(Le,S0[p])]:
   for col,c in b.items():
    for row,x in A[col].items():v[row]=plus(v.get(row,ZERO),mul(ex(c),ex(x)))
  E[p]=v
 nativebounds={v['case619']:[tuple(map(int,p)) for p in v['four_native_partial_positive_bounds']] for v in nd['all70_native_partial_positive_bounds']}
 casebounds=[];global4=[ZERO]*4
 for rec in recs:
  old={i:Q(c) for i,c in rec['original768_constant_combination'] if i<640};a1=sp(rec['actual_a1_27_times9']);cb=[]
  for channel,row in enumerate(Gamma):
   oldbound=ZERO;carry=ZERO
   for slot,c in row.items():
    source,pair=divmod(slot,8);u,v=PAIRS[pair]
    for j,c0 in old.items():
     root=prows[j].get(source)
     if root is not None:oldbound=plus(oldbound,mul(ex(c),mul(ex(c0),jets[root][10+pair])))
    for p in needed:
     au=a1.get(9*p+u,Q(0));av=a1.get(9*p+v,Q(0))
     if au:
      xp=ex(Xp[p][v].get(source,Q(0))) if v<8 else E[p].get(source,ZERO);carry=plus(carry,mul(ex(c),mul(ex(au),xp)))
     if u!=v and av:
      xp=ex(Xp[p][u].get(source,Q(0))) if u<8 else E[p].get(source,ZERO);carry=plus(carry,mul(ex(c),mul(ex(av),xp)))
   native=mark(nativebounds[rec['case619']][channel]);total=plus(plus(oldbound,native),carry);cb.append({'old':[str(x) for x in oldbound],'native':[str(x) for x in native],'individual_pivot_carry':[str(x) for x in carry],'combined':[str(x) for x in total]})
   D,U=global4[channel];d,z=total;dd=math.lcm(D,d);global4[channel]=mark((dd,max(U*(dd//D),z*(dd//d))))
  casebounds.append({'case619':rec['case619'],'four_component_positive_envelopes':cb})
 maximum=max(maximum,nr['source_committed_bits']);temp=3*maximum+1;assert maximum<=384 and temp<=1153
 # Double decoded input and explicit scalar working blocks. No dense70 S/N
 # arrays, full9 jet, whole105 solve or rank matrix is needed.
 scratch=18*len(closure)+2*(28*72+8*72+8*72+4*72)+27*9*len(rows)+70*20+10000
 future=2*(base+1000)+scratch
 rawbound=2_000_000+future*(2*math.ceil(maximum*math.log10(2))+45)
 sourceowned=count([data,jets,E,frame,kernel,mug,raw,allfirst,casebounds,nd])+10000
 assert sourceowned<150000 and future<150000 and rawbound<32*1024*1024
 data.update({'all70_four_component_positive_envelopes':casebounds,'all70_four_images_global_D_U':[[str(D),str(U)] for D,U in global4],'whole_source_committed_bound':maximum,'whole_source_temporary_bound':temp,'future_whole_numeric_bound':future,'future_whole_RAW_bound':rawbound,'future_scratch_explicit':scratch,'source_owned_bound':sourceowned,'future_scope':'Witness-firstscanonlyof70actualLOWERfourclasses;completeold+native−individualXp1a1. Oneboundnonzerovectorgivesactualquadraticrelation. AllzeroONLYlower619familynegative,notupper2/HF2/length.'})
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);ep=out/'complete70-cached-input-source-and-positive-whole-image-bounds.json';ep.write_text(json.dumps(data,separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'whole_native70_source_only':True,'actual_native70_values_not_evaluated':True,'old_typed_nodes':len(closure),'old_cached_first_leaves':len(leaves),'all13_honest_second_callbacks_bound_or_eligible':True,'seven_primitive_second_raw33_overlap_verified':True,'source_committed_bits':maximum,'source_temporary_bits':temp,'source_owned_bound':sourceowned,'future_numeric_bound':future,'future_RAW_bound':rawbound,'input_sha256':pins,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'prior':maximum,'future':future,'RAW':rawbound,'owned':sourceowned,'seconds':report['seconds']}))
if __name__=='__main__':main()
