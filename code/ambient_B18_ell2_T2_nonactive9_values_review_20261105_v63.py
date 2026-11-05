#!/usr/bin/env python3
"""Independent complete9 cached component/sign/fixed-rank postcheck."""
import json,resource,signal,time
from fractions import Fraction as Q
from pathlib import Path
from ambient_B18_ell2_T2_nonactive171_source_review_20261105_v63 import sha,sp,count,add,image
ROOT=Path(__file__).resolve().parents[1];BASE='AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUES_20261105_v63';STEM=BASE+'_REVIEW'
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024));inputs={}
 def read(p):inputs[p]={'path':p,'sha256':sha(ROOT/p)};return json.loads((ROOT/p).read_text())
 rp=f'reports/{BASE}/report.json';r=read(rp);assert r['all_passed'] and r['no_new_physical_gauge_implicit_or_jet_coefficients']
 assert sha(ROOT/'code/ambient_B18_ell2_T2_nonactive9_values_20261105_v63.py')==r['code_sha256']
 for p,h in r['input_sha256'].items():assert sha(ROOT/p)==h
 pre='reports/AMBIENT_B18_ELL2_T2_NONACTIVE9_VALUE_INPUT_GATE_20261105_v63';source=read(pre+'_REVIEW/report.json');assert source['all_passed']
 d=read(pre+'/whole-old-native-pivot-nine-classes-positive-source-input.json')
 phase1=read(f'reports/{BASE}/phase1-cached34-pivot-chart-firsts.json');phase2=read(f'reports/{BASE}/phase2-actual18-psi-and-nine-quadratic-classes.json');rank=read(f'reports/{BASE}/actual-nine-quadratic-fixed-integer-rank-and-all-factorizations.json')
 for p,h in r['evidence_sha256'].items():assert sha(ROOT/f'reports/{BASE}'/p)==h
 L0,Lr,Le,S0,Sr,Se,S2=[list(map(sp,d[k])) for k in ['native_L0','native_Lr','native_Le','native_S0','native_Sr','native_Se','native_S2rell']];Gam=list(map(sp,d['Gamma_two']));QQ=list(map(sp,d['Q_two6_columns']));CV={i:Q(x) for i,x in d['literal_old104_second_values']}
 def dot(v,w):return sum((x*w.get(i,Q(0)) for i,x in v.items()),Q(0))
 pR=list(map(sp,d['pivot27_r_columns']));pE=[]
 for mode,value in d['pivot27_ell_recipes']:
  if mode=='literal':pE.append(sp(value))
  else:assert mode=='native';pE.append(add(image(L0,Se[value]),image(Le,S0[value])))
 f=d['carry_first_inputs'];C={j:sp(v) for j,v in f['selected99_original768']};R={j:sp(v) for j,v in f['original_r_first_columns']};E={j:sp(v) for j,v in f['original640_ell_first_columns']};inv=list(map(sp,f['fixed_B30_inverse_columns']));free=f['free30_V_indices'];carry={}
 for k in {k for c in C.values() for k in c if k>=640}:E[k]=image(L0,Se[k-640])
 for j,c in C.items():
  assert not image(S0,{k-640:x for k,x in c.items() if k>=640});pair=[]
  for raw in [R,E]:
   v=image(raw,c);a=image(inv,{i:v[row] for i,row in enumerate(free) if row in v});pair.append({i:x for i,x in a.items() if i<27})
  carry[j]=pair
 assert carry=={j:list(map(sp,pair)) for j,pair in phase1['all34_a1_r_ell']}
 def envelope(x,b):
  D,U=map(int,b);assert D%x.denominator==0 and abs(x.numerator)*(D//x.denominator)<=U
 psi=[]
 for n,(rec,wrec,expected,bounds) in enumerate(zip(d['old18_formal_results'],f['formal_eighteen_weights'],phase2['all18_exact_component_identities'],d['all18_component_positive_bounds'],strict=True)):
  assert (rec['original418_source'],rec['contraction'])==(expected['source418'],expected['contraction'])==(wrec['original418_source'],wrec['contraction'])
  aa=[[{},{}],[{},{}]]
  for flat,c in wrec['weights']:
   j,ch=divmod(flat,2)
   for direction in range(2):aa[ch][direction]=add(aa[ch][direction],carry[j][direction],Q(c))
  old=sum((Q(c)*CV[i] for i,c in rec['exact_independent_original_C0_second_sensitivity']),Q(0));envelope(old,bounds['old'])
  native=[];pivot=[]
  for ch,record in enumerate(rec['exact_original_native128_combinations_two_channels']):
   cn=sp(record);assert not image(S0,cn)
   vec=add(add(image(L0,image(S2,cn)),image(Lr,image(Se,cn))),image(Le,image(Sr,cn)))
   x=dot(Gam[ch],vec);native.append(x);envelope(x,bounds['native_channels'][ch])
   ar,ae=aa[ch];p=dot(Gam[ch],add(image(pR,ae),image(pE,ar)));pivot.append(p);envelope(p,bounds['pivot_channels'][ch])
  value=old+sum(native,Q(0))-sum(pivot,Q(0));envelope(value,bounds['psi']);assert not value
  assert Q(expected['old_cached'])==old and Q(expected['native_coherent'])==sum(native,Q(0)) and Q(expected['pivot_to_subtract'])==sum(pivot,Q(0)) and Q(expected['psi'])==value
  assert list(map(Q,expected['native_two_channels']))==native and list(map(Q,expected['pivot_two_channels']))==pivot
  assert expected['M1_identically_zero'] and expected['G2quadratic_identically_zero'];psi.append(value)
 assert psi==[Q(0)]*18
 assert list(map(sp,phase2['actual_nine_coordinates2']))==[{}]*9 and list(map(sp,phase2['actual_nine_vectors6']))==[{}]*9
 assert rank['actual_integer_matrix2_by9']==[[0,0]]*9 and rank['rank']==0 and rank['selected_original9_columns']==[] and rank['fixed_minor_rows']==[] and rank['fixed_integer_minor']=='0'
 assert list(map(sp,rank['all9_selected_basis_factorizations']))==[{}]*9 and list(map(sp,rank['cleared_integer_kernel_vectors']))==[{i:Q(1)} for i in range(9)]
 assert rank['remaining619_mixed_cycle_family_open'] and r['complete_G3_projected393_quadratic_rank']==0
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);e=out/'independent18-component-identities-nine-zero-classes-and-fixed-rank.json';e.write_text(json.dumps({'all34_first_pivot_carries_rebuilt':True,'all18_old_native_pivot_component_values_exact':True,'all3_native_products_and_individual_LeS0_retained':True,'all18_psi_zero':True,'all9_six_coordinate_vectors_zero':True,'fixed_rank0_and_identity9_kernel_verified':True,'complete393_G3_constant_projected_kernel_quadratic_image_zero':True,'remaining619_cycle_family_open':True,'actual_HF2_or_total_not_decided':True,'scope':'Active222zeroand162structuralzeroarefrozenpremises.The9completevectorsareindependentlyzero;thisclosesG3projected393only.NoactualM1responseornewphysical/gaugecoefficientwascomputed.'},separators=(',',':'))+'\n')
 report={'id':STEM,'all_passed':True,'all18_psi_and9_vectors_zero':True,'complete_G3_projected393_quadratic_rank':0,'no_HF2_equality_or_total_claim':True,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{e.name:sha(e)},'seconds':time.monotonic()-t};(out/'report.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({'all_passed':True,'report_sha256':sha(out/'report.json'),'rank':0}))
if __name__=='__main__':main()
