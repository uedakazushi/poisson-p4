#!/usr/bin/env python3
"""Independent27+99 formal unit Schur/source proof; no new coefficients."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib,json,resource,signal,time
from ambient_B18_ell2_raw33_quotient_review_20261105_v63 import sp,enc,add,scale,mv,mp,eye
from ambient_B18_ell2_T2_mu_review_20261105_v63 import fields
ROOT=Path(__file__).resolve().parents[1]
A='AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63';STEM=A+'_REVIEW'
P={
 'primary':f'reports/{A}/report.json','constant':f'reports/{A}/actual27-plus99-unit-complements-and-source-combinations.json',
 'op':f'reports/{A}/actual-honest-T2-minimal6-by1060-factored-source-operator.json',
 'mu_peer':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63_REVIEW/report.json',
 'rank':'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json',
 'gate':'reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json',
 'frame':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-active-unit-bases-and26-frame-transport.json',
 'source':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-honest-T2-132-by1186-factored-source-operator.json',
 'source_peer':'reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63_REVIEW/report.json',
 'origin':'reports/AMBIENT_B18_COMPARISON_X_ORIGIN_20261104_v63/actual-X0-and-unit-cancellation.json',
 'origin_peer':'reports/AMBIENT_B18_ELL2_T2_ORIGIN_BRIDGE_20261105_v63_REVIEW/report.json',
}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(o,n,x):(o/n).write_text(json.dumps(x,separators=(',',':'))+'\n')
def det(M):
 a=[list(r) for r in M];out=Q(1);n=len(a)
 for j in range(n):
  k=next(i for i in range(j,n) if a[i][j])
  if k!=j:a[k],a[j]=a[j],a[k];out=-out
  d=a[j][j];out*=d
  for i in range(j+1,n):
   c=a[i][j]/d
   for k in range(j+1,n):a[i][k]-=c*a[j][k]
   a[i][j]=Q(0)
 return out
def main():
 t=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000))
 o=ROOT/'reports'/STEM;o.mkdir(parents=True,exist_ok=False)
 inputs={k:{'path':p,'sha256':sha(ROOT/p)} for k,p in P.items()};dump(o,'phase.json',{'phase':'ACTUAL27_PLUS99_UNIT_SOURCE_SCHUR','inputs':inputs})
 D={k:json.loads((ROOT/p).read_text()) for k,p in P.items() if k!='origin'}
 assert D['primary']['all_passed'] and all(D[k]['all_passed'] for k in P if k.endswith('_peer'))
 for p,h in D['primary']['input_sha256'].items():assert sha(ROOT/p)==h
 for p,h in D['primary']['evidence_sha256'].items():assert sha(ROOT/'reports'/A/p)==h
 O=fields(ROOT/P['origin'],['selected_independent_columns','full_X0_56_by768'])
 C=D['constant'];piv=O['selected_independent_columns'];assert len(piv)==27 and C['original27_X_source_columns']==piv
 freeV=[i for i in range(56) if i not in D['frame']['same_selected26_V_columns']];assert len(freeV)==30 and C['free30_V_indices']==freeV
 Bp=[{i:Q(O['full_X0_56_by768'][row][j]) for i,row in enumerate(freeV) if Q(O['full_X0_56_by768'][row][j])} for j in piv];del O
 basis={}
 for j,r in enumerate(C['fixed27_receipts']):
  assert r['column']==j;v=dict(Bp[j])
  for p,c in r['subtract']:
   assert v.get(p)==Q(c);v=add(v,scale(basis[p],-Q(c)))
  p=r['pivot'];c=Q(r['divide_by']);assert c and min(v)==p and v[p]==c and p not in basis
  basis[p]=scale(v,1/c);assert enc(basis[p])==r['normalized']
 free3=[i for i in range(30) if i not in basis];assert len(free3)==3 and C['remaining_free3_bottom_coordinates']==free3
 B30=Bp+[{i:Q(1)} for i in free3];Bi=list(map(sp,C['fixed_bottom30_inverse_columns']))
 assert list(map(enc,B30))==C['fixed_bottom30_basis_Xp27_then_Efree3'] and mp(B30,Bi)==mp(Bi,B30)==eye(30)
 determinant=det([[B30[j].get(i,Q(0)) for j in range(30)] for i in range(30)])
 assert str(determinant)==C['fixed_bottom30_determinant'] and determinant
 assert C['remaining_free3_originalV_rows']==[freeV[i] for i in free3]==[23,29,30]
 R=D['rank'];selected=R['actual_selected_Y0_source_columns'];Y=list(map(sp,R['actual_selected_Y0_columns102']));factors=list(map(sp,R['all741_factorizations_in_original_selected_Y0_basis']))
 freeY=[i for i in range(102) if i not in R['fixed_original_minor_rows']];assert freeY==C['remaining_free3_Y_rows']==[1,36,38]
 assert selected==C['selected99_kernel741_columns'] and len(selected)==99
 rest=[i for i in range(741) if i not in selected];assert rest==C['remaining642_kernel741_columns'] and len(rest)==642
 kernels=D['gate']['all741_zero_origin_source_combinations'];zero=[]
 for k in rest:
  coeff=add({k:Q(1)},{selected[j]:-a for j,a in factors[k].items()});assert not mv(factors,coeff)
  raw={}
  for j,a in coeff.items():raw=add(raw,scale(sp(kernels[j]['original768']),a))
  zero.append({'kernel741_index':k,'constant741_combination':enc(coeff),'original768_combination':enc(raw)})
 assert zero==C['all642_zero_constant_combinations'] and sum(len(r['original768_combination']) for r in zero)==935
 assert C['final_minimal6_target_order']==['Y_source_rows',freeY,'free30_coordinates',free3]
 E=D['op'];N=E['operator_nodes'];S=D['source'];sourcepath=P['source']
 def one(tag):
  matches=[i for i,n in enumerate(N) if n[0]==tag];assert len(matches)==1;return matches[0]
 def check(i,tag,*args):assert N[i]==[tag,*args]
 raw=E['original132_operator_node'];check(raw,'external_complete_matrix_operator',sourcepath,S['complete_T2_presentation_node'],[132,1186])
 b1=[i for i,n in enumerate(N) if n[0]=='horizontal_unit_basis' and n[-1].get('column_blocks')==[27,102,3]];assert len(b1)==1;b1=b1[0]
 for n,j in zip(N[b1][1],piv):check(n,'matrix_column',raw,418+j,[132])
 assert len(N[b1][1])==27 and N[b1][2:4]==[['standard_top102_in_B132'],['standard_bottom_coordinates_in_B132',free3]]
 assert N[b1][-1]['actual_unknown_top_Yp_retained'] and N[b1][-1]['unit_constant_determinant']==str(determinant)
 b1i=next(i for i,n in enumerate(N) if n[0]=='formal_regular_unit_inverse' and n[1]==b1)
 pi1=E['first27_quotient_node'];check(pi1,'compose',['coordinate_projection',132,list(range(27,132))],b1i)
 p1=one('matrix_of_columns');stage1=N[p1][1];assert len(stage1)==1159 and N[p1][2]==[105,1159]
 for j,v in enumerate(stage1[:418]):
  assert N[v][:2]==['apply_vector',pi1];check(N[v][2],'matrix_column',raw,j,[132])
 for r,v in zip(kernels,stage1[418:]):
  assert N[v][:2]==['apply_vector',pi1];check(N[v][2],'constant_linear_combination_of_matrix_columns',raw,[[418+j,c] for j,c in r['original768']],[132])
 b2=next(i for i,n in enumerate(N) if n[0]=='horizontal_unit_basis' and n[-1].get('column_blocks')==[99,3,3])
 assert N[b2][1:4]==[[stage1[418+j] for j in selected],['standard_top102_free_coordinates',freeY],['standard_bottom3_coordinates']]
 assert N[b2][-1]['unknown_higher_coefficients_retained']
 b2i=next(i for i,n in enumerate(N) if n[0]=='formal_regular_unit_inverse' and n[1]==b2)
 pi2=E['next99_quotient_node'];check(pi2,'compose',['coordinate_projection',105,list(range(99,105))],b2i)
 final=E['complete_minimal_T2_node'];info=N[final][1];assert N[final][0]=='complete_minimal_honest_T2_presentation' and info['shape']==[6,1060]
 remaining=list(range(418))+[418+j for j in rest]
 assert len(info['matrix_columns'])==len(remaining)==1060
 for v,j in zip(info['matrix_columns'],remaining):check(v,'apply_vector',pi2,stage1[j])
 assert info['target_order']==C['final_minimal6_target_order']
 ep=o/'independent-actual27-plus99-unit-complete-minimal6-proof.json'
 dump(o,ep.name,{'inputs':inputs,'all27_fixed_bottom30_unit_columns_and_two_sided_inverse_verified':True,'bottom30_unit_determinant':str(determinant),'unknown_Yp_actual_columns_retained_not_set_zero':True,'all99_actual_Y_origin_units_and_free6_completion_source_bound':True,'free3_Y_rows':freeY,'free3_originalV_bottom_rows':[23,29,30],'all642_zero_constant_source_combinations_and935coefficients_verified':True,'all418_G3_and768_original_X_sources_retained':True,'complete_minimal_T2_source_shape':[6,1060],'all6_by1060_constant_entries_zero_by_unit_projection':True,'conditional_first_Hilbert_economy':'on642 zeroY0 combinations, bottom3 firsts known; ontheirconstantkernel fullX1=Xp0*a1, subtractactualXp*a1 so(MXeff)2=M0Xeff2; no cachecompleteness or topfirstvalues claimed','new_physical_higher_or_T2_length_values_computed':False})
 report={'id':STEM,'all_passed':True,'complete_actual_minimal_T2_typed6_by1060_accepted':True,'all27_and99_actual_unit_pairs_with_unknown_Yp_retained_accepted':True,'all642_source_combinations_and_conditional_first_Hilbert_method_accepted':True,'new_source_higher_coefficient_T2_length_or_H0_computed':False,'input_sha256':inputs,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/f'{STEM}.md'),'evidence_sha256':{ep.name:sha(ep)},'seconds':time.monotonic()-t}
 dump(o,'report.json',report);dump(o,'phase.json',{'phase':'FINAL','report_sha256':sha(o/'report.json')});print(json.dumps(report),flush=True)
if __name__=='__main__':main()
