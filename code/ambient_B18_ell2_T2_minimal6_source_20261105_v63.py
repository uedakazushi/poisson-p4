#!/usr/bin/env python3
"""Actual allorder6x1060 source presentation, by frozen27+99 unit pairs.

No new comparison, source, or implicit coefficient is evaluated.
"""
import argparse,json,resource,signal,time
from fractions import Fraction as F
from pathlib import Path
import ambient_B18_ell2_T2_mu_20261105_v63 as H
from ambient_B18_ell2_raw33_quotient_20261105_v63 import add,scale,mv,srec
ROOT=H.ROOT;STEM='AMBIENT_B18_ELL2_T2_MINIMAL6_SOURCE_OPERATOR_20261105_v63'
def parse(v):return H.parse(v)
def fixed_basis(cols,n):
 piv={};receipts=[]
 for j,v0 in enumerate(cols):
  v=dict(v0);steps=[]
  while v and min(v) in piv:
   p=min(v);c=v[p];steps.append([p,str(c)]);v=H.guard(add(v,piv[p],-c))
  assert v;p=min(v);c=v[p];piv[p]=H.guard(scale(v,1/c));receipts.append({'column':j,'pivot':p,'divide_by':str(c),'subtract':steps,'normalized':srec(piv[p])})
 free=[i for i in range(n) if i not in piv]
 B=cols+[{i:F(1)} for i in free];assert len(B)==n
 # Fixed unit basis inverse, not an independent rank search.
 aug=[[B[j].get(i,F(0)) for j in range(n)]+[F(i==j) for j in range(n)] for i in range(n)];det=F(1);ops=[]
 for j in range(n):
  k=next(i for i in range(j,n) if aug[i][j]);swap=k
  if k!=j:aug[k],aug[j]=aug[j],aug[k];det=-det
  c=aug[j][j];det*=c;aug[j]=[a/c for a in aug[j]];clear=[]
  for i in range(n):
   if i!=j and aug[i][j]:
    a=aug[i][j];clear.append([i,str(a)]);aug[i]=[x-a*y for x,y in zip(aug[i],aug[j])]
  ops.append({'column':j,'swap_with':swap,'divide_by':str(c),'clear':clear})
  for r in aug:H.guard({i:c for i,c in enumerate(r) if c})
 inv=[{i:aug[i][n+j] for i in range(n) if aug[i][n+j]} for j in range(n)]
 assert all(mv(B,v)=={j:F(1)} for j,v in enumerate(inv)) and all(mv(inv,v)=={j:F(1)} for j,v in enumerate(B))
 return free,B,inv,det,receipts,ops
def main(out,mupeer):
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2_000_000_000,2_000_000_000));resource.setrlimit(resource.RLIMIT_FSIZE,(32*1024*1024,32*1024*1024))
 out.mkdir(parents=True,exist_ok=False);H.dump(out/'phase.json',{'phase':'source_unit_operator_started_no_newcoefficients'})
 hs={}
 def read(p):hs[p]=H.sha(ROOT/p);return json.loads((ROOT/p).read_text())
 pp=read(mupeer);assert pp['all_passed']
 mu=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/report.json');assert mu['actual_Y0_rank']==99 and mu['actual_mu_B_T2']==6
 assert any((v.get('sha256') if isinstance(v,dict) else v)==H.sha(ROOT/'reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/report.json') for v in pp['input_sha256'].values()),'peer must bind primary'
 addendum=read('reports/AMBIENT_B18_ELL2_T2_MU_RESOURCE_ADDENDUM_20261105_v63/report.json');assert addendum['corrected_conservative_peak_QQ_source_records']<150000
 gate=read('reports/AMBIENT_B18_ELL2_T2_MU_INPUT_GATE_20261105_v63/actual-T2-mu-input-source-and-cache-gate.json')
 rank=read('reports/AMBIENT_B18_ELL2_T2_MU_20261105_v63/actual741-Y0-fixed-rank-and-T2-minimal-generators.json')
 t=read(H.PATHS['transport']);sourcepath='reports/AMBIENT_B18_ELL2_HONEST_T2_ACTIVE_UNIT_TRANSPORT_20261105_v63/actual-honest-T2-132-by1186-factored-source-operator.json';source=read(sourcepath)
 hs[H.PATHS['origin']]=H.sha(ROOT/H.PATHS['origin']);piv=H.selected(H.PATHS['origin'],'selected_independent_columns');X0=H.selected(H.PATHS['origin'],'full_X0_56_by768')
 C=[i for i in range(56) if i not in t['same_selected26_V_columns']];assert len(C)==30
 Bp=[{i:F(X0[row][j]) for i,row in enumerate(C) if F(X0[row][j])} for j in piv];del X0
 free3,B30,B30inv,det30,rec27,invops=fixed_basis(Bp,30);assert len(free3)==3
 Y=list(map(parse,rank['actual_selected_Y0_columns102']));fixedrows=rank['fixed_original_minor_rows'];freeY=[i for i in range(102) if i not in fixedrows];assert len(freeY)==3
 selected=rank['actual_selected_Y0_source_columns'];rest=[j for j in range(741) if j not in selected];assert len(rest)==642
 factors=list(map(parse,rank['all741_factorizations_in_original_selected_Y0_basis']));kernels=gate['all741_zero_origin_source_combinations'];constzero=[]
 for k in rest:
  c=add({k:F(1)},{selected[j]:a for j,a in factors[k].items()},-1)
  assert not mv(factors,c)
  raw={}
  for j,a in c.items():raw=add(raw,parse(kernels[j]['original768']),a)
  constzero.append({'kernel741_index':k,'constant741_combination':srec(c),'original768_combination':srec(raw)})
 constant={'original27_X_source_columns':piv,'free30_V_indices':C,
  'fixed_bottom30_basis_Xp27_then_Efree3':[srec(v) for v in B30],
  'fixed_bottom30_inverse_columns':[srec(v) for v in B30inv],
  'fixed_bottom30_determinant':str(det30),'fixed27_receipts':rec27,'fixed30_inverse_receipts':invops,
  'remaining_free3_bottom_coordinates':free3,'remaining_free3_originalV_rows':[C[i] for i in free3],
  'selected99_kernel741_columns':selected,'remaining_free3_Y_rows':freeY,
  'remaining642_kernel741_columns':rest,'all642_zero_constant_combinations':constzero,
  'final_minimal6_target_order':['Y_source_rows',freeY,'free30_coordinates',free3]}
 cp=out/'actual27-plus99-unit-complements-and-source-combinations.json';H.dump(cp,constant)
 nodes=[]
 def put(op):nodes.append(op);return len(nodes)-1
 P=put(['external_complete_matrix_operator',sourcepath,source['complete_T2_presentation_node'],[132,1186]])
 first27=[put(['matrix_column',P,418+j,[132]]) for j in piv]
 B1=put(['horizontal_unit_basis',first27,['standard_top102_in_B132'],['standard_bottom_coordinates_in_B132',free3],
  {'column_blocks':[27,102,3],'actual_unknown_top_Yp_retained':True,'unit_constant_bottom30_certificate':cp.name,'unit_constant_determinant':str(det30)}])
 B1inv=put(['formal_regular_unit_inverse',B1,{'shape':[132,132],'constant_unit_proof':'bottom30certificate; top102 identity; no Yp0 query'}])
 pi1=put(['compose',['coordinate_projection',132,list(range(27,132))],B1inv])
 stage1=[]
 for j in range(418):stage1.append(put(['apply_vector',pi1,put(['matrix_column',P,j,[132]])]))
 for rec in kernels:
  v=put(['constant_linear_combination_of_matrix_columns',P,[[418+j,c] for j,c in rec['original768']],[132]])
  stage1.append(put(['apply_vector',pi1,v]))
 P1=put(['matrix_of_columns',stage1,[105,1159],{'exact_cokernel_preserved':'first27pivotcolumns sent tofirst27basisvectors, then quotient'}])
 B2=put(['horizontal_unit_basis',[stage1[418+j] for j in selected],['standard_top102_free_coordinates',freeY],['standard_bottom3_coordinates'],
  {'shape':[105,105],'column_blocks':[99,3,3],'constant_certificate':'acceptedselected99Yminor; originalrowsandpivotproduct inmuprimary','unknown_higher_coefficients_retained':True}])
 B2inv=put(['formal_regular_unit_inverse',B2,{'shape':[105,105],'unit_proof':'nonzero fixed99originalYminor with identityfree6completion'}])
 pi2=put(['compose',['coordinate_projection',105,list(range(99,105))],B2inv])
 remaining=list(range(418))+[418+j for j in rest]
 final=[put(['apply_vector',pi2,stage1[j]]) for j in remaining]
 pres=put(['complete_minimal_honest_T2_presentation',{'shape':[6,1060],'matrix_columns':final,
  'source_order':'418completeG3columns,642remainingXkernelcolumns','target_order':constant['final_minimal6_target_order']}])
 # The constant zero family is used only for first-order source profiling. It
 # is an invertible triangular reparameterization of the642 remainingcolumns.
 proof={'base':'B=(k[ell]/ell²)[[zA0,...,zA6,rA]]',
  'all_original_sources':'all original418G3 and768fullR0X columns are retained throughinvertibleconstant768sourcechange andformalunitquotients',
  'first27_basis':'actual132 sourcepivotcolumns,102topstandardvectors,3bottomcomplementvectors; unknownYp0 is retained andnotsetzero',
  'first27_quotient':'B1^-1 then drop27; no matrix-source scalarhomotopy query',
  'remaining741_origins':'bottom0 andtopactualY0 frommuprimary; unknownYp0 multiplieszero constantbottomcoefficients',
  'next99_basis':'actualstage1selected99columns plus6standardfreecoordinates; nonzerofixedminor givesformalunit',
  'minimality':'G3(0)=0; all741Y0 are inselected99span; all final1060entries vanishatorigin',
  'exact_complete_shape':[6,1060],'T2_length_and_K_exponents_not_computed':True,
  'first_hilbert_source_economy':'on642constantY0kernel, inspectbottomfree3linear fromactualX1first. Onitskernel fullX1isXp0*a1; subtractactualXp*a1, thenXeff startsdegree2 and(MXeff)2=M0Xeff2, soM1drops. OnlyprojectedX2neededforremainingtop3linearpart; pairs/cacheavailabilitynotyetclaimed',
  'no_new_physical_implicit_gauge_N1_or_higher_source_coefficients':True}
 op=out/'actual-honest-T2-minimal6-by1060-factored-source-operator.json';H.dump(op,{'operator_nodes':nodes,'complete_minimal_T2_node':pres,'original132_operator_node':P,'first27_quotient_node':pi1,'next99_quotient_node':pi2,'source_proof':proof})
 slots=sum(map(len,B30+B30inv+Y))+sum(len(r['subtract']) for r in rec27)+sum(len(q['constant741_combination'])+len(q['original768_combination']) for q in constzero)+len(nodes)+10000
 assert slots<150000 and sum(p.stat().st_size for p in out.iterdir())<32*1024*1024
 report={'id':STEM,'all_passed':True,'complete_actual_T2_minimal_presentation_shape':[6,1060],
  'unit_pairs':[27,99],'all_original768_source_columns_retained':True,'all1060_constant_entries_zero':True,
  'new_physical_or_higher_coefficients':False,'actual_T2_length_or_K_exponents_computed':False,
  'input_sha256':hs,'code_sha256':H.sha(Path(__file__)),'note_sha256':H.sha(ROOT/'research/notes'/f'{STEM}.md'),
  'evidence_sha256':{p.name:H.sha(p) for p in [cp,op]},'owned_source_slot_bound':slots,'seconds':time.monotonic()-start}
 H.dump(out/'report.json',report);H.dump(out/'phase.json',{'phase':'FINAL_complete_minimal6_source_operator_no_newcoefficients'});print(json.dumps({k:report[k] for k in ['all_passed','complete_actual_T2_minimal_presentation_shape','owned_source_slot_bound','seconds']}))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);p.add_argument('--mu-peer',required=True);a=p.parse_args();main(a.output.resolve(),a.mu_peer)
