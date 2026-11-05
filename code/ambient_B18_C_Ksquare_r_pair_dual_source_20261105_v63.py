#!/usr/bin/env python3
"""Unpaired SOURCE r-pair quartic-dual tensor; no pair input substituted."""
import hashlib,json,signal,time
from fractions import Fraction as Q
from pathlib import Path
import ambient_B18_C_genuine_projection7_source_v2_20261105_v63 as S
import ambient_B18_C_Ksquare_leading_pair_source_20261105_v63 as L
ROOT=Path(__file__).resolve().parents[1]
STEM='AMBIENT_B18_C_KSQUARE_R_PAIR_DUAL_SOURCE_20261105_v63'
PATHS={
 'leading':'reports/AMBIENT_B18_C_KSQUARE_LEADING_PAIR_SOURCE_20261105_v63/unpaired-leading-three-pair-extraction-tensor-and-rank2-sufficient-criterion.json',
 'linear_zero':'reports/AMBIENT_B18_C_KSQUARE_LEADING_LINEAR_SPAN_SOURCE_20261105_v63/unpaired-known174-membership-and-pure-r-H2-source.json',
 'point':L.PATHS['point'],'Pi':L.PATHS['Pi'],'zero3':L.PATHS['zero3'],
 'quadbasis':'reports/AMBIENT_B18_HONEST_NIL_T_QUADRATIC_CLOSURE_20261104_v63/known-quadratic-closure-and-complete814-source-projection.json',
 'rightsource':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_SOURCE_PREFLIGHT_20261105_v63/right106-full-first-span-and-cache-footprint.json',
 'rightvalues':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_VALUES_20261105_v63/phase1-all35-nine-convolution-projections.json',
 'rightpeer':'reports/AMBIENT_B18_HONEST_NIL_T_RIGHT106_20261105_v63_REVIEW/report.json',
 'mu_source':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63/complete-twelve-generator-kernel-and-two-cache-callback-gate.json',
 'mu_source_peer':'reports/AMBIENT_B18_NIL_R0_MINIMAL_GENERATORS_SOURCE_20261105_v63_REVIEW/report.json',
 'quartic':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63/actual-eight-relations-and-complete-quartic-dual.json',
 'quartic_peer':'reports/AMBIENT_B18_NBAR_CUBE_QUARTIC_DUAL_VALUES_20261105_v63_REVIEW/report.json',
 'Loewy':'reports/AMBIENT_B18_TNIL_EXACT_LOEWY16_20261105_v63/report.json',
}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 t=time.monotonic();signal.alarm(60);pins={p:sha(ROOT/p)for p in PATHS.values()};receipts=[]
 fields={
 'leading':['Wq44_sparse','Wbeta44_sparse','linear_B_h0','fraction_majorant_before_arithmetic'],
 'linear_zero':['all324_quotient6_zero'],
 'point':['source_m0_24_by44'],
 'Pi':['constant_minimal_projection5_by128_sparse'],
 'zero3':['actual_three_minimal_first_columns40'],
 'quadbasis':['fixed_rank_pivot_columns','fixed_rank_pivot_rows','fixed_minor_inverse','all280_source_labels_relation_and_variable'],
 'rightsource':['selected35_actual_first_relations','all_nonzero_convolution_duals','all_actual_effective_C1_H2_callbacks'],
 'rightvalues':['all35_projected_M2_lower_representatives','all35_source_input'],
 'mu_source':['callback_means','individual_callback_U_sensitivity_zero_reason','individual_callback_V_sensitivity_zero_reason'],
 'quartic':['actual_complete_quartic_dual1650_sparse'],
 }
 # The full174 inverse is decoded BEFORE other arithmetic buffers, then
 # only the four needed rows survive. No dense174 Fraction copy is made.
 d={};q=S.selective(ROOT/PATHS['quadbasis'],fields['quadbasis'],receipts)
 parse_peak=S.count(q)+10000;assert parse_peak<150000
 cp=q['fixed_rank_pivot_columns'];rr=q['fixed_rank_pivot_rows'];labels=q['all280_source_labels_relation_and_variable'];assert len(rr)==174
 source_rel=[labels[8*i][0]for i in range(35)]
 selected_ids=[623,627,673,675]
 selected_positions=[source_rel.index(x)for x in selected_ids]
 needed={i:cp.index(8*i+7)for i in selected_positions}
 inverse_rows={i:q['fixed_minor_inverse'][j]for i,j in needed.items()}
 # Copy only selected rows to immutable sparse QQ lists, releasing the
 # original30276 rational strings before other fields are decoded.
 inverse_rows={i:[[rr[k],c]for k,c in enumerate(row)if Q(c)]for i,row in inverse_rows.items()}
 del q
 d['quadbasis']={'r_coefficient_inverse_rows35_by180_sparse':inverse_rows,'selected_source_relations35':source_rel}
 for k,p in PATHS.items():
  if k=='quadbasis':continue
  d[k]=S.selective(ROOT/p,fields[k],receipts)if k in fields else json.loads((ROOT/p).read_text())
  if k.endswith('peer'):assert d[k]['all_passed']
 assert d['linear_zero']['all324_quotient6_zero']
 assert d['quartic']['actual_complete_quartic_dual1650_sparse']==[[329,'1']]
 assert d['Loewy']['exact_actual_HF']==[5,5,3,2,1]
 assert all(not Q(col.get('7','0'))for col in d['zero3']['actual_three_minimal_first_columns40'])
 assert d['rightsource']['selected35_actual_first_relations']==source_rel
 assert d['rightsource']['all_nonzero_convolution_duals'][0]==[0,7,[[35,'1']]]
 eff=next(v[2]for v in d['rightsource']['all_actual_effective_C1_H2_callbacks']if v[:2]==[0,7]);assert eff==[[27,7,'-8/3']]
 q35=[Q(row[0])for row in d['rightvalues']['all35_projected_M2_lower_representatives']]
 assert [source_rel[i]for i,x in enumerate(q35)if x]==selected_ids
 assert [r['source_relation']for r in d['rightvalues']['all35_source_input']]==source_rel
 # Bound input-leaf/common-denominator operations before SOURCE folds.
 max_input=max([L.ht(Q(c))for i,row in inverse_rows.items()for j,c in row]+[L.ht(c)for c in q35])
 sourceH=d['leading']['fraction_majorant_before_arithmetic']['H']+2*max_input+30
 assert sourceH<384 and 2*sourceH+1<1153
 m=[[Q(c)for c in row]for row in d['point']['source_m0_24_by44']]
 P=[dict((i,Q(c))for i,c in row)for row in d['Pi']['constant_minimal_projection5_by128_sparse']]
 Wq=dict((i,Q(c))for i,c in d['leading']['Wq44_sparse']);Wb=dict((i,Q(c))for i,c in d['leading']['Wbeta44_sparse']);beta=dict((tuple(e),Q(c))for e,c in d['leading']['linear_B_h0'])
 Drows={i:{j:Q(c)for j,c in row}for i,row in inverse_rows.items()}
 tensor=[];witness_hash=hashlib.sha256()
 # Only zz monomials occur in L_U L_V: no unprovided relation removes a
 # pure-r coefficient, the actual input-cache first r coefficients are0.
 for mi,e in enumerate(L.M2):
  if e[7]:continue
  for qvar in range(9):
   raw={}
   if qvar<8:
    for j,c in Wq.items():L.put(raw,(j,L.shift(e,qvar)),c)
   else:
    for f,b in beta.items():
     ex=tuple(x+y for x,y in zip(e,f))
     for j,c in Wb.items():L.put(raw,(j,ex),L.mul(b,c))
   vv={}
   for (j,ex),c in raw.items():
    for i in range(24):
     if m[i][j]:L.put(vv,(i,ex),L.mul(m[i][j],c))
   hh=L.H3(vv);mini={}
   for (j,ex),c in hh.items():
    for target in range(5):
     if 44+j in P[target]:L.put(mini,36*target+L.I2[ex],L.mul(P[target][44+j],c))
   carry=L.mul(Q(-8,3),hh.get((27,L.mon([7,7])),Q(0)))
   linear_clear=Q(0);cofs=[]
   for i,row in Drows.items():
    u=Q(0)
    for j,c in mini.items():
     if j in row:u=L.add(u,L.mul(row[j],c))
    cofs.append([source_rel[i],str(u)])
    linear_clear=L.add(linear_clear,L.mul(q35[i],u))
   total=L.add(carry,-linear_clear)
   tensor.append([mi,qvar,str(carry),str(linear_clear),str(total)])
   witness_hash.update((json.dumps(cofs,separators=(',',':'))+'\n').encode())
 allzero=all(row[4]=='0'for row in tensor)
 evidence={
 'scope':'UNPAIRED SOURCE. No LU/LS or actualq0..7 substituted. CachedG2first-direction valuesenterasoperatorcoefficients, notnewpairresponses.',
 'input_sha256':pins,'all252_unpaired_r_pair_dual_coefficients':tensor,'all252_r_dual_coefficients_zero':allzero,
 'fixed174_four_r_lift_rows_sparse':inverse_rows,'four_lift_coefficient_stream_sha256':witness_hash.hexdigest(),
 'actual35_cached_M2_head0_r2_values':[str(x)for x in q35],
 'literal_C1_on_this_raw_input':eff,
 'domain_certificate_needed':'raw128 degree2=(free44zero,H2). First24and18 originunit cancellations preservezero free complement andH2; theaccepted effectiveC1 callback thenapplies toH2 withnooriginfreecleanup. No genericinputdomain identificationisassumed.',
 'exact_actual_relation_clearing':'w2=M1*u1 inknown174; subtractTRUEallorderM*u1. Atpure-r cubic head0, onlyu1_r*M2_head0r2 enters. r timescorrectedrepresentativestartsdegree4, hencecompletephi4readsthisexactcoefficient. M3/Gamma2/section1/m1/Pi1-otherdirections doNOTenterthisspecificdual; actualremainingcubicclass maystilldependonthem.',
 'raw42_cubic_pure_r_zero':'L_U,L_S havezero rcoefficients; Product44c3_42=2L_U L_V q*v42_b6, soitsr3coefficientzero. Thee43cross termneverentersraw42.',
 'conditional_if_allzero':'Allthreeactualr(pair)classes vanish. Acceptedq2ann+Jclosure thenr*Kbar²=0. SinceKbar²⊂m3Tnil andr:m3→m4rank1, imageKbar²hasQQdimension≤2; generationofH0bar3byK²isimpossible. NooriginalH0type/nilindexinference.',
 'conditional_if_nonzero':'A nonzerooperatorcoefficientdoesnotmeanactualproductnonzero. MustbindactualqfirstandLU/LS,freshSOURCE/valuegatebeforepairresponses.',
 'source_fraction_bound':{'max_input_leaf':max_input,'conservative_stored':sourceH,'pre':2*sourceH+1},
 'source_max_stored_H':L.STORED,'source_max_pre_H':L.PRE,'all_new_product_values_UNCOMPUTED':True,'selective_receipts':receipts,
 }
 owned=max(parse_peak,S.count(d)+S.count(inverse_rows)+S.count(tensor)+S.count(evidence)+25000);assert owned<150000;evidence['whole_SOURCE_owned_upper']=owned
 out=ROOT/'reports'/STEM;out.mkdir(exist_ok=False);p=out/'unpaired-r-pair-quartic-dual-with-coherent-linear-clearing-source.json';p.write_text(json.dumps(evidence,indent=2)+'\n')
 rep={'id':STEM,'all_passed':True,'status':'SOURCE_UNPAIRED_R_PAIR_DUAL_DOMAIN_PENDING_INDEPENDENT_CHECK','all252_r_dual_coefficients_zero':allzero,'whole_SOURCE_owned_upper':owned,'code_sha256':sha(Path(__file__)),'note_sha256':sha(ROOT/'research/notes'/(STEM+'.md')),'input_sha256':pins,'evidence_sha256':{p.name:sha(p)},'seconds':time.monotonic()-t}
 (out/'report.json').write_text(json.dumps(rep,indent=2)+'\n');print(json.dumps(rep));signal.alarm(0)
if __name__=='__main__':main()
