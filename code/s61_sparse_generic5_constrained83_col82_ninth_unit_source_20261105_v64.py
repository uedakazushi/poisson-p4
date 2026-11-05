#!/usr/bin/env python3
"""SOURCE-only ninth pivot, full168x15/RHS and next fixed-column costs."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time,signal,resource,gc
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_NINTH_UNIT_SOURCE_20261105_v64'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def main(a):
 assert a.build
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def ev(p,r,name):
  f=p.parent/name;assert sha(f)==r['evidence_sha256'][name];inputs[f]=sha(f);return read(f)
 p,r=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL62_EIGHTH_UNIT_SOURCE_20261105_v64','0b1d02dc5f3d638c7a201f81ff4cea5cfd997691eb2cb7553b76b8d59a17ba0c');op=ev(p,r,'actual-row52-eighth-unit-and-full169x16-operator.json');oldcost=ev(p,r,'remaining16-literal-source-cost-comparison.json')['all16_fixed_source_tuple_majorants']
 sp,sr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL82_ROW35_EXISTING_UNIT_STRIP_SOURCE_V2_20261105_v64','dd03cd595b076abd7a1098dad61f95be2cb8c0ad4d11b52fa50dbd164213e8e6');strip=ev(sp,sr,'fixed-row35-existing-unit-divisions.json');assert strip['only_existing_unit_factors']=={'Q':1,'R=q0+1':1,'E':0}and strip['actual_Laurent_monomial_shift5']==[0,0,-1,-1,-2];F=strip['reduced_nonzero_polynomial'];assert F['QQ_content']=='1/62208000000'and len(F['integer_terms'])==75;terms=[(tuple(e),int(c))for e,c in F['integer_terms']];assert len({e for e,c in terms})==75 and all(c for e,c in terms)
 cp,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_COEFFICIENTS_20261105_v64','f7f9bc898a54b9a457a59ba4146f38b0e5a8ce9f2cbda906379a71db20e5f9cb');api=ev(cp,cr,'actual-cleared-col82-and-commonQ3RE-representation.json');vals={z['row']:z for z in api['actual16_col82_numerator_values']}
 rp,rr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_COL82_READBACK_FORMULA_REVIEW_V2_20261105_v64','7eb3518c7812c2ef74b03d4d8eb11b4f9912dd66f66d6d3f2eef47139fe61363');assert rr['all4160_QQ_coefficients_staged_formula_equal']
 rows=[i for i in op['selected16_row_ids']if i!=35];cols=[j for j in op['residual16_column_ids']if j!=82];full=[i for i in op['full169_row_ids']if i!=35];top=[i for i in op['top103_row_ids']if i!=35];assert len(rows)==len(cols)==15 and len(full)==168 and len(top)==102
 selected=[]
 for i in rows:
  v=vals[i];f=cp.parent/v['file'];assert sha(f)==v['sha256'];inputs[f]=sha(f);selected.append({'row':i,'actual_N8_i82_path':str(f.relative_to(ROOT)),'sha256':v['sha256'],'typed_multiplier':'62208000000*s0*q0*t0^2*N8_i82/(Q*(q0+1)*F9)'})
 assert sum(vals[i]['terms']for i in rows)==4049
 evidence={'all_passed':True,'source_only':True,'pivot_row':35,'pivot_column':82,'primitive_F9_integer_terms75':F['integer_terms'],'F9_nonzero_witness':{'exponent5':list(terms[0][0]),'coefficient':str(terms[0][1]),'distinct75_monomials':True},'actual_p9':'F9/(62208000000*s0*q0*t0^2*Q^2*E)','typed_inverse':'62208000000*s0*q0*t0^2*Q^2*E/F9','additional_chart_factor':'F9 (the specified nonzero75-term polynomial ONLY)','new_and_all_old_chart_complements_OPEN':True,'no_irreducibility_or_component_coverage_assertion':True,'full168_row_ids':full,'top102_row_ids':top,'selected15_row_ids':rows,'residual15_column_ids':cols,'all87_other_top_RHS':op['all87_other_top_RHS'],'all66_zero_constraints':op['all66_zero_constraints'],'all168x15_RF_roots':'S9_i,j=S8_i,j-lambda_i*S8_35,j, lambda_i=S8_i82/p9; all original full callbacks retained.','RHS':'h9_i=h8_i-lambda_i*h8_35; ALL169 incoming RHS retained by inverse substitution and all168 residual RHS.','pivot_row_NOT_raw_or_row58':'The WHOLE row35 and h8_35 remain the complete eighth-stage callbacks; only col82 was computed. No identification with row58 from one equal cell.','actual15_selected_lambda':selected,'other87_lambda_UNEVALUATED':'Use actual formal S8_i82/p9 callbacks; no binary value, zero, or missing coefficient is invented.','zero66_lambda':'Formal zero coefficient rows preserved; RHS consistency still retained.','substitution':'x82=(h8_35-sum_j S8_35,j*xj)/p9, then restore previous eight pivots, one P290 scale, ten constraint units and49 old units.','full_cokernel':'full169x16 unit-equivalent tofull168x15; full168coker=top102coker direct_sumK5^66. Selected15square is NOT samecokernel.','future_common_representation':'D9*S9_i,j=P7*F9*N8_i,j-62208000000*s0*q0*t0^2*N8_i82*N8_35,j; P7=Q*(q0+1), D9=Q^4*(q0+1)^2*E*F9. Representation ONLY.','all_other14_columns_stay_typed_no_updates':True,'P290onlyonce':True,'no_new_matrix_cell_rank_inverse_or_exhaustion_claim':True}
 (out/'actual-row35-ninth-unit-and-full168x15-operator.json').write_bytes(enc(evidence))
 # Fixed-column source costs ONLY, using the already frozen scalar supports.
 ep,er=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');ME={(z['row'],z['column']):z for z in ev(ep,er,'actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 dp,dr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_ROWS11_44_LITERAL_SOURCE_20261105_v64','0ca442b4eb828e535bdb527357826c6f4cd7edb48f23735f08742b798ada9371');dv=ev(dp,dr,'whole-cached-rows11-44-exact-comparison.json');parse=25000+numeric(ME)+numeric(dv)+numeric(evidence)+numeric(api);assert parse<150000;nd={z['column']:len(z['difference_terms'])for z in dv['all24_complete_actual_row_comparisons']};del dv
 tp,tr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_EIGHTH_UNIT_COL82_SOURCE_20261105_v64','6925f8f7072565b98804d3c87262a5a20b5bdc05cc5f3171884f6f73d97bfdfb');idx=ev(tp,tr,'row-stream-input-and-plan-index16.json');ref=next(z for z in idx['row_records16']if z['row']==35);f=ROOT/ref['row_operand_file'];assert sha(f)==ref['row_operand_sha256'];inputs[f]=sha(f);b=read(f)['row_operands'];n={k:len(v['original_literal']['terms'])for k,v in b.items()}
 def m2(i,j):return ME[i,j]['terms']+ME[i,49]['terms']*(ME[44,j]['terms']+nd[j])+ME[i,51]['terms']*nd[j]
 def c6(j):return 3*m2(35,j)+2*n['U35']*m2(31,j)+3*n['U35']*m2(83,j)+3*n['V35']*m2(83,j)+2*n['C35']*(m2(14,j)+m2(35,j))+n['F35']*m2(104,j)
 def c8(j):return 3*(3*c6(j)+n['N35']*(ME[51,j]['terms']+2*(ME[11,j]['terms']+ME[44,j]['terms'])))+n['M35']*ME[52,j]['terms']
 assert c8(82)==ref['source_tuple_majorant'];costs=[]
 for z in oldcost:
  j=z['column']
  if j==82:continue
  row35=c8(j);other=z['new16_fixed_source_tuple_majorant']-row35;assert other>=0;co=225*other+4049*row35;costs.append({'column':j,'eighth_source_total_excluding35':other,'unevaluated_row35_eighth_source_cost':row35,'next15_literal_source_tuple_majorant':co})
 costs.sort(key=lambda z:(z['next15_literal_source_tuple_majorant'],z['column']));assert len(costs)==15
 (out/'remaining15-literal-source-cost-comparison.json').write_bytes(enc({'all_passed':True,'source_only':True,'all15_fixed_source_costs':costs,'minimum_in_this_exact_source_recipe':costs[0]['column'],'no_matrix_cell_computation':True,'new_support_height_lifetime_RAW_gate_REQUIRED_before_future_coefficients':True}));report={'id':ID,'all_passed':True,'source_only':True,'ninth_actual_unit_on_F9_chart':True,'F9_terms':75,'full_residual_shape':[168,15],'selected_subblock_shape':[15,15],'next_cheapest_source_column':costs[0]['column'],'next_literal_tuple_majorant':costs[0]['next15_literal_source_tuple_majorant'],'source_numeric_upper':parse,'actual_READ_RAW_bytes':sum(p.stat().st_size for p in inputs),'all169_RHS_87plus66_other14_oneP290_and_complements_kept':True,'no_new_coefficient_or_matrix_update':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
