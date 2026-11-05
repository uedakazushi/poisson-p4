#!/usr/bin/env python3
"""SOURCE-only eighth typed pivot and next-column literal cost."""
from pathlib import Path
from fractions import Fraction as Q
import argparse,json,hashlib,time,signal,resource
from s61_sparse_generic5_primitive_binary_codec_20261105_v64 import decode
ROOT=Path(__file__).resolve().parents[1];ID='S61_SPARSE_GENERIC5_CONSTRAINED83_COL62_EIGHTH_UNIT_SOURCE_20261105_v64'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(65536),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def enc(v):return(json.dumps(v,separators=(',',':'))+'\n').encode()
def numeric(v):
 if isinstance(v,bool):return 0
 if isinstance(v,(int,float)):return 1
 if isinstance(v,list):return sum(map(numeric,v))
 if isinstance(v,dict):return sum(map(numeric,v.values()))
 return 0
def poly(p):return{tuple(p['shift'][j]+t[j+1]for j in range(5)):Q(p['content'])*int(t[0])for t in p['terms']}
def mul(A,B):
 R={}
 for e,q in A.items():
  for f,z in B.items():
   g=tuple(x+y for x,y in zip(e,f));R[g]=R.get(g,Q(0))+q*z
 return{e:q for e,q in R.items()if q}
def scale(A,q,shift=(0,0,0,0,0)):return{tuple(e[j]+shift[j]for j in range(5)):q*z for e,z in A.items()}
def main(a):
 assert a.build
 start=time.monotonic();signal.alarm(60);resource.setrlimit(resource.RLIMIT_AS,(2000000000,2000000000));out=ROOT/'reports'/ID;out.mkdir(exist_ok=False);inputs={Path(__file__):sha(Path(__file__)),ROOT/f'research/notes/{ID}.md':sha(ROOT/f'research/notes/{ID}.md'),ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py':sha(ROOT/'code/s61_sparse_generic5_primitive_binary_codec_20261105_v64.py')}
 def gate(stem,h):
  p=ROOT/f'reports/{stem}/report.json';assert sha(p)==h;r=read(p);assert r['all_passed'];inputs[p]=h;return p,r
 def ev(p,r,name):
  v=p.parent/name;assert sha(v)==r['evidence_sha256'][name];inputs[v]=sha(v);return read(v)
 p,s7=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_COL74_SEVENTH_UNIT_SOURCE_20261105_v64','23dc07db98639cc1a418192064e922d6c467ea8003b11b4b4dc55113dc1c5319');op=ev(p,s7,'actual-row51-seventh-unit-and-full170x17-operator.json');cost7=ev(p,s7,'remaining17-literal-source-cost-comparison.json')['all17_fixed_source_tuple_majorants']
 p,s=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_UNIT_COL62_SOURCE_20261105_v64','f76a6096796dc0917773d963c9a1df89e7826880f785126b1411399796e6a99c');common=ev(p,s,'common-source-operands.json')['common_operands'];row52=ev(p,s,'rows/row-52-source-operands.json')['row_operands']
 ap,ar=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_IMPLEMENTATION_GATE_20261105_v64','08741e886fec1f973698f201efa13371289a555157f5dda9fab717f30ac9bd87');ac=ev(ap,ar,'implementation-custody-storage-bit-addendum.json');assert ac['source_custody_counts']=={'raw':66,'explicit':3,'delta':1,'reuse':68,'N':17};del ac
 cp,cr=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_SEVENTH_COL62_COEFFICIENTS_20261105_v64','5d1637a8014757fe187a2457159112409106adde3795d2360fce7e0a03042b5e');api=ev(cp,cr,'actual-cleared-col62-and-commonQ3R-representation.json');values={z['row']:z for z in api['actual17_col62_numerator_values']};v=values[52];file=cp.parent/v['file'];assert sha(file)==v['sha256'];inputs[file]=sha(file);blob=file.read_bytes();assert blob[:9]==b'Q5CELL1\n\x01';literal=decode(blob[9:]);assert len(literal['terms'])==11
 Qp=poly(common['Q']['original_literal']);Q2=poly(common['Q2']['original_literal']);P7=poly(common['P7']['original_literal']);D7=mul(Q2,P7);E={(0,0,1,0,0):Q(32),(0,0,0,1,0):Q(15),(0,0,0,0,0):Q(35)};pivot=scale(E,Q(32768,273375),(3,0,0,0,-2));rawpivot=poly(row52['A52_62']['original_literal']);assert rawpivot==pivot and poly(literal)==mul(D7,pivot)
 for label in('A52_49','A52_51','U52','V52','C52','F52','N52'):assert not row52[label]['original_literal']['terms']
 ep,er=gate('S61_SPARSE_GENERIC5_CONSTRAINED83_FIXED24_INVERSE_SOURCE_PLAN_20261105_v64','a904b75b0c72e4831b2f4cb5bf99e85df30b2a4a15e8ffdc334d291da0421b4a');ME={(z['row'],z['column']):z for z in ev(ep,er,'actual-fixed24-entry-source-and-exponent-plan.json')['actual_entry_sources576']}
 assert row52['A52_62']['source']==ME[52,62]['source']and ME[52,62]['terms']==3
 rows=[i for i in op['selected17_row_ids']if i!=52];cols=[j for j in op['residual17_column_ids']if j!=62];full=[i for i in op['full170_row_ids']if i!=52];top=[i for i in op['top104_row_ids']if i!=52];assert len(rows)==len(cols)==16 and len(full)==169 and len(top)==103
 nf=sum(values[i]['terms']for i in rows);assert nf==8784;costs=[]
 for z in cost7:
  j=z['column']
  if j==62:continue
  n=ME[52,j]['terms'];previous=z['new17_fixed_source_tuple_majorant']-9*n;c=3*previous+nf*n;assert previous>=0;costs.append({'column':j,'previous17_source_total_excluding52':previous,'raw_row52_terms':n,'new16_fixed_source_tuple_majorant':c})
 costs.sort(key=lambda z:(z['new16_fixed_source_tuple_majorant'],z['column']));assert len(costs)==16 and costs[0]['column']==82
 evidence={'all_passed':True,'source_only':True,'pivot_row':52,'pivot_column':62,'actual11_term_numerator':literal,'literal_D7_times_raw3term_pivot_identity':True,'actual_S7_52_62':'32768*H^3*(32*s0+15*q0+35)/(273375*t0^2)','inverse':'273375*t0^2/(32768*H^3*(32*s0+15*q0+35))','additional_chart_factor':'E=32*s0+15*q0+35','new_and_all_previous_chart_complements_OPEN':True,'accepted_zero_A52_49_A52_51_U52_V52_C52_F52_N52_proves_fullrow_RHS_alias':'S7_52,j=A52,j and h7_52=h0_52. Here A is the old full177 source with P290 already applied once to its four columns.','full169_row_ids':full,'top103_row_ids':top,'selected16_row_ids':rows,'residual16_column_ids':cols,'all87_other_top_RHS':op['all87_other_top_RHS'],'all66_zero_constraints':op['all66_zero_constraints'],'all169x16_RF_roots':'S8_i,j=S7_i,j-(S7_i62/p8)*A52,j; all source callbacks retained.','RHS':'h8_i=h7_i-(S7_i62/p8)*h0_52, all87+66 retained.','substitution':'x62=(h0_52-sum_j A52,j*xj)/p8, then restore prior seven pivots, one P290 scale, ten constraints and49 units.','cokernel':'full170x17 unit-equivalent tofull169x16; full169coker=top103coker direct_sum K5^66. Selected16 square is NOT samecokernel.','future_common_representation':'D7*E*S8_i,j=E*N7_i,j-J8*M7_i*A52,j; D7=Q^3*(q0+1), M7_i=accepted numerator of S7_i62, J8=273375*t0^2/(32768*H^3). Representation ONLY.','all_other15_columns_unchanged_typed':True,'P290onlyonce':True,'no_new_next_column_values_rank_inverse_or_exhaustion_claim':True}
 (out/'actual-row52-eighth-unit-and-full169x16-operator.json').write_bytes(enc(evidence));(out/'remaining16-literal-source-cost-comparison.json').write_bytes(enc({'all_passed':True,'source_only':True,'all16_fixed_source_tuple_majorants':costs,'selected_minimum_in_this_frozen_source_model':82,'only_metadata_no_new_QQ_matrix_cells':True,'future_support_height_custody_storage_gate_REQUIRED_before_values':True}));bound=25000+numeric(ME)+numeric(common)+numeric(row52)+numeric(op)+numeric(api)+numeric(cost7)+numeric(evidence)+numeric(costs)+4096;assert bound<150000
 report={'id':ID,'all_passed':True,'source_only':True,'eighth_actual_unit_factor_verified':True,'additional_chart_factor':'32*s0+15*q0+35','full_residual_shape':[169,16],'selected_subblock_shape':[16,16],'next_source_cheapest_column':82,'next_source_tuple_majorant':costs[0]['new16_fixed_source_tuple_majorant'],'source_numeric_upper':bound,'actual_READ_RAW_bytes':sum(p.stat().st_size for p in inputs),'no_new_column_values':True,'input_sha256':{str(p.relative_to(ROOT)):h for p,h in inputs.items()},'evidence_sha256':{p.name:sha(p)for p in out.iterdir()},'seconds':time.monotonic()-start};(out/'report.json').write_bytes(enc(report));signal.alarm(0);print(json.dumps({k:v for k,v in report.items()if k not in('input_sha256','evidence_sha256')}));print('report_sha256',sha(out/'report.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',action='store_true');main(p.parse_args())
